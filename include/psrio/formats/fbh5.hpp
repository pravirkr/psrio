#pragma once

#include "psrio/astro.hpp"
#include "psrio/common/types.hpp"
#include "psrio/detail/exceptions.hpp"
#include "psrio/detail/skip.hpp"
#include "psrio/packed.hpp"

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <format>
#include <highfive/highfive.hpp>
#include <limits>
#include <memory>
#include <span>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace psrio::formats::fbh5 {

namespace detail {

inline bool equals_ignore_case(std::string_view left, std::string_view right) {
    if (left.size() != right.size()) {
        return false;
    }
    for (std::size_t index = 0; index < left.size(); ++index) {
        const auto a = static_cast<unsigned char>(left[index]);
        const auto b = static_cast<unsigned char>(right[index]);
        if (std::tolower(a) != std::tolower(b)) {
            return false;
        }
    }
    return true;
}

inline std::string trim_trailing(std::string value) {
    while (!value.empty() && (value.back() == '\0' || value.back() == ' ')) {
        value.pop_back();
    }
    return value;
}

inline bool has_attr(const HighFive::DataSet& dataset,
                     const HighFive::File& file,
                     const std::string& name) {
    return dataset.hasAttribute(name) || file.hasAttribute(name);
}

inline HighFive::Attribute get_attr(const HighFive::DataSet& dataset,
                                    const HighFive::File& file,
                                    const std::string& name) {
    if (dataset.hasAttribute(name)) {
        return dataset.getAttribute(name);
    }
    if (file.hasAttribute(name)) {
        return file.getAttribute(name);
    }
    throw FormatError(std::format(
        "psrio: HDF5 filterbank missing required attribute '{}'", name));
}

inline double read_number(const HighFive::Attribute& attr,
                          std::string_view name) {
    if (attr.getSpace().getElementCount() != 1) {
        throw FormatError(std::format(
            "psrio: HDF5 attribute '{}' must be a numeric scalar", name));
    }
    const auto type = attr.getDataType();
    const auto cls  = type.getClass();
    if (cls != HighFive::DataTypeClass::Float &&
        cls != HighFive::DataTypeClass::Integer) {
        throw FormatError(std::format(
            "psrio: HDF5 attribute '{}' is not a numeric scalar", name));
    }
    double value = 0.0;
    try {
        attr.read(value);
    } catch (...) {
        throw FormatError(
            std::format("psrio: HDF5 attribute '{}' could not be read", name));
    }
    return value;
}

inline int read_int(const HighFive::Attribute& attr, std::string_view name) {
    const double value = read_number(attr, name);
    if (!std::isfinite(value) || std::trunc(value) != value ||
        value > static_cast<double>(std::numeric_limits<int>::max()) ||
        value < static_cast<double>(std::numeric_limits<int>::min())) {
        throw FormatError(
            std::format("psrio: HDF5 attribute '{}' is not an integer", name));
    }
    return static_cast<int>(value);
}

inline std::string read_text(const HighFive::Attribute& attr,
                             std::string_view name) {
    const auto type = attr.getDataType();
    if (type.getClass() != HighFive::DataTypeClass::String) {
        throw FormatError(
            std::format("psrio: HDF5 attribute '{}' is not a string", name));
    }
    std::string text;
    try {
        attr.read(text);
    } catch (...) {
        throw FormatError(std::format(
            "psrio: HDF5 attribute '{}' could not be read as a string", name));
    }
    return trim_trailing(std::move(text));
}

inline void require_class(const HighFive::File& file) {
    if (!file.hasAttribute("CLASS")) {
        throw FormatError(
            "psrio: HDF5 filterbank missing required attribute 'CLASS'");
    }
    const auto text = read_text(file.getAttribute("CLASS"), "CLASS");
    if (text != "FILTERBANK") {
        throw FormatError(std::format("psrio: expected HDF5 CLASS attribute to "
                                      "be 'FILTERBANK' but saw '{}'",
                                      text));
    }
    if (!file.hasAttribute("VERSION")) {
        throw FormatError(
            "psrio: HDF5 filterbank missing required attribute 'VERSION'");
    }
}

inline std::string read_version(const HighFive::File& file) {
    if (!file.hasAttribute("VERSION")) {
        throw FormatError(
            "psrio: HDF5 filterbank missing required attribute 'VERSION'");
    }
    const auto attr = file.getAttribute("VERSION");
    const auto cls  = attr.getDataType().getClass();
    if (cls == HighFive::DataTypeClass::String) {
        const auto text = read_text(attr, "VERSION");
        if (text.empty()) {
            throw FormatError(
                "psrio: HDF5 filterbank VERSION attribute is empty");
        }
        return text;
    }
    if (cls == HighFive::DataTypeClass::Float ||
        cls == HighFive::DataTypeClass::Integer) {
        return std::format("{}", read_number(attr, "VERSION"));
    }
    throw FormatError(
        "psrio: HDF5 filterbank VERSION attribute is not a string or number");
}

inline SampleType sample_type_from_dataset(const HighFive::DataSet& dataset) {
    const auto type = dataset.getDataType();
    const auto cls  = type.getClass();
    const auto size = type.getSize();
    if (cls == HighFive::DataTypeClass::Integer) {
        if (H5Tget_sign(type.getId()) != H5T_SGN_NONE) {
            throw FormatError(
                "psrio: unsupported HDF5 sample type: signed integers are not "
                "read (expected uint8, uint16, uint32, or float32)");
        }
        if (size == 1U) {
            return SampleType::kUInt8;
        }
        if (size == 2U) {
            return SampleType::kUInt16;
        }
        if (size == 4U) {
            return SampleType::kUInt32;
        }
    } else if (cls == HighFive::DataTypeClass::Float && size == sizeof(float)) {
        return SampleType::kFloat32;
    }
    throw FormatError(
        "psrio: unsupported HDF5 sample type: expected uint8, uint16, uint32, "
        "or float32");
}

inline bool nbits_matches(SampleType type, int nbits) {
    switch (type) {
    case SampleType::kUInt8:
        return nbits == 8;
    case SampleType::kUInt16:
        return nbits == 16;
    case SampleType::kUInt32:
        return nbits == 32;
    case SampleType::kFloat32:
        return nbits == 32 || nbits == -32;
    case SampleType::kUInt1:
    case SampleType::kUInt2:
    case SampleType::kUInt4:
    case SampleType::kInt8:
        break;
    }
    return false;
}

inline void ensure_filters_available(const HighFive::DataSet& dataset) {
    const auto plist   = dataset.getCreatePropertyList();
    const auto pid     = plist.getId();
    const int nfilters = H5Pget_nfilters(pid);
    if (nfilters < 0) {
        throw FormatError("psrio: HDF5 filterbank could not read the dataset "
                          "filter pipeline");
    }
    for (int index = 0; index < nfilters; ++index) {
        unsigned int flags         = 0;
        std::size_t cd_nelmts      = 0;
        unsigned int filter_config = 0;
        std::array<char, 128> name{};
        const auto id = H5Pget_filter2(
            pid, static_cast<unsigned>(index), &flags, &cd_nelmts, nullptr,
            name.size() - 1U, name.data(), &filter_config);
        if (id < 0) {
            throw FormatError(
                "psrio: HDF5 filterbank could not identify a dataset filter");
        }
        if (H5Zfilter_avail(id) <= 0) {
            const std::string label =
                name[0] != '\0' ? std::string{name.data()} : std::to_string(id);
            throw FormatError(std::format(
                "psrio: unsupported compression: HDF5 filter '{}' ({}) is not "
                "registered (bitshuffle is not linked)",
                label, static_cast<int>(id)));
        }
    }
}

} // namespace detail

/**
 * @brief Header fields blimpy stores on a Breakthrough Listen filterbank.
 *
 * Frequencies are MHz, `tsamp` is seconds, and `tstart` is MJD. `foff` keeps
 * its sign. `nifs` and `nbeams` describe the file. A reader streams one beam.
 */
struct FilterbankHeader {
    std::string source_name;
    std::string version;
    std::string src_raj_text;
    std::string src_dej_text;
    double src_raj         = 0.0;
    double src_dej         = 0.0;
    double az_start        = 0.0;
    double za_start        = 0.0;
    double tstart          = 0.0;
    double tsamp           = 0.0;
    double fch1            = 0.0;
    double foff            = 0.0;
    std::uint64_t nchans   = 0;
    std::uint64_t nifs     = 0;
    std::uint64_t nbeams   = 0;
    std::uint64_t nsamples = 0;
    int nbits              = 0;
    int telescope_id       = 0;
    int machine_id         = 0;
    int data_type          = 0;
    int ibeam              = 0;
    int barycentric        = 0;
    int pulsarcentric      = 0;
    bool has_src_raj       = false;
    bool has_src_dej       = false;
    bool has_telescope_id  = false;
};

/// True when @p path ends in `.h5` or `.hdf5`, or its bytes are an HDF5 file.
[[nodiscard]] inline bool
is_hdf5_filterbank(const std::filesystem::path& path) {
    const auto extension = path.extension().string();
    if (detail::equals_ignore_case(extension, ".h5") ||
        detail::equals_ignore_case(extension, ".hdf5")) {
        return true;
    }
    if (!std::filesystem::exists(path)) {
        return false;
    }
    try {
        HighFive::SilenceHDF5 silence;
        return H5Fis_hdf5(path.c_str()) > 0;
    } catch (...) {
        return false;
    }
}

/**
 * @brief One beam of a Breakthrough Listen HDF5 filterbank.
 *
 * `data` is rank 3: time, beam, frequency. `read_block` copies one hyperslab
 * for the selected beam. Channel 0 is `fch1`, stepping by signed `foff`.
 * Negative `foff` is left in that order. The `mask` dataset is not read.
 * This type models `concepts::BlockReader`.
 */
class FilterbankFile {
public:
    explicit FilterbankFile(const std::filesystem::path& path,
                            std::uint64_t beam_index = 0);

    FilterbankFile(const FilterbankFile&)                = delete;
    FilterbankFile& operator=(const FilterbankFile&)     = delete;
    FilterbankFile(FilterbankFile&&) noexcept            = default;
    FilterbankFile& operator=(FilterbankFile&&) noexcept = default;
    ~FilterbankFile()                                    = default;

    [[nodiscard]] const FilterbankHeader& header() const noexcept {
        return m_header;
    }

    [[nodiscard]] std::uint64_t nchans() const noexcept {
        return m_header.nchans;
    }

    /// The stream is one selected beam, so the block has a single IF.
    [[nodiscard]] std::uint64_t nifs() const noexcept { return 1; }

    [[nodiscard]] int nbits() const noexcept { return m_header.nbits; }

    [[nodiscard]] SampleType sample_type() const noexcept {
        return m_sample_type;
    }

    [[nodiscard]] std::uint64_t bytes_per_sample() const noexcept {
        return m_header.nchans * element_bytes();
    }

    [[nodiscard]] std::uint64_t nsamples() const noexcept {
        return m_header.nsamples;
    }

    [[nodiscard]] bool has_nsamples() const noexcept { return true; }

    [[nodiscard]] double tsamp() const noexcept { return m_header.tsamp; }
    [[nodiscard]] double tstart() const noexcept { return m_header.tstart; }
    [[nodiscard]] std::string_view source_name() const noexcept {
        return m_header.source_name;
    }

    [[nodiscard]] std::string_view telescope() const noexcept {
        return "Unknown";
    }

    [[nodiscard]] double raj() const noexcept { return m_header.src_raj; }
    [[nodiscard]] double dej() const noexcept { return m_header.src_dej; }

    [[nodiscard]] double fch1() const noexcept { return m_header.fch1; }
    [[nodiscard]] double foff() const noexcept { return m_header.foff; }

    [[nodiscard]] double bandwidth() const noexcept {
        return std::abs(m_header.foff) * static_cast<double>(m_header.nchans);
    }

    [[nodiscard]] double center_frequency() const noexcept {
        return m_header.fch1 +
               (m_header.foff * (static_cast<double>(m_header.nchans) - 1.0) /
                2.0);
    }

    [[nodiscard]] int beam() const noexcept { return m_beam; }

    [[nodiscard]] double spectra_rate() const noexcept {
        return m_header.tsamp > 0.0 ? 1.0 / m_header.tsamp : 0.0;
    }

    [[nodiscard]] std::time_t utc_start() const noexcept {
        return astro::mjd_to_time(m_header.tstart);
    }

    [[nodiscard]] std::uint64_t tell() const noexcept { return m_sample; }

    void seek(std::uint64_t sample) {
        if (sample > m_header.nsamples) {
            throw ValidationError("psrio: seek is past the readable samples");
        }
        m_sample = sample;
    }

    void rewind() { seek(0); }

    void skip(std::int64_t delta) {
        m_sample =
            ::psrio::detail::apply_skip(m_sample, delta, m_header.nsamples);
    }

    void set_fswap(bool enable) noexcept { m_apply_fswap = enable; }
    [[nodiscard]] bool fswap_enabled() const noexcept { return m_apply_fswap; }

    /// Copy the next packed samples. The destination must hold the requested
    /// count. The return count is short at the end of the time axis.
    std::uint64_t read_block(std::uint64_t count, std::span<std::byte> dest);

    /// Read the next samples as floats. Integer datasets are converted to
    /// their numeric values. The destination length must equal the request.
    std::uint64_t read_samples(std::uint64_t count, std::span<float> dest);

    /// Read next samples as uint8_t values.
    std::uint64_t read_samples(std::uint64_t count,
                               std::span<std::uint8_t> dest);

    /// Read next samples as uint16_t values.
    std::uint64_t read_samples(std::uint64_t count,
                               std::span<std::uint16_t> dest);

    /// Read raw packed bytes directly into caller storage.
    std::uint64_t read_bytes(std::uint64_t nbytes, std::span<std::byte> dest);

    /// Convenience allocating read: copy next @p count packed time samples.
    [[nodiscard]] std::vector<std::byte> read_block(std::uint64_t count);

    /// Convenience allocating read: unpack next @p count time samples into a
    /// new vector.
    template <typename T = float>
    [[nodiscard]] std::vector<T> read_samples(std::uint64_t count);

    /// Convenience allocating read: copy next @p nbytes payload bytes.
    [[nodiscard]] std::vector<std::byte> read_bytes(std::uint64_t nbytes);

private:
    [[nodiscard]] std::uint64_t element_bytes() const noexcept {
        switch (m_sample_type) {
        case SampleType::kUInt8:
            return 1;
        case SampleType::kUInt16:
            return 2;
        case SampleType::kUInt32:
        case SampleType::kFloat32:
            return 4;
        case SampleType::kUInt1:
        case SampleType::kUInt2:
        case SampleType::kUInt4:
        case SampleType::kInt8:
            break;
        }
        return 0;
    }

    void load(const std::filesystem::path& path);
    void load_header(const std::vector<size_t>& dims,
                     const std::filesystem::path& path);
    void read_optional();
    void read_angle(const HighFive::Attribute& attr,
                    const char* name,
                    double& numeric,
                    std::string& text,
                    bool& has_numeric);

    std::unique_ptr<HighFive::File> m_file;
    std::unique_ptr<HighFive::DataSet> m_dataset;
    FilterbankHeader m_header{};
    SampleType m_sample_type{SampleType::kUInt8};
    std::uint64_t m_beam_index{0};
    std::uint64_t m_sample{0};
    int m_beam{0};
    bool m_apply_fswap{false};
};

inline FilterbankFile::FilterbankFile(const std::filesystem::path& path,
                                      std::uint64_t beam_index)
    : m_beam_index(beam_index) {
    HighFive::SilenceHDF5 silence;
    const auto native = path.string();
    try {
        m_file =
            std::make_unique<HighFive::File>(native, HighFive::File::ReadOnly);
    } catch (...) {
        throw IoError(
            std::format("psrio: failed to open HDF5 filterbank {}", native));
    }
    load(path);
}

inline void FilterbankFile::load(const std::filesystem::path& path) {
    detail::require_class(*m_file);
    m_header.version = detail::read_version(*m_file);
    if (!m_file->exist("data")) {
        throw FormatError(
            "psrio: HDF5 filterbank missing required dataset 'data'");
    }
    try {
        m_dataset =
            std::make_unique<HighFive::DataSet>(m_file->getDataSet("data"));
    } catch (...) {
        throw FormatError(
            "psrio: HDF5 filterbank missing required dataset 'data'");
    }
    detail::ensure_filters_available(*m_dataset);

    const auto dims = m_dataset->getDimensions();
    if (dims.size() != 3) {
        throw FormatError(std::format(
            "psrio: expected HDF5 data.ndim to be 3 (time, beam, frequency) "
            "but saw {}",
            dims.size()));
    }
    if (dims[0] == 0 || dims[1] == 0 || dims[2] == 0) {
        throw FormatError(
            "psrio: HDF5 data shape must be non-empty (time, beam, frequency)");
    }
    load_header(dims, path);
    if (m_beam_index >= dims[1]) {
        throw ValidationError(std::format(
            "psrio: HDF5 beam index {} is out of range for {} beam(s) in {}",
            m_beam_index, dims[1], path.string()));
    }
}

inline void FilterbankFile::load_header(const std::vector<size_t>& dims,
                                        const std::filesystem::path& /*path*/) {
    if (!detail::has_attr(*m_dataset, *m_file, "source_name")) {
        throw FormatError(
            "psrio: HDF5 filterbank missing required attribute 'source_name'");
    }
    m_header.source_name = detail::read_text(
        detail::get_attr(*m_dataset, *m_file, "source_name"), "source_name");

    const auto require_number = [&](const char* name) {
        if (!detail::has_attr(*m_dataset, *m_file, name)) {
            throw FormatError(std::format(
                "psrio: HDF5 filterbank missing required attribute '{}'",
                name));
        }
        return detail::read_number(detail::get_attr(*m_dataset, *m_file, name),
                                   name);
    };

    const auto require_count = [&](const char* name) {
        if (!detail::has_attr(*m_dataset, *m_file, name)) {
            throw FormatError(std::format(
                "psrio: HDF5 filterbank missing required attribute '{}'",
                name));
        }
        const int value =
            detail::read_int(detail::get_attr(*m_dataset, *m_file, name), name);
        if (value <= 0) {
            throw ValidationError(
                std::format("psrio: HDF5 {} must be positive", name));
        }
        return static_cast<std::uint64_t>(value);
    };

    m_header.tstart = require_number("tstart");
    m_header.tsamp  = require_number("tsamp");
    m_header.fch1   = require_number("fch1");
    m_header.foff   = require_number("foff");
    m_header.nchans = require_count("nchans");
    m_header.nbits  = [&] {
        if (!detail::has_attr(*m_dataset, *m_file, "nbits")) {
            throw FormatError(
                "psrio: HDF5 filterbank missing required attribute 'nbits'");
        }
        return detail::read_int(detail::get_attr(*m_dataset, *m_file, "nbits"),
                                "nbits");
    }();

    if (!std::isfinite(m_header.tstart)) {
        throw ValidationError("psrio: HDF5 tstart must be finite");
    }
    if (!(m_header.tsamp > 0.0) || !std::isfinite(m_header.tsamp)) {
        throw ValidationError("psrio: HDF5 tsamp must be positive");
    }
    if (!std::isfinite(m_header.fch1)) {
        throw ValidationError("psrio: HDF5 fch1 must be finite");
    }
    if (m_header.foff == 0.0 || !std::isfinite(m_header.foff)) {
        throw ValidationError(
            "psrio: HDF5 foff must be a non-zero finite channel width");
    }
    if (m_header.nchans != dims[2]) {
        throw FormatError(std::format(
            "psrio: HDF5 nchans ({}) does not match the frequency axis ({})",
            m_header.nchans, dims[2]));
    }

    const bool has_nifs   = detail::has_attr(*m_dataset, *m_file, "nifs");
    const bool has_nbeams = detail::has_attr(*m_dataset, *m_file, "nbeams");
    if (!has_nifs && !has_nbeams) {
        throw FormatError(
            "psrio: HDF5 filterbank missing required attribute 'nifs' or "
            "'nbeams'");
    }
    if (has_nifs) {
        m_header.nifs = require_count("nifs");
        if (m_header.nifs != dims[1]) {
            throw FormatError(std::format(
                "psrio: HDF5 nifs ({}) does not match the beam axis ({})",
                m_header.nifs, dims[1]));
        }
    }
    if (has_nbeams) {
        m_header.nbeams = require_count("nbeams");
        if (!has_nifs && m_header.nbeams != dims[1]) {
            throw FormatError(std::format(
                "psrio: HDF5 nbeams ({}) does not match the beam axis ({})",
                m_header.nbeams, dims[1]));
        }
    }
    if (!has_nifs) {
        m_header.nifs = dims[1];
    }
    if (!has_nbeams) {
        m_header.nbeams = m_header.nifs;
    }

    if (detail::has_attr(*m_dataset, *m_file, "nsamples")) {
        const int declared = detail::read_int(
            detail::get_attr(*m_dataset, *m_file, "nsamples"), "nsamples");
        if (declared < 0 || static_cast<std::uint64_t>(declared) != dims[0]) {
            throw FormatError(std::format(
                "psrio: HDF5 nsamples ({}) does not match the time axis ({})",
                declared, dims[0]));
        }
    }
    m_header.nsamples = dims[0];

    m_sample_type = detail::sample_type_from_dataset(*m_dataset);
    if (!detail::nbits_matches(m_sample_type, m_header.nbits)) {
        throw FormatError(std::format(
            "psrio: HDF5 nbits ({}) does not match the data dataset type",
            m_header.nbits));
    }
    if (m_sample_type == SampleType::kUInt1 ||
        m_sample_type == SampleType::kUInt2 ||
        m_sample_type == SampleType::kUInt4) {
        throw FormatError(
            "psrio: unsupported HDF5 sample type: packed bit depths are not "
            "written by blimpy");
    }

    read_optional();

    if (dims[1] == 1 && m_header.ibeam > 0) {
        m_beam = m_header.ibeam;
    } else if (m_beam_index >=
               static_cast<std::uint64_t>(std::numeric_limits<int>::max())) {
        throw ValidationError("psrio: HDF5 beam index does not fit in int");
    } else {
        m_beam = static_cast<int>(m_beam_index + 1U);
    }
}

inline void FilterbankFile::read_optional() {
    if (detail::has_attr(*m_dataset, *m_file, "telescope_id")) {
        m_header.telescope_id = detail::read_int(
            detail::get_attr(*m_dataset, *m_file, "telescope_id"),
            "telescope_id");
        m_header.has_telescope_id = true;
    }
    if (detail::has_attr(*m_dataset, *m_file, "machine_id")) {
        m_header.machine_id = detail::read_int(
            detail::get_attr(*m_dataset, *m_file, "machine_id"), "machine_id");
    }
    if (detail::has_attr(*m_dataset, *m_file, "data_type")) {
        m_header.data_type = detail::read_int(
            detail::get_attr(*m_dataset, *m_file, "data_type"), "data_type");
    }
    if (detail::has_attr(*m_dataset, *m_file, "ibeam")) {
        m_header.ibeam = detail::read_int(
            detail::get_attr(*m_dataset, *m_file, "ibeam"), "ibeam");
    }
    if (detail::has_attr(*m_dataset, *m_file, "barycentric")) {
        m_header.barycentric = detail::read_int(
            detail::get_attr(*m_dataset, *m_file, "barycentric"),
            "barycentric");
    }
    if (detail::has_attr(*m_dataset, *m_file, "pulsarcentric")) {
        m_header.pulsarcentric = detail::read_int(
            detail::get_attr(*m_dataset, *m_file, "pulsarcentric"),
            "pulsarcentric");
    }
    if (detail::has_attr(*m_dataset, *m_file, "az_start")) {
        m_header.az_start = detail::read_number(
            detail::get_attr(*m_dataset, *m_file, "az_start"), "az_start");
    }
    if (detail::has_attr(*m_dataset, *m_file, "za_start")) {
        m_header.za_start = detail::read_number(
            detail::get_attr(*m_dataset, *m_file, "za_start"), "za_start");
    }
    if (detail::has_attr(*m_dataset, *m_file, "src_raj")) {
        read_angle(detail::get_attr(*m_dataset, *m_file, "src_raj"), "src_raj",
                   m_header.src_raj, m_header.src_raj_text,
                   m_header.has_src_raj);
    }
    if (detail::has_attr(*m_dataset, *m_file, "src_dej")) {
        read_angle(detail::get_attr(*m_dataset, *m_file, "src_dej"), "src_dej",
                   m_header.src_dej, m_header.src_dej_text,
                   m_header.has_src_dej);
    }
}

inline void FilterbankFile::read_angle(const HighFive::Attribute& attr,
                                       const char* name,
                                       double& numeric,
                                       std::string& text,
                                       bool& has_numeric) {
    const auto cls = attr.getDataType().getClass();
    if (cls == HighFive::DataTypeClass::String) {
        text        = detail::read_text(attr, name);
        has_numeric = false;
        return;
    }
    numeric     = detail::read_number(attr, name);
    has_numeric = true;
}

inline std::uint64_t FilterbankFile::read_block(std::uint64_t count,
                                                std::span<std::byte> dest) {
    const auto stride = bytes_per_sample();
    if (count > 0U &&
        stride > (std::numeric_limits<std::uint64_t>::max() / count)) {
        throw ValidationError("psrio: requested sample block is too large");
    }
    const auto bytes_needed = count * stride;
    if (dest.size() < bytes_needed) {
        throw ValidationError(std::format(
            "psrio: destination has {} bytes but the block needs {}",
            dest.size(), bytes_needed));
    }
    if (m_sample > m_header.nsamples) {
        throw ValidationError(
            "psrio: reader cursor is past the readable samples");
    }
    const auto to_read = std::min(count, m_header.nsamples - m_sample);
    if (to_read == 0U) {
        return 0;
    }
    const auto element = element_bytes();
    if (element > 1U &&
        reinterpret_cast<std::uintptr_t>(dest.data()) % element != 0U) {
        throw ValidationError(
            "psrio: HDF5 destination is not aligned for the sample type");
    }
    try {
        HighFive::SilenceHDF5 silence;
        auto selection = m_dataset->select(
            HighFive::HyperSlab(HighFive::RegularHyperSlab(
                {static_cast<size_t>(m_sample),
                 static_cast<size_t>(m_beam_index), 0},
                {static_cast<size_t>(to_read), 1,
                 static_cast<size_t>(m_header.nchans)})),
            HighFive::DataSpace(
                {static_cast<size_t>(to_read * m_header.nchans)}));

        switch (m_sample_type) {
        case SampleType::kUInt8:
            selection.read_raw(reinterpret_cast<std::uint8_t*>(dest.data()));
            break;
        case SampleType::kUInt16:
            selection.read_raw(reinterpret_cast<std::uint16_t*>(dest.data()));
            break;
        case SampleType::kUInt32:
            selection.read_raw(reinterpret_cast<std::uint32_t*>(dest.data()));
            break;
        case SampleType::kFloat32:
            selection.read_raw(reinterpret_cast<float*>(dest.data()));
            break;
        default:
            throw ValidationError(
                "psrio: unsupported sample type for HDF5 read_block");
        }
    } catch (const ValidationError&) {
        throw;
    } catch (const std::exception& e) {
        throw IoError(
            std::format("psrio: HDF5 filterbank read failed: {}", e.what()));
    }

    if (m_apply_fswap && to_read > 0U) {
        reverse_channels(dest.first(static_cast<std::size_t>(to_read * stride)),
                         to_read, m_header.nchans, m_header.nbits);
    }
    m_sample += to_read;
    return to_read;
}

inline std::uint64_t FilterbankFile::read_samples(std::uint64_t count,
                                                  std::span<float> dest) {
    if (count > 0U &&
        m_header.nchans > (std::numeric_limits<std::uint64_t>::max() / count)) {
        throw ValidationError("psrio: requested sample block is too large");
    }
    const auto expected = count * m_header.nchans;
    if (expected != dest.size()) {
        throw ValidationError(std::format(
            "psrio: destination has {} values but the request needs {}",
            dest.size(), expected));
    }
    if (m_sample > m_header.nsamples) {
        throw ValidationError(
            "psrio: reader cursor is past the readable samples");
    }
    const auto to_read = std::min(count, m_header.nsamples - m_sample);
    if (to_read == 0U) {
        return 0;
    }
    try {
        HighFive::SilenceHDF5 silence;
        auto selection = m_dataset->select(
            HighFive::HyperSlab(HighFive::RegularHyperSlab(
                {static_cast<size_t>(m_sample),
                 static_cast<size_t>(m_beam_index), 0},
                {static_cast<size_t>(to_read), 1,
                 static_cast<size_t>(m_header.nchans)})),
            HighFive::DataSpace(
                {static_cast<size_t>(to_read * m_header.nchans)}));
        selection.read_raw(dest.data());
    } catch (const std::exception& e) {
        throw IoError(
            std::format("psrio: HDF5 filterbank read failed: {}", e.what()));
    }
    m_sample += to_read;
    return to_read;
}

inline std::uint64_t
FilterbankFile::read_samples(std::uint64_t count,
                             std::span<std::uint8_t> dest) {
    if (m_sample_type != SampleType::kUInt8) {
        throw ValidationError(
            "psrio: uint8_t output requires uint8 sample type");
    }
    if (count > 0U &&
        m_header.nchans > (std::numeric_limits<std::uint64_t>::max() / count)) {
        throw ValidationError("psrio: requested sample block is too large");
    }
    const auto expected = count * m_header.nchans;
    if (expected != dest.size()) {
        throw ValidationError(std::format(
            "psrio: destination has {} values but the request needs {}",
            dest.size(), expected));
    }
    if (m_sample > m_header.nsamples) {
        throw ValidationError(
            "psrio: reader cursor is past the readable samples");
    }
    const auto to_read = std::min(count, m_header.nsamples - m_sample);
    if (to_read == 0U) {
        return 0;
    }
    try {
        HighFive::SilenceHDF5 silence;
        auto selection = m_dataset->select(
            HighFive::HyperSlab(HighFive::RegularHyperSlab(
                {static_cast<size_t>(m_sample),
                 static_cast<size_t>(m_beam_index), 0},
                {static_cast<size_t>(to_read), 1,
                 static_cast<size_t>(m_header.nchans)})),
            HighFive::DataSpace(
                {static_cast<size_t>(to_read * m_header.nchans)}));
        selection.read_raw(dest.data());
    } catch (const std::exception& e) {
        throw IoError(
            std::format("psrio: HDF5 filterbank read failed: {}", e.what()));
    }
    m_sample += to_read;
    return to_read;
}

inline std::uint64_t
FilterbankFile::read_samples(std::uint64_t count,
                             std::span<std::uint16_t> dest) {
    if (m_sample_type != SampleType::kUInt16) {
        throw ValidationError(
            "psrio: uint16_t output requires uint16 sample type");
    }
    if (count > 0U &&
        m_header.nchans > (std::numeric_limits<std::uint64_t>::max() / count)) {
        throw ValidationError("psrio: requested sample block is too large");
    }
    const auto expected = count * m_header.nchans;
    if (expected != dest.size()) {
        throw ValidationError(std::format(
            "psrio: destination has {} values but the request needs {}",
            dest.size(), expected));
    }
    if (m_sample > m_header.nsamples) {
        throw ValidationError(
            "psrio: reader cursor is past the readable samples");
    }
    const auto to_read = std::min(count, m_header.nsamples - m_sample);
    if (to_read == 0U) {
        return 0;
    }
    try {
        HighFive::SilenceHDF5 silence;
        auto selection = m_dataset->select(
            HighFive::HyperSlab(HighFive::RegularHyperSlab(
                {static_cast<size_t>(m_sample),
                 static_cast<size_t>(m_beam_index), 0},
                {static_cast<size_t>(to_read), 1,
                 static_cast<size_t>(m_header.nchans)})),
            HighFive::DataSpace(
                {static_cast<size_t>(to_read * m_header.nchans)}));
        selection.read_raw(dest.data());
    } catch (const std::exception& e) {
        throw IoError(
            std::format("psrio: HDF5 filterbank read failed: {}", e.what()));
    }
    m_sample += to_read;
    return to_read;
}

inline std::uint64_t FilterbankFile::read_bytes(std::uint64_t nbytes,
                                                std::span<std::byte> dest) {
    const auto stride = bytes_per_sample();
    if (stride == 0U || nbytes % stride != 0U) {
        throw ValidationError(
            "psrio: byte request must be a multiple of the sample stride");
    }
    const auto samples = nbytes / stride;
    return read_block(samples, dest);
}

inline std::vector<std::byte> FilterbankFile::read_block(std::uint64_t count) {
    const auto stride = bytes_per_sample();
    std::vector<std::byte> out(static_cast<std::size_t>(count * stride));
    const auto actual = read_block(count, std::span<std::byte>(out));
    out.resize(static_cast<std::size_t>(actual * stride));
    return out;
}

template <typename T>
inline std::vector<T> FilterbankFile::read_samples(std::uint64_t count) {
    const auto values_needed = count * m_header.nchans;
    std::vector<T> out(static_cast<std::size_t>(values_needed));
    const auto actual = read_samples(count, std::span<T>(out));
    out.resize(static_cast<std::size_t>(actual * m_header.nchans));
    return out;
}

inline std::vector<std::byte> FilterbankFile::read_bytes(std::uint64_t nbytes) {
    std::vector<std::byte> out(static_cast<std::size_t>(nbytes));
    const auto actual = read_bytes(nbytes, std::span<std::byte>(out));
    out.resize(static_cast<std::size_t>(actual));
    return out;
}

} // namespace psrio::formats::fbh5

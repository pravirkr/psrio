#pragma once

#include "psrio/astro.hpp"
#include "psrio/common/types.hpp"
#include "psrio/detail/endian.hpp"
#include "psrio/detail/exceptions.hpp"
#include "psrio/detail/parallel.hpp"
#include "psrio/detail/readahead.hpp"
#include "psrio/detail/rfile.hpp"
#include "psrio/detail/skip.hpp"
#include "psrio/formats/psrfits/table.hpp"
#include "psrio/formats/sigproc/header.hpp"
#include "psrio/packed.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <ctime>
#include <exception>
#include <filesystem>
#include <fmt/format.h>
#include <limits>
#include <memory>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace psrio::formats::psrfits::detail {

inline constexpr std::array<float, 4> kTwoBit{-2.0F, -0.5F, 0.5F, 2.0F};
inline constexpr double kFreqAtolMhz = 1.0e-4;

struct DataLayout {
    std::array<int, 3> shape{};
    int time_axis{0};
    int pol_axis{1};
    int chan_axis{2};
    int elem_size{1};
    int group{1};
    int stride{0};
    bool fast{false};
};

struct PolPlan {
    bool sum{false};
};

[[nodiscard]] inline PolPlan plan_polarisation(std::string_view type,
                                               int npol) {
    PolPlan plan;
    const auto fail = [&](std::string_view message) {
        throw FormatError(std::string{message});
    };
    int allowed_a = 0;
    int allowed_b = 0;
    bool sum      = false;
    if (type == "AABB" || type == "XXYY" || type == "LLRR") {
        allowed_a = 2;
        sum       = true;
    } else if (type == "AABBCRCI" || type == "XXYYCRCI" || type == "LLRRCRCI") {
        allowed_a = 4;
        sum       = true;
    } else if (type == "INTEN" || type == "AA+BB") {
        allowed_a = 1;
    } else if (type == "STOKE" || type == "STOKES" || type == "IQUV") {
        allowed_a = 1;
        allowed_b = 4;
    } else {
        fail(
            fmt::format("psrio: unsupported PSRFITS POL_TYPE '{}' with NPOL {}",
                        type, npol));
    }
    if (npol != allowed_a && npol != allowed_b) {
        fail(fmt::format(
            "psrio: PSRFITS POL_TYPE '{}' is inconsistent with NPOL {}", type,
            npol));
    }
    plan.sum = sum;
    return plan;
}

[[nodiscard]] inline DataLayout make_layout(const Cards& header,
                                            const Column& data,
                                            int field,
                                            int nchan,
                                            int npol,
                                            int nbits,
                                            int nsblk) {
    if (nbits != 1 && nbits != 2 && nbits != 4 && nbits != 8 && nbits != 16) {
        throw FormatError(fmt::format(
            "psrio: PSRFITS NBITS {} is not 1, 2, 4, 8, or 16", nbits));
    }
    const int samples_per_byte = nbits < 8 ? 8 / nbits : 1;
    if (nbits < 8 && (nchan * nbits) % 8 != 0) {
        throw FormatError(
            fmt::format("psrio: NCHAN {} is not a whole number of {}-bit bytes",
                        nchan, nbits));
    }
    if (nbits < 8 && nsblk % samples_per_byte != 0) {
        throw FormatError(
            fmt::format("psrio: NSBLK {} is not a whole number of {}-bit bytes",
                        nsblk, nbits));
    }
    const int time_len = nbits < 8 ? nsblk / samples_per_byte : nsblk;
    std::vector<int> tdim;
    const auto tdim_key = fmt::format("TDIM{}", field);
    if (header.contains(tdim_key)) {
        tdim = parse_tdim(header.value(tdim_key));
    } else {
        tdim = {nchan, npol, time_len};
    }
    if (tdim.size() != 3U) {
        throw FormatError(fmt::format(
            "psrio: DATA TDIM has {} axes, expected 3", tdim.size()));
    }
    std::uint64_t product = 1;
    for (const int axis : tdim) {
        product =
            checked_mul(product, static_cast<std::uint64_t>(axis), "TDIM");
    }
    if (!data.heap && product != data.repeat) {
        throw FormatError(fmt::format(
            "psrio: DATA TDIM product {} does not match TFORM repeat {}",
            product, data.repeat));
    }
    DataLayout layout;
    layout.shape     = {tdim[2], tdim[1], tdim[0]};
    layout.elem_size = static_cast<int>(data.element_bytes);
    if ((nbits == 16 && layout.elem_size != 2) ||
        (nbits != 16 && layout.elem_size != 1)) {
        throw FormatError(
            fmt::format("psrio: {}-bit DATA element size is {} byte(s)", nbits,
                        layout.elem_size));
    }
    const std::array<int, 3> wanted{time_len, npol, nchan};
    constexpr auto kPerms = std::to_array<std::array<int, 3>>({
        std::array<int, 3>{0, 1, 2},
        std::array<int, 3>{0, 2, 1},
        std::array<int, 3>{1, 0, 2},
        std::array<int, 3>{1, 2, 0},
        std::array<int, 3>{2, 0, 1},
        std::array<int, 3>{2, 1, 0},
    });
    int chosen            = -1;
    if (layout.shape == wanted) {
        chosen = 0;
    } else {
        for (std::size_t index = 0; index < kPerms.size(); ++index) {
            const auto& perm = kPerms[index];
            if (layout.shape[static_cast<std::size_t>(perm[0])] == time_len &&
                layout.shape[static_cast<std::size_t>(perm[1])] == npol &&
                layout.shape[static_cast<std::size_t>(perm[2])] == nchan) {
                if (chosen >= 0) {
                    throw FormatError(fmt::format(
                        "psrio: DATA shape ({}, {}, {}) is ambiguous for "
                        "time {}, pol {}, chan {}",
                        layout.shape[0], layout.shape[1], layout.shape[2],
                        time_len, npol, nchan));
                }
                chosen = static_cast<int>(index);
            }
        }
    }
    if (chosen < 0) {
        throw FormatError(fmt::format(
            "psrio: DATA shape ({}, {}, {}) does not match time {}, pol {}, "
            "chan {}",
            layout.shape[0], layout.shape[1], layout.shape[2], time_len, npol,
            nchan));
    }
    const auto& perm = kPerms[static_cast<std::size_t>(chosen)];
    layout.time_axis = perm[0];
    layout.pol_axis  = perm[1];
    layout.chan_axis = perm[2];
    const int chan_len =
        layout.shape[static_cast<std::size_t>(layout.chan_axis)];
    layout.stride = nbits < 8 ? (nchan * nbits) / 8 : nchan * layout.elem_size;
    const int byte_len = chan_len * layout.elem_size;
    if (layout.stride <= 0 || byte_len % layout.stride != 0) {
        throw FormatError(
            "psrio: DATA channel axis does not pack into samples");
    }
    layout.group = byte_len / layout.stride;
    if (layout.group <= 0 ||
        checked_mul(static_cast<std::uint64_t>(time_len),
                    static_cast<std::uint64_t>(layout.group),
                    "samples") != static_cast<std::uint64_t>(nsblk)) {
        throw FormatError(
            fmt::format("psrio: DATA packing does not cover NSBLK {}", nsblk));
    }
    layout.fast = nbits >= 8 && layout.group == 1 && layout.time_axis == 0 &&
                  layout.pol_axis == 1 && layout.chan_axis == 2 &&
                  layout.shape[0] == time_len && layout.shape[1] == npol &&
                  layout.shape[2] == nchan;
    return layout;
}

[[nodiscard]] inline int narrow_positive(std::int64_t value,
                                         std::string_view name) {
    if (value <= 0 ||
        value > static_cast<std::int64_t>(std::numeric_limits<int>::max())) {
        throw FormatError(
            fmt::format("psrio: PSRFITS {} {} is not usable", name, value));
    }
    return static_cast<int>(value);
}

/// Optional numeric card that ignores PSRFITS template unset placeholders
/// (empty string or "*").
[[nodiscard]] inline std::optional<double>
find_optional_number(const Cards& cards, std::string_view key) {
    const auto text = cards.find(key);
    if (!text.has_value()) {
        return std::nullopt;
    }
    const auto cleaned = trim(*text);
    if (cleaned.empty() || cleaned == "*") {
        return std::nullopt;
    }
    return cards.find_number(key);
}

/// Optional integer card that ignores PSRFITS template unset placeholders
/// (empty string or "*").
[[nodiscard]] inline std::optional<std::int64_t>
find_optional_integer(const Cards& cards, std::string_view key) {
    const auto text = cards.find(key);
    if (!text.has_value()) {
        return std::nullopt;
    }
    const auto cleaned = trim(*text);
    if (cleaned.empty() || cleaned == "*") {
        return std::nullopt;
    }
    return cards.find_integer(key);
}

} // namespace psrio::formats::psrfits::detail

namespace psrio::formats::psrfits {

/// How SearchFile turns stored DATA into float samples.
struct DecodeOptions {
    /// Apply `DAT_SCL`, `DAT_OFFS`, and `ZERO_OFF` (dspsr `--scloffs`).
    /// Without them samples are the stored digits (2-bit data keeps the dspsr
    /// levels).
    bool scloffs{true};
    /// Multiply by `DAT_WTS`.
    bool weights{true};
    /// Stored polarisation product to read (0 .. NPOL-1, for example 0 for
    /// AA and 1 for BB), or -1 for Stokes I (the POL_TYPE rule below).
    int product{-1};
};

/**
 * @brief Search-mode PSRFITS file, streamed as calibrated Stokes I.
 *
 * The file is read with `pread`, not mapped: headers and the per-row
 * frequency, scale, offset, and weight cells are read as needed, and DATA is
 * staged a few MiB of time samples at a time. A reader holds one open file
 * descriptor and buffers of about kStageBytes.
 *
 * One uncompressed SUBINT binary table. Samples are float32 in `DAT_FREQ`
 * order, one IF. Coherence and two-product feeds (`AABBCRCI`, `AABB`, and the
 * XX/LL equivalents) sum the first two products. Intensity and Stokes labels
 * copy product 0. `DAT_SCL`, `DAT_OFFS`, `ZERO_OFF`, and `DAT_WTS` follow the
 * dspsr rules used by sigpyproc. Fold mode and tile compression are rejected.
 */
class SearchFile {
public:
    /// Bytes of DATA staged per read from a time-major SUBINT row.
    static constexpr std::uint64_t kStageBytes = std::uint64_t{16} << 20U;

    /// Blocks of ReadAhead::kBlockBytes read ahead of the caller.
    static constexpr std::size_t kReadAheadBlocks = 4;

    /// @throws IoError if @p path cannot be opened.
    /// @throws FormatError if @p path is not an uncompressed search-mode file.
    explicit SearchFile(const std::filesystem::path& path)
        : SearchFile(path, DecodeOptions{}) {}

    /// @throws ValidationError if @p options.product is not -1 or a stored
    ///         product.
    SearchFile(const std::filesystem::path& path, DecodeOptions options);

    SearchFile(const SearchFile&)                = delete;
    SearchFile& operator=(const SearchFile&)     = delete;
    SearchFile(SearchFile&&) noexcept            = default;
    SearchFile& operator=(SearchFile&&) noexcept = default;
    ~SearchFile()                                = default;

    /// SIGPROC header for the calibrated float stream (nbits 32, one IF).
    [[nodiscard]] sigproc::FilterbankHeader sigproc_header() const;

    /// Bits stored in the SUBINT DATA column, before calibration.
    [[nodiscard]] int input_nbits() const noexcept { return m_input_nbits; }

    /// Polarisation products stored in DATA, before the Stokes I sum.
    [[nodiscard]] std::uint64_t input_npol() const noexcept { return m_npol; }

    /// Samples in one SUBINT row (`NSBLK`).
    [[nodiscard]] std::uint64_t subint_samples() const noexcept {
        return m_nsblk;
    }

    /// True when samples are the sum of the first two stored products.
    [[nodiscard]] bool sums_polarisations() const noexcept { return m_sum; }

    [[nodiscard]] const DecodeOptions& decode_options() const noexcept {
        return m_options;
    }

    /// SUBINT rows in this file.
    [[nodiscard]] std::uint64_t nrows() const noexcept { return m_table.nrows; }

    /// `NSUBOFFS` from the sub-integration header, or the primary header.
    [[nodiscard]] std::optional<std::int64_t> nsuboffs() const {
        if (const auto value =
                detail::find_optional_integer(m_table.subint, "NSUBOFFS")) {
            return value;
        }
        return detail::find_optional_integer(m_table.primary, "NSUBOFFS");
    }

    /// Exact split start time (`STT_IMJD`, `STT_SMJD`, `STT_OFFS`).
    [[nodiscard]] std::int64_t stt_imjd() const noexcept { return m_stt_imjd; }
    [[nodiscard]] double stt_smjd() const noexcept { return m_stt_smjd; }
    [[nodiscard]] double stt_offs() const noexcept { return m_stt_offs; }

    /// Path and size in bytes of this file on disk.
    [[nodiscard]] const std::filesystem::path& path() const noexcept {
        return m_file.path();
    }
    [[nodiscard]] std::uint64_t file_size() const noexcept {
        return m_file.size();
    }
    /// True while the underlying OS file descriptor is open.
    [[nodiscard]] bool is_open() const noexcept { return m_file.is_open(); }

    [[nodiscard]] std::string_view pol_type() const noexcept {
        return m_pol_type;
    }

    /// Primary header `FRONTEND` (receiver), or empty.
    [[nodiscard]] std::string_view frontend() const noexcept {
        return m_frontend;
    }
    /// Primary header `BACKEND`, or empty.
    [[nodiscard]] std::string_view backend() const noexcept {
        return m_backend;
    }
    /// Primary header `OBS_MODE` (always SEARCH for a readable file).
    [[nodiscard]] std::string_view obs_mode() const noexcept {
        return m_obs_mode;
    }

    /// `CHAN_DM`: the DM the data were coherently dedispersed to, when the
    /// primary header records one.
    [[nodiscard]] std::optional<double> refdm() const noexcept {
        return m_has_dm ? std::optional<double>(m_dm) : std::nullopt;
    }

    [[nodiscard]] std::uint64_t nchans() const noexcept { return m_nchan; }
    [[nodiscard]] static std::uint64_t nifs() noexcept { return 1; }
    [[nodiscard]] static int nbits() noexcept { return 32; }
    [[nodiscard]] static SampleType sample_type() noexcept {
        return SampleType::kFloat32;
    }
    [[nodiscard]] std::uint64_t bytes_per_sample() const noexcept {
        return m_bytes_per_sample;
    }
    [[nodiscard]] std::uint64_t nsamples() const noexcept { return m_nsamples; }
    [[nodiscard]] static bool has_nsamples() noexcept { return true; }
    [[nodiscard]] double tsamp() const noexcept { return m_tsamp; }
    [[nodiscard]] double tstart() const noexcept { return m_tstart; }
    [[nodiscard]] std::string_view source_name() const noexcept {
        return m_source;
    }
    [[nodiscard]] std::string_view telescope() const noexcept {
        return m_telescope;
    }
    [[nodiscard]] double raj() const noexcept { return m_raj; }
    [[nodiscard]] double dej() const noexcept { return m_dej; }
    [[nodiscard]] double fch1() const noexcept { return m_fch1; }
    [[nodiscard]] double foff() const noexcept { return m_foff; }
    [[nodiscard]] double bandwidth() const noexcept { return m_bandwidth; }
    [[nodiscard]] double center_frequency() const noexcept { return m_centre; }
    [[nodiscard]] int beam() const noexcept { return m_beam; }
    [[nodiscard]] double spectra_rate() const noexcept {
        return m_tsamp > 0.0 ? 1.0 / m_tsamp : 0.0;
    }
    [[nodiscard]] std::time_t utc_start() const noexcept {
        return astro::mjd_to_time(m_tstart);
    }
    [[nodiscard]] std::uint64_t tell() const noexcept { return m_sample; }

    void set_fswap(bool enable) noexcept { m_fswap = enable; }
    [[nodiscard]] bool fswap_enabled() const noexcept { return m_fswap; }

    /// Drop cached file pages once their samples are read (Linux). Off by
    /// default.
    void set_drop_cache(bool enable) noexcept { m_file.set_drop_cache(enable); }

    /// Once reads run forward, stream the rest of the DATA on a background
    /// thread, read_ahead_blocks() blocks ahead of the caller. On by default.
    void set_prefetch(bool enable);

    /// Number of ReadAhead::kBlockBytes blocks buffered ahead of the caller.
    void set_read_ahead_blocks(std::size_t blocks) noexcept {
        m_read_ahead_blocks = std::max<std::size_t>(blocks, 1U);
    }

    /// Start streaming DATA from the beginning of this file on a background
    /// thread before the first read_samples() call.
    void prefetch_start();

    /// Free the staging buffers, for example once a stream has moved past
    /// this file. A later read allocates them again.
    void release_buffers();

    /// Free staging and calibration buffers and close the OS file descriptor.
    /// Metadata stays available, and a later read reopens the file on demand.
    void close_handle();

    void seek(std::uint64_t sample);
    void rewind() { seek(0); }
    void skip(std::int64_t delta);

    std::uint64_t read_block(std::uint64_t count, std::span<std::byte> dest);
    std::uint64_t read_samples(std::uint64_t count, std::span<float> dest);
    static std::uint64_t read_samples(std::uint64_t count,
                                      std::span<std::uint8_t> dest);
    static std::uint64_t read_samples(std::uint64_t count,
                                      std::span<std::uint16_t> dest);
    std::uint64_t read_bytes(std::uint64_t nbytes, std::span<std::byte> dest);

    template <typename T>
        requires std::same_as<T, float> || std::same_as<T, std::uint8_t> ||
                 std::same_as<T, std::uint16_t>
    std::uint64_t read(std::uint64_t count, std::span<T> dest) {
        return read_samples(count, dest);
    }

    [[nodiscard]] std::vector<std::byte> read_block(std::uint64_t count);
    [[nodiscard]] std::vector<std::byte> read_bytes(std::uint64_t nbytes);

    template <typename T = float>
        requires std::same_as<T, float> || std::same_as<T, std::uint8_t> ||
                 std::same_as<T, std::uint16_t>
    [[nodiscard]] std::vector<T> read_samples(std::uint64_t count) {
        const auto per_sample = m_nchan;
        std::vector<T> buffer(static_cast<std::size_t>(count * per_sample));
        const auto actual = read_samples(count, std::span<T>(buffer));
        buffer.resize(static_cast<std::size_t>(actual * per_sample));
        return buffer;
    }

private:
    void load(const std::filesystem::path& path);
    void load_geometry();
    void load_frequencies();
    void load_timing();
    void load_calibration(std::uint64_t row);

    /// DATA bytes holding samples [first, last) of one row. @p cell points at
    /// sample @p first for a time-major layout, or at the start of the row.
    struct Staged {
        const std::byte* cell{nullptr};
        std::uint64_t first{0};
        std::uint64_t last{0};
    };

    /// Make sample @p isamp of row @p isub, and those after it, available.
    Staged stage(std::uint64_t isub, std::uint64_t isamp);

    /// DATA range of row @p isub, checked against the expected cell size.
    [[nodiscard]] detail::CellRange data_range(std::uint64_t isub);

    /// Stream the file from @p offset to the end of the last DATA cell on a
    /// background thread.
    void start_read_ahead(std::uint64_t offset);

    /// Decode @p count samples from @p cell, starting at sample @p isamp
    /// relative to the cell, into @p dest.
    void decode_run(const std::byte* cell,
                    std::uint64_t isamp,
                    std::uint64_t count,
                    std::span<float> dest) const;

    [[nodiscard]] detail::Source source() const noexcept {
        return detail::Source(m_file);
    }
    void check_sample_request(std::uint64_t count,
                              std::size_t dest_values) const;

    template <bool Signed, bool Sum, bool Calibrated>
    void decode_u8(const std::byte* cell,
                   std::uint64_t isamp,
                   std::span<float> dest) const;

    template <bool Sum>
    void decode_i16(const std::byte* cell,
                    std::uint64_t isamp,
                    std::span<float> dest) const;

    template <bool Sum>
    void decode_packed(const std::byte* cell,
                       std::uint64_t isamp,
                       std::span<float> dest) const;

    template <bool Sum>
    void decode_indexed(const std::byte* cell,
                        std::uint64_t isamp,
                        std::span<float> dest) const;

    [[nodiscard]] std::uint8_t packed_byte(const std::byte* cell,
                                           std::uint64_t isamp,
                                           int pol,
                                           int byte_index) const;

    [[nodiscard]] std::size_t
    element_index(std::uint64_t isamp, int pol, int chan) const;

    [[nodiscard]] double scalar_at(std::uint64_t row,
                                   const detail::Column& column);

    /// Bytes of a fixed-width metadata cell. The metadata columns of a row
    /// are read in one request and kept until another row is asked for.
    [[nodiscard]] std::span<const std::byte>
    meta_cell(std::uint64_t row, const detail::Column& column);

    ::psrio::detail::ReadFile m_file;
    detail::ParsedTable m_table{};
    detail::Column m_data{};
    detail::Column m_scales_col{};
    detail::Column m_offsets_col{};
    detail::Column m_weights_col{};
    detail::DataLayout m_layout{};
    std::uint64_t m_nchan{0};
    std::uint64_t m_npol{0};
    std::uint64_t m_nsblk{0};
    std::uint64_t m_data_bytes{0};
    std::uint64_t m_nsamples{0};
    std::uint64_t m_sample{0};
    std::uint64_t m_bytes_per_sample{0};
    std::uint64_t m_cal_row{std::numeric_limits<std::uint64_t>::max()};
    int m_input_nbits{0};
    int m_beam{0};
    bool m_sum{false};
    bool m_signint{false};
    bool m_fswap{false};
    bool m_uniform{false};
    bool m_has_raj{false};
    bool m_has_dej{false};
    bool m_has_dm{false};
    float m_zero{0.0F};
    /// First stored product read by the decoders.
    std::uint64_t m_pol0{0};
    DecodeOptions m_options{};
    double m_tsamp{0.0};
    double m_tstart{0.0};
    std::int64_t m_stt_imjd{0};
    double m_stt_smjd{0.0};
    double m_stt_offs{0.0};
    double m_fch1{0.0};
    double m_foff{0.0};
    double m_bandwidth{0.0};
    double m_centre{0.0};
    double m_raj{0.0};
    double m_dej{0.0};
    double m_dm{0.0};
    std::optional<double> m_az_start;
    std::optional<double> m_za_start;
    std::string m_source;
    std::string m_telescope;
    std::string m_pol_type;
    std::string m_frontend;
    std::string m_backend;
    std::string m_obs_mode;
    std::vector<double> m_freqs;
    std::vector<double> m_scale_raw;
    std::vector<double> m_offset_raw;
    std::vector<double> m_weight_raw;
    std::vector<float> m_scales;
    std::vector<float> m_offsets;
    std::vector<float> m_weights;
    std::vector<float> m_floats;
    std::vector<std::byte> m_cell;
    std::vector<std::byte> m_meta;
    std::uint64_t m_meta_row{std::numeric_limits<std::uint64_t>::max()};
    std::uint64_t m_meta_first{0};
    std::uint64_t m_meta_bytes{0};
    std::vector<double> m_row_tsub;
    std::vector<double> m_row_offs;
    std::vector<std::byte> m_stage;
    detail::CellRange m_data_range{};
    std::uint64_t m_data_row{std::numeric_limits<std::uint64_t>::max()};
    std::uint64_t m_stage_row{std::numeric_limits<std::uint64_t>::max()};
    std::uint64_t m_stage_first{0};
    std::uint64_t m_stage_last{0};
    detail::CellRange m_stage_range{};

    /// DATA streamed ahead on another thread while the caller decodes.
    /// Destroyed before m_file.
    std::unique_ptr<::psrio::detail::ReadAhead> m_ahead;
    bool m_prefetch{true};
    std::size_t m_read_ahead_blocks{kReadAheadBlocks};
};

inline SearchFile::SearchFile(const std::filesystem::path& path,
                              DecodeOptions options)
    : m_options(options) {
    load(path);
}

inline void SearchFile::load(const std::filesystem::path& path) {
    try {
        // Opening reads a little metadata from every row, so no sequential
        // read-ahead; DATA is read in large requests regardless.
        m_file =
            ::psrio::detail::ReadFile(path, ::psrio::detail::Access::kNormal);
    } catch (const std::exception& ex) {
        throw IoError(ex.what());
    }
    m_table             = detail::parse_search(source(), path.string());
    const auto& primary = m_table.primary;
    if (!detail::equals_ignore_case(primary.value("FITSTYPE"), "PSRFITS")) {
        throw FormatError(
            fmt::format("{}: psrio: FITSTYPE is not PSRFITS", path.string()));
    }
    const auto mode = primary.value("OBS_MODE");
    if (!detail::equals_ignore_case(mode, "SEARCH")) {
        throw FormatError(
            fmt::format("{}: psrio: OBS_MODE is '{}'; only SEARCH is supported",
                        path.string(), mode));
    }
    m_source    = primary.value("SRC_NAME");
    m_telescope = primary.value("TELESCOP");
    m_frontend  = primary.find("FRONTEND").value_or("");
    m_backend   = primary.find("BACKEND").value_or("");
    m_obs_mode  = mode;
    if (const auto ra = primary.find("RA")) {
        try {
            m_raj     = astro::string_to_packed(*ra, "right ascension");
            m_has_raj = true;
        } catch (const std::invalid_argument&) {
            m_has_raj = false;
        }
    }
    if (const auto dec = primary.find("DEC")) {
        try {
            m_dej     = astro::string_to_packed(*dec, "declination");
            m_has_dej = true;
        } catch (const std::invalid_argument&) {
            m_has_dej = false;
        }
    }
    if (const auto beam = detail::find_optional_integer(primary, "IBEAM")) {
        if (*beam < std::numeric_limits<int>::min() ||
            *beam > std::numeric_limits<int>::max()) {
            throw FormatError(
                fmt::format("{}: psrio: IBEAM {} does not fit in an int",
                            path.string(), *beam));
        }
        m_beam = static_cast<int>(*beam);
    }
    if (const auto dm = detail::find_optional_number(primary, "CHAN_DM")) {
        m_dm     = *dm;
        m_has_dm = true;
    }
    m_stt_imjd = primary.integer("STT_IMJD");
    m_stt_smjd = primary.number("STT_SMJD");
    m_stt_offs = primary.number("STT_OFFS");
    m_tstart =
        static_cast<double>(m_stt_imjd) + ((m_stt_smjd + m_stt_offs) / 86400.0);
    load_geometry();
    load_frequencies();
    load_timing();
    // Free per-row timing scratch once validation is done.
    m_row_tsub = {};
    m_row_offs = {};
}

inline void SearchFile::load_geometry() {
    const auto& sub = m_table.subint;
    const int npol  = detail::narrow_positive(sub.integer("NPOL"), "NPOL");
    const int nchan = detail::narrow_positive(sub.integer("NCHAN"), "NCHAN");
    const int nbits = detail::narrow_positive(sub.integer("NBITS"), "NBITS");
    const int nsblk = detail::narrow_positive(sub.integer("NSBLK"), "NSBLK");
    if (nchan < 2) {
        throw FormatError(fmt::format(
            "psrio: search files need at least two channels, NCHAN is {}",
            nchan));
    }
    m_tsamp = sub.number("TBIN");
    if (!(m_tsamp > 0.0)) {
        throw FormatError(fmt::format(
            "psrio: TBIN {} is not a positive sample interval", m_tsamp));
    }
    m_pol_type      = sub.value("POL_TYPE");
    const auto plan = detail::plan_polarisation(m_pol_type, npol);
    m_sum           = plan.sum;
    m_input_nbits   = nbits;
    m_npol          = static_cast<std::uint64_t>(npol);
    m_nchan         = static_cast<std::uint64_t>(nchan);
    m_nsblk         = static_cast<std::uint64_t>(nsblk);
    auto data_bits  = detail::checked_mul(m_nchan, m_npol, "DATA");
    data_bits       = detail::checked_mul(data_bits, m_nsblk, "DATA");
    data_bits = detail::checked_mul(data_bits,
                                    static_cast<std::uint64_t>(nbits), "DATA");
    if (data_bits % 8U != 0U) {
        throw FormatError("psrio: DATA does not fill a whole number of bytes");
    }
    m_data_bytes = data_bits / 8U;
    if (m_options.product >= 0) {
        if (std::cmp_greater_equal(m_options.product, npol)) {
            throw ValidationError(fmt::format(
                "psrio: polarisation product {} is out of range for NPOL {}",
                m_options.product, npol));
        }
        m_sum  = false;
        m_pol0 = static_cast<std::uint64_t>(m_options.product);
    } else if (m_options.product != -1) {
        throw ValidationError(fmt::format(
            "psrio: polarisation product must be -1 or 0 .. {}, got {}",
            npol - 1, m_options.product));
    }
    // dspsr flips negative ZERO_OFF to positive before subtracting it.
    m_zero =
        m_options.scloffs
            ? static_cast<float>(std::abs(
                  detail::find_optional_number(sub, "ZERO_OFF").value_or(0.0)))
            : 0.0F;
    m_signint = detail::find_optional_integer(sub, "SIGNINT").value_or(0) != 0;
    m_data    = detail::require_column(m_table, "DATA");
    const auto count = m_npol * m_nchan;
    if (m_options.scloffs) {
        m_scales_col  = detail::require_column(m_table, "DAT_SCL");
        m_offsets_col = detail::require_column(m_table, "DAT_OFFS");
        if ((!m_scales_col.heap && m_scales_col.repeat < count) ||
            (!m_offsets_col.heap && m_offsets_col.repeat < count)) {
            throw FormatError("psrio: DAT_SCL and DAT_OFFS must cover every "
                              "polarisation and channel");
        }
    }
    if (m_options.weights) {
        m_weights_col = detail::require_column(m_table, "DAT_WTS");
        if (!m_weights_col.heap && m_weights_col.repeat < m_nchan) {
            throw FormatError("psrio: DAT_WTS must cover every channel");
        }
    }
    // The span of a row holding the fixed-width metadata columns this reader
    // actually touches. Skipping unused DAT_SCL/DAT_OFFS/DAT_WTS avoids
    // reading hundreds of KiB per row when scloffs and weights are off.
    {
        std::uint64_t first    = std::numeric_limits<std::uint64_t>::max();
        std::uint64_t last     = 0;
        const auto include_col = [&](std::string_view name) {
            const auto* column = detail::find_column(m_table, name);
            if (column == nullptr || column->heap) {
                return;
            }
            first = std::min(first, column->offset);
            last  = std::max(last, column->offset + column->nbytes);
        };
        for (const auto* name : {
                 "TSUBINT",
                 "OFFS_SUB",
                 "TEL_AZ",
                 "TEL_ZEN",
                 "DAT_FREQ",
             }) {
            include_col(name);
        }
        if (m_options.weights) {
            include_col("DAT_WTS");
        }
        if (m_options.scloffs) {
            include_col("DAT_OFFS");
            include_col("DAT_SCL");
        }
        constexpr std::uint64_t kMetaMaxBytes = std::uint64_t{8} << 20U;
        if (last > first && last - first <= kMetaMaxBytes &&
            last <= m_table.row_bytes) {
            m_meta_first = first;
            m_meta_bytes = last - first;
        }
    }
    int field = 0;
    for (std::size_t index = 0; index < m_table.columns.size(); ++index) {
        if (m_table.columns[index].name == "DATA") {
            field = static_cast<int>(index) + 1;
            break;
        }
    }
    m_layout =
        detail::make_layout(sub, m_data, field, nchan, npol, nbits, nsblk);
    m_bytes_per_sample = detail::checked_mul(m_nchan, 4U, "float sample");
}

inline void SearchFile::load_frequencies() {
    const auto& freq_col = detail::require_column(m_table, "DAT_FREQ");
    if (!freq_col.heap && freq_col.repeat < m_nchan) {
        throw FormatError("psrio: DAT_FREQ is shorter than NCHAN");
    }
    m_freqs.resize(static_cast<std::size_t>(m_nchan));
    std::vector<double> row(static_cast<std::size_t>(m_nchan));
    // One pass over the rows also collects TSUBINT and OFFS_SUB for
    // load_timing(), so each row's metadata is read once.
    const auto* durations = detail::find_column(m_table, "TSUBINT");
    const auto* centres   = detail::find_column(m_table, "OFFS_SUB");
    const auto* az_col    = detail::find_column(m_table, "TEL_AZ");
    const auto* zen_col   = detail::find_column(m_table, "TEL_ZEN");
    m_row_tsub.assign(durations != nullptr ? m_table.nrows : 0U, 0.0);
    // OFFS_SUB is only checked between rows.
    m_row_offs.assign(
        centres != nullptr && m_table.nrows >= 2U ? m_table.nrows : 0U, 0.0);
    for (std::uint64_t isub = 0; isub < m_table.nrows; ++isub) {
        if (durations != nullptr) {
            m_row_tsub[static_cast<std::size_t>(isub)] =
                scalar_at(isub, *durations);
        }
        if (!m_row_offs.empty()) {
            m_row_offs[static_cast<std::size_t>(isub)] =
                scalar_at(isub, *centres);
        }
        if (isub == 0U) {
            if (az_col != nullptr) {
                try {
                    m_az_start = scalar_at(0U, *az_col);
                } catch (const FormatError&) {
                    m_az_start.reset();
                }
            }
            if (zen_col != nullptr) {
                try {
                    m_za_start = scalar_at(0U, *zen_col);
                } catch (const FormatError&) {
                    m_za_start.reset();
                }
            }
        }
        detail::load_reals(meta_cell(isub, freq_col), freq_col, row);
        if (isub == 0U) {
            m_freqs = row;
            continue;
        }
        for (std::uint64_t chan = 0; chan < m_nchan; ++chan) {
            if (std::abs(row[static_cast<std::size_t>(chan)] -
                         m_freqs[static_cast<std::size_t>(chan)]) >
                detail::kFreqAtolMhz) {
                throw FormatError(fmt::format(
                    "psrio: DAT_FREQ of subint {} differs from subint 0",
                    isub));
            }
        }
    }
    m_fch1            = m_freqs.front();
    m_uniform         = true;
    const double step = m_freqs[1] - m_freqs[0];
    for (std::size_t chan = 1; chan < m_freqs.size(); ++chan) {
        const double delta = m_freqs[chan] - m_freqs[chan - 1];
        if (std::abs(delta - step) > detail::kFreqAtolMhz) {
            m_uniform = false;
            break;
        }
    }
    m_foff = m_uniform ? step : 0.0;
    if (m_uniform) {
        m_bandwidth = std::abs(m_foff) * static_cast<double>(m_nchan);
        m_centre =
            m_fch1 + (m_foff * (static_cast<double>(m_nchan) - 1.0) / 2.0);
    } else {
        m_bandwidth = std::abs(m_freqs.back() - m_freqs.front());
        double sum  = 0.0;
        for (const double freq : m_freqs) {
            sum += freq;
        }
        m_centre = sum / static_cast<double>(m_freqs.size());
    }
}

inline double SearchFile::scalar_at(std::uint64_t row,
                                    const detail::Column& column) {
    double value = 0.0;
    detail::load_reals(meta_cell(row, column), column,
                       std::span<double>(&value, 1));
    return value;
}

inline std::span<const std::byte>
SearchFile::meta_cell(std::uint64_t row, const detail::Column& column) {
    m_file.ensure_open();
    const bool together =
        !column.heap && m_meta_bytes > 0U && column.offset >= m_meta_first &&
        column.offset + column.nbytes <= m_meta_first + m_meta_bytes &&
        row < m_table.nrows;
    if (!together) {
        return detail::read_cell(source(), m_table, row, column, m_cell);
    }
    if (m_meta_row != row) {
        const auto start =
            m_table.data_offset +
            detail::checked_mul(row, m_table.row_bytes, "SUBINT row") +
            m_meta_first;
        m_meta.resize(static_cast<std::size_t>(m_meta_bytes));
        m_file.read_at(start, m_meta);
        m_meta_row = row;
    }
    return std::span<const std::byte>(m_meta).subspan(
        static_cast<std::size_t>(column.offset - m_meta_first),
        static_cast<std::size_t>(column.nbytes));
}

inline void SearchFile::load_timing() {
    const auto* durations = detail::find_column(m_table, "TSUBINT");
    const auto* centres   = detail::find_column(m_table, "OFFS_SUB");
    const double full     = static_cast<double>(m_nsblk) * m_tsamp;
    std::vector<double> tsub(static_cast<std::size_t>(m_table.nrows), full);
    if (durations != nullptr) {
        for (std::uint64_t row = 0; row < m_table.nrows; ++row) {
            tsub[static_cast<std::size_t>(row)] =
                m_row_tsub[static_cast<std::size_t>(row)];
            if (!(tsub[static_cast<std::size_t>(row)] > 0.0)) {
                throw FormatError(
                    "psrio: TSUBINT must be positive for every row");
            }
        }
        for (std::uint64_t row = 0; row + 1U < m_table.nrows; ++row) {
            if (std::abs(tsub[static_cast<std::size_t>(row)] - full) >
                (0.5 * m_tsamp)) {
                throw FormatError(
                    "psrio: only the last sub-integration may be short");
            }
        }
    }
    std::uint64_t trimmed =
        detail::checked_mul(m_nsblk, m_table.nrows, "sample count");
    if (durations != nullptr) {
        const double last = tsub.back();
        if (std::abs(last - full) > (0.5 * m_tsamp)) {
            const auto samples = std::llround(last / m_tsamp);
            if (samples < 1 ||
                std::cmp_greater(static_cast<std::uint64_t>(samples),
                                 m_nsblk)) {
                throw FormatError(fmt::format(
                    "psrio: last TSUBINT {} s is not between one sample and "
                    "NSBLK",
                    last));
            }
            trimmed = ((m_table.nrows - 1U) * m_nsblk) +
                      static_cast<std::uint64_t>(samples);
        }
    }
    if (m_table.nrows >= 2U) {
        if (centres == nullptr) {
            throw FormatError("psrio: split files and row gaps are not "
                              "supported (no OFFS_SUB)");
        }
        std::vector<double> offs(static_cast<std::size_t>(m_table.nrows));
        for (std::uint64_t row = 0; row < m_table.nrows; ++row) {
            offs[static_cast<std::size_t>(row)] =
                m_row_offs[static_cast<std::size_t>(row)];
        }
        for (std::uint64_t row = 0; row + 1U < m_table.nrows; ++row) {
            const double expected =
                0.5 * (tsub[static_cast<std::size_t>(row)] +
                       tsub[static_cast<std::size_t>(row + 1U)]);
            const double delta = offs[static_cast<std::size_t>(row + 1U)] -
                                 offs[static_cast<std::size_t>(row)];
            if (std::abs(delta - expected) > (0.5 * m_tsamp)) {
                throw FormatError("psrio: sub-integrations are not contiguous");
            }
        }
    }
    const auto capacity =
        detail::checked_mul(m_nsblk, m_table.nrows, "sample count");
    if (const auto nstot =
            detail::find_optional_integer(m_table.subint, "NSTOT")) {
        if (*nstot < 0 ||
            std::cmp_greater(static_cast<std::uint64_t>(*nstot), capacity)) {
            throw FormatError(fmt::format("psrio: NSTOT {} is outside 0..{}",
                                          *nstot, capacity));
        }
        // NSTOT == 0 is a common "unset" sentinel in PSRFITS templates.
        if (*nstot > 0) {
            if (std::cmp_greater(static_cast<std::uint64_t>(*nstot), trimmed)) {
                throw FormatError(fmt::format(
                    "psrio: NSTOT {} extends past the last sub-integration",
                    *nstot));
            }
            m_nsamples = static_cast<std::uint64_t>(*nstot);
        } else {
            m_nsamples = trimmed;
        }
    } else {
        m_nsamples = trimmed;
    }
}

inline sigproc::FilterbankHeader SearchFile::sigproc_header() const {
    sigproc::FilterbankHeader header;
    if (m_nchan >
        static_cast<std::uint64_t>(std::numeric_limits<std::int32_t>::max())) {
        throw ValidationError("psrio: NCHAN does not fit in a SIGPROC header");
    }
    header.nchans    = static_cast<std::int32_t>(m_nchan);
    header.nifs      = 1;
    header.nbits     = 32;
    header.tsamp     = m_tsamp;
    header.tstart    = m_tstart;
    header.data_type = 1;
    header.fch1      = m_fch1;
    if (!m_source.empty()) {
        header.source_name = m_source;
    }
    if (const auto filename = m_file.path().filename().string();
        !filename.empty()) {
        header.rawdatafile = filename;
    }
    if (m_has_raj) {
        header.src_raj = m_raj;
    }
    if (m_has_dej) {
        header.src_dej = m_dej;
    }
    if (m_az_start.has_value()) {
        header.az_start = m_az_start;
    }
    if (m_za_start.has_value()) {
        header.za_start = m_za_start;
    }
    if (m_has_dm) {
        header.refdm = m_dm;
    }
    header.ibeam = m_beam;
    if (const auto id = sigproc::telescope_id_from_name(m_telescope)) {
        header.telescope_id = id;
    }
    if (const auto id = sigproc::machine_id_from_name(m_backend)) {
        header.machine_id = id;
    }
    if (m_uniform) {
        header.foff = m_foff;
    } else {
        header.frequency_table = m_freqs;
    }
    return header;
}

inline void SearchFile::seek(std::uint64_t sample) {
    if (sample > m_nsamples) {
        throw ValidationError("psrio: seek is past the readable samples");
    }
    m_sample = sample;
}

inline void SearchFile::skip(std::int64_t delta) {
    m_sample = ::psrio::detail::apply_skip(m_sample, delta, m_nsamples);
}

inline void SearchFile::check_sample_request(std::uint64_t count,
                                             std::size_t dest_values) const {
    if (count > 0U &&
        m_nchan > (std::numeric_limits<std::uint64_t>::max() / count)) {
        throw ValidationError("psrio: requested sample block is too large");
    }
    if (dest_values != static_cast<std::size_t>(count * m_nchan)) {
        throw ValidationError(fmt::format(
            "psrio: destination has {} values but the request needs {}",
            dest_values, count * m_nchan));
    }
    if (m_sample > m_nsamples) {
        throw ValidationError(
            "psrio: reader cursor is past the readable samples");
    }
}

inline void SearchFile::load_calibration(std::uint64_t row) {
    if (m_cal_row == row) {
        return;
    }
    if (!m_options.scloffs && !m_options.weights) {
        m_cal_row = row;
        return;
    }
    const auto count = static_cast<std::size_t>(m_npol * m_nchan);
    const auto nchan = static_cast<std::size_t>(m_nchan);
    m_scales.resize(count);
    m_offsets.resize(count);
    m_weights.resize(nchan);
    if (m_options.scloffs) {
        m_scale_raw.resize(count);
        m_offset_raw.resize(count);
        detail::load_reals(meta_cell(row, m_scales_col), m_scales_col,
                           m_scale_raw);
        detail::load_reals(meta_cell(row, m_offsets_col), m_offsets_col,
                           m_offset_raw);
        for (std::size_t index = 0; index < count; ++index) {
            m_scales[index]  = static_cast<float>(m_scale_raw[index]);
            m_offsets[index] = static_cast<float>(m_offset_raw[index]);
        }
    } else {
        std::ranges::fill(m_scales, 1.0F);
        std::ranges::fill(m_offsets, 0.0F);
    }
    if (m_options.weights) {
        m_weight_raw.resize(nchan);
        detail::load_reals(meta_cell(row, m_weights_col), m_weights_col,
                           m_weight_raw);
        for (std::size_t index = 0; index < nchan; ++index) {
            m_weights[index] = static_cast<float>(m_weight_raw[index]);
        }
    } else {
        std::ranges::fill(m_weights, 1.0F);
    }
    m_cal_row = row;
}

inline std::size_t
SearchFile::element_index(std::uint64_t isamp, int pol, int chan) const {
    std::array<int, 3> coord{};
    coord[static_cast<std::size_t>(m_layout.time_axis)] =
        static_cast<int>(isamp);
    coord[static_cast<std::size_t>(m_layout.pol_axis)]  = pol;
    coord[static_cast<std::size_t>(m_layout.chan_axis)] = chan;
    const auto row = (static_cast<std::size_t>(coord[0]) *
                      static_cast<std::size_t>(m_layout.shape[1])) +
                     static_cast<std::size_t>(coord[1]);
    return (row * static_cast<std::size_t>(m_layout.shape[2])) +
           static_cast<std::size_t>(coord[2]);
}

inline std::uint8_t SearchFile::packed_byte(const std::byte* cell,
                                            std::uint64_t isamp,
                                            int pol,
                                            int byte_index) const {
    const auto slot =
        static_cast<int>(isamp / static_cast<std::uint64_t>(m_layout.group));
    const auto within =
        static_cast<int>(isamp % static_cast<std::uint64_t>(m_layout.group));
    const int element = (within * m_layout.stride) + byte_index;
    std::array<int, 3> coord{};
    coord[static_cast<std::size_t>(m_layout.time_axis)] = slot;
    coord[static_cast<std::size_t>(m_layout.pol_axis)]  = pol;
    coord[static_cast<std::size_t>(m_layout.chan_axis)] = element;
    const auto row    = (static_cast<std::size_t>(coord[0]) *
                         static_cast<std::size_t>(m_layout.shape[1])) +
                        static_cast<std::size_t>(coord[1]);
    const auto linear = (row * static_cast<std::size_t>(m_layout.shape[2])) +
                        static_cast<std::size_t>(coord[2]);
    return static_cast<std::uint8_t>(cell[linear]);
}

template <bool Signed, bool Sum, bool Calibrated>
inline void SearchFile::decode_u8(const std::byte* cell,
                                  std::uint64_t isamp,
                                  std::span<float> dest) const {
    const auto nchan = static_cast<std::size_t>(m_nchan);
    const auto pol0  = static_cast<std::size_t>(m_pol0);
    // NOLINTNEXTLINE(cppcoreguidelines-pro-type-reinterpret-cast)
    const auto* raw0 = reinterpret_cast<const std::uint8_t*>(
        cell + (((isamp * m_npol) + m_pol0) * m_nchan));
    const auto* raw1    = raw0 + nchan;
    float* out          = dest.data();
    const auto as_float = [](std::uint8_t raw) {
        if constexpr (Signed) {
            return static_cast<float>(static_cast<std::int8_t>(raw));
        } else {
            return static_cast<float>(raw);
        }
    };
    if constexpr (!Calibrated) {
        for (std::size_t chan = 0; chan < nchan; ++chan) {
            // NOLINTNEXTLINE(misc-const-correctness): summed when Sum
            float value = as_float(raw0[chan]);
            if constexpr (Sum) {
                value += as_float(raw1[chan]);
            }
            out[chan] = value;
        }
        return;
    }
    const float* scl0 = m_scales.data() + (pol0 * nchan);
    const float* off0 = m_offsets.data() + (pol0 * nchan);
    const float* scl1 = scl0 + nchan;
    const float* off1 = off0 + nchan;
    const float* wts  = m_weights.data();
    const float zero  = m_zero;
    PSRIO_OMP(omp simd)
    for (std::size_t chan = 0; chan < nchan; ++chan) {
        // NOLINTNEXTLINE(misc-const-correctness): summed when Sum
        float value = ((as_float(raw0[chan]) - zero) * scl0[chan]) + off0[chan];
        if constexpr (Sum) {
            value += ((as_float(raw1[chan]) - zero) * scl1[chan]) + off1[chan];
        }
        out[chan] = value * wts[chan];
    }
}

template <bool Sum>
inline void SearchFile::decode_i16(const std::byte* cell,
                                   std::uint64_t isamp,
                                   std::span<float> dest) const {
    const auto nchan  = static_cast<std::size_t>(m_nchan);
    const auto pol0   = static_cast<std::size_t>(m_pol0);
    const auto* base0 = cell + ((((isamp * m_npol) + m_pol0) * m_nchan) * 2U);
    const auto* base1 = base0 + (nchan * 2U);
    const bool calibrated = m_options.scloffs || m_options.weights;
    float* out            = dest.data();
    if (!calibrated) {
        for (std::size_t chan = 0; chan < nchan; ++chan) {
            // NOLINTNEXTLINE(misc-const-correctness): summed when Sum
            auto value = static_cast<float>(
                detail::load_be<std::int16_t>(base0 + (chan * 2U)));
            if constexpr (Sum) {
                value += static_cast<float>(
                    detail::load_be<std::int16_t>(base1 + (chan * 2U)));
            }
            out[chan] = value;
        }
        return;
    }
    const float* scl0 = m_scales.data() + (pol0 * nchan);
    const float* off0 = m_offsets.data() + (pol0 * nchan);
    const float* scl1 = scl0 + nchan;
    const float* off1 = off0 + nchan;
    const float* wts  = m_weights.data();
    for (std::size_t chan = 0; chan < nchan; ++chan) {
        const auto first = static_cast<float>(
            detail::load_be<std::int16_t>(base0 + (chan * 2U)));
        // NOLINTNEXTLINE(misc-const-correctness): summed when Sum
        float value = (first * scl0[chan]) + off0[chan];
        if constexpr (Sum) {
            const auto second = static_cast<float>(
                detail::load_be<std::int16_t>(base1 + (chan * 2U)));
            value += (second * scl1[chan]) + off1[chan];
        }
        out[chan] = value * wts[chan];
    }
}

template <bool Sum>
inline void SearchFile::decode_packed(const std::byte* cell,
                                      std::uint64_t isamp,
                                      std::span<float> dest) const {
    const unsigned samples_per_byte = 8U / static_cast<unsigned>(m_input_nbits);
    const unsigned mask = (1U << static_cast<unsigned>(m_input_nbits)) - 1U;
    const auto pol0     = static_cast<int>(m_pol0);
    const auto value_of = [&](int pol, std::uint64_t chan) {
        const auto byte_index = static_cast<int>(chan / samples_per_byte);
        const auto field      = static_cast<unsigned>(chan % samples_per_byte);
        const auto shift      = (samples_per_byte - 1U - field) *
                                static_cast<unsigned>(m_input_nbits);
        const unsigned digit =
            (static_cast<unsigned>(packed_byte(cell, isamp, pol, byte_index)) >>
             shift) &
            mask;
        float raw = 0.0F;
        if (m_input_nbits == 2) {
            raw = detail::kTwoBit[static_cast<std::size_t>(digit)];
        } else {
            raw = static_cast<float>(digit) - m_zero;
        }
        return raw;
    };
    const bool calibrated = m_options.scloffs || m_options.weights;
    const auto first0     = m_pol0 * m_nchan;
    const auto first1     = first0 + m_nchan;
    for (std::uint64_t chan = 0; chan < m_nchan; ++chan) {
        float value = value_of(pol0, chan);
        if (calibrated) {
            value =
                (value * m_scales[first0 + chan]) + m_offsets[first0 + chan];
        }
        if constexpr (Sum) {
            float second = value_of(pol0 + 1, chan);
            if (calibrated) {
                second = (second * m_scales[first1 + chan]) +
                         m_offsets[first1 + chan];
            }
            value += second;
        }
        dest[static_cast<std::size_t>(chan)] =
            calibrated ? (value * m_weights[chan]) : value;
    }
}

template <bool Sum>
inline void SearchFile::decode_indexed(const std::byte* cell,
                                       std::uint64_t isamp,
                                       std::span<float> dest) const {
    const auto raw_at = [&](int pol, std::uint64_t chan) {
        const auto linear = element_index(isamp, pol, static_cast<int>(chan));
        const auto* source =
            cell + (linear * static_cast<std::size_t>(m_layout.elem_size));
        if (m_input_nbits == 16) {
            return static_cast<float>(detail::load_be<std::int16_t>(source));
        }
        const auto byte   = static_cast<std::uint8_t>(source[0]);
        const float value = (m_signint && byte >= 128U)
                                ? static_cast<float>(byte) - 256.0F
                                : static_cast<float>(byte);
        return value - m_zero;
    };
    const bool calibrated = m_options.scloffs || m_options.weights;
    const auto pol0       = static_cast<int>(m_pol0);
    const auto first0     = m_pol0 * m_nchan;
    const auto first1     = first0 + m_nchan;
    for (std::uint64_t chan = 0; chan < m_nchan; ++chan) {
        float value = raw_at(pol0, chan);
        if (calibrated) {
            value =
                (value * m_scales[first0 + chan]) + m_offsets[first0 + chan];
        }
        if constexpr (Sum) {
            float second = raw_at(pol0 + 1, chan);
            if (calibrated) {
                second = (second * m_scales[first1 + chan]) +
                         m_offsets[first1 + chan];
            }
            value += second;
        }
        dest[static_cast<std::size_t>(chan)] =
            calibrated ? (value * m_weights[chan]) : value;
    }
}

inline detail::CellRange SearchFile::data_range(std::uint64_t isub) {
    m_file.ensure_open();
    const auto range = detail::locate_cell(source(), m_table, isub, m_data);
    if (range.nbytes != m_data_bytes) {
        throw FormatError(
            fmt::format("psrio: SUBINT row {} DATA has {} bytes, expected {}",
                        isub, range.nbytes, m_data_bytes));
    }
    return range;
}

inline SearchFile::Staged SearchFile::stage(std::uint64_t isub,
                                            std::uint64_t isamp) {
    if (m_stage_row == isub && isamp >= m_stage_first && isamp < m_stage_last) {
        return Staged{
            .cell  = m_stage.data(),
            .first = m_stage_first,
            .last  = m_stage_last,
        };
    }
    m_file.ensure_open();
    // A read that starts where the staged run ends is part of a forward
    // stream; only then does the rest of the file stream ahead.
    const bool forward =
        (m_stage_row == isub && isamp == m_stage_last) ||
        (m_stage_last == m_nsblk && m_stage_row + 1U == isub && isamp == 0U) ||
        (m_stage_row == std::numeric_limits<std::uint64_t>::max() &&
         isub == 0U && isamp == 0U);
    m_file.release(m_stage_range.offset, m_stage_range.nbytes);
    if (m_data_row != isub) {
        m_data_range = data_range(isub);
        m_data_row   = isub;
    }
    std::uint64_t first = 0;
    std::uint64_t last  = m_nsblk;
    detail::CellRange range{m_data_range};
    if (m_layout.fast) {
        // Time-major rows: read a run of whole samples.
        const auto stride = m_data_bytes / m_nsblk;
        const auto run    = std::max<std::uint64_t>(kStageBytes / stride, 1U);
        first             = isamp;
        last              = std::min(m_nsblk, isamp + run);
        range.offset      = m_data_range.offset + (first * stride);
        range.nbytes      = (last - first) * stride;
    }
    // Other layouts interleave time with the other axes: read the row.
    m_stage.resize(static_cast<std::size_t>(range.nbytes));
    if (m_ahead != nullptr && !m_ahead->covers(range.offset, range.nbytes)) {
        m_ahead.reset();
    }
    // Time-major rows stored in the table stream through m_ahead.
    if (m_ahead == nullptr && forward && m_prefetch && m_layout.fast &&
        !m_data.heap) {
        start_read_ahead(range.offset);
    }
    if (m_ahead != nullptr) {
        m_ahead->read(range.offset, m_stage);
    } else {
        m_file.read_at(range.offset, m_stage);
    }
    m_stage_row   = isub;
    m_stage_first = first;
    m_stage_last  = last;
    m_stage_range = range;
    return Staged{
        .cell  = m_stage.data(),
        .first = m_stage_first,
        .last  = m_stage_last,
    };
}

inline void SearchFile::start_read_ahead(std::uint64_t offset) {
    m_file.ensure_open();
    const auto last = data_range(m_table.nrows - 1U);
    m_ahead         = std::make_unique<::psrio::detail::ReadAhead>(
        m_file.duplicate(), offset, last.offset + last.nbytes,
        m_read_ahead_blocks);
}

inline void SearchFile::prefetch_start() {
    if (!m_prefetch || !m_layout.fast || m_data.heap || m_ahead != nullptr ||
        m_table.nrows == 0U) {
        return;
    }
    const auto first = data_range(0U);
    start_read_ahead(first.offset);
}

inline void SearchFile::set_prefetch(bool enable) {
    m_prefetch = enable;
    if (!enable) {
        m_ahead.reset();
    }
}

inline void SearchFile::release_buffers() {
    m_ahead.reset();
    m_file.release(m_stage_range.offset, m_stage_range.nbytes);
    m_stage       = {};
    m_floats      = {};
    m_cell        = {};
    m_meta        = {};
    m_meta_row    = std::numeric_limits<std::uint64_t>::max();
    m_stage_row   = std::numeric_limits<std::uint64_t>::max();
    m_stage_first = 0;
    m_stage_last  = 0;
    m_stage_range = {};
}

inline void SearchFile::close_handle() {
    release_buffers();
    m_scale_raw  = {};
    m_offset_raw = {};
    m_weight_raw = {};
    m_scales     = {};
    m_offsets    = {};
    m_weights    = {};
    m_cal_row    = std::numeric_limits<std::uint64_t>::max();
    m_data_row   = std::numeric_limits<std::uint64_t>::max();
    m_file.close_fd();
}

inline void SearchFile::decode_run(const std::byte* cell,
                                   std::uint64_t isamp,
                                   std::uint64_t count,
                                   std::span<float> dest) const {
    // Samples decode independently into disjoint rows of dest.
    [[maybe_unused]] const bool wide =
        ::psrio::detail::parallel_worthwhile(count * m_nchan);
    const auto nchan = static_cast<std::size_t>(m_nchan);
    const auto run   = [&](auto&& decode) {
        PSRIO_OMP(omp parallel for default(none) schedule(static) if (wide)
                      shared(count, decode, cell, isamp, dest, nchan))
        for (std::uint64_t index = 0; index < count; ++index) {
            decode(
                cell, isamp + index,
                dest.subspan(static_cast<std::size_t>(index) * nchan, nchan));
        }
    };
    if (m_input_nbits == 8 && m_layout.fast) {
        const bool calibrated = m_options.scloffs || m_options.weights;
        const auto run_u8     = [&]<bool Signed, bool Sum> {
            if (calibrated) {
                run([this](const std::byte* bytes, std::uint64_t sample,
                           std::span<float> out) {
                    decode_u8<Signed, Sum, true>(bytes, sample, out);
                });
            } else {
                run([this](const std::byte* bytes, std::uint64_t sample,
                           std::span<float> out) {
                    decode_u8<Signed, Sum, false>(bytes, sample, out);
                });
            }
        };
        if (m_sum) {
            if (m_signint) {
                run_u8.template operator()<true, true>();
            } else {
                run_u8.template operator()<false, true>();
            }
        } else if (m_signint) {
            run_u8.template operator()<true, false>();
        } else {
            run_u8.template operator()<false, false>();
        }
        return;
    }
    if (m_input_nbits == 16 && m_layout.fast) {
        if (m_sum) {
            run([this](const std::byte* bytes, std::uint64_t sample,
                       std::span<float> out) {
                decode_i16<true>(bytes, sample, out);
            });
        } else {
            run([this](const std::byte* bytes, std::uint64_t sample,
                       std::span<float> out) {
                decode_i16<false>(bytes, sample, out);
            });
        }
        return;
    }
    if (m_input_nbits < 8) {
        if (m_sum) {
            run([this](const std::byte* bytes, std::uint64_t sample,
                       std::span<float> out) {
                decode_packed<true>(bytes, sample, out);
            });
        } else {
            run([this](const std::byte* bytes, std::uint64_t sample,
                       std::span<float> out) {
                decode_packed<false>(bytes, sample, out);
            });
        }
        return;
    }
    if (m_sum) {
        run([this](const std::byte* bytes, std::uint64_t sample,
                   std::span<float> out) {
            decode_indexed<true>(bytes, sample, out);
        });
    } else {
        run([this](const std::byte* bytes, std::uint64_t sample,
                   std::span<float> out) {
            decode_indexed<false>(bytes, sample, out);
        });
    }
}

inline std::uint64_t SearchFile::read_samples(std::uint64_t count,
                                              std::span<float> dest) {
    check_sample_request(count, dest.size());
    const auto to_read  = std::min(count, m_nsamples - m_sample);
    std::uint64_t index = 0;
    while (index < to_read) {
        const auto global = m_sample + index;
        const auto isub   = global / m_nsblk;
        const auto isamp  = global % m_nsblk;
        load_calibration(isub);
        const auto staged = stage(isub, isamp);
        const auto run    = std::min(to_read - index, staged.last - isamp);
        decode_run(staged.cell, isamp - staged.first, run,
                   dest.subspan(static_cast<std::size_t>(index * m_nchan),
                                static_cast<std::size_t>(run * m_nchan)));
        index += run;
    }
    if (m_fswap && to_read > 0U) {
        reverse_unpacked_channels(
            dest.first(static_cast<std::size_t>(to_read * m_nchan)), to_read,
            m_nchan, 1U);
    }
    m_sample += to_read;
    return to_read;
}

inline std::uint64_t
SearchFile::read_samples(std::uint64_t /*count*/,
                         std::span<std::uint8_t> /*dest*/) {
    throw ValidationError("psrio: calibrated PSRFITS samples unpack to float");
}

inline std::uint64_t
SearchFile::read_samples(std::uint64_t /*count*/,
                         std::span<std::uint16_t> /*dest*/) {
    throw ValidationError("psrio: calibrated PSRFITS samples unpack to float");
}

inline std::uint64_t SearchFile::read_block(std::uint64_t count,
                                            std::span<std::byte> dest) {
    if (count > 0U && m_bytes_per_sample >
                          (std::numeric_limits<std::uint64_t>::max() / count)) {
        throw ValidationError("psrio: requested sample block is too large");
    }
    const auto bytes_needed = count * m_bytes_per_sample;
    if (dest.size() < bytes_needed) {
        throw ValidationError(fmt::format(
            "psrio: destination has {} bytes but the block needs {}",
            dest.size(), bytes_needed));
    }
    if (m_sample > m_nsamples) {
        throw ValidationError(
            "psrio: reader cursor is past the readable samples");
    }
    const auto to_read = std::min(count, m_nsamples - m_sample);
    m_floats.resize(static_cast<std::size_t>(to_read * m_nchan));
    const auto got = read_samples(to_read, std::span<float>(m_floats));
    auto* out      = dest.data();
    if constexpr (std::endian::native == std::endian::little) {
        std::memcpy(out, m_floats.data(), m_floats.size() * sizeof(float));
    } else {
        for (std::size_t index = 0; index < m_floats.size(); ++index) {
            auto bits = std::bit_cast<std::uint32_t>(m_floats[index]);
            bits      = ::psrio::detail::to_little_endian(bits);
            std::memcpy(out + (index * 4U), &bits, sizeof(bits));
        }
    }
    return got;
}

inline std::uint64_t SearchFile::read_bytes(std::uint64_t nbytes,
                                            std::span<std::byte> dest) {
    if (nbytes == 0U || m_bytes_per_sample == 0U ||
        nbytes % m_bytes_per_sample != 0U) {
        throw ValidationError("psrio: byte request must be a positive multiple "
                              "of the sample stride");
    }
    if (dest.size() != nbytes) {
        throw ValidationError(fmt::format(
            "psrio: destination has {} bytes but the request needs {}",
            dest.size(), nbytes));
    }
    const auto samples = nbytes / m_bytes_per_sample;
    if (m_sample > m_nsamples || samples > m_nsamples - m_sample) {
        throw ValidationError(
            "psrio: byte request extends past the readable samples");
    }
    const auto got = read_block(samples, dest);
    return got * m_bytes_per_sample;
}

inline std::vector<std::byte> SearchFile::read_block(std::uint64_t count) {
    std::vector<std::byte> out(
        static_cast<std::size_t>(count * m_bytes_per_sample));
    const auto actual = read_block(count, out);
    out.resize(static_cast<std::size_t>(actual * m_bytes_per_sample));
    return out;
}

inline std::vector<std::byte> SearchFile::read_bytes(std::uint64_t nbytes) {
    std::vector<std::byte> out(static_cast<std::size_t>(nbytes));
    const auto actual = read_bytes(nbytes, out);
    out.resize(static_cast<std::size_t>(actual));
    return out;
}

} // namespace psrio::formats::psrfits

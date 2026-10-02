#include "psrio/astro.hpp"
#include "psrio/block_source.hpp"
#include "psrio/detail/exceptions.hpp"
#include "psrio/formats/fbh5.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <hdf5.h>
#include <span>
#include <string>
#include <type_traits>
#include <vector>

using Catch::Matchers::ContainsSubstring;
using psrio::formats::fbh5::FilterbankFile;
using psrio::formats::fbh5::is_hdf5_filterbank;

namespace {

class Hid {
public:
    explicit Hid(hid_t id) noexcept : m_id(id) {}
    ~Hid() { reset(); }
    Hid(const Hid&)            = delete;
    Hid& operator=(const Hid&) = delete;
    [[nodiscard]] hid_t get() const noexcept { return m_id; }
    [[nodiscard]] bool valid() const noexcept { return m_id >= 0; }

private:
    void reset() noexcept {
        if (m_id < 0) {
            return;
        }
        switch (H5Iget_type(m_id)) {
        case H5I_FILE:
            H5Fclose(m_id);
            break;
        case H5I_DATASET:
            H5Dclose(m_id);
            break;
        case H5I_DATASPACE:
            H5Sclose(m_id);
            break;
        case H5I_ATTR:
            H5Aclose(m_id);
            break;
        case H5I_DATATYPE:
            H5Tclose(m_id);
            break;
        case H5I_GENPROP_LST:
            H5Pclose(m_id);
            break;
        default:
            break;
        }
        m_id = H5I_INVALID_HID;
    }

    hid_t m_id;
};

class TempFile {
public:
    explicit TempFile(std::string suffix) {
        static int sequence = 0;
        m_path =
            std::filesystem::temp_directory_path() /
            ("psrio-fbh5-" + std::to_string(++sequence) + "-" +
             std::to_string(
                 std::chrono::steady_clock::now().time_since_epoch().count()) +
             suffix);
    }

    ~TempFile() { std::filesystem::remove(m_path); }

    TempFile(const TempFile&)            = delete;
    TempFile& operator=(const TempFile&) = delete;

    [[nodiscard]] const std::filesystem::path& path() const { return m_path; }

private:
    std::filesystem::path m_path;
};

struct HeaderSpec {
    std::string source_name = "B0329+54";
    double tstart           = 58000.0;
    double tsamp            = 1.0e-3;
    double fch1             = 1400.0;
    double foff             = -1.0;
    int nbits               = 8;
    int telescope_id        = 6;
    int ibeam               = 3;
    double src_raj          = 123456.7;
    bool write_tsamp        = true;
    bool write_data         = true;
    bool gzip               = false;
    bool foff_on_file       = false;
};

void write_string_attr(hid_t location,
                       const char* name,
                       const std::string& text) {
    Hid space{H5Screate(H5S_SCALAR)};
    Hid type{H5Tcopy(H5T_C_S1)};
    REQUIRE(space.valid());
    REQUIRE(type.valid());
    REQUIRE(H5Tset_size(type.get(), H5T_VARIABLE) >= 0);
    REQUIRE(H5Tset_cset(type.get(), H5T_CSET_UTF8) >= 0);
    Hid attr{H5Acreate2(location, name, type.get(), space.get(), H5P_DEFAULT,
                        H5P_DEFAULT)};
    REQUIRE(attr.valid());
    const char* value = text.c_str();
    REQUIRE(H5Awrite(attr.get(), type.get(), &value) >= 0);
}

void write_double_attr(hid_t location, const char* name, double value) {
    Hid space{H5Screate(H5S_SCALAR)};
    Hid attr{H5Acreate2(location, name, H5T_NATIVE_DOUBLE, space.get(),
                        H5P_DEFAULT, H5P_DEFAULT)};
    REQUIRE(attr.valid());
    REQUIRE(H5Awrite(attr.get(), H5T_NATIVE_DOUBLE, &value) >= 0);
}

void write_int_attr(hid_t location, const char* name, int value) {
    Hid space{H5Screate(H5S_SCALAR)};
    Hid attr{H5Acreate2(location, name, H5T_NATIVE_INT, space.get(),
                        H5P_DEFAULT, H5P_DEFAULT)};
    REQUIRE(attr.valid());
    REQUIRE(H5Awrite(attr.get(), H5T_NATIVE_INT, &value) >= 0);
}

template <typename T> hid_t native_type() {
    if constexpr (std::same_as<T, std::uint8_t>) {
        return H5T_NATIVE_UINT8;
    } else if constexpr (std::same_as<T, std::uint16_t>) {
        return H5T_NATIVE_UINT16;
    } else if constexpr (std::same_as<T, std::uint32_t>) {
        return H5T_NATIVE_UINT32;
    } else if constexpr (std::same_as<T, std::int32_t>) {
        return H5T_NATIVE_INT32;
    } else if constexpr (std::same_as<T, float>) {
        return H5T_NATIVE_FLOAT;
    } else {
        static_assert(sizeof(T) == 0, "unsupported HDF5 test type");
        return H5I_INVALID_HID;
    }
}

template <typename T>
void write_filterbank(const std::filesystem::path& path,
                      std::size_t ntime,
                      std::size_t nbeams,
                      std::size_t nchans,
                      const HeaderSpec& spec,
                      const std::vector<T>& cube) {
    if (spec.write_data) {
        REQUIRE(cube.size() == ntime * nbeams * nchans);
    }
    Hid file{H5Fcreate(path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT)};
    REQUIRE(file.valid());
    write_string_attr(file.get(), "CLASS", "FILTERBANK");
    write_string_attr(file.get(), "VERSION", "1.0");
    if (!spec.write_data) {
        return;
    }

    const std::array<hsize_t, 3> dims{ntime, nbeams, nchans};
    Hid space{H5Screate_simple(3, dims.data(), nullptr)};
    REQUIRE(space.valid());
    Hid props{H5Pcreate(H5P_DATASET_CREATE)};
    REQUIRE(props.valid());
    if (spec.gzip) {
        REQUIRE(H5Zfilter_avail(H5Z_FILTER_DEFLATE) > 0);
        const std::array<hsize_t, 3> chunk{std::min<hsize_t>(ntime, 2), 1,
                                           nchans};
        REQUIRE(H5Pset_chunk(props.get(), 3, chunk.data()) >= 0);
        REQUIRE(H5Pset_deflate(props.get(), 4) >= 0);
    }
    Hid data{H5Dcreate2(file.get(), "data", native_type<T>(), space.get(),
                        H5P_DEFAULT, props.get(), H5P_DEFAULT)};
    REQUIRE(data.valid());
    REQUIRE(H5Dwrite(data.get(), native_type<T>(), H5S_ALL, H5S_ALL,
                     H5P_DEFAULT, cube.data()) >= 0);

    write_string_attr(data.get(), "source_name", spec.source_name);
    write_double_attr(data.get(), "tstart", spec.tstart);
    if (spec.write_tsamp) {
        write_double_attr(data.get(), "tsamp", spec.tsamp);
    }
    write_double_attr(data.get(), "fch1", spec.fch1);
    if (spec.foff_on_file) {
        write_double_attr(file.get(), "foff", spec.foff);
    } else {
        write_double_attr(data.get(), "foff", spec.foff);
    }
    write_int_attr(data.get(), "nchans", static_cast<int>(nchans));
    write_int_attr(data.get(), "nifs", static_cast<int>(nbeams));
    write_int_attr(data.get(), "nbits", spec.nbits);
    write_int_attr(data.get(), "telescope_id", spec.telescope_id);
    write_int_attr(data.get(), "ibeam", spec.ibeam);
    write_double_attr(data.get(), "src_raj", spec.src_raj);
}

} // namespace

static_assert(psrio::concepts::BlockReader<FilterbankFile>);

TEST_CASE("HDF5 filterbank header round-trips", "[hdf5]") {
    constexpr std::size_t kNTime  = 2;
    constexpr std::size_t kNBeams = 1;
    constexpr std::size_t kNChans = 4;
    const TempFile file(".h5");
    HeaderSpec spec;
    spec.foff = -0.5;
    const std::vector<std::uint8_t> cube(kNTime * kNBeams * kNChans, 0);
    write_filterbank(file.path(), kNTime, kNBeams, kNChans, spec, cube);

    const FilterbankFile source(file.path());
    const auto& header = source.header();
    CHECK(header.source_name == "B0329+54");
    CHECK(header.version == "1.0");
    CHECK(header.tstart == 58000.0);
    CHECK(header.tsamp == 1.0e-3);
    CHECK(header.fch1 == 1400.0);
    CHECK(header.foff == -0.5);
    CHECK(header.nchans == 4U);
    CHECK(header.nifs == 1U);
    CHECK(header.nbeams == 1U);
    CHECK(header.nbits == 8);
    CHECK(header.nsamples == 2U);
    CHECK(header.telescope_id == 6);
    CHECK(header.has_telescope_id);
    CHECK(header.ibeam == 3);
    CHECK(header.has_src_raj);
    CHECK(header.src_raj == 123456.7);
    CHECK(source.nchans() == kNChans);
    CHECK(source.nifs() == 1U);
    CHECK(source.sample_type() == psrio::SampleType::kUInt8);
    CHECK(source.bytes_per_sample() == kNChans);
    CHECK(source.tsamp() == header.tsamp);
    CHECK(source.beam() == 3);
    CHECK(source.fch1() == 1400.0);
    CHECK(source.foff() == -0.5);
    CHECK(source.has_nsamples());
    CHECK(source.spectra_rate() == 1000.0);
    CHECK(source.utc_start() == psrio::astro::mjd_to_time(58000.0));
}

TEST_CASE("HDF5 filterbank reads a channel ramp in file order", "[hdf5]") {
    constexpr std::size_t kNTime  = 3;
    constexpr std::size_t kNBeams = 1;
    constexpr std::size_t kNChans = 4;
    const TempFile file(".h5");
    std::vector<std::uint8_t> cube(kNTime * kNChans);
    for (std::size_t t = 0; t < kNTime; ++t) {
        for (std::size_t c = 0; c < kNChans; ++c) {
            cube[(t * kNChans) + c] =
                static_cast<std::uint8_t>((t * 10) + c + 1);
        }
    }
    write_filterbank(file.path(), kNTime, kNBeams, kNChans, HeaderSpec{}, cube);

    FilterbankFile source(file.path());
    std::vector<std::byte> gulp(kNTime * kNChans);
    CHECK(source.read_block(kNTime, gulp) == kNTime);
    CHECK(std::memcmp(gulp.data(), cube.data(), cube.size()) == 0);
    std::vector<std::byte> extra(kNChans, std::byte{0xFF});
    CHECK(source.read_block(1, extra) == 0U);

    source.rewind();
    CHECK(source.read_block(2, gulp) == 2U);
    source.skip(-1);
    REQUIRE(source.tell() == 1U);
    CHECK(source.read_block(1, std::span<std::byte>{gulp}.first(kNChans)) ==
          1U);
    CHECK(gulp[0] == std::byte{11});
    CHECK(gulp[3] == std::byte{14});

    source.rewind();
    std::vector<float> floats(kNChans);
    CHECK(source.read_samples(1, floats) == 1U);
    CHECK(floats[0] == 1.0F);
    CHECK(floats[3] == 4.0F);

    psrio::BlockSource block{FilterbankFile{file.path()}};
    CHECK(block.nchans() == kNChans);
    CHECK(block.read_block(1, std::span<std::byte>{gulp}.first(kNChans)) == 1U);
    CHECK(gulp[0] == std::byte{1});
}

TEST_CASE("HDF5 filterbank gulps concatenate without overlap or gap",
          "[hdf5]") {
    constexpr std::size_t kNTime  = 5;
    constexpr std::size_t kNBeams = 1;
    constexpr std::size_t kNChans = 3;
    const TempFile file(".h5");
    std::vector<std::uint8_t> cube(kNTime * kNChans);
    for (std::size_t index = 0; index < cube.size(); ++index) {
        cube[index] = static_cast<std::uint8_t>(index + 1);
    }
    HeaderSpec spec;
    spec.gzip = true;
    write_filterbank(file.path(), kNTime, kNBeams, kNChans, spec, cube);

    FilterbankFile source(file.path());
    std::vector<std::uint8_t> got;
    got.reserve(cube.size());
    std::vector<std::byte> gulp(2 * kNChans);
    std::size_t reads = 0;
    while (true) {
        const auto nread = source.read_block(2, gulp);
        if (nread == 0U) {
            break;
        }
        REQUIRE(nread <= 2U);
        for (std::uint64_t index = 0; index < nread * kNChans; ++index) {
            got.push_back(static_cast<std::uint8_t>(
                gulp[static_cast<std::size_t>(index)]));
        }
        ++reads;
    }
    CHECK(reads == 3U);
    REQUIRE(got.size() == cube.size());
    CHECK(std::memcmp(got.data(), cube.data(), cube.size()) == 0);
}

TEST_CASE("negative foff keeps sigproc channel order", "[hdf5]") {
    constexpr std::size_t kNTime  = 1;
    constexpr std::size_t kNBeams = 1;
    constexpr std::size_t kNChans = 4;
    const TempFile file(".h5");
    const std::vector<std::uint8_t> cube{10, 20, 30, 40};
    HeaderSpec spec;
    spec.fch1 = 1600.0;
    spec.foff = -2.5;
    write_filterbank(file.path(), kNTime, kNBeams, kNChans, spec, cube);

    FilterbankFile source(file.path());
    CHECK(source.fch1() == 1600.0);
    CHECK(source.foff() == -2.5);
    CHECK(source.header().fch1 == 1600.0);
    CHECK(source.header().foff == -2.5);
    std::vector<std::byte> gulp(kNChans);
    CHECK(source.read_block(1, gulp) == 1U);
    CHECK(gulp[0] == std::byte{10});
    CHECK(gulp[1] == std::byte{20});
    CHECK(gulp[3] == std::byte{40});
}

TEST_CASE("HDF5 filterbank reads the selected beam", "[hdf5]") {
    constexpr std::size_t kNTime  = 2;
    constexpr std::size_t kNBeams = 2;
    constexpr std::size_t kNChans = 3;
    const TempFile file(".h5");
    std::vector<std::uint8_t> cube(kNTime * kNBeams * kNChans, 0);
    for (std::size_t t = 0; t < kNTime; ++t) {
        for (std::size_t c = 0; c < kNChans; ++c) {
            cube[(t * kNBeams * kNChans) + c]           = 1;
            cube[(t * kNBeams * kNChans) + kNChans + c] = 2;
        }
    }
    write_filterbank(file.path(), kNTime, kNBeams, kNChans, HeaderSpec{}, cube);

    FilterbankFile beam0(file.path());
    CHECK(beam0.header().nifs == 2U);
    CHECK(beam0.nifs() == 1U);
    CHECK(beam0.beam() == 1);
    std::vector<std::byte> gulp(kNTime * kNChans);
    CHECK(beam0.read_block(kNTime, gulp) == kNTime);
    CHECK(gulp[0] == std::byte{1});
    CHECK(gulp[kNChans] == std::byte{1});

    FilterbankFile beam1(file.path(), 1);
    CHECK(beam1.beam() == 2);
    std::fill(gulp.begin(), gulp.end(), std::byte{0});
    CHECK(beam1.read_block(kNTime, gulp) == kNTime);
    CHECK(gulp[0] == std::byte{2});
    CHECK(gulp[(kNChans * 2) - 1] == std::byte{2});

    CHECK_THROWS_WITH(FilterbankFile(file.path(), 2),
                      ContainsSubstring("out of range"));
    CHECK_THROWS_AS(FilterbankFile(file.path(), 2), psrio::ValidationError);
}

TEST_CASE("HDF5 filterbank reads float32 samples", "[hdf5]") {
    constexpr std::size_t kNTime  = 2;
    constexpr std::size_t kNBeams = 1;
    constexpr std::size_t kNChans = 2;
    const TempFile file(".h5");
    const std::vector<float> cube{1.0F, 2.0F, 3.0F, 4.0F};
    HeaderSpec spec;
    spec.nbits = 32;
    write_filterbank(file.path(), kNTime, kNBeams, kNChans, spec, cube);

    FilterbankFile source(file.path());
    CHECK(source.sample_type() == psrio::SampleType::kFloat32);
    CHECK(source.bytes_per_sample() == kNChans * sizeof(float));
    std::vector<std::byte> bytes(cube.size() * sizeof(float));
    CHECK(source.read_block(kNTime, bytes) == kNTime);
    std::vector<float> got(cube.size());
    std::memcpy(got.data(), bytes.data(), bytes.size());
    CHECK(got == cube);

    const TempFile negative(".h5");
    spec.nbits = -32;
    write_filterbank(negative.path(), kNTime, kNBeams, kNChans, spec, cube);
    FilterbankFile signed_bits(negative.path());
    CHECK(signed_bits.sample_type() == psrio::SampleType::kFloat32);
    CHECK(signed_bits.nbits() == -32);
}

TEST_CASE("HDF5 filterbank reads uint16 and uint32 samples", "[hdf5]") {
    constexpr std::size_t kNTime  = 1;
    constexpr std::size_t kNBeams = 1;
    constexpr std::size_t kNChans = 3;
    const TempFile narrow(".h5");
    const std::vector<std::uint16_t> cube16{1000, 2000, 3000};
    HeaderSpec spec;
    spec.nbits = 16;
    write_filterbank(narrow.path(), kNTime, kNBeams, kNChans, spec, cube16);

    FilterbankFile source(narrow.path());
    CHECK(source.sample_type() == psrio::SampleType::kUInt16);
    std::vector<std::byte> bytes(cube16.size() * sizeof(std::uint16_t));
    CHECK(source.read_block(1, bytes) == 1U);
    std::vector<std::uint16_t> got(cube16.size());
    std::memcpy(got.data(), bytes.data(), bytes.size());
    CHECK(got == cube16);
    source.rewind();
    std::vector<float> floats(kNChans);
    CHECK(source.read_samples(1, floats) == 1U);
    CHECK(floats[0] == 1000.0F);
    CHECK(floats[2] == 3000.0F);

    const TempFile wide(".h5");
    const std::vector<std::uint32_t> cube32{10, 20, 30};
    spec.nbits = 32;
    write_filterbank(wide.path(), kNTime, kNBeams, kNChans, spec, cube32);
    FilterbankFile wide_source(wide.path());
    CHECK(wide_source.sample_type() == psrio::SampleType::kUInt32);
    std::vector<float> wide_floats(kNChans);
    CHECK(wide_source.read_samples(1, wide_floats) == 1U);
    CHECK(wide_floats[1] == 20.0F);
}

TEST_CASE("HDF5 filterbank missing tsamp throws", "[hdf5]") {
    const TempFile file(".h5");
    HeaderSpec spec;
    spec.write_tsamp = false;
    const std::vector<std::uint8_t> cube(4, 1);
    write_filterbank(file.path(), 1, 1, 4, spec, cube);
    CHECK_THROWS_AS(FilterbankFile(file.path()), psrio::FormatError);
    CHECK_THROWS_WITH(FilterbankFile(file.path()), ContainsSubstring("tsamp"));
}

TEST_CASE("HDF5 filterbank missing data dataset throws", "[hdf5]") {
    const TempFile file(".h5");
    HeaderSpec spec;
    spec.write_data = false;
    write_filterbank(file.path(), 0, 0, 0, spec, std::vector<std::uint8_t>{});
    CHECK_THROWS_AS(FilterbankFile(file.path()), psrio::FormatError);
    CHECK_THROWS_WITH(FilterbankFile(file.path()), ContainsSubstring("data"));
}

TEST_CASE("HDF5 filterbank rejects an unsupported dtype", "[hdf5]") {
    const TempFile file(".h5");
    HeaderSpec spec;
    spec.nbits = 32;
    const std::vector<std::int32_t> cube{1, 2, 3, 4};
    write_filterbank(file.path(), 1, 1, 4, spec, cube);
    CHECK_THROWS_AS(FilterbankFile(file.path()), psrio::FormatError);
    CHECK_THROWS_WITH(FilterbankFile(file.path()),
                      ContainsSubstring("unsupported HDF5 sample type"));
}

TEST_CASE("HDF5 magic selects the filterbank reader without a .h5 suffix",
          "[hdf5]") {
    const TempFile file(".fil");
    const std::vector<std::uint8_t> cube{4, 5, 6, 7};
    write_filterbank(file.path(), 1, 1, 4, HeaderSpec{}, cube);
    CHECK(is_hdf5_filterbank(file.path()));
    CHECK(is_hdf5_filterbank("anything.h5"));
    CHECK(is_hdf5_filterbank("anything.HDF5"));
    CHECK_FALSE(is_hdf5_filterbank("missing-not-hdf5.fil"));

    const TempFile plain(".fil");
    {
        std::ofstream stream(plain.path(), std::ios::binary);
        stream << "not hdf5";
    }
    CHECK_FALSE(is_hdf5_filterbank(plain.path()));

    FilterbankFile source(file.path());
    std::vector<std::byte> gulp(4);
    CHECK(source.read_block(1, gulp) == 1U);
    CHECK(gulp[0] == std::byte{4});
    CHECK(gulp[3] == std::byte{7});
}

TEST_CASE("HDF5 foff may live on the file instead of the dataset", "[hdf5]") {
    const TempFile file(".h5");
    HeaderSpec spec;
    spec.foff_on_file = true;
    spec.foff         = 0.25;
    const std::vector<std::uint8_t> cube(4, 1);
    write_filterbank(file.path(), 1, 1, 4, spec, cube);
    FilterbankFile source(file.path());
    CHECK(source.foff() == 0.25);
}

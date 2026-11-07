#pragma once

#include "psrio/block_source.hpp"
#include "psrio/common/types.hpp"
#include "psrio/detail/exceptions.hpp"
#include "psrio/formats/sigproc/reader.hpp"

#include <algorithm>
#include <array>
#include <cctype>
#include <cstdint>
#include <filesystem>
#include <fmt/format.h>
#include <fstream>
#include <string>

namespace psrio {

namespace open_detail {

inline constexpr std::array<char, 8> kHdf5Signature{
    '\x89', 'H', 'D', 'F', '\r', '\n', '\x1a', '\n',
};

/// True when @p path has an HDF5 extension or starts with the HDF5 signature.
[[nodiscard]] inline bool looks_like_hdf5(const std::filesystem::path& path) {
    std::string ext = path.extension().string();
    std::ranges::transform(ext, ext.begin(), [](unsigned char c) {
        return static_cast<char>(std::tolower(c));
    });
    if (ext == ".h5" || ext == ".hdf5") {
        return true;
    }
    std::array<char, kHdf5Signature.size()> head{};
    std::ifstream in(path, std::ios::binary);
    in.read(head.data(), static_cast<std::streamsize>(head.size()));
    return in.gcount() == static_cast<std::streamsize>(head.size()) &&
           head == kHdf5Signature;
}

} // namespace open_detail

/// Open a filterbank file as a type-erased BlockSource.
///
/// SIGPROC `.fil` files open as FilterbankReader. An HDF5 filterbank needs
/// the opt-in FBH5 reader: include `<psrio/formats/fbh5.hpp>` and call
/// `formats::fbh5::open_block_source`.
/// @param beam Beam index. SIGPROC files hold one beam, so it must be 0.
/// @throws ValidationError if @p beam is not 0.
/// @throws FormatError if @p path is an HDF5 file or not SIGPROC.
/// @throws IoError if @p path cannot be opened.
[[nodiscard]] inline BlockSource
open_block_source(const std::filesystem::path& path,
                  std::uint64_t beam = 0,
                  BitOrder bit_order = BitOrder::kLsbFirst) {
    if (open_detail::looks_like_hdf5(path)) {
        throw FormatError(fmt::format(
            "psrio: {} is HDF5; include <psrio/formats/fbh5.hpp> and use "
            "formats::fbh5::open_block_source",
            path.string()));
    }
    if (beam != 0U) {
        throw ValidationError(fmt::format(
            "psrio: {} is a single-beam SIGPROC file; beam {} does not exist",
            path.string(), beam));
    }
    return BlockSource(formats::sigproc::FilterbankReader(path, bit_order));
}

} // namespace psrio

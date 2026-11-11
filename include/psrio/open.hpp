#pragma once

#include "psrio/block_source.hpp"
#include "psrio/common/types.hpp"
#include "psrio/detail/exceptions.hpp"
#include "psrio/formats/psrfits.hpp"
#include "psrio/formats/sigproc/reader.hpp"
#include "psrio/formats/sigproc/stream.hpp"

#include <algorithm>
#include <array>
#include <cctype>
#include <cstdint>
#include <filesystem>
#include <fmt/format.h>
#include <fstream>
#include <ios>
#include <span>
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

/// Open a filterbank or search-mode PSRFITS file as a type-erased BlockSource.
///
/// SIGPROC `.fil` files open as FilterbankReader; search-mode PSRFITS `.sf`,
/// `.fits`, and `.sfits` files open as `formats::psrfits::SearchFile`. An HDF5
/// filterbank needs the opt-in FBH5 reader: include `<psrio/formats/fbh5.hpp>`
/// and call `formats::fbh5::open_block_source`.
/// @param beam Beam index. SIGPROC and PSRFITS files hold one beam, so it must
///        be 0.
/// @throws ValidationError if @p beam is not 0.
/// @throws FormatError if @p path is an HDF5 file or not a supported format.
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
            "psrio: {} is a single-beam file; beam {} does not exist",
            path.string(), beam));
    }
    if (formats::psrfits::is_search_path(path)) {
        return BlockSource(formats::psrfits::SearchFile(path));
    }
    return BlockSource(formats::sigproc::FilterbankReader(path, bit_order));
}

/// Open one or more contiguous SIGPROC or search-mode PSRFITS files in time
/// order as a single type-erased BlockSource.
[[nodiscard]] inline BlockSource
open_block_source(std::span<const std::filesystem::path> paths,
                  BitOrder bit_order = BitOrder::kLsbFirst) {
    if (paths.empty()) {
        throw ValidationError(
            "psrio: open_block_source needs at least one file");
    }
    if (paths.size() == 1U) {
        return open_block_source(paths.front(), 0U, bit_order);
    }
    const bool search = formats::psrfits::is_search_path(paths.front());
    for (const auto& path : paths) {
        if (open_detail::looks_like_hdf5(path)) {
            throw ValidationError(
                "psrio: an HDF5 filterbank does not join with other files");
        }
        if (formats::psrfits::is_search_path(path) != search) {
            throw ValidationError(
                "psrio: every file in a stream must be the same format");
        }
    }
    if (search) {
        return BlockSource(formats::psrfits::SearchStream(paths));
    }
    return BlockSource(formats::sigproc::FilterbankStream(paths, bit_order));
}

} // namespace psrio

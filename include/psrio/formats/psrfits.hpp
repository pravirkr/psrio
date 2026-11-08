#pragma once

#include "psrio/block_source.hpp"
#include "psrio/formats/psrfits/reader.hpp"
#include "psrio/formats/psrfits/stream.hpp"

#include <cctype>
#include <filesystem>
#include <string>

namespace psrio::formats::psrfits {

/// True for a search-mode PSRFITS extension: `.sf`, `.fits`, or `.sfits`.
[[nodiscard]] inline bool is_search_path(const std::filesystem::path& path) {
    std::string ext = path.extension().string();
    for (char& character : ext) {
        character = static_cast<char>(
            std::tolower(static_cast<unsigned char>(character)));
    }
    return ext == ".sf" || ext == ".fits" || ext == ".sfits";
}

/// Open a search-mode PSRFITS file as a calibrated Stokes I block.
[[nodiscard]] inline BlockSource
open_block_source(const std::filesystem::path& path) {
    return BlockSource(SearchFile(path));
}

} // namespace psrio::formats::psrfits

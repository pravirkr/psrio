#pragma once

#include <ostream>
#include <span>
#include <string>

namespace psrio::cli {

/// Run the `psrio` command line with @p args (without the program name).
///
/// Normal output goes to @p out; progress, warnings, and errors to @p err.
/// @return The process exit code: 0 on success.
int run(std::span<const std::string> args,
        std::ostream& out,
        std::ostream& err);

} // namespace psrio::cli

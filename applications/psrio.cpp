#include "cli.hpp"

#include <cstddef>
#include <exception>
#include <iostream>
#include <span>
#include <string>
#include <vector>

int main(int argc, char** argv) {
    try {
        const std::span<char*> raw(argv, static_cast<std::size_t>(argc));
        const std::vector<std::string> args(raw.begin() + 1, raw.end());
        return psrio::cli::run(args, std::cout, std::cerr);
    } catch (const std::exception& error) {
        std::cerr << "psrio: fatal: " << error.what() << '\n';
    } catch (...) {
        std::cerr << "psrio: fatal: unknown error\n";
    }
    return 1;
}

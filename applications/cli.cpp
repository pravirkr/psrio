#include "cli.hpp"

#include "common.hpp"

#include <CLI/CLI.hpp>
#include <exception>
#include <ostream>
#include <span>
#include <string>
#include <string_view>
#include <vector>

#ifndef PSRIO_CLI_VERSION
#define PSRIO_CLI_VERSION "unknown"
#endif

namespace psrio::cli {

namespace {

std::string_view strip_prefix(std::string_view message) {
    constexpr std::string_view kPrefix = "psrio: ";
    if (message.starts_with(kPrefix)) {
        message.remove_prefix(kPrefix.size());
    }
    return message;
}

} // namespace

int run(std::span<const std::string> args,
        std::ostream& out,
        std::ostream& err) {
    Context ctx{.out = out, .err = err};

    CLI::App app{"Read, inspect, and transform pulsar and FRB time-domain data",
                 "psrio"};
    app.set_version_flag("-V,--version", PSRIO_CLI_VERSION);
    app.require_subcommand(1);
    app.fallthrough();
    app.get_formatter()->column_width(32);
    app.add_flag("-v,--verbose", ctx.verbosity, "Print more detail");
    app.add_flag("-q,--quiet", ctx.quiet, "Print only warnings and errors");
    app.footer("Run 'psrio <command> --help' for the options of a command.");

    add_inspect(app, ctx);
    add_slice(app, ctx);
    add_splice(app, ctx);
    add_decimate(app, ctx);
    add_convert(app, ctx);
    add_digifil(app, ctx);
    add_dedisperse(app, ctx);

    std::vector<std::string> reversed(args.rbegin(), args.rend());
    try {
        app.parse(reversed);
    } catch (const CLI::ParseError& e) {
        return app.exit(e, out, err);
    } catch (const std::exception& e) {
        const auto* sub = app.get_subcommands().empty()
                              ? nullptr
                              : app.get_subcommands().front();
        err << "psrio";
        if (sub != nullptr) {
            err << ' ' << sub->get_name();
        }
        err << ": error: " << strip_prefix(e.what()) << '\n';
        return 1;
    }
    return 0;
}

} // namespace psrio::cli

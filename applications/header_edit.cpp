#include "common.hpp"

#include <algorithm>
#include <charconv>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fmt/format.h>
#include <fstream>
#include <ios>
#include <optional>
#include <psrio/detail/exceptions.hpp>
#include <psrio/detail/rfile.hpp>
#include <psrio/formats/sigproc/header.hpp>
#include <span>
#include <string>
#include <string_view>
#include <system_error>
#include <type_traits>
#include <utility>
#include <vector>

namespace psrio::cli {

namespace {

[[noreturn]] void
bad_value(std::string_view key, std::string_view value, std::string_view why) {
    throw ValidationError(fmt::format("psrio: {}={}: {}", key, value, why));
}

double parse_double(std::string_view key, std::string_view value) {
    double out        = 0.0;
    const auto* first = value.data();
    const auto* last  = value.data() + value.size();
    const auto result = std::from_chars(first, last, out);
    if (result.ec != std::errc{} || result.ptr != last || !std::isfinite(out)) {
        bad_value(key, value, "not a finite number");
    }
    return out;
}

std::int64_t parse_int(std::string_view key, std::string_view value) {
    std::int64_t out  = 0;
    const auto* first = value.data();
    const auto* last  = value.data() + value.size();
    const auto result = std::from_chars(first, last, out);
    if (result.ec != std::errc{} || result.ptr != last) {
        bad_value(key, value, "not an integer");
    }
    return out;
}

template <typename T> std::string show(const std::optional<T>& value) {
    if (!value.has_value()) {
        return "";
    }
    if constexpr (std::is_same_v<T, std::string>) {
        return *value;
    } else {
        return fmt::format("{}", +*value);
    }
}

/// Set @p field from @p value; the key must already be in the file.
template <typename T>
void set_existing(std::optional<T>& field,
                  std::string_view key,
                  T value,
                  HeaderChange& change) {
    if (!field.has_value()) {
        throw ValidationError(fmt::format(
            "psrio: '{}' is not in this header; an in-place edit cannot add "
            "keys, because that changes the header size",
            key));
    }
    change.before = show(field);
    field         = std::move(value);
    change.after  = show(field);
}

void apply(sigproc::FilterbankHeader& header,
           std::string_view key,
           std::string_view value,
           HeaderChange& change) {
    const auto non_negative = [&](double number) {
        if (number < 0.0) {
            bad_value(key, value, "must not be negative");
        }
    };
    const auto int32 = [&](bool flag, bool nonneg) {
        const auto number = parse_int(key, value);
        if (flag && number != 0 && number != 1) {
            bad_value(key, value, "must be 0 or 1");
        }
        if (nonneg && number < 0) {
            bad_value(key, value, "must not be negative");
        }
        if (number < INT32_MIN || number > INT32_MAX) {
            bad_value(key, value, "does not fit a 32-bit integer");
        }
        return static_cast<std::int32_t>(number);
    };
    if (key == "nbits" || key == "nchans" || key == "nifs") {
        throw ValidationError(fmt::format(
            "psrio: '{}' cannot be changed in place: it sets how the samples "
            "are laid out (use convert or digifil)",
            key));
    }
    if (key == "source_name" || key == "rawdatafile") {
        auto& field =
            key == "source_name" ? header.source_name : header.rawdatafile;
        const auto old_size = field.has_value() ? field->size() : 0U;
        if (value.size() > old_size) {
            bad_value(key, value,
                      fmt::format("longer than the stored {} characters, "
                                  "which would change the header size",
                                  old_size));
        }
        std::string padded(value);
        padded.resize(old_size, ' ');
        set_existing(field, key, std::move(padded), change);
        return;
    }
    if (key == "telescope_id" || key == "machine_id" || key == "ibeam" ||
        key == "nbeams" || key == "nbins" || key == "npuls" ||
        key == "data_type") {
        std::optional<std::int32_t>* field = &header.data_type;
        for (const auto& [name, target] : {
                 std::pair{"telescope_id", &header.telescope_id},
                 std::pair{"machine_id", &header.machine_id},
                 std::pair{"ibeam", &header.ibeam},
                 std::pair{"nbeams", &header.nbeams},
                 std::pair{"nbins", &header.nbins},
                 std::pair{"npuls", &header.npuls},
             }) {
            if (key == name) {
                field = target;
            }
        }
        set_existing(*field, key, int32(false, true), change);
        return;
    }
    if (key == "barycentric" || key == "pulsarcentric") {
        auto& field =
            key == "barycentric" ? header.barycentric : header.pulsarcentric;
        set_existing(field, key, int32(true, false), change);
        return;
    }
    if (key == "signed") {
        set_existing(header.signed_data, key,
                     static_cast<std::int8_t>(int32(true, false)), change);
        return;
    }
    if (key == "nsamples") {
        const auto number = parse_int(key, value);
        if (number <= 0 || std::cmp_greater(number, UINT32_MAX)) {
            bad_value(key, value, "must be a positive 32-bit count");
        }
        set_existing(header.declared_nsamples, key,
                     static_cast<std::uint32_t>(number), change);
        return;
    }
    const auto number = parse_double(key, value);
    if ((key == "tsamp" || key == "tstart" || key == "fch1") && number <= 0.0) {
        bad_value(key, value, "must be positive");
    }
    if (key == "refdm" || key == "period") {
        non_negative(number);
    }
    if (key == "az_start" && (number < 0.0 || number > 360.0)) {
        bad_value(key, value, "must be within [0, 360]");
    }
    if (key == "za_start" && (number < 0.0 || number > 180.0)) {
        bad_value(key, value, "must be within [0, 180]");
    }
    std::optional<double>* field = nullptr;
    for (const auto& [name, target] : {
             std::pair{"tsamp", &header.tsamp},
             std::pair{"tstart", &header.tstart},
             std::pair{"fch1", &header.fch1},
             std::pair{"foff", &header.foff},
             std::pair{"src_raj", &header.src_raj},
             std::pair{"src_dej", &header.src_dej},
             std::pair{"az_start", &header.az_start},
             std::pair{"za_start", &header.za_start},
             std::pair{"refdm", &header.refdm},
             std::pair{"period", &header.period},
         }) {
        if (key == name) {
            field = target;
        }
    }
    if (field == nullptr) {
        throw ValidationError(
            fmt::format("psrio: '{}' is not a SIGPROC header key", key));
    }
    set_existing(*field, key, number, change);
}

} // namespace

std::vector<HeaderChange>
edit_sigproc_header(const std::filesystem::path& path,
                    std::span<const std::string> assignments,
                    bool dry_run) {
    if (assignments.empty()) {
        throw ValidationError("psrio: no header keys to update");
    }
    auto header = [&] {
        const detail::ReadFile file(path, detail::Access::kNormal);
        return sigproc::FilterbankHeader::read(file);
    }();
    const auto original_bytes = header.header_bytes;
    std::vector<HeaderChange> changes;
    std::vector<std::string> seen;
    for (const auto& assignment : assignments) {
        const auto equals = assignment.find('=');
        if (equals == std::string::npos || equals == 0U) {
            throw ValidationError(
                fmt::format("psrio: '{}' is not key=value", assignment));
        }
        const std::string key = assignment.substr(0, equals);
        const std::string_view value =
            std::string_view(assignment).substr(equals + 1U);
        if (std::ranges::find(seen, key) != seen.end()) {
            throw ValidationError(
                fmt::format("psrio: '{}' is given more than once", key));
        }
        seen.push_back(key);
        HeaderChange change{.key = key, .before = {}, .after = {}};
        apply(header, key, value, change);
        changes.push_back(std::move(change));
    }
    const auto bytes = header.to_bytes();
    if (bytes.size() != original_bytes) {
        throw ValidationError(fmt::format(
            "psrio: the edited header has {} bytes, not {}; this file's "
            "header cannot be rewritten in place",
            bytes.size(), original_bytes));
    }
    // The new header must read back with the requested values.
    const auto check =
        sigproc::FilterbankHeader::parse_prefix(bytes, header.file_bytes);
    if (check.to_bytes() != bytes) {
        throw ValidationError(
            "psrio: the edited header did not read back as written; file not "
            "changed");
    }
    if (dry_run) {
        return changes;
    }
    std::fstream file(path, std::ios::in | std::ios::out | std::ios::binary);
    if (!file) {
        throw IoError(
            fmt::format("psrio: cannot open {} for writing", path.string()));
    }
    file.write(reinterpret_cast<const char*>(bytes.data()), // NOLINT
               static_cast<std::streamsize>(bytes.size()));
    file.flush();
    if (!file) {
        throw IoError(fmt::format("psrio: failed to write the header of {}",
                                  path.string()));
    }
    return changes;
}

} // namespace psrio::cli

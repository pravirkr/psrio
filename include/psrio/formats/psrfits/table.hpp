#pragma once

#include "psrio/detail/endian.hpp"
#include "psrio/detail/exceptions.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <cctype>
#include <charconv>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <fmt/format.h>
#include <limits>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <system_error>
#include <utility>
#include <vector>

namespace psrio::formats::psrfits::detail {

inline constexpr std::size_t kRecordBytes    = 2880;
inline constexpr std::size_t kCardBytes      = 80;
inline constexpr std::size_t kCardsPerRecord = kRecordBytes / kCardBytes;

[[nodiscard]] inline std::string_view trim(std::string_view text) noexcept {
    while (!text.empty() && (text.front() == ' ' || text.front() == '\0')) {
        text.remove_prefix(1);
    }
    while (!text.empty() && (text.back() == ' ' || text.back() == '\0')) {
        text.remove_suffix(1);
    }
    return text;
}

[[nodiscard]] inline bool equals_ignore_case(std::string_view left,
                                             std::string_view right) noexcept {
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

[[nodiscard]] inline std::uint64_t
checked_mul(std::uint64_t left, std::uint64_t right, std::string_view what) {
    if (left != 0U &&
        right > std::numeric_limits<std::uint64_t>::max() / left) {
        throw FormatError(
            fmt::format("psrio: FITS {} overflows a 64-bit size", what));
    }
    return left * right;
}

[[nodiscard]] inline std::uint64_t pad_record(std::uint64_t bytes) {
    const auto rem = bytes % kRecordBytes;
    if (rem == 0U) {
        return bytes;
    }
    const auto gap = kRecordBytes - rem;
    if (bytes > std::numeric_limits<std::uint64_t>::max() - gap) {
        throw FormatError("psrio: FITS record padding overflows");
    }
    return bytes + gap;
}

/// Keyword/value cards from one header unit. Values are unquoted and trimmed.
class Cards {
public:
    void add(std::string key, std::string value) {
        m_entries.emplace_back(std::move(key), std::move(value));
    }

    [[nodiscard]] bool contains(std::string_view key) const {
        return find_entry(key) != nullptr;
    }

    [[nodiscard]] const std::string& value(std::string_view key) const {
        const auto* entry = find_entry(key);
        if (entry == nullptr) {
            throw FormatError(
                fmt::format("psrio: FITS header missing keyword '{}'", key));
        }
        return entry->second;
    }

    [[nodiscard]] std::optional<std::string> find(std::string_view key) const {
        const auto* entry = find_entry(key);
        if (entry == nullptr) {
            return std::nullopt;
        }
        return entry->second;
    }

    [[nodiscard]] double number(std::string_view key) const {
        return parse_number(value(key), key);
    }

    [[nodiscard]] std::optional<double>
    find_number(std::string_view key) const {
        const auto text = find(key);
        if (!text.has_value()) {
            return std::nullopt;
        }
        return parse_number(*text, key);
    }

    [[nodiscard]] std::int64_t integer(std::string_view key) const {
        return parse_integer(number(key), key);
    }

    [[nodiscard]] std::optional<std::int64_t>
    find_integer(std::string_view key) const {
        const auto parsed = find_number(key);
        if (!parsed.has_value()) {
            return std::nullopt;
        }
        return parse_integer(*parsed, key);
    }

    [[nodiscard]] bool flag(std::string_view key) const {
        return parse_flag(value(key), key);
    }

private:
    using Entry = std::pair<std::string, std::string>;

    [[nodiscard]] const Entry* find_entry(std::string_view key) const {
        for (const auto& entry : m_entries) {
            if (entry.first == key) {
                return &entry;
            }
        }
        return nullptr;
    }

    [[nodiscard]] static double parse_number(std::string_view text,
                                             std::string_view key) {
        text = trim(text);
        if (text.empty()) {
            throw FormatError(fmt::format(
                "psrio: FITS keyword '{}' has an empty numeric value", key));
        }
        if (!text.empty() && text.front() == '+') {
            text.remove_prefix(1);
        }
        std::array<char, 64> buffer{};
        if (text.empty() || text.size() >= buffer.size()) {
            throw FormatError(fmt::format(
                "psrio: FITS keyword '{}' is not a numeric value", key));
        }
        std::ranges::transform(text, buffer.begin(), [](char c) {
            return (c == 'D' || c == 'd') ? 'E' : c;
        });
        double parsed        = 0.0;
        const auto* end      = buffer.data() + text.size();
        const auto [ptr, ec] = std::from_chars(buffer.data(), end, parsed,
                                               std::chars_format::general);
        if (ec != std::errc{} || ptr != end || !std::isfinite(parsed)) {
            throw FormatError(fmt::format(
                "psrio: FITS keyword '{}' is not a numeric value", key));
        }
        return parsed;
    }

    [[nodiscard]] static std::int64_t parse_integer(double parsed,
                                                    std::string_view key) {
        if (std::trunc(parsed) != parsed ||
            parsed >
                static_cast<double>(std::numeric_limits<std::int64_t>::max()) ||
            parsed <
                static_cast<double>(std::numeric_limits<std::int64_t>::min())) {
            throw FormatError(
                fmt::format("psrio: FITS keyword '{}' is not an integer", key));
        }
        return static_cast<std::int64_t>(parsed);
    }

    [[nodiscard]] static bool parse_flag(std::string_view text,
                                         std::string_view key) {
        text = trim(text);
        if (text == "T") {
            return true;
        }
        if (text == "F") {
            return false;
        }
        throw FormatError(
            fmt::format("psrio: FITS keyword '{}' is not a T/F flag", key));
    }

    std::vector<Entry> m_entries;
};

struct Column {
    std::string name;
    /// Element type. For a heap column this is the array element, not `P`/`Q`.
    char code{'?'};
    /// Elements in a fixed-width cell. A heap column stores one descriptor.
    std::uint64_t repeat{0};
    std::uint64_t element_bytes{0};
    std::uint64_t offset{0};
    /// Bytes occupied in the row. A heap column occupies its descriptor.
    std::uint64_t nbytes{0};
    bool heap{false};
    /// 8 for a `P` descriptor, 16 for `Q`. Zero for a fixed-width column.
    std::uint64_t descriptor_bytes{0};
};

struct ParsedTable {
    Cards primary;
    Cards subint;
    std::uint64_t data_offset{0};
    std::uint64_t row_bytes{0};
    std::uint64_t nrows{0};
    /// Byte offset of the heap from the start of the table data (`THEAP`).
    std::uint64_t heap_offset{0};
    /// Heap size in bytes (`PCOUNT`).
    std::uint64_t heap_bytes{0};
    std::vector<Column> columns;
};

[[nodiscard]] inline std::string unquote(std::string_view field,
                                         std::string_view key) {
    field = trim(field);
    if (field.empty() || field.front() != '\'') {
        const auto slash = field.find('/');
        if (slash != std::string_view::npos) {
            field = trim(field.substr(0, slash));
        }
        return std::string{field};
    }
    std::string out;
    for (std::size_t index = 1; index < field.size(); ++index) {
        if (field[index] != '\'') {
            out.push_back(field[index]);
            continue;
        }
        if (index + 1 < field.size() && field[index + 1] == '\'') {
            out.push_back('\'');
            ++index;
            continue;
        }
        return std::string{trim(out)};
    }
    throw FormatError(fmt::format(
        "psrio: FITS keyword '{}' has an unterminated string", key));
}

[[nodiscard]] inline Cards parse_header(std::span<const std::byte> file,
                                        std::uint64_t start,
                                        std::uint64_t& data_offset) {
    Cards cards;
    std::uint64_t record = start;
    bool ended           = false;
    for (int block = 0; block < 64 && !ended; ++block) {
        if (record > file.size() || file.size() - record < kRecordBytes) {
            throw FormatError("psrio: FITS header is truncated");
        }
        for (std::size_t card = 0; card < kCardsPerRecord; ++card) {
            const auto* raw = file.data() + (static_cast<std::size_t>(record) +
                                             (card * kCardBytes));
            const std::string_view text(reinterpret_cast<const char*>(raw),
                                        kCardBytes);
            const auto key = std::string{trim(text.substr(0, 8))};
            if (key == "END") {
                ended = true;
                break;
            }
            if (key.empty() || text[8] != '=') {
                continue;
            }
            cards.add(key, unquote(text.substr(10), key));
        }
        record += kRecordBytes;
    }
    if (!ended) {
        throw FormatError("psrio: FITS header has no END card");
    }
    data_offset = record;
    return cards;
}

[[nodiscard]] inline std::uint64_t data_bytes(const Cards& header) {
    const auto naxis = header.integer("NAXIS");
    if (naxis < 0 || naxis > 999) {
        throw FormatError(
            fmt::format("psrio: FITS NAXIS {} is not usable", naxis));
    }
    if (naxis == 0) {
        return 0;
    }
    std::uint64_t product = 1;
    for (std::int64_t axis = 1; axis <= naxis; ++axis) {
        const auto length = header.integer(fmt::format("NAXIS{}", axis));
        if (length < 0) {
            throw FormatError(
                fmt::format("psrio: FITS NAXIS{} is negative", axis));
        }
        product = checked_mul(product, static_cast<std::uint64_t>(length),
                              "axis length");
    }
    const auto bitpix = header.integer("BITPIX");
    if (bitpix == 0 || bitpix % 8 != 0) {
        throw FormatError(fmt::format(
            "psrio: FITS BITPIX {} is not a whole number of bytes", bitpix));
    }
    const auto bytes_per = static_cast<std::uint64_t>(std::llabs(bitpix) / 8);
    const auto pcount =
        static_cast<std::uint64_t>(header.find_integer("PCOUNT").value_or(0));
    const auto gcount = header.find_integer("GCOUNT").value_or(1);
    if (gcount < 1) {
        throw FormatError("psrio: FITS GCOUNT must be positive");
    }
    if (pcount > std::numeric_limits<std::uint64_t>::max() - product) {
        throw FormatError("psrio: FITS data size overflows");
    }
    auto bytes = checked_mul(product + pcount,
                             static_cast<std::uint64_t>(gcount), "GCOUNT");
    bytes      = checked_mul(bytes, bytes_per, "BITPIX");
    return bytes;
}

struct TForm {
    std::uint64_t repeat{0};
    char code{'?'};
    std::uint64_t element_bytes{0};
    bool heap{false};
    std::uint64_t descriptor_bytes{0};
};

[[nodiscard]] inline std::uint64_t element_width(char code,
                                                 std::string_view text) {
    switch (code) {
    case 'L':
    case 'B':
    case 'A':
        return 1;
    case 'I':
        return 2;
    case 'J':
    case 'E':
        return 4;
    case 'K':
    case 'D':
        return 8;
    default:
        throw FormatError(fmt::format(
            "psrio: FITS TFORM '{}' uses an unsupported type", text));
    }
}

[[nodiscard]] inline TForm parse_tform(std::string_view text) {
    text = trim(text);
    if (text.empty()) {
        throw FormatError("psrio: FITS column has an empty TFORM");
    }
    std::size_t index    = 0;
    std::uint64_t repeat = 0;
    bool saw_digit       = false;
    while (index < text.size() && text[index] >= '0' && text[index] <= '9') {
        saw_digit        = true;
        const auto digit = static_cast<std::uint64_t>(text[index] - '0');
        repeat           = checked_mul(repeat, 10U, "TFORM repeat");
        if (repeat > std::numeric_limits<std::uint64_t>::max() - digit) {
            throw FormatError("psrio: FITS TFORM repeat overflows");
        }
        repeat += digit;
        ++index;
    }
    if (!saw_digit) {
        repeat = 1;
    }
    if (index >= text.size()) {
        throw FormatError(
            fmt::format("psrio: FITS TFORM '{}' has no type code", text));
    }
    const auto code = static_cast<char>(
        std::toupper(static_cast<unsigned char>(text[index])));
    const auto rest = trim(text.substr(index + 1));
    TForm parsed;
    parsed.repeat = repeat;
    if (code == 'P' || code == 'Q') {
        if (rest.size() != 1U) {
            throw FormatError(fmt::format(
                "psrio: FITS heap TFORM '{}' needs one element type", text));
        }
        const auto element = static_cast<char>(
            std::toupper(static_cast<unsigned char>(rest.front())));
        if (repeat != 1U) {
            throw FormatError(fmt::format(
                "psrio: FITS heap TFORM '{}' must describe one array", text));
        }
        parsed.heap             = true;
        parsed.code             = element;
        parsed.element_bytes    = element_width(element, text);
        parsed.descriptor_bytes = code == 'P' ? 8U : 16U;
        return parsed;
    }
    if (!rest.empty()) {
        throw FormatError(fmt::format(
            "psrio: FITS TFORM '{}' is not a fixed-width binary column", text));
    }
    parsed.code          = code;
    parsed.element_bytes = element_width(code, text);
    return parsed;
}

inline void reject_compressed(const Cards& header) {
    if (header.contains("ZCMPTYPE") || header.contains("ZBITPIX") ||
        header.contains("ZNAXIS")) {
        throw FormatError("psrio: compressed PSRFITS tables are not supported");
    }
    if (header.contains("ZIMAGE") && header.flag("ZIMAGE")) {
        throw FormatError("psrio: compressed PSRFITS tables are not supported");
    }
}

[[nodiscard]] inline std::vector<Column> columns_of(const Cards& header) {
    const auto nfields = header.integer("TFIELDS");
    if (nfields < 1 || nfields > 999) {
        throw FormatError(
            fmt::format("psrio: FITS TFIELDS {} is not usable", nfields));
    }
    std::vector<Column> columns;
    columns.reserve(static_cast<std::size_t>(nfields));
    std::uint64_t offset = 0;
    for (std::int64_t index = 1; index <= nfields; ++index) {
        Column column;
        column.name = header.value(fmt::format("TTYPE{}", index));
        const auto form =
            parse_tform(header.value(fmt::format("TFORM{}", index)));
        column.code             = form.code;
        column.repeat           = form.repeat;
        column.element_bytes    = form.element_bytes;
        column.heap             = form.heap;
        column.descriptor_bytes = form.descriptor_bytes;
        column.offset           = offset;
        const auto width =
            form.heap ? form.descriptor_bytes : form.element_bytes;
        column.nbytes = checked_mul(form.repeat, width, "TFORM");
        if (offset >
            std::numeric_limits<std::uint64_t>::max() - column.nbytes) {
            throw FormatError("psrio: FITS row width overflows");
        }
        offset += column.nbytes;
        for (const auto& earlier : columns) {
            if (earlier.name == column.name) {
                throw FormatError(fmt::format(
                    "psrio: FITS table repeats column '{}'", column.name));
            }
        }
        columns.push_back(std::move(column));
    }
    const auto naxis1 = header.integer("NAXIS1");
    if (naxis1 < 0 ||
        std::cmp_not_equal(static_cast<std::uint64_t>(naxis1), offset)) {
        throw FormatError(fmt::format(
            "psrio: SUBINT column widths sum to {} bytes, NAXIS1 is {}", offset,
            naxis1));
    }
    return columns;
}

[[nodiscard]] inline std::string extension_name(const Cards& header) {
    const auto name = header.find("EXTNAME");
    if (!name.has_value()) {
        return {};
    }
    return *name;
}

[[nodiscard]] inline ParsedTable parse_search(std::span<const std::byte> file,
                                              std::string_view path) {
    try {
        if (file.size() < kRecordBytes) {
            throw FormatError("psrio: file is too small to be FITS");
        }
        const std::string_view simple(
            reinterpret_cast<const char*>(file.data()), 9);
        if (simple != "SIMPLE  =") {
            throw FormatError("psrio: file does not start with a FITS header");
        }
        ParsedTable parsed;
        std::uint64_t cursor = 0;
        bool saw_subint      = false;
        for (int hdu = 0; hdu < 16 && !saw_subint; ++hdu) {
            std::uint64_t data_at = 0;
            Cards header          = parse_header(file, cursor, data_at);
            reject_compressed(header);
            const auto bytes = data_bytes(header);
            if (data_at > file.size() || bytes > file.size() - data_at) {
                throw FormatError(
                    "psrio: FITS data extends past the end of the file");
            }
            const auto padded = pad_record(bytes);
            // The final HDU may omit padding after its last byte. The bytes
            // themselves were checked above.
            if ((data_at > std::numeric_limits<std::uint64_t>::max() - padded ||
                 data_at + padded > file.size()) &&
                bytes != 0U && data_at + bytes > file.size()) {
                throw FormatError(
                    "psrio: FITS data extends past the end of the file");
            }
            const bool primary = hdu == 0;
            if (primary) {
                parsed.primary = std::move(header);
            } else if (equals_ignore_case(extension_name(header), "SUBINT")) {
                const auto xtension = header.find("XTENSION").value_or("");
                if (!equals_ignore_case(xtension, "BINTABLE")) {
                    throw FormatError(
                        "psrio: SUBINT extension is not a binary table");
                }
                if (header.integer("BITPIX") != 8) {
                    throw FormatError("psrio: SUBINT BITPIX must be 8");
                }
                if (header.integer("NAXIS") != 2) {
                    throw FormatError("psrio: SUBINT NAXIS must be 2");
                }
                const auto gcount = header.find_integer("GCOUNT").value_or(1);
                if (gcount != 1) {
                    throw FormatError(
                        "psrio: SUBINT GCOUNT other than 1 is not supported");
                }
                const auto nrows = header.integer("NAXIS2");
                if (nrows < 1) {
                    throw FormatError("psrio: SUBINT table has no rows");
                }
                parsed.subint      = std::move(header);
                parsed.data_offset = data_at;
                parsed.row_bytes =
                    static_cast<std::uint64_t>(parsed.subint.integer("NAXIS1"));
                parsed.nrows   = static_cast<std::uint64_t>(nrows);
                parsed.columns = columns_of(parsed.subint);
                const auto main =
                    checked_mul(parsed.row_bytes, parsed.nrows, "SUBINT table");
                const auto theap = parsed.subint.find_integer("THEAP");
                if (theap.has_value() && *theap < 0) {
                    throw FormatError("psrio: SUBINT THEAP is negative");
                }
                parsed.heap_offset = theap.has_value()
                                         ? static_cast<std::uint64_t>(*theap)
                                         : main;
                const auto pcount  = parsed.subint.find_integer("PCOUNT");
                parsed.heap_bytes  = pcount.has_value() && *pcount > 0
                                         ? static_cast<std::uint64_t>(*pcount)
                                         : 0U;
                saw_subint         = true;
            }
            if (saw_subint) {
                break;
            }
            if (padded > file.size() - data_at) {
                cursor = data_at + bytes;
            } else {
                cursor = data_at + padded;
            }
            if (cursor == file.size()) {
                break;
            }
            if (cursor > file.size()) {
                throw FormatError("psrio: FITS header walk ran past the file");
            }
        }
        if (!saw_subint) {
            throw FormatError("psrio: PSRFITS file has no SUBINT table");
        }
        const auto table_bytes =
            checked_mul(parsed.row_bytes, parsed.nrows, "SUBINT table");
        if (parsed.data_offset > file.size() ||
            table_bytes > file.size() - parsed.data_offset) {
            throw FormatError(
                "psrio: SUBINT table extends past the end of the file");
        }
        return parsed;
    } catch (const FormatError& error) {
        throw FormatError(fmt::format("{}: {}", path, error.what()));
    }
}

[[nodiscard]] inline const Column* find_column(const ParsedTable& table,
                                               std::string_view name) {
    for (const auto& column : table.columns) {
        if (column.name == name) {
            return &column;
        }
    }
    return nullptr;
}

[[nodiscard]] inline const Column& require_column(const ParsedTable& table,
                                                  std::string_view name) {
    const auto* column = find_column(table, name);
    if (column == nullptr) {
        throw FormatError(
            fmt::format("psrio: SUBINT table has no {} column", name));
    }
    return *column;
}

template <typename T>
[[nodiscard]] inline T load_be(const std::byte* source) noexcept {
    if constexpr (std::same_as<T, float>) {
        std::uint32_t bits = 0;
        std::memcpy(&bits, source, sizeof(bits));
        if (::psrio::detail::host_is_little_endian()) {
            bits = ::psrio::detail::swap_endian(bits);
        }
        return std::bit_cast<float>(bits);
    } else if constexpr (std::same_as<T, double>) {
        std::uint64_t bits = 0;
        std::memcpy(&bits, source, sizeof(bits));
        if (::psrio::detail::host_is_little_endian()) {
            bits = ::psrio::detail::swap_endian(bits);
        }
        return std::bit_cast<double>(bits);
    } else {
        T value{};
        std::memcpy(&value, source, sizeof(T));
        if (::psrio::detail::host_is_little_endian()) {
            value = ::psrio::detail::swap_endian(value);
        }
        return value;
    }
}

[[nodiscard]] inline double load_real(const std::byte* source, char code) {
    switch (code) {
    case 'B':
        return static_cast<double>(std::to_integer<unsigned char>(source[0]));
    case 'I':
        return static_cast<double>(load_be<std::int16_t>(source));
    case 'J':
        return static_cast<double>(load_be<std::int32_t>(source));
    case 'K':
        return static_cast<double>(load_be<std::int64_t>(source));
    case 'E':
        return static_cast<double>(load_be<float>(source));
    case 'D':
        return load_be<double>(source);
    default:
        throw FormatError(
            fmt::format("psrio: FITS column type '{}' is not numeric", code));
    }
}

inline void load_reals(std::span<const std::byte> cell,
                       const Column& column,
                       std::span<double> dest) {
    if (column.code != 'B' && column.code != 'I' && column.code != 'J' &&
        column.code != 'K' && column.code != 'E' && column.code != 'D') {
        throw FormatError(fmt::format(
            "psrio: SUBINT column '{}' is not numeric", column.name));
    }
    if (column.element_bytes == 0U ||
        cell.size() / static_cast<std::size_t>(column.element_bytes) <
            dest.size()) {
        throw FormatError(fmt::format(
            "psrio: SUBINT column '{}' has {} bytes, needed {}", column.name,
            cell.size(),
            dest.size() * static_cast<std::size_t>(column.element_bytes)));
    }
    for (std::size_t index = 0; index < dest.size(); ++index) {
        const auto* source =
            cell.data() +
            (index * static_cast<std::size_t>(column.element_bytes));
        const double value = load_real(source, column.code);
        if (!std::isfinite(value)) {
            throw FormatError(fmt::format(
                "psrio: SUBINT column '{}' contains a non-finite value",
                column.name));
        }
        dest[index] = value;
    }
}

[[nodiscard]] inline std::span<const std::byte>
column_cell(std::span<const std::byte> file,
            const ParsedTable& table,
            std::uint64_t row,
            const Column& column) {
    if (row >= table.nrows) {
        throw FormatError(fmt::format("psrio: SUBINT row {} is outside 0..{}",
                                      row, table.nrows - 1U));
    }
    const auto start = table.data_offset +
                       checked_mul(row, table.row_bytes, "SUBINT row") +
                       column.offset;
    if (start > file.size() || column.nbytes > file.size() - start) {
        throw FormatError(fmt::format(
            "psrio: SUBINT row {} column '{}' extends past the file", row,
            column.name));
    }
    if (!column.heap) {
        return file.subspan(static_cast<std::size_t>(start),
                            static_cast<std::size_t>(column.nbytes));
    }
    if (column.descriptor_bytes != 8U && column.descriptor_bytes != 16U) {
        throw FormatError(
            fmt::format("psrio: SUBINT column '{}' has a bad heap descriptor",
                        column.name));
    }
    if (column.nbytes > file.size() - start) {
        throw FormatError(fmt::format(
            "psrio: SUBINT row {} column '{}' extends past the file", row,
            column.name));
    }
    const auto* descriptor = file.data() + static_cast<std::size_t>(start);
    std::int64_t count     = 0;
    std::int64_t heap_at   = 0;
    if (column.descriptor_bytes == 8U) {
        count   = load_be<std::int32_t>(descriptor);
        heap_at = load_be<std::int32_t>(descriptor + 4);
    } else {
        count   = load_be<std::int64_t>(descriptor);
        heap_at = load_be<std::int64_t>(descriptor + 8);
    }
    if (count < 0 || heap_at < 0) {
        throw FormatError(fmt::format(
            "psrio: SUBINT row {} column '{}' has a negative heap descriptor",
            row, column.name));
    }
    const auto nbytes = checked_mul(static_cast<std::uint64_t>(count),
                                    column.element_bytes, "heap array");
    const auto main_bytes =
        checked_mul(table.row_bytes, table.nrows, "SUBINT table");
    if (table.data_offset > file.size() ||
        main_bytes > file.size() - table.data_offset ||
        table.heap_bytes > file.size() - table.data_offset - main_bytes) {
        throw FormatError(fmt::format(
            "psrio: SUBINT row {} column '{}' heap array extends past the file",
            row, column.name));
    }
    const auto hdu_bytes = main_bytes + table.heap_bytes;
    const auto heap_pos  = static_cast<std::uint64_t>(heap_at);
    if (table.heap_offset > hdu_bytes ||
        heap_pos > hdu_bytes - table.heap_offset ||
        nbytes > hdu_bytes - table.heap_offset - heap_pos) {
        throw FormatError(fmt::format(
            "psrio: SUBINT row {} column '{}' heap array extends past the file",
            row, column.name));
    }
    const auto origin = table.data_offset + table.heap_offset + heap_pos;
    return file.subspan(static_cast<std::size_t>(origin),
                        static_cast<std::size_t>(nbytes));
}

[[nodiscard]] inline std::vector<int> parse_tdim(std::string_view text) {
    text = trim(text);
    if (!text.empty() && text.front() == '(' && text.back() == ')') {
        text = text.substr(1, text.size() - 2);
    }
    text = trim(text);
    if (text.empty()) {
        throw FormatError("psrio: FITS TDIM is empty");
    }
    std::vector<int> axes;
    std::size_t cursor = 0;
    while (cursor < text.size()) {
        const auto comma = text.find(',', cursor);
        const auto token = trim(text.substr(
            cursor, comma == std::string_view::npos ? text.size() - cursor
                                                    : comma - cursor));
        int length       = 0;
        const auto [ptr, ec] =
            std::from_chars(token.data(), token.data() + token.size(), length);
        if (token.empty() || ec != std::errc{} ||
            ptr != token.data() + token.size() || length <= 0) {
            throw FormatError(
                fmt::format("psrio: FITS TDIM '{}' has a bad axis", text));
        }
        axes.push_back(length);
        if (comma == std::string_view::npos) {
            break;
        }
        cursor = comma + 1;
    }
    return axes;
}

} // namespace psrio::formats::psrfits::detail

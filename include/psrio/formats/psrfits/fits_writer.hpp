#pragma once

#include "psrio/detail/endian.hpp"
#include "psrio/detail/exceptions.hpp"
#include "psrio/formats/psrfits/table.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fmt/format.h>
#include <fstream>
#include <ios>
#include <span>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace psrio::formats::psrfits {

/**
 * @brief One FITS header card: keyword, formatted value, and comment.
 *
 * Values are formatted the way CFITSIO writes them, so files read back the
 * same in CFITSIO, astropy, and psrchive: strings quoted and padded to at
 * least eight characters, integers right-justified to column 30, and
 * doubles as `%.15G` with a decimal point always present.
 */
struct Card {
    std::string key;
    /// Value field exactly as written (quoted for strings). Empty for a
    /// commentary card.
    std::string value;
    std::string comment;
    /// True for a COMMENT/HISTORY card with no value indicator.
    bool commentary{false};

    [[nodiscard]] static Card
    string(std::string key, std::string_view text, std::string comment = {}) {
        std::string quoted{"'"};
        for (const char character : text) {
            quoted.push_back(character);
            if (character == '\'') {
                quoted.push_back('\'');
            }
        }
        while (quoted.size() < 9U) {
            quoted.push_back(' ');
        }
        quoted.push_back('\'');
        return Card{
            .key        = std::move(key),
            .value      = std::move(quoted),
            .comment    = std::move(comment),
            .commentary = false,
        };
    }

    [[nodiscard]] static Card
    integer(std::string key, std::int64_t number, std::string comment = {}) {
        return Card{
            .key        = std::move(key),
            .value      = fmt::format("{:>20}", number),
            .comment    = std::move(comment),
            .commentary = false,
        };
    }

    [[nodiscard]] static Card
    real(std::string key, double number, std::string comment = {}) {
        return Card{
            .key        = std::move(key),
            .value      = fmt::format("{:>20}", format_double(number)),
            .comment    = std::move(comment),
            .commentary = false,
        };
    }

    [[nodiscard]] static Card
    logical(std::string key, bool flag, std::string comment = {}) {
        return Card{
            .key        = std::move(key),
            .value      = fmt::format("{:>20}", flag ? "T" : "F"),
            .comment    = std::move(comment),
            .commentary = false,
        };
    }

    [[nodiscard]] static Card note(std::string key, std::string text) {
        return Card{
            .key        = std::move(key),
            .value      = {},
            .comment    = std::move(text),
            .commentary = true,
        };
    }

    /// CFITSIO `ffd2e(value, -15)`: `%.15G`, with a decimal point added when
    /// the result has none.
    [[nodiscard]] static std::string format_double(double number) {
        std::string text = fmt::format("{:.15G}", number);
        if (text.find('.') == std::string::npos &&
            text.find('N') == std::string::npos &&
            text.find('I') == std::string::npos) {
            const auto exponent = text.find('E');
            if (exponent == std::string::npos) {
                text.push_back('.');
            } else {
                text.insert(exponent, ".");
            }
        }
        return text;
    }

    /// The 80-character card image.
    [[nodiscard]] std::string image() const {
        if (key.size() > 8U) {
            throw ValidationError(
                fmt::format("psrio: FITS keyword '{}' is longer than 8", key));
        }
        std::string out(80, ' ');
        out.replace(0, key.size(), key);
        if (commentary) {
            const auto text = std::string_view(comment).substr(
                0, std::min<std::size_t>(comment.size(), 72U));
            out.replace(8, text.size(), text);
            return out;
        }
        if (key == "END") {
            return out;
        }
        out[8] = '=';
        if (value.size() > 70U) {
            throw ValidationError(
                fmt::format("psrio: FITS value of '{}' does not fit", key));
        }
        out.replace(10, value.size(), value);
        if (!comment.empty()) {
            const std::size_t at =
                std::max<std::size_t>(10U + value.size(), 30U);
            if (at + 3U < 80U) {
                out.replace(at, 3, " / ");
                const auto room = 80U - (at + 3U);
                const auto text = std::string_view(comment).substr(
                    0, std::min(room, comment.size()));
                out.replace(at + 3U, text.size(), text);
            }
        }
        return out;
    }
};

/// One binary-table column definition.
struct TableColumn {
    std::string name;
    std::string tform;
    std::string tunit;
    std::string tdim;
    std::string comment;
};

/// A column definition (designated so the field order cannot slip).
[[nodiscard]] inline TableColumn column(std::string name,
                                        std::string tform,
                                        std::string tunit,
                                        std::string tdim,
                                        std::string comment) {
    return TableColumn{
        .name    = std::move(name),
        .tform   = std::move(tform),
        .tunit   = std::move(tunit),
        .tdim    = std::move(tdim),
        .comment = std::move(comment),
    };
}

/// Bytes one row of @p columns occupies.
[[nodiscard]] inline std::uint64_t
row_width(std::span<const TableColumn> columns) {
    std::uint64_t width = 0;
    for (const auto& column : columns) {
        const auto form = detail::parse_tform(column.tform);
        if (form.heap) {
            throw ValidationError(
                "psrio: the FITS writer does not write heap columns");
        }
        width += detail::checked_mul(form.repeat, form.element_bytes, "TFORM");
    }
    return width;
}

/**
 * @brief Header cards of a binary-table extension.
 *
 * Writes the mandatory keywords, the column definitions in order, the
 * @p extra cards, and `EXTNAME`.
 */
[[nodiscard]] inline std::vector<Card>
bintable_header(std::string_view extname,
                std::span<const TableColumn> columns,
                std::uint64_t nrows,
                std::span<const Card> extra,
                std::string_view title = {}) {
    std::vector<Card> cards;
    cards.push_back(Card::string("XTENSION", "BINTABLE", std::string{title}));
    cards.push_back(Card::integer("BITPIX", 8, "Bytes"));
    cards.push_back(Card::integer("NAXIS", 2, "Rows and columns"));
    cards.push_back(Card::integer("NAXIS1",
                                  static_cast<std::int64_t>(row_width(columns)),
                                  "Bytes per row"));
    cards.push_back(
        Card::integer("NAXIS2", static_cast<std::int64_t>(nrows), "Rows"));
    cards.push_back(Card::integer("PCOUNT", 0, "No heap"));
    cards.push_back(Card::integer("GCOUNT", 1, "Always 1"));
    cards.push_back(Card::integer(
        "TFIELDS", static_cast<std::int64_t>(columns.size()), "Columns"));
    for (std::size_t index = 0; index < columns.size(); ++index) {
        const auto& column = columns[index];
        const auto number  = index + 1U;
        cards.push_back(Card::string(fmt::format("TTYPE{}", number),
                                     column.name, column.comment));
        cards.push_back(
            Card::string(fmt::format("TFORM{}", number), column.tform));
        if (!column.tunit.empty()) {
            cards.push_back(Card::string(fmt::format("TUNIT{}", number),
                                         column.tunit, "Units"));
        }
        if (!column.tdim.empty()) {
            cards.push_back(
                Card::string(fmt::format("TDIM{}", number), column.tdim));
        }
    }
    cards.insert(cards.end(), extra.begin(), extra.end());
    cards.push_back(Card::string("EXTNAME", extname, "Table name"));
    cards.push_back(Card::integer("EXTVER", 1, "Table version"));
    return cards;
}

/// Store @p value big-endian at @p out.
template <typename T>
    requires std::integral<T> || std::same_as<T, float> ||
             std::same_as<T, double>
inline void store_be(std::byte* out, T value) noexcept {
    if constexpr (std::same_as<T, float>) {
        store_be(out, std::bit_cast<std::uint32_t>(value));
    } else if constexpr (std::same_as<T, double>) {
        store_be(out, std::bit_cast<std::uint64_t>(value));
    } else {
        if constexpr (sizeof(T) > 1U) {
            if (::psrio::detail::host_is_little_endian()) {
                value = ::psrio::detail::swap_endian(value);
            }
        }
        std::memcpy(out, &value, sizeof(T));
    }
}

/// Store @p values big-endian into @p out, which must hold them all.
template <typename T>
inline void store_be(std::span<std::byte> out, std::span<const T> values) {
    if (out.size() < values.size() * sizeof(T)) {
        throw ValidationError("psrio: FITS cell buffer is too small");
    }
    for (std::size_t index = 0; index < values.size(); ++index) {
        store_be(out.data() + (index * sizeof(T)), values[index]);
    }
}

/**
 * @brief Sequential FITS file writer.
 *
 * Write a header, then that HDU's data, then the next header. Headers and
 * data are padded to 2880-byte records (spaces and zeros respectively).
 */
class FitsWriter {
public:
    /// @throws IoError if @p path cannot be created.
    explicit FitsWriter(const std::filesystem::path& path)
        : m_path(path),
          m_out(path, std::ios::binary | std::ios::trunc) {
        if (!m_out) {
            throw IoError(
                fmt::format("psrio: cannot create '{}'", path.string()));
        }
    }

    FitsWriter(const FitsWriter&)                = delete;
    FitsWriter& operator=(const FitsWriter&)     = delete;
    FitsWriter(FitsWriter&&) noexcept            = default;
    FitsWriter& operator=(FitsWriter&&) noexcept = default;
    ~FitsWriter()                                = default;

    /// Write @p cards, the END card, and padding.
    void write_header(std::span<const Card> cards) {
        finish_data();
        std::string block;
        block.reserve((cards.size() + 1U) * 80U);
        for (const auto& card : cards) {
            block += card.image();
        }
        block +=
            Card{
                .key        = "END",
                .value      = {},
                .comment    = {},
                .commentary = false,
            }
                .image();
        const auto rem = block.size() % detail::kRecordBytes;
        if (rem != 0U) {
            block.append(detail::kRecordBytes - rem, ' ');
        }
        put(block.data(), block.size());
        m_in_data = true;
        m_data    = 0;
    }

    /// Append bytes to the current HDU's data.
    void write_data(std::span<const std::byte> bytes) {
        if (!m_in_data) {
            throw ValidationError("psrio: FITS data written before a header");
        }
        // NOLINTNEXTLINE(cppcoreguidelines-pro-type-reinterpret-cast)
        put(reinterpret_cast<const char*>(bytes.data()), bytes.size());
        m_data += bytes.size();
    }

    /// Pad the last data unit and flush.
    /// @throws IoError if writing fails.
    void close() {
        finish_data();
        m_out.flush();
        if (!m_out) {
            throw IoError(
                fmt::format("psrio: writing '{}' failed", m_path.string()));
        }
        m_out.close();
    }

private:
    void finish_data() {
        if (!m_in_data) {
            return;
        }
        const auto rem = m_data % detail::kRecordBytes;
        if (rem != 0U) {
            const std::string zeros(detail::kRecordBytes - rem, '\0');
            put(zeros.data(), zeros.size());
        }
        m_in_data = false;
        m_data    = 0;
    }

    void put(const char* bytes, std::size_t count) {
        m_out.write(bytes, static_cast<std::streamsize>(count));
        if (!m_out) {
            throw IoError(
                fmt::format("psrio: writing '{}' failed", m_path.string()));
        }
    }

    std::filesystem::path m_path;
    std::ofstream m_out;
    std::uint64_t m_data{0};
    bool m_in_data{false};
};

} // namespace psrio::formats::psrfits

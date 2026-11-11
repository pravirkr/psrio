#pragma once

#include "psrio/detail/exceptions.hpp"
#include "psrio/detail/write.hpp"
#include "psrio/header.hpp"

#include <cstdint>
#include <filesystem>
#include <fmt/format.h>
#include <fstream>
#include <ios>
#include <span>
#include <utility>

namespace psrio::formats::presto {

/// Streaming writer for a PRESTO `.dat` time series and its `.inf` file.
///
/// The constructor truncates `<base>.dat`. Each write appends little-endian
/// float32 samples. close() writes `<base>.inf` from the header, whose
/// nsamples is the number of samples written, so the `.inf` file exists only
/// once the series is complete.
class PrestoTimeSeriesWriter {
public:
    /// @param base_or_dat_path Output base name. A `.dat` or `.inf` suffix is
    ///        removed.
    /// @throws ValidationError if tsamp is not positive.
    /// @throws IoError if the `.dat` file cannot be created.
    explicit PrestoTimeSeriesWriter(
        const std::filesystem::path& base_or_dat_path, Header header);

    PrestoTimeSeriesWriter(const PrestoTimeSeriesWriter&)            = delete;
    PrestoTimeSeriesWriter& operator=(const PrestoTimeSeriesWriter&) = delete;
    PrestoTimeSeriesWriter(PrestoTimeSeriesWriter&&) noexcept        = default;
    PrestoTimeSeriesWriter&
    operator=(PrestoTimeSeriesWriter&&) noexcept = default;

    /// Closes the writer. A failure is ignored; call close() to observe it.
    ~PrestoTimeSeriesWriter() noexcept;

    /// Append @p samples to the `.dat` file.
    /// @throws IoError if the writer is closed or the write fails.
    void write_samples(std::span<const float> samples);

    /// Alias for write_samples().
    void write(std::span<const float> samples) { write_samples(samples); }

    [[nodiscard]] std::uint64_t nsamples_written() const noexcept {
        return m_samples_written;
    }

    /// Header for the `.inf` file. nsamples counts the samples written.
    [[nodiscard]] const Header& header() const noexcept { return m_header; }

    [[nodiscard]] const std::filesystem::path& dat_path() const noexcept {
        return m_dat_path;
    }

    [[nodiscard]] const std::filesystem::path& inf_path() const noexcept {
        return m_inf_path;
    }

    [[nodiscard]] bool is_open() const noexcept { return m_dat.is_open(); }

    /// Close the `.dat` file and write the `.inf` file. Does nothing when
    /// already closed.
    /// @throws IoError if either file cannot be written.
    void close();

private:
    std::ofstream m_dat;
    Header m_header;
    std::filesystem::path m_dat_path;
    std::filesystem::path m_inf_path;
    std::uint64_t m_samples_written{0};
};

inline PrestoTimeSeriesWriter::PrestoTimeSeriesWriter(
    const std::filesystem::path& base_or_dat_path, Header header)
    : m_header(std::move(header)) {
    if (!(m_header.tsamp > 0.0)) {
        throw ValidationError("psrio: time series tsamp must be positive");
    }
    std::filesystem::path base = base_or_dat_path;
    const auto ext             = base.extension();
    if (ext == ".dat" || ext == ".inf") {
        base.replace_extension();
    }
    m_dat_path        = std::filesystem::path(base).concat(".dat");
    m_inf_path        = std::filesystem::path(base).concat(".inf");
    m_header.filename = m_dat_path.string();
    m_header.nbits    = 32;
    m_header.nsamples = 0;

    m_dat.open(m_dat_path, std::ios::binary | std::ios::trunc);
    if (!m_dat) {
        throw IoError(fmt::format("psrio: cannot create .dat file '{}'",
                                  m_dat_path.string()));
    }
}

inline PrestoTimeSeriesWriter::~PrestoTimeSeriesWriter() noexcept {
    try {
        close();
    } catch (...) {
        // Destructors must not throw; close() reports the failure.
    }
}

inline void
PrestoTimeSeriesWriter::write_samples(std::span<const float> samples) {
    if (!m_dat.is_open()) {
        throw IoError(fmt::format("psrio: .dat writer for '{}' is closed",
                                  m_dat_path.string()));
    }
    if (!::psrio::detail::write_float32_le(m_dat, samples)) {
        throw IoError(fmt::format("psrio: failed to write .dat file '{}'",
                                  m_dat_path.string()));
    }
    m_samples_written += samples.size();
    m_header.nsamples = m_samples_written;
}

inline void PrestoTimeSeriesWriter::close() {
    if (!m_dat.is_open()) {
        return;
    }
    m_dat.close();
    if (!m_dat) {
        throw IoError(fmt::format("psrio: failed to close .dat file '{}'",
                                  m_dat_path.string()));
    }
    m_header.make_inf(m_inf_path);
}

} // namespace psrio::formats::presto

#pragma once

#include "common.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <fmt/format.h>
#include <iostream>
#include <string>
#include <string_view>
#include <unistd.h>

namespace psrio::cli {

enum class ProgressMode : std::uint8_t {
    kAuto,
    kBar,
    kLog,
    kNone,
};

[[nodiscard]] inline ProgressMode parse_progress_mode(std::string_view text) {
    if (text == "auto") {
        return ProgressMode::kAuto;
    }
    if (text == "bar") {
        return ProgressMode::kBar;
    }
    if (text == "log") {
        return ProgressMode::kLog;
    }
    if (text == "none") {
        return ProgressMode::kNone;
    }
    throw ValidationError(fmt::format(
        "psrio: --progress must be auto, bar, log, or none, got '{}'", text));
}

[[nodiscard]] inline std::string format_bytes(std::uint64_t bytes) {
    constexpr double kKiB = 1024.0;
    constexpr double kMiB = kKiB * 1024.0;
    constexpr double kGiB = kMiB * 1024.0;
    constexpr double kTiB = kGiB * 1024.0;
    const auto value      = static_cast<double>(bytes);
    if (value >= kTiB) {
        return fmt::format("{:.2f} TiB", value / kTiB);
    }
    if (value >= kGiB) {
        return fmt::format("{:.2f} GiB", value / kGiB);
    }
    if (value >= kMiB) {
        return fmt::format("{:.1f} MiB", value / kMiB);
    }
    if (value >= kKiB) {
        return fmt::format("{:.1f} KiB", value / kKiB);
    }
    return fmt::format("{} B", bytes);
}

[[nodiscard]] inline std::string format_duration(double seconds) {
    if (!std::isfinite(seconds) || seconds < 0.0) {
        return "--:--";
    }
    const auto total = static_cast<std::uint64_t>(std::llround(seconds));
    const auto hrs   = total / 3600U;
    const auto mins  = (total % 3600U) / 60U;
    const auto secs  = total % 60U;
    if (hrs > 0U) {
        return fmt::format("{}:{:02d}:{:02d}", hrs, mins, secs);
    }
    return fmt::format("{:02d}:{:02d}", mins, secs);
}

[[nodiscard]] inline std::string format_rate(double bytes_per_sec) {
    if (!std::isfinite(bytes_per_sec) || bytes_per_sec <= 0.0) {
        return "-- MB/s";
    }
    constexpr double kGB = 1.0e9;
    constexpr double kMB = 1.0e6;
    if (bytes_per_sec >= kGB) {
        return fmt::format("{:.2f} GB/s", bytes_per_sec / kGB);
    }
    return fmt::format("{:.1f} MB/s", bytes_per_sec / kMB);
}

/// Adaptive progress reporter for long-running CLI jobs.
///
/// Renders an in-place `\r` progress bar when stderr is an interactive TTY,
/// and clean newline-delimited status lines when redirected to a batch log
/// (Slurm / HTCondor). Respects `--quiet` and `--progress`.
class ProgressReporter {
public:
    using Clock = std::chrono::steady_clock;

    ProgressReporter(const Context& ctx,
                     ProgressMode mode,
                     std::string_view label)
        : m_ctx(&ctx),
          m_label(label),
          m_start(Clock::now()),
          m_last_emit(m_start) {
        if (ctx.quiet || mode == ProgressMode::kNone) {
            m_mode = ProgressMode::kNone;
        } else if (mode == ProgressMode::kAuto) {
            const bool tty =
                (&ctx.err == &std::cerr) && (::isatty(STDERR_FILENO) != 0);
            m_mode = tty ? ProgressMode::kBar : ProgressMode::kLog;
        } else {
            m_mode = mode;
        }
    }

    ~ProgressReporter() noexcept {
        try {
            clear_line();
        } catch (...) {
            // Destructors must not throw.
        }
    }

    ProgressReporter(const ProgressReporter&)            = delete;
    ProgressReporter& operator=(const ProgressReporter&) = delete;
    ProgressReporter(ProgressReporter&&)                 = delete;
    ProgressReporter& operator=(ProgressReporter&&)      = delete;

    [[nodiscard]] ProgressMode mode() const noexcept { return m_mode; }

    /// Report header-opening progress when scanning many input files.
    void update_scan(std::size_t done, std::size_t total) {
        if (m_mode == ProgressMode::kNone || total <= 8U) {
            return;
        }
        const auto now = Clock::now();
        const double since_last =
            std::chrono::duration<double>(now - m_last_emit).count();
        const double elapsed =
            std::chrono::duration<double>(now - m_start).count();
        if (m_mode == ProgressMode::kBar) {
            if (m_scan_emitted && done < total && since_last < 0.2) {
                return;
            }
            m_scan_emitted = true;
            m_last_emit    = now;
            const double pct =
                100.0 * static_cast<double>(done) / static_cast<double>(total);
            emit_bar_line(fmt::format("{}: opening headers {}/{} ({:5.1f}%)",
                                      m_label, done, total, pct));
            return;
        }
        // Log mode: emit once opening exceeds 2 s, then every 5 s or 10%,
        // plus a final 100% line if any scan line was printed.
        const auto bucket    = (done * 10U) / total;
        const bool first_due = !m_scan_emitted && elapsed >= 2.0;
        const bool next_due =
            m_scan_emitted &&
            (done >= total || since_last >= 5.0 || bucket > m_last_scan_bucket);
        if (!first_due && !next_due) {
            return;
        }
        m_scan_emitted     = true;
        m_last_emit        = now;
        m_last_scan_bucket = std::max(m_last_scan_bucket, bucket);
        const double pct =
            100.0 * static_cast<double>(done) / static_cast<double>(total);
        m_ctx->err << fmt::format("{}: opened {}/{} headers ({:.0f}%, {})\n",
                                  m_label, done, total, pct,
                                  format_duration(elapsed));
        m_ctx->err.flush();
    }

    void finish_scan() {
        clear_line();
        m_start     = Clock::now();
        m_last_emit = m_start;
    }

    /// Configure the processing phase before streaming blocks.
    void begin_stream(std::uint64_t total_samples,
                      double tsamp,
                      std::uint64_t total_input_bytes,
                      std::size_t nfiles) {
        clear_line();
        m_total_samples     = total_samples;
        m_tsamp             = tsamp;
        m_total_input_bytes = total_input_bytes;
        m_nfiles            = nfiles;
        m_start             = Clock::now();
        m_last_emit         = m_start;
        m_last_log_bucket   = 0;
        m_last_file_logged  = 0;
        m_stream_started    = true;
        m_stream_emitted    = false;
    }

    /// Report streaming progress after one block when total samples are known.
    void update_stream(std::uint64_t samples_done, std::size_t file_index) {
        if (m_mode == ProgressMode::kNone || m_total_samples == 0U) {
            return;
        }
        const double frac = std::clamp(static_cast<double>(samples_done) /
                                           static_cast<double>(m_total_samples),
                                       0.0, 1.0);
        const bool finished   = samples_done >= m_total_samples;
        const auto bytes_done = static_cast<std::uint64_t>(
            frac * static_cast<double>(m_total_input_bytes));
        emit_stream_progress(samples_done, m_total_samples, file_index, frac,
                             bytes_done, finished);
    }

    /// Report streaming progress for lazy `--cont` streams where total samples
    /// are not pre-scanned up front.
    void update_stream_lazy(std::uint64_t samples_done,
                            std::size_t file_index,
                            double file_frac,
                            std::uint64_t bytes_done,
                            bool finished) {
        if (m_mode == ProgressMode::kNone || m_nfiles == 0U) {
            return;
        }
        const double clamped_file = std::clamp(file_frac, 0.0, 1.0);
        const double frac =
            finished
                ? 1.0
                : std::clamp((static_cast<double>(file_index) + clamped_file) /
                                 static_cast<double>(m_nfiles),
                             0.0, 1.0);
        emit_stream_progress(samples_done, 0U, file_index, frac, bytes_done,
                             finished);
    }

private:
    void emit_stream_progress(std::uint64_t samples_done,
                              std::uint64_t total_samples,
                              std::size_t file_index,
                              double frac,
                              std::uint64_t bytes_done,
                              bool finished) {
        const auto now = Clock::now();
        const double since_last =
            std::chrono::duration<double>(now - m_last_emit).count();
        const double elapsed =
            std::chrono::duration<double>(now - m_start).count();

        if (m_mode == ProgressMode::kBar) {
            if (m_stream_emitted && !finished && since_last < 0.25) {
                return;
            }
            m_stream_emitted  = true;
            m_last_emit       = now;
            const double rate = elapsed > 1.0e-3
                                    ? static_cast<double>(bytes_done) / elapsed
                                    : 0.0;
            const double eta  = (frac > 1.0e-6 && elapsed > 0.25)
                                    ? elapsed * ((1.0 - frac) / frac)
                                    : -1.0;
            constexpr std::size_t kWidth = 24;
            const auto filled            = static_cast<std::size_t>(
                std::llround(frac * static_cast<double>(kWidth)));
            std::string bar(kWidth, ' ');
            for (std::size_t i = 0; i < std::min(filled, kWidth); ++i) {
                bar[i] = '=';
            }
            if (filled < kWidth) {
                bar[filled] = '>';
            }
            const double obs_done = static_cast<double>(samples_done) * m_tsamp;
            std::string line;
            if (total_samples > 0U) {
                const double obs_tot =
                    static_cast<double>(total_samples) * m_tsamp;
                if (m_nfiles > 1U) {
                    line = fmt::format(
                        "[{}] {:5.1f}% | file {}/{} | {}/{} | {} | ETA {}", bar,
                        frac * 100.0, file_index + 1U, m_nfiles,
                        format_duration(obs_done), format_duration(obs_tot),
                        format_rate(rate), format_duration(eta));
                } else {
                    line =
                        fmt::format("[{}] {:5.1f}% | {}/{} | {} | ETA {}", bar,
                                    frac * 100.0, format_duration(obs_done),
                                    format_duration(obs_tot), format_rate(rate),
                                    format_duration(eta));
                }
            } else {
                line = fmt::format(
                    "[{}] {:5.1f}% | file {}/{} | {} obs | {} | ETA {}", bar,
                    frac * 100.0, file_index + 1U, m_nfiles,
                    format_duration(obs_done), format_rate(rate),
                    format_duration(eta));
            }
            emit_bar_line(line);
            return;
        }

        // Log mode (non-TTY / batch logs): emit first update after 5 s (or 5%),
        // then at 5% steps or every 15 s.
        const auto bucket    = static_cast<std::size_t>(frac * 20.0);
        const bool first_due = !m_stream_emitted && elapsed >= 5.0;
        const bool due = (finished && m_last_log_bucket < 20U) || first_due ||
                         (bucket > m_last_log_bucket) || (since_last >= 15.0);
        if (!due) {
            return;
        }
        m_stream_emitted  = true;
        m_last_emit       = now;
        m_last_log_bucket = std::max(m_last_log_bucket, bucket);
        const double rate =
            elapsed > 1.0e-3 ? static_cast<double>(bytes_done) / elapsed : 0.0;
        const double eta = (frac > 1.0e-6 && elapsed > 0.25)
                               ? elapsed * ((1.0 - frac) / frac)
                               : -1.0;
        if (total_samples > 0U) {
            if (m_nfiles > 1U) {
                m_ctx->err << fmt::format(
                    "{}: {:5.1f}% (file {}/{}, {}/{} samples, {}, elapsed {}, "
                    "ETA {})\n",
                    m_label, frac * 100.0, file_index + 1U, m_nfiles,
                    samples_done, total_samples, format_rate(rate),
                    format_duration(elapsed), format_duration(eta));
            } else {
                m_ctx->err << fmt::format(
                    "{}: {:5.1f}% ({}/{} samples, {}, elapsed {}, ETA {})\n",
                    m_label, frac * 100.0, samples_done, total_samples,
                    format_rate(rate), format_duration(elapsed),
                    format_duration(eta));
            }
        } else {
            m_ctx->err << fmt::format(
                "{}: {:5.1f}% (file {}/{}, {} samples, {}, elapsed {}, "
                "ETA {})\n",
                m_label, frac * 100.0, file_index + 1U, m_nfiles, samples_done,
                format_rate(rate), format_duration(elapsed),
                format_duration(eta));
        }
        m_ctx->err.flush();
    }

public:
    /// Clear any active `\r` line and return total elapsed seconds since
    /// begin_stream().
    [[nodiscard]] double finish() {
        clear_line();
        if (!m_stream_started) {
            return 0.0;
        }
        return std::chrono::duration<double>(Clock::now() - m_start).count();
    }

    void clear_line() {
        if (m_mode == ProgressMode::kBar && m_last_bar_len > 0U) {
            m_ctx->err << '\r' << std::string(m_last_bar_len, ' ') << '\r';
            m_ctx->err.flush();
            m_last_bar_len = 0U;
        }
    }

private:
    void emit_bar_line(const std::string& line) {
        const std::size_t pad =
            line.size() < m_last_bar_len ? (m_last_bar_len - line.size()) : 0U;
        m_ctx->err << '\r' << line;
        if (pad > 0U) {
            m_ctx->err << std::string(pad, ' ');
        }
        m_ctx->err.flush();
        m_last_bar_len = line.size();
    }

    const Context* m_ctx{nullptr};
    std::string m_label;
    ProgressMode m_mode{ProgressMode::kNone};
    Clock::time_point m_start;
    Clock::time_point m_last_emit;
    std::size_t m_last_bar_len{0};
    std::size_t m_last_scan_bucket{0};
    std::size_t m_last_log_bucket{0};
    std::size_t m_last_file_logged{0};
    std::uint64_t m_total_samples{0};
    double m_tsamp{0.0};
    std::uint64_t m_total_input_bytes{0};
    std::size_t m_nfiles{1};
    bool m_scan_emitted{false};
    bool m_stream_started{false};
    bool m_stream_emitted{false};
};

} // namespace psrio::cli

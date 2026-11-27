#pragma once

#include "psrfits_bytes.hpp"
#include "psrio/formats/psrfits/extract.hpp"

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fmt/format.h>
#include <string>
#include <utility>
#include <vector>

namespace psrio::test {

/// A synthetic search-mode PSRFITS observation with one dispersed pulse.
///
/// AABBCRCI, 8-bit, time-major rows. Every product holds `kBaseline`
/// except the pulse: at output sample `pulse` (band-centre time) channel c
/// of AA and BB holds `kPulse` at sample `pulse - D_c`, where D_c is dspsr's
/// delay of channel c relative to OBSFREQ. CR and CI hold
/// `kCross + chan`. DAT_SCL is 2, DAT_OFFS is 1, and ZERO_OFF is 127.5, so
/// --scloffs changes every value.
struct BurstFixture {
    static constexpr int kNchan             = 8;
    static constexpr int kNpol              = 4;
    static constexpr int kNsblk             = 64;
    static constexpr double kTbin           = 1.0e-3;
    static constexpr double kObsFreq        = 1400.0;
    static constexpr double kObsBw          = 400.0;
    static constexpr std::uint8_t kBaseline = 100;
    static constexpr std::uint8_t kPulse    = 200;
    static constexpr std::uint8_t kCross    = 50;

    int nrows{16};
    double dm{50.0};
    std::uint64_t pulse{512};

    [[nodiscard]] std::vector<std::int64_t> delays() const {
        const auto freqs = formats::psrfits::dspsr_channel_frequencies(
            kObsFreq, kObsBw, kNchan);
        std::vector<std::int64_t> out;
        out.reserve(freqs.size());
        for (const double freq : freqs) {
            out.push_back(formats::psrfits::dspsr_sample_delay(
                dm, kObsFreq, freq, 1.0 / kTbin));
        }
        return out;
    }

    [[nodiscard]] std::uint8_t
    raw(std::uint64_t sample, int pol, int chan) const {
        if (pol >= 2) {
            return static_cast<std::uint8_t>(kCross + chan);
        }
        const auto at = static_cast<std::int64_t>(pulse) -
                        delays()[static_cast<std::size_t>(chan)];
        return std::cmp_equal(sample, at) ? kPulse : kBaseline;
    }

    /// Write rows [row0, row0 + count) as one file, with STT and NSUBOFFS
    /// for that position in the observation.
    void write(const std::filesystem::path& path, int row0, int count) const {
        FitsFileBytes file;
        file.begin_hdu();
        file.key("SIMPLE", "T");
        file.key("BITPIX", "8");
        file.key("NAXIS", "0");
        file.key("EXTEND", "T");
        file.key("FITSTYPE", "'PSRFITS'");
        file.key("OBS_MODE", "'SEARCH'");
        file.key("TELESCOP", "'Parkes'");
        file.key("FRONTEND", "'UWL'");
        file.key("BACKEND", "'Medusa'");
        file.key("SRC_NAME", "'FRBTest'");
        file.key("RA", "'05:08:03.5080'");
        file.key("DEC", "'+26:03:38.5040'");
        file.key("IBEAM", "1");
        file.key("NRCVR", "2");
        file.key("FD_HAND", "-1");
        file.key("OBSFREQ", fmt::format("{}", kObsFreq));
        file.key("OBSBW", fmt::format("{}", kObsBw));
        file.key("OBSNCHAN", fmt::format("{}", kNchan));
        file.key("STT_IMJD", "59630");
        file.key("STT_SMJD", "3600");
        file.key("STT_OFFS",
                 fmt::format("{:.17g}", 0.5 + (row0 * kNsblk * kTbin)));
        file.end_header();

        const auto per_row = kNchan * kNpol * kNsblk;
        const auto row_bytes =
            16 + (8 * kNchan) + (4 * kNchan) + (8 * kNchan * kNpol) + per_row;
        file.begin_hdu();
        file.key("XTENSION", "'BINTABLE'");
        file.key("BITPIX", "8");
        file.key("NAXIS", "2");
        file.key("NAXIS1", std::to_string(row_bytes));
        file.key("NAXIS2", std::to_string(count));
        file.key("PCOUNT", "0");
        file.key("GCOUNT", "1");
        file.key("TFIELDS", "7");
        file.key("EXTNAME", "'SUBINT'");
        file.key("TTYPE1", "'TSUBINT'");
        file.key("TFORM1", "'1D'");
        file.key("TTYPE2", "'OFFS_SUB'");
        file.key("TFORM2", "'1D'");
        file.key("TTYPE3", "'DAT_FREQ'");
        file.key("TFORM3", fmt::format("'{}D'", kNchan));
        file.key("TTYPE4", "'DAT_WTS'");
        file.key("TFORM4", fmt::format("'{}E'", kNchan));
        file.key("TTYPE5", "'DAT_OFFS'");
        file.key("TFORM5", fmt::format("'{}E'", kNchan * kNpol));
        file.key("TTYPE6", "'DAT_SCL'");
        file.key("TFORM6", fmt::format("'{}E'", kNchan * kNpol));
        file.key("TTYPE7", "'DATA'");
        file.key("TFORM7", fmt::format("'{}B'", per_row));
        file.key("TDIM7", fmt::format("'({},{},{})'", kNchan, kNpol, kNsblk));
        file.key("NPOL", std::to_string(kNpol));
        file.key("NCHAN", std::to_string(kNchan));
        file.key("NBITS", "8");
        file.key("NSBLK", std::to_string(kNsblk));
        file.key("TBIN", fmt::format("{}", kTbin));
        file.key("POL_TYPE", "'AABBCRCI'");
        file.key("ZERO_OFF", "127.5");
        file.key("SIGNINT", "0");
        file.key("NSUBOFFS", std::to_string(row0));
        file.end_header();

        const auto freqs = formats::psrfits::dspsr_channel_frequencies(
            kObsFreq, kObsBw, kNchan);
        std::vector<std::byte> table;
        for (int row = row0; row < row0 + count; ++row) {
            const double tsub = kNsblk * kTbin;
            append_be(table, tsub);
            append_be(table, (row - row0 + 0.5) * tsub);
            for (const double freq : freqs) {
                append_be(table, freq);
            }
            for (int chan = 0; chan < kNchan; ++chan) {
                append_be(table, chan == 0 ? 0.0F : 1.0F);
            }
            for (int index = 0; index < kNchan * kNpol; ++index) {
                append_be(table, 1.0F);
            }
            for (int index = 0; index < kNchan * kNpol; ++index) {
                append_be(table, 2.0F);
            }
            for (int isamp = 0; isamp < kNsblk; ++isamp) {
                const auto sample = (static_cast<std::uint64_t>(row) * kNsblk) +
                                    static_cast<std::uint64_t>(isamp);
                for (int pol = 0; pol < kNpol; ++pol) {
                    for (int chan = 0; chan < kNchan; ++chan) {
                        append_be(table, raw(sample, pol, chan));
                    }
                }
            }
        }
        file.raw(table);
        file.pad_data();
        file.write(path);
    }
};

} // namespace psrio::test

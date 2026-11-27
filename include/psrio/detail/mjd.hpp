#pragma once

#include <cmath>
#include <cstdint>

namespace psrio::detail {

/**
 * @brief An MJD time stamp split into whole days, whole seconds of the day,
 *        and a fraction of a second.
 *
 * A double MJD resolves only about a microsecond, so PSRFITS stores start
 * times as STT_IMJD, STT_SMJD and STT_OFFS. This type keeps the same three
 * parts and does its arithmetic on them: whole seconds are exact integers,
 * and only the sub-second fraction is a double.
 *
 * **Normal form.** After every operation the three parts share the sign of
 * the time they represent: for a time at or after MJD 0, `second()` is in
 * [0, 86400) and `fraction()` in [0, 1); for a negative interval both are
 * zero or negative. Signed values appear only as differences.
 *
 * **Rounding.** The only inexact steps are double operations on the
 * fraction: splitting an added interval into whole and fractional seconds,
 * adding fractions, moving one second between the fraction and the whole
 * part, and halving for a midpoint. Keeping to these steps, and to this
 * order, makes STT_OFFS and OFFS_SUB come out exactly as PSRCHIVE and dspsr
 * store them (checked against dspsr output in tests/test_cutout.cpp).
 */
class Mjd {
public:
    static constexpr std::int64_t kSecondsPerDay = 86400;

    Mjd() = default;

    /// Time from its three parts, which may be out of range or of mixed sign.
    Mjd(std::int64_t day, std::int64_t second, double fraction) {
        assign((day * kSecondsPerDay) + second, fraction);
    }

    [[nodiscard]] std::int64_t day() const noexcept { return m_day; }
    [[nodiscard]] std::int64_t second() const noexcept { return m_second; }
    [[nodiscard]] double fraction() const noexcept { return m_fraction; }

    /// Seconds since MJD 0, as one double.
    [[nodiscard]] double in_seconds() const noexcept {
        return (static_cast<double>(m_day) *
                static_cast<double>(kSecondsPerDay)) +
               static_cast<double>(m_second) + m_fraction;
    }

    /// Seconds from @p earlier to this time (negative if this is earlier).
    [[nodiscard]] double seconds_since(const Mjd& earlier) const {
        return Mjd(m_day - earlier.m_day, m_second - earlier.m_second,
                   m_fraction - earlier.m_fraction)
            .in_seconds();
    }

    /// The time halfway between @p first and @p second.
    [[nodiscard]] static Mjd midpoint(const Mjd& first, const Mjd& second) {
        const Mjd sum(first.m_day + second.m_day,
                      first.m_second + second.m_second,
                      first.m_fraction + second.m_fraction);
        const auto whole = sum.whole_seconds();
        const double odd = (whole % 2 != 0) ? 0.5 : 0.0;
        return {0, whole / 2, (sum.m_fraction * 0.5) + odd};
    }

    /// This time moved by @p seconds.
    friend Mjd operator+(const Mjd& time, double seconds) {
        const double whole = std::trunc(seconds);
        return {time.m_day, time.m_second + static_cast<std::int64_t>(whole),
                time.m_fraction + (seconds - whole)};
    }

    friend Mjd operator-(const Mjd& time, double seconds) {
        return time + -seconds;
    }

    Mjd& operator+=(double seconds) { return *this = *this + seconds; }
    Mjd& operator-=(double seconds) { return *this = *this - seconds; }

    friend bool operator==(const Mjd&, const Mjd&) = default;

private:
    [[nodiscard]] std::int64_t whole_seconds() const noexcept {
        return (m_day * kSecondsPerDay) + m_second;
    }

    /// Store @p whole seconds plus @p fraction in normal form.
    void assign(std::int64_t whole, double fraction) {
        // Whole seconds that built up in the fraction.
        const double carry = std::trunc(fraction);
        whole += static_cast<std::int64_t>(carry);
        fraction -= carry;
        // Give the fraction the sign of the whole time.
        const bool negative = whole < 0 || (whole == 0 && fraction < 0.0);
        if (!negative && fraction < 0.0) {
            whole -= 1;
            fraction += 1.0;
        } else if (negative && fraction > 0.0) {
            whole += 1;
            fraction -= 1.0;
        }
        // Integer division truncates, so day and second share that sign.
        m_day      = whole / kSecondsPerDay;
        m_second   = whole % kSecondsPerDay;
        m_fraction = fraction;
    }

    std::int64_t m_day{0};
    std::int64_t m_second{0};
    double m_fraction{0.0};
};

} // namespace psrio::detail

#pragma once

#include "psrio/common/types.hpp"
#include "psrio/detail/exceptions.hpp"

#include <concepts>
#include <cstdint>
#include <ctime>
#include <memory>
#include <span>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

namespace psrio::concepts {

/**
 * @brief Packed intensity stream read one block of time samples at a time.
 *
 * Modeled by filterbank, in-memory, FBH5, and PSRDADA readers. GUPPI RAW is
 * baseband and does not model this concept.
 */
template <typename R>
concept BlockReader = requires(R& reader,
                               const R& reader_const,
                               std::uint64_t count,
                               std::int64_t delta,
                               std::span<std::byte> bytes,
                               std::span<float> floats,
                               std::span<std::uint8_t> u8s,
                               std::span<std::uint16_t> u16s) {
    { reader_const.nchans() } -> std::same_as<std::uint64_t>;
    { reader_const.nifs() } -> std::same_as<std::uint64_t>;
    { reader_const.nbits() } -> std::same_as<int>;
    { reader_const.sample_type() } -> std::same_as<SampleType>;
    { reader_const.bytes_per_sample() } -> std::same_as<std::uint64_t>;
    { reader_const.nsamples() } -> std::same_as<std::uint64_t>;
    { reader_const.has_nsamples() } -> std::same_as<bool>;
    { reader_const.tsamp() } -> std::same_as<double>;
    { reader_const.tstart() } -> std::same_as<double>;
    { reader_const.source_name() } -> std::convertible_to<std::string_view>;
    { reader_const.telescope() } -> std::convertible_to<std::string_view>;
    { reader_const.raj() } -> std::same_as<double>;
    { reader_const.dej() } -> std::same_as<double>;
    { reader_const.fch1() } -> std::same_as<double>;
    { reader_const.foff() } -> std::same_as<double>;
    { reader_const.bandwidth() } -> std::same_as<double>;
    { reader_const.center_frequency() } -> std::same_as<double>;
    { reader_const.beam() } -> std::same_as<int>;
    { reader_const.spectra_rate() } -> std::same_as<double>;
    { reader_const.utc_start() } -> std::same_as<std::time_t>;
    { reader_const.tell() } -> std::same_as<std::uint64_t>;
    { reader.set_fswap(bool{}) };
    { reader_const.fswap_enabled() } -> std::same_as<bool>;
    { reader.seek(count) };
    { reader.rewind() };
    { reader.skip(delta) };
    { reader.read_block(count, bytes) } -> std::same_as<std::uint64_t>;
    { reader.read_samples(count, floats) } -> std::same_as<std::uint64_t>;
    { reader.read_samples(count, u8s) } -> std::same_as<std::uint64_t>;
    { reader.read_samples(count, u16s) } -> std::same_as<std::uint64_t>;
    { reader.read_bytes(count, bytes) } -> std::same_as<std::uint64_t>;
};

} // namespace psrio::concepts

namespace psrio {

/// Move-only type-erased intensity block reader.
///
/// Concrete readers stay free of a base class. This handle is the runtime
/// front a pipeline stores when the file type is chosen at run time. It
/// accepts only types that model `concepts::BlockReader`.
class BlockSource {
public:
    template <concepts::BlockReader R>
        requires(!std::same_as<std::remove_cvref_t<R>, BlockSource>)
    explicit BlockSource(R reader)
        : m_impl(std::make_unique<Model<R>>(std::move(reader))) {}

    BlockSource(const BlockSource&)                = delete;
    BlockSource& operator=(const BlockSource&)     = delete;
    BlockSource(BlockSource&&) noexcept            = default;
    BlockSource& operator=(BlockSource&&) noexcept = default;
    ~BlockSource()                                 = default;

    [[nodiscard]] std::uint64_t nchans() const { return impl().nchans(); }
    [[nodiscard]] std::uint64_t nifs() const { return impl().nifs(); }
    [[nodiscard]] int nbits() const { return impl().nbits(); }
    [[nodiscard]] SampleType sample_type() const {
        return impl().sample_type();
    }
    [[nodiscard]] std::uint64_t bytes_per_sample() const {
        return impl().bytes_per_sample();
    }
    [[nodiscard]] std::uint64_t nsamples() const { return impl().nsamples(); }
    [[nodiscard]] bool has_nsamples() const { return impl().has_nsamples(); }
    [[nodiscard]] double tsamp() const { return impl().tsamp(); }
    [[nodiscard]] double tstart() const { return impl().tstart(); }
    [[nodiscard]] std::string_view source_name() const {
        return impl().source_name();
    }
    [[nodiscard]] std::string_view telescope() const {
        return impl().telescope();
    }
    [[nodiscard]] double raj() const { return impl().raj(); }
    [[nodiscard]] double dej() const { return impl().dej(); }
    [[nodiscard]] double fch1() const { return impl().fch1(); }
    [[nodiscard]] double foff() const { return impl().foff(); }
    [[nodiscard]] double bandwidth() const { return impl().bandwidth(); }
    [[nodiscard]] double center_frequency() const {
        return impl().center_frequency();
    }
    [[nodiscard]] int beam() const { return impl().beam(); }
    [[nodiscard]] double spectra_rate() const { return impl().spectra_rate(); }
    [[nodiscard]] std::time_t utc_start() const { return impl().utc_start(); }
    [[nodiscard]] std::uint64_t tell() const { return impl().tell(); }

    void set_fswap(bool enable) { impl().set_fswap(enable); }
    [[nodiscard]] bool fswap_enabled() const { return impl().fswap_enabled(); }

    void seek(std::uint64_t sample) { impl().seek(sample); }
    void rewind() { impl().rewind(); }
    void skip(std::int64_t delta) { impl().skip(delta); }

    std::uint64_t read_block(std::uint64_t count, std::span<std::byte> dest) {
        return impl().read_block(count, dest);
    }

    std::uint64_t read_samples(std::uint64_t count, std::span<float> dest) {
        return impl().read_samples(count, dest);
    }

    std::uint64_t read_samples(std::uint64_t count,
                               std::span<std::uint8_t> dest) {
        return impl().read_samples(count, dest);
    }

    std::uint64_t read_samples(std::uint64_t count,
                               std::span<std::uint16_t> dest) {
        return impl().read_samples(count, dest);
    }

    std::uint64_t read_bytes(std::uint64_t nbytes, std::span<std::byte> dest) {
        return impl().read_bytes(nbytes, dest);
    }

    /// Convenience allocating read: copy next @p count packed time samples.
    [[nodiscard]] std::vector<std::byte> read_block(std::uint64_t count) {
        const auto bytes_needed = count * bytes_per_sample();
        std::vector<std::byte> out(static_cast<std::size_t>(bytes_needed));
        const auto actual = read_block(count, std::span<std::byte>(out));
        out.resize(static_cast<std::size_t>(actual * bytes_per_sample()));
        return out;
    }

    /// Convenience allocating read: copy next @p nbytes raw payload bytes.
    [[nodiscard]] std::vector<std::byte> read_bytes(std::uint64_t nbytes) {
        std::vector<std::byte> out(static_cast<std::size_t>(nbytes));
        const auto actual = read_bytes(nbytes, std::span<std::byte>(out));
        out.resize(static_cast<std::size_t>(actual));
        return out;
    }

    /// Convenience allocating read: unpack next @p count time samples into a
    /// new vector.
    template <typename T = float>
    [[nodiscard]] std::vector<T> read_samples(std::uint64_t count) {
        const auto values_needed = count * nifs() * nchans();
        std::vector<T> out(static_cast<std::size_t>(values_needed));
        const auto actual = read_samples(count, std::span<T>(out));
        out.resize(static_cast<std::size_t>(actual * nifs() * nchans()));
        return out;
    }

private:
    struct Concept {
        Concept(const Concept&)            = delete;
        Concept& operator=(const Concept&) = delete;
        Concept(Concept&&)                 = delete;
        Concept& operator=(Concept&&)      = delete;

        virtual ~Concept() = default;

        [[nodiscard]] virtual std::uint64_t nchans() const           = 0;
        [[nodiscard]] virtual std::uint64_t nifs() const             = 0;
        [[nodiscard]] virtual int nbits() const                      = 0;
        [[nodiscard]] virtual SampleType sample_type() const         = 0;
        [[nodiscard]] virtual std::uint64_t bytes_per_sample() const = 0;
        [[nodiscard]] virtual std::uint64_t nsamples() const         = 0;
        [[nodiscard]] virtual bool has_nsamples() const              = 0;
        [[nodiscard]] virtual double tsamp() const                   = 0;
        [[nodiscard]] virtual double tstart() const                  = 0;
        [[nodiscard]] virtual std::string_view source_name() const   = 0;
        [[nodiscard]] virtual std::string_view telescope() const     = 0;
        [[nodiscard]] virtual double raj() const                     = 0;
        [[nodiscard]] virtual double dej() const                     = 0;
        [[nodiscard]] virtual double fch1() const                    = 0;
        [[nodiscard]] virtual double foff() const                    = 0;
        [[nodiscard]] virtual double bandwidth() const               = 0;
        [[nodiscard]] virtual double center_frequency() const        = 0;
        [[nodiscard]] virtual int beam() const                       = 0;
        [[nodiscard]] virtual double spectra_rate() const            = 0;
        [[nodiscard]] virtual std::time_t utc_start() const          = 0;
        [[nodiscard]] virtual std::uint64_t tell() const             = 0;
        virtual void set_fswap(bool enable)                          = 0;
        [[nodiscard]] virtual bool fswap_enabled() const             = 0;
        virtual void seek(std::uint64_t sample)                      = 0;
        virtual void rewind()                                        = 0;
        virtual void skip(std::int64_t delta)                        = 0;
        [[nodiscard]] virtual std::uint64_t
        read_block(std::uint64_t count, std::span<std::byte> dest) = 0;
        [[nodiscard]] virtual std::uint64_t
        read_samples(std::uint64_t count, std::span<float> dest) = 0;
        [[nodiscard]] virtual std::uint64_t
        read_samples(std::uint64_t count, std::span<std::uint8_t> dest) = 0;
        [[nodiscard]] virtual std::uint64_t
        read_samples(std::uint64_t count, std::span<std::uint16_t> dest) = 0;
        [[nodiscard]] virtual std::uint64_t
        read_bytes(std::uint64_t nbytes, std::span<std::byte> dest) = 0;

    protected:
        Concept() = default;
    };

    template <concepts::BlockReader R> struct Model final : Concept {
        explicit Model(R reader) : m_reader(std::move(reader)) {}

        [[nodiscard]] std::uint64_t nchans() const override {
            return m_reader.nchans();
        }
        [[nodiscard]] std::uint64_t nifs() const override {
            return m_reader.nifs();
        }
        [[nodiscard]] int nbits() const override { return m_reader.nbits(); }
        [[nodiscard]] SampleType sample_type() const override {
            return m_reader.sample_type();
        }
        [[nodiscard]] std::uint64_t bytes_per_sample() const override {
            return m_reader.bytes_per_sample();
        }
        [[nodiscard]] std::uint64_t nsamples() const override {
            return m_reader.nsamples();
        }
        [[nodiscard]] bool has_nsamples() const override {
            return m_reader.has_nsamples();
        }
        [[nodiscard]] double tsamp() const override { return m_reader.tsamp(); }
        [[nodiscard]] double tstart() const override {
            return m_reader.tstart();
        }
        [[nodiscard]] std::string_view source_name() const override {
            return m_reader.source_name();
        }
        [[nodiscard]] std::string_view telescope() const override {
            return m_reader.telescope();
        }
        [[nodiscard]] double raj() const override { return m_reader.raj(); }
        [[nodiscard]] double dej() const override { return m_reader.dej(); }
        [[nodiscard]] double fch1() const override { return m_reader.fch1(); }
        [[nodiscard]] double foff() const override { return m_reader.foff(); }
        [[nodiscard]] double bandwidth() const override {
            return m_reader.bandwidth();
        }
        [[nodiscard]] double center_frequency() const override {
            return m_reader.center_frequency();
        }
        [[nodiscard]] int beam() const override { return m_reader.beam(); }
        [[nodiscard]] double spectra_rate() const override {
            return m_reader.spectra_rate();
        }
        [[nodiscard]] std::time_t utc_start() const override {
            return m_reader.utc_start();
        }
        [[nodiscard]] std::uint64_t tell() const override {
            return m_reader.tell();
        }
        void set_fswap(bool enable) override { m_reader.set_fswap(enable); }
        [[nodiscard]] bool fswap_enabled() const override {
            return m_reader.fswap_enabled();
        }
        void seek(std::uint64_t sample) override { m_reader.seek(sample); }
        void rewind() override { m_reader.rewind(); }
        void skip(std::int64_t delta) override { m_reader.skip(delta); }
        [[nodiscard]] std::uint64_t
        read_block(std::uint64_t count, std::span<std::byte> dest) override {
            return m_reader.read_block(count, dest);
        }
        [[nodiscard]] std::uint64_t
        read_samples(std::uint64_t count, std::span<float> dest) override {
            return m_reader.read_samples(count, dest);
        }
        [[nodiscard]] std::uint64_t
        read_samples(std::uint64_t count,
                     std::span<std::uint8_t> dest) override {
            return m_reader.read_samples(count, dest);
        }
        [[nodiscard]] std::uint64_t
        read_samples(std::uint64_t count,
                     std::span<std::uint16_t> dest) override {
            return m_reader.read_samples(count, dest);
        }
        [[nodiscard]] std::uint64_t
        read_bytes(std::uint64_t nbytes, std::span<std::byte> dest) override {
            return m_reader.read_bytes(nbytes, dest);
        }

        R m_reader;
    };

    [[nodiscard]] const Concept& impl() const {
        if (m_impl == nullptr) {
            throw ValidationError("psrio: block source is empty");
        }
        return *m_impl;
    }

    [[nodiscard]] Concept& impl() {
        if (m_impl == nullptr) {
            throw ValidationError("psrio: block source is empty");
        }
        return *m_impl;
    }

    std::unique_ptr<Concept> m_impl;
};

} // namespace psrio

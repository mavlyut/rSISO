#ifndef GRAY_CODE_H
#define GRAY_CODE_H

#include <cstdint>
#include <iterator>
#include <vector>

#include "defines.h"

namespace long_domain {
    struct gray_code {
    private:
        using _digits_t = std::size_t;

        static_assert(chunk_size <= std::numeric_limits<_digits_t>::digits);
        static inline constexpr _digits_t __min(_digits_t const& a, _digits_t const& b) {
            return (a < b ? a : b);
        }
        static constexpr const _digits_t MAX_DIGIT
            = (chunk_size < 64)
            ? __min(std::numeric_limits<_digits_t>::max(), (1ull << chunk_size))
            : std::numeric_limits<_digits_t>::max();

    public:
        gray_code(unsigned, bool with_zero = true);

        class gray_code_iterator : std::iterator<
                        std::input_iterator_tag,
                        unsigned,
                        std::ptrdiff_t,
                        unsigned const*,
                        unsigned const&
                    > {
        friend struct gray_code;

        private:
            explicit gray_code_iterator(unsigned, bool with_zero, bool is_final);

        public:
            gray_code_iterator& operator++();
            gray_code_iterator operator++(signed);

            bool operator==(gray_code_iterator const&) const;
            bool operator!=(gray_code_iterator const&) const;

            reference operator*() const;

        private:
            unsigned n, change_bit;
            std::vector<_digits_t> ind;

            void next_vector();
        };

        gray_code_iterator begin() const;
        gray_code_iterator end() const;

    private:
        bool with_zero;
        unsigned n;
    };
}

namespace short_domain {
    struct gray_code {
        gray_code(unsigned);

        class gray_code_iterator : std::iterator<
                        std::input_iterator_tag,
                        unsigned,
                        std::ptrdiff_t,
                        unsigned const*,
                        unsigned const&
                    > {
        friend struct gray_code;

        private:
            explicit gray_code_iterator(unsigned, std::size_t);

        public:
            gray_code_iterator& operator++();
            gray_code_iterator operator++(signed);

            bool operator==(gray_code_iterator const&) const;
            bool operator!=(gray_code_iterator const&) const;

            reference operator*() const;

        private:
            unsigned n, change_bit;
            std::size_t tmp_num;
        };

        gray_code_iterator begin() const;
        gray_code_iterator end() const;

    private:
        unsigned n;
    };
}

#endif // GRAY_CODE_H

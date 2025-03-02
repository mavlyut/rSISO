#include <cstdint>
#include <iterator>
#include <limits>
#include <vector>

static const unsigned UNINIT = -1;
using _digits_t = std::size_t;
static const unsigned chunk_size = 64;
static_assert(chunk_size <= std::numeric_limits<_digits_t>::digits);
static inline constexpr _digits_t __min(_digits_t const& a, _digits_t const& b) {
    return (a < b ? a : b);
}
static constexpr const _digits_t MAX_DIGIT
    = (chunk_size < 64)
    ? __min(std::numeric_limits<_digits_t>::max(), (1ull << chunk_size))
    : std::numeric_limits<_digits_t>::max();


struct GrayCode {
public:
    GrayCode(unsigned n) : n(n) {}

    class iterator : std::iterator<
                    std::input_iterator_tag,
                    unsigned,
                    std::ptrdiff_t,
                    unsigned const*,
                    unsigned const&
                > {
    public:
        explicit iterator(unsigned n, bool is_final)
            : n(n), change_bit(UNINIT), ind((n + chunk_size) / chunk_size) {
            if (is_final) {
                ind[n / chunk_size] = (1 << (n % chunk_size));
            }
        }

        iterator& operator++() {
            next_vector();
            return *this;
        }

        iterator operator++(signed) {
            iterator ret = *this;
            ++*this;
            return ret;
        }

        bool operator==(iterator const& other) const {
            return ind == other.ind;
        }
        
        bool operator!=(iterator const& other) const {
            return !(*this == other);
        }
        
        reference operator*() const {
            return change_bit;
        }

    private:
        unsigned n, change_bit;
        std::vector<_digits_t> ind;

        void next_vector() {
            for (unsigned i = 0; i < ind.size(); i++) {
                _digits_t& c = ind[i];
                if (c == MAX_DIGIT) {
                    c = 0;
                } else {
                    _digits_t d = c + 1;
                    for (unsigned j = chunk_size; j-- > 0; ) {
                        if (((c >> j) & 1) != ((d >> j) & 1)) {
                            change_bit = i * chunk_size + j;
                            break;
                        }
                    }
                    c = d;
                    break;
                }
            }
        }
    };

    iterator begin() const {
        return iterator(n, false);
    }

    iterator end() const {
        return iterator(n, true);
    }

private:
    unsigned n;
};

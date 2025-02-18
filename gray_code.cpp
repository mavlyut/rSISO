#include <iterator>
#include <cstdint>

struct GrayCode {
    GrayCode(unsigned n) : n(n) {}

    class iterator : std::iterator<
                    std::input_iterator_tag,
                    unsigned,
                    std::ptrdiff_t,
                    unsigned const*,
                    unsigned const&
                > {
    public:
        explicit iterator(unsigned n, std::size_t ind) : n(n), i(-1), ind(ind) {}

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
            return i;
        }

    private:
        unsigned n, i;
        std::size_t ind;

        void next_vector() {
            std::size_t a = ind, b = ind + 1;
            for (i = n; i-- > 0; ) {
                if (((a >> i) & 1) != ((b >> i) & 1)) {
                    break;
                }
            }
            ind++;
        }
    };

    iterator begin() const {
        return iterator(n, 0);
    }

    iterator end() const {
        return iterator(n, (1 << n));
    }

private:
    unsigned n;
};

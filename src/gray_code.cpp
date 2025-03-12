#include <cstdint>
#include <iterator>
#include <limits>
#include <vector>

#include "../include/defines.h"
#include "../include/gray_code.h"

gray_code::gray_code(unsigned n) : n(n) {}

gray_code::gray_code_iterator::gray_code_iterator(unsigned n, bool is_final)
    : n(n), change_bit(UNINIT), ind((n + chunk_size) / chunk_size) {
    if (is_final) {
        ind[n / chunk_size] = (1 << (n % chunk_size));
    }
}

gray_code::gray_code_iterator& gray_code::gray_code_iterator::operator++() {
    next_vector();
    return *this;
}

gray_code::gray_code_iterator gray_code::gray_code_iterator::operator++(signed) {
    gray_code_iterator ret = *this;
    ++*this;
    return ret;
}

bool gray_code::gray_code_iterator::operator==(gray_code_iterator const& other) const {
    return ind == other.ind;
}

bool gray_code::gray_code_iterator::operator!=(gray_code_iterator const& other) const {
    return !(*this == other);
}

gray_code::gray_code_iterator::reference gray_code::gray_code_iterator::operator*() const {
    return change_bit;
}

void gray_code::gray_code_iterator::next_vector() {
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

gray_code::gray_code_iterator gray_code::begin() const {
    return gray_code_iterator(n, false);
}

gray_code::gray_code_iterator gray_code::end() const {
    return gray_code_iterator(n, true);
}

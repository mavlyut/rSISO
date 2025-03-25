#include "../include/turbo_decoder.h"

turbo_decoder::turbo_decoder(soft_decoder const& dec_row, soft_decoder const& dec_col)
    : n1(dec_row.length()), n2(dec_col.length()), k1(dec_row.dim()), k2(dec_col.dim())
    , dec_row(dec_row), dec_col(dec_col) {}

unsigned turbo_decoder::length() const {
    return n1 * n2;
}

unsigned turbo_decoder::dim() const {
    return k1 * k2;
}

binvector turbo_decoder::encode(binvector const& c) const {
    fail(c.size() == dim(), "turbo, encode: incorrect dim");
    matrix enc1t(n1, binvector(k2));
    for (unsigned i = 0; i < k2; i++) {
        binvector row = dec_row.encode(c.subvector(i * k1, i * (k1 + 1)));
        for (unsigned j = 0; j < n1; j++) {
            enc1t[j].set(i, row[j]);
        }
    }
    binvector ans(length());
    for (unsigned j = 0; j < n1; j++) {
        binvector row = dec_col.encode(enc1t[j]);
        for (unsigned i = 0; i < n2; i++) {
            ans.set(j * n2 + i, row[i]);
        }
    }
    return ans;
}

std::vector<double> turbo_decoder::soft_decode(std::vector<double> const& L0) {
    fail(L0.size() == length(), "turbo, decode: incorrect dim");
    std::vector<std::vector<double>> dec2t(n2, std::vector<double>(n1));
    std::vector<double> L(L0);
    for (unsigned cnt = 0; cnt < ITER_CNT; cnt++) {
        for (unsigned i = 0; i < n1; i++) {
            auto col = dec_col.decode_soft(std::vector(L.begin() + i * n2, L.begin() + (i + 1) * n2));
            for (unsigned j = 0; j < n2; j++) {
                dec2t[j][i] = col[j];
            }
        }
        for (unsigned j = 0; j < n2; j++) {
            auto row = dec_row.decode_soft(dec2t[j]);
            for (unsigned i = 0; i < n1; i++) {
                L[i * n1 + j] = row[i];
            }
        }
    }
    return L;
}

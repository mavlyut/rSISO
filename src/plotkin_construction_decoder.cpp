#include <cmath>

#include "../include/plotkin_construction_decoder.h"

plotkin_construction_decoder::plotkin_construction_decoder(soft_decoder* dec1, soft_decoder* dec2)
        : soft_decoder(dec1->length() * 2, dec1->dim() + dec2->dim()), dec1(dec1), dec2(dec2) {
    fail(dec1->length() == dec2->length(), "Pt, ctor: incompatible sizes");
    __log("Pt, created" << std::endl);
}

binvector plotkin_construction_decoder::encode(binvector const& c) const {
    fail(c.size() == dim(), "Pt, encode: incorrect dim");
    __log("Pt, encode: " << c << std::endl);
    binvector u1 = dec1->encode(c.subvector(0, dec1->dim()));
    binvector u2 = (dec2->encode(c.subvector(dec1->dim(), dec1->dim() + dec2->dim())) ^ u1);
    binvector ans(length());
    unsigned i = 0;
    for (; i < u1.size(); i++) {
        ans.set(i, u1[i]);
    }
    for (unsigned j = 0; j < u2.size(); i++, j++) {
        ans.set(i, u2[j]);
    }
    return ans;
}

double sign(double x) {
    if (x < 0) {
        return -1;
    }
    if (x > 0) {
        return 1;
    }
    return 0;
}

std::vector<double> plotkin_construction_decoder::decode_soft(std::vector<double> const& L0) {
    fail(L0.size() == length(), "Pt, decode: incorrect size");
    __log("Pt, decode: " << L0 << std::endl);
    unsigned n = length(), m = n / 2;
    std::vector<std::array<double, 3>> Lr(m);
    std::vector<double> Lq0(m, 0), Lq1(m), Lq2(m);
    for (unsigned i = 0; i < m; i++) {
        Lq1[i] = truncate(L0[i]);
        Lq2[i] = truncate(L0[i + m]);
    }
    for (unsigned _t = 0; _t < MAX_ITER_COUNT; _t++) {
        Lq0 = dec1->decode_soft(Lq0);
        Lq2 = dec2->decode_soft(Lq2);

        for (unsigned j = 0; j < m; j++) {
            Lr[j][0] = truncate((sign(Lq1[j]) * sign(Lq2[j])) * phi(phi(std::abs(Lq1[j])) + phi(std::abs(Lq2[j]))));
            Lr[j][1] = truncate((sign(Lq2[j]) * sign(Lq0[j])) * phi(phi(std::abs(Lq2[j])) + phi(std::abs(Lq0[j]))));
            Lr[j][2] = truncate((sign(Lq0[j]) * sign(Lq1[j])) * phi(phi(std::abs(Lq0[j])) + phi(std::abs(Lq1[j]))));
        }

        for (unsigned i = 0; i < m; i++) {
            Lq0[i] = truncate(Lq0[i] + Lr[i][0]);
            Lq1[i] = truncate(Lq1[i] + Lr[i][1]);
            Lq2[i] = truncate(Lq2[i] + Lr[i][2]);
        }

        __log("After " << _t + 1 << " iteration: " << Lq1 << " " << Lq2 << std::endl);

        bool is_codeword = true;
        for (unsigned j = 0; j < m; j++) {
            if (!((Lq0[j] < 0) ^ (Lq1[j] < 0) ^ (Lq2[j] < 0))) {
                is_codeword = false;
                break;
            }
        }
        if (is_codeword) {
            break;
        }
    }
    std::vector<double> L(n);
    for (unsigned i = 0; i < m; i++) {
        L[i] = Lq1[i];
        L[i + m] = Lq2[i];
    }
    return L;
}

double plotkin_construction_decoder::phi(double x) {
    if (x >= 7) {
        return 2.0 / (exp(x) + 1);
    }
    if (x <= 5e-2) {
        return log(2.0 / x);
    }
    return -log(tanh(x / 2.0));  
}

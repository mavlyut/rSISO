#include <cmath>

#include "../include/plotkin_construction_decoder.h"

namespace short_domain {
    plotkin_construction_decoder::plotkin_construction_decoder(soft_decoder* dec1, short_domain::soft_decoder* dec2)
            : soft_decoder(dec1->length() * 2, dec1->dim() + dec2->dim()), m(dec1->length())
            , Lq0(m, 0), Lq1(m), Lq2(m), L_out(n), Lr0(m), Lr1(m), Lr2(m), L_ext_0(m), L_ext_2(m)
            , dec1(dec1), dec2(dec2) {
        fail(dec1->length() == dec2->length(), "Pt, ctor: incompatible sizes");
        __log("Pt, created" << std::endl);
    }

    binvector plotkin_construction_decoder::encode(binvector const& c) {
        _log_bv("Pt, encode: ", k, c); __log(std::endl);
        _log_bv("Pt, c1: ", dec1->dim(), subvector(c, 0, dec1->dim())); __log(std::endl);
        _log_bv("Pt, c2: ", dec2->dim(), subvector(c, dec1->dim(), dec1->dim() + dec2->dim())); __log(std::endl);
        auto u1 = dec1->encode(subvector(c, 0, dec1->dim()));
        auto u2 = dec2->encode(subvector(c, dec1->dim(), dec1->dim() + dec2->dim()));
        _log_bv("Pt, part1: ", m, u1); __log(std::endl);
        _log_bv("Pt, part2: ", m, u2); __log(std::endl);
        _log_bv("Pt, encoded: ", n, (((u1 ^ u2) << m) | u1)); __log(std::endl);
        return (((u1 ^ u2) << m) | u1);
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

    const unsigned plotkin_construction_decoder::MAX_ITER_COUNT = 20;
    std::vector<double> plotkin_construction_decoder::decode_soft(std::vector<double> const& L0) {
        fail(L0.size() == length(), "Pt, decode: incorrect size");
        __log("Pt, decode: " << L0 << std::endl);

        for (unsigned i = 0; i < m; ++i) {
            Lq0[i] = truncate(L0[i]);
            Lq1[i] = truncate(L0[i + m]);
            // double ea = exp(L0[i]), eb = exp(L0[i + m]);
            // Lq2[i] = truncate(log((ea * eb + 1) / (ea + eb)));
            Lq2[i] = 0.0;
            Lr0[i] = Lr1[i] = Lr2[i] = 0.0;
        }

        for (unsigned _t = 1; _t <= MAX_ITER_COUNT; ++_t) {
            for (unsigned j = 0; j < m; ++j) {
                Lr0[j] = truncate((sign(Lq1[j]) * sign(Lq2[j])) * phi(phi(std::abs(Lq1[j])) + phi(std::abs(Lq2[j]))));
                Lr1[j] = truncate((sign(Lq2[j]) * sign(Lq0[j])) * phi(phi(std::abs(Lq2[j])) + phi(std::abs(Lq0[j]))));
                Lr2[j] = truncate((sign(Lq0[j]) * sign(Lq1[j])) * phi(phi(std::abs(Lq0[j])) + phi(std::abs(Lq1[j]))));
            }

            L_ext_0 = dec1->decode_soft(Lr0);
            L_ext_2 = dec2->decode_soft(Lr2);
            for (unsigned i = 0; i < m; i++) {
                Lq0[i] = truncate(Lq0[i] + L_ext_0[i]);
                Lq1[i] = truncate(Lq1[i] + Lr1[i]);
                Lq2[i] = truncate(Lq2[i] + L_ext_2[i]);
            }
    
            __log("After " << _t + 1 << " iteration: " << Lq0 << " " << Lq1 << std::endl);
        }

        for (unsigned i = 0; i < m; ++i) {
            L_out[i] = Lq0[i];
            L_out[i + m] = Lq1[i];
        }
        return L_out;
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

	matrix plotkin_construction_decoder::generate_matrix() const {
        matrix G1 = dec1->generate_matrix();
        matrix G2 = dec2->generate_matrix();
        matrix G;
        for (auto& i : G1) {
            G.push_back(i | (i << m));
        }
        for (auto& i : G2) {
            G.push_back(i << m);
        }
        return G;
    }
}

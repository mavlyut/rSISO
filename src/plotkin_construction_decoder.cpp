#include <cmath>

#include "../include/plotkin_construction_decoder.h"

plotkin_construction_decoder::plotkin_construction_decoder(recursive_decoder* dec1, recursive_decoder* dec2)
        : n(dec1->length() * 2), m(dec1->length()), k(dec1->dim() + dec2->dim())
        , Lq0(m, 0), Lq1(m), Lq2(m), L_out(n), Lr0(m), Lr1(m), Lr2(m), L_ext_0(m), L_ext_2(m)
        , dec1(dec1), dec2(dec2) {
    fail(dec1->length() == dec2->length(), "Pt, ctor: incompatible sizes");
    __log("Pt, created" << std::endl);
}

#define getbit(v, i) (((v) >> (i)) & 1)
#define subvector(v, x, y) (((v) >> (x)) & ((1ull << ((y) - (x))) - 1))
#define setbit(v, i, b) if (b) {v |= (1ull << (i));} else {v &= ~(1ull << (i));}
#define printbv(msg, n, v) 					\
	__log(msg);								\
	for (unsigned j = 0; j < n; ++j) { 		\
		if (j != 0) __log(" ");				\
		__log((((v >> j) & 1) ? 1 : 0));	\
	}                                       \
    __log(std::endl);

std::size_t plotkin_construction_decoder::encode(std::size_t const& c) const {
    printbv("Pt, encode: ", k, c);
    auto u1 = dec1->encode(subvector(c, 0, dec1->dim()));
    auto u2 = dec2->encode(subvector(c, dec1->dim(), dec1->dim() + dec2->dim()));
    printbv("Pt, encode, half1-encoded: ", m, u1);
    printbv("Pt, encode, half2-encoded: ", m, u2);
    printbv("Pt, encode, half1^half2  : ", m, (u1 ^ u2));
    std::size_t ans = (((u1 ^ u2) << m) | u1);
    printbv("Pt, encoded: ", n, ans);
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

const double MAX_TRUNC = 32;
double truncate(double const &a) {
    if (a > MAX_TRUNC) {
        return MAX_TRUNC;
    }
    if (a < -MAX_TRUNC) {
        return -MAX_TRUNC;
    }
    return a;
}

const unsigned plotkin_construction_decoder::MAX_ITER_COUNT = 20;
std::vector<double> plotkin_construction_decoder::decode_soft(std::vector<double> const& L0) {
    fail(L0.size() == length(), "Pt, decode: incorrect size");
    __log("Pt, decode: " << L0 << std::endl);

    for (unsigned i = 0; i < m; ++i) {
        Lq0[i] = truncate(L0[i]);
        Lq1[i] = truncate(L0[i + m]);
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

std::pair<double, double> plotkin_construction_decoder::simulate(double snr, unsigned iter_cnt, unsigned max_error) {
    _Float64 sigma = sqrt(0.5 * pow(10.0, -snr / 10.0) * n / k);
    std::normal_distribution<_Float64> norm(0.0, sigma);
    unsigned errr = 0, berr = 0, sim_cnt = 0;
    while (errr < max_error && sim_cnt < iter_cnt) {
        std::size_t x = 0;
        for (unsigned t = 0; t < k; ++t) {
			if (rand() % 2) {
				x |= (1ull << t);
			}
        }
        std::size_t enc = encode(x);
        _Float64 coef = 2.0f / (sigma * sigma);
        std::vector<_Float64> L_in(n);
        for (unsigned t = 0; t < n; ++t) {
            L_in[t] = coef * ((((enc >> t) & 1) ? -1 : 1) + norm(gen));
        }
        std::vector<double> soft_dec = decode_soft(L_in);
		std::size_t dec = 0;
		for (unsigned t = 0; t < n; ++t) {
			if (soft_dec[t] < 0) {
				dec |= (1ull << t);
			}
			if ((soft_dec[t] < 0) != ((enc >> t) & 1)) {
				++berr;
			}
		}
        errr += (dec != enc);
        ++sim_cnt;
    }
    return {(errr + 0.0) / sim_cnt, (berr + 0.0) / (sim_cnt * n)};
}

unsigned plotkin_construction_decoder::length() const {
    return n;
}

unsigned plotkin_construction_decoder::dim() const {
    return k;
}

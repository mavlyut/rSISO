#include "../include/decoder.h"
#include "../include/defines.h"

// static double truncate(double const& a) {
//     if (a > MAX_TRUNC) {
//         return MAX_TRUNC;
//     }
//     if (a < -MAX_TRUNC) {
//         return -MAX_TRUNC;
//     }
//     return a;
// }

namespace short_domain {
    soft_decoder::soft_decoder(unsigned n, unsigned k) : n(n), k(k) {}

    soft_decoder::~soft_decoder() = default;

    binvector soft_decoder::encode(binvector const&) {
        throw std::runtime_error("Method encode is not implemented");
    }

    std::vector<double> soft_decoder::decode_soft(std::vector<double> const&) {
        throw std::runtime_error("Method decode is not implemented");
    }

    binvector soft_decoder::decode(std::vector<double> const& L_in) {
        std::vector<double> L_out = decode_soft(L_in);
        binvector ans = get_empty(n);
        for (unsigned i = 0; i < n; ++i) {
            setbit(ans, i, L_out[i] < 0);
        }
        return ans;
    }

    std::pair<double, double> soft_decoder::simulate(_Float64 snr, unsigned iter_cnt, unsigned max_error) {
        _Float64 sigma = sqrt(0.5 * pow(10.0, -snr / 10.0) * length() / dim());
        std::normal_distribution<_Float64> norm(0.0, sigma);
        unsigned errr = 0, berr = 0, sim_cnt = 0;
        while (errr < max_error && sim_cnt < iter_cnt) {
            binvector x = get_empty(k);
            for (unsigned t = 0; t < dim(); ++t) {
                setbit(x, t, rand() % 2);
            }
            _log_bv("Simulate, x=", k, x); __log(" '" << x << "'" << std::endl);
            binvector enc = encode(x);
            _log_bv("Simulate, enc=", n, enc); __log(std::endl);
            _Float64 coef = 2.0f / (sigma * sigma);
            std::vector<_Float64> L_in(length());
            for (unsigned t = 0; t < length(); ++t) {
                L_in[t] = coef * ((getbit(enc, t) ? -1 : 1) + norm(gen));
            }
            binvector dec = decode(L_in);
            errr += (dec != enc);
            berr += wt(dec ^ enc);
            ++sim_cnt;
        }
        return {(errr + 0.0) / sim_cnt, (berr + 0.0) / (sim_cnt * n)};
    }

    matrix soft_decoder::generate_matrix() const {
        throw std::runtime_error("Method generate_matrix is not implemented");
    }

    unsigned soft_decoder::length() const {
        return n;
    }

    unsigned soft_decoder::dim() const {
        return k;
    }
}

namespace long_domain {
    soft_decoder::soft_decoder(unsigned n, unsigned k) : n(n), k(k) {}

    soft_decoder::~soft_decoder() = default;

    binvector soft_decoder::encode(binvector const&) {
        throw std::runtime_error("Method encode is not implemented");
    }

    std::vector<double> soft_decoder::decode_soft(std::vector<double> const&) {
        throw std::runtime_error("Method decode is not implemented");
    }

    binvector soft_decoder::decode(std::vector<double> const& L_in) {
        std::vector<double> L_out = decode_soft(L_in);
        binvector ans = get_empty(n);
        for (unsigned i = 0; i < n; ++i) {
            setbit(ans, i, L_out[i] < 0);
        }
        return ans;
    }

    std::pair<double, double> soft_decoder::simulate(_Float64 snr, unsigned iter_cnt, unsigned max_error) {
        _Float64 sigma = sqrt(0.5 * pow(10.0, -snr / 10.0) * length() / dim());
        std::normal_distribution<_Float64> norm(0.0, sigma);
        unsigned errr = 0, berr = 0, sim_cnt = 0;
        while (errr < max_error && sim_cnt < iter_cnt) {
            binvector x = get_empty(k);
            for (unsigned t = 0; t < dim(); ++t) {
                setbit(x, t, rand() % 2);
            }
            binvector enc = encode(x);
            _Float64 coef = 2.0f / (sigma * sigma);
            std::vector<_Float64> L_in(length());
            for (unsigned t = 0; t < length(); ++t) {
                L_in[t] = coef * ((getbit(enc, t) ? -1 : 1) + norm(gen));
            }
            binvector dec = decode(L_in);
            errr += (dec != enc);
            berr += wt(dec ^ enc);
            ++sim_cnt;
        }
        return {(errr + 0.0) / sim_cnt, (berr + 0.0) / (sim_cnt * n)};
    }

    unsigned soft_decoder::length() const {
        return n;
    }

    unsigned soft_decoder::dim() const {
        return k;
    }
}

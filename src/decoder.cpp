#include "../include/decoder.h"
#include "../include/defines.h"

decoder::decoder(unsigned n, unsigned k) : n(n), k(k) {}

decoder::~decoder() = default;

binvector decoder::encode(binvector const&) const {
    throw std::runtime_error("Method encode is not implemented");
}

binvector decoder::decode(std::vector<double> const&) {
    throw std::runtime_error("Method decode is not implemented");
}

std::pair<double, double> decoder::simulate(_Float64 snr, unsigned iter_cnt, unsigned max_error) {
    _Float64 sigma = sqrt(0.5 * pow(10.0, -snr / 10.0) * length() / dim());
    std::normal_distribution<_Float64> norm(0.0, sigma);
    unsigned errr = 0, berr = 0, sim_cnt = 0;
    while (errr < max_error && sim_cnt < iter_cnt) {
        binvector x(dim());
        for (unsigned t = 0; t < dim(); t++) {
            x.set(t, rand() % 2);
        }
        binvector enc = encode(x);
        _Float64 coef = 2.0f / (sigma * sigma);
        std::vector<_Float64> L_in(length());
        for (unsigned t = 0; t < length(); t++) {
            L_in[t] = coef * ((enc[t] ? -1 : 1) + norm(gen));
        }
        binvector dec = decode(L_in);
        errr += (dec != enc);
        berr += (dec ^ enc).wt();
        sim_cnt++;
    }
    return {(errr + 0.0) / sim_cnt, (berr + 0.0) / (sim_cnt * n)};
}

unsigned decoder::length() const {
    return n;
}

unsigned decoder::dim() const {
    return k;
}

soft_decoder::soft_decoder(unsigned n, unsigned k) : decoder(n, k) {}

soft_decoder::~soft_decoder() = default;

binvector soft_decoder::decode(std::vector<double> const& L_in) {
    std::vector<double> L_out = decode_soft(L_in);
    binvector ans(length());
    for (unsigned i = 0; i < ans.size(); i++) {
        ans.set(i, L_out[i] < 0);
    }
    return ans;
}

const double soft_decoder::MAX_TRUNC = 12;
double soft_decoder::truncate(double const &a) {
    if (a > MAX_TRUNC) {
        return MAX_TRUNC;
    }
    if (a < -MAX_TRUNC) {
        return -MAX_TRUNC;
    }
    return a;
}

linear_soft_decoder::linear_soft_decoder(matrix const& G) : soft_decoder(G.front().size(), G.size()), G(G) {}

linear_soft_decoder::~linear_soft_decoder() = default;

std::vector<double> soft_decoder::decode_soft(std::vector<double> const&) {
    throw std::runtime_error("Method decode_soft is not implemented");
}

binvector linear_soft_decoder::encode(binvector const& c) const {
    fail(c.size() == dim(), "encode: incorrect input dim");
    binvector ans(length());
    for (unsigned i = 0; i < dim(); i++) {
        if (c[i]) {
            ans ^= G[i];
        }
    }
    return ans;
}

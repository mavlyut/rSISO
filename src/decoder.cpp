#include "../include/decoder.h"
#include "../include/defines.h"

binvector decoder::encode(binvector const&) const {
    throw std::runtime_error("Not implemented");
}

binvector decoder::decode(std::vector<double> const&) {
    throw std::runtime_error("Not implemented");
}

std::pair<unsigned, unsigned> decoder::simulate(_Float64 snr, unsigned iter_cnt, unsigned max_error) {
    _Float64 sigma = sqrt(0.5 * pow(10.0, -snr / 10.0) * length() / dim());
    std::normal_distribution<_Float64> norm(0.0, sigma);
    int errr = 0, sim_cnt = 0;
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
        sim_cnt++;
    }
    return {errr, sim_cnt};
}

unsigned decoder::length() const {
    throw std::runtime_error("Not implemented");
}

unsigned decoder::dim() const {
    throw std::runtime_error("Not implemented");
}

linear_decoder::linear_decoder(matrix const& G) : n(G.front().size()), k(G.size()), G(G) {}

binvector linear_decoder::encode(binvector const& c) const {
    fail(c.size() == dim(), "encode: incorrect input dim");
    binvector ans(length());
    for (unsigned i = 0; i < dim(); i++) {
        if (c[i]) {
            ans ^= G[i];
        }
    }
    return ans;
}

unsigned linear_decoder::length() const {
    return n;
}

unsigned linear_decoder::dim() const {
    return k;
}

linear_soft_decoder::linear_soft_decoder(matrix const& G) : linear_decoder(G) {}

std::vector<double> soft_decoder::decode_soft(std::vector<double> const&) {
    throw std::runtime_error("Not implemented");
}

binvector soft_decoder::decode(std::vector<double> const& L_in) {
    std::vector<double> L_out = decode_soft(L_in);
    binvector ans(length());
    for (unsigned i = 0; i < ans.size(); i++) {
        ans.set(i, L_out[i] < 0);
    }
    return ans;
}

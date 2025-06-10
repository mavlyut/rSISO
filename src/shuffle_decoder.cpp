#include "../include/shuffle_decoder.h"

namespace short_domain {
    bool is_permutation(std::vector<unsigned> const& p) {
        unsigned n = p.size();
        std::vector<unsigned> q(n, UNINIT);
        for (unsigned i = 0; i < n; i++) {
            if (q[p[i]] != UNINIT) {
                return false;
            }
            q[p[i]] = i;
        }
        return true;
    }

    shuffle_decoder::shuffle_decoder(std::vector<unsigned> const& perm, soft_decoder* inner)
            : soft_decoder(inner->length(), inner->dim()), perm(perm), inner_decoder(inner)
            , encoded_perm(get_empty(n)), input_perm(n), output_perm(n) {   
        fail(perm.size() == n, "shuffle_dec, ctor: incorrect size of permutation");
        for (unsigned i = 0; i < n; ++i) {
            fail(0 <= perm[i] && perm[i] < n, "shuffle_dec, ctor: perm[i] not in [0,n)");
            fail(is_permutation(perm), "shuffle_dec: is not permutation");
        }
    }

    binvector shuffle_decoder::encode(binvector const& x) {
        _log_bv("Shuffle, encode: ", k, x); __log(std::endl);
        auto enc = inner_decoder->encode(x);
        for (unsigned i = 0; i < n; ++i) {
            setbit(encoded_perm, i, getbit(enc, perm[i]));
        }
        _log_bv("Shuffle, enc, before perm: ", n, enc); __log(std::endl);
        _log_bv("Shuffle, encoded: ", n, encoded_perm); __log(std::endl);
        return encoded_perm;
    }

    std::vector<double> shuffle_decoder::decode_soft(std::vector<double> const& input) {
        __log("Shuffle, decode: " << input << std::endl);
        for (unsigned i = 0; i < n; ++i) {
            input_perm[perm[i]] = input[i];
        }
        auto output = inner_decoder->decode_soft(input_perm);
        for (unsigned i = 0; i < n; ++i) {
            output_perm[i] = output[perm[i]];
        }
        __log("Shuffle, decoded: " << output_perm << std::endl);
        return output_perm;
    }

    matrix shuffle_decoder::generate_matrix() const {
        matrix G = inner_decoder->generate_matrix(), Ans(k, get_empty(n));
        // printmatrix(std::cout, "inner gen:", n, G); std::cout << std::endl;
        for (unsigned i = 0; i < k; ++i) {
            for (unsigned j = 0; j < n; ++j) {
                setbit(Ans[i], j, getbit(G[i], perm[j]));
            }
        }
        return Ans;
    }
}

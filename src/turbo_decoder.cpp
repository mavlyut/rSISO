#include "../include/turbo_decoder.h"
#include "../include/recursive_decoder.h"

namespace long_domain {
    turbo_decoder::turbo_decoder(short_domain::soft_decoder* dec_row, short_domain::soft_decoder* dec_col)
        : soft_decoder(dec_row->length() * dec_col->length(), dec_row->dim() * dec_col->dim())
        , n1(dec_row->length()), n2(dec_col->length()), k1(dec_row->dim()), k2(dec_col->dim())
        , dec_row(dec_row), dec_col(dec_col)
        , ans(n1 * n2)
        , enc1t(n1, binvector(k2))
        , dec2t(n2, std::vector<double>(n1))
        , colL(n2)
        , L_out(n1 * n2)
    {
        __log("Turbo.ctor" << std::endl);
    }

    binvector turbo_decoder::encode(binvector const& c) {
        fail(c.size() == dim(), "turbo, encode: incorrect dim");
        __log("Turbo.encode: " << c << std::endl);
        for (unsigned i = 0; i < k2; ++i) {
            auto row = dec_row->encode(c.subvector(k1 * i, k1 * (i + 1)));
            for (unsigned j = 0; j < n1; ++j) {
                enc1t[j].set(i, ((row >> j) & 1));
            }
        }
        __log("Turbo.row_encoded:\n" << enc1t);
        for (unsigned j = 0; j < n1; ++j) {
            auto row = dec_col->encode(enc1t[j]);
            for (unsigned i = 0; i < n2; ++i) {
                ans.set(i * n1 + j, ((row >> i) & 1));
            }
        }
        __log("Turbo.encoded: " << ans << std::endl);
        return ans;
    }

    // C=0.35
    const unsigned turbo_decoder::ITER_CNT = 4;
    std::vector<double> turbo_decoder::decode_soft(std::vector<double> const& L0) {
        fail(L0.size() == length(), "turbo, decode: incorrect dim");
        __log("Turbo.decode: " << L0 << std::endl);
        L_out = L0;
        for (unsigned cnt = 0; cnt < ITER_CNT; ++cnt) {
            for (unsigned i = 0; i < n1; ++i) {
                for (unsigned j = 0; j < n2; ++j) {
                    colL[j] = L_out[j * n1 + i];
                }
                auto colDec = dec_col->decode_soft(colL);
                for (unsigned j = 0; j < n2; ++j) {
                    dec2t[j][i] = colL[j] + colDec[j];
                    // dec2t[j][i] = colL[j] + 0.35 * colDec[j];
                    // dec2t[j][i] = colL[j] + 0.35 * (colDec[j] - colL[j]);
                    // dec2t[j][i] = colDec[j];
                }
            }
            for (unsigned j = 0; j < n2; ++j) {
                auto row = dec_row->decode_soft(dec2t[j]);
                for (unsigned i = 0; i < n1; ++i) {
                    L_out[j * n1 + i] = dec2t[j][i] + row[i];
                    // L_out[j * n1 + i] = dec2t[j][i] + 0.35 * row[i];
                    // L_out[j * n1 + i] = dec2t[j][i] + 0.35 * (row[i] - dec2t[j][i]);
                    // L_out[j * n1 + i] = row[i];
                }
            }
            __log("Turbo.decoded, iter=" << cnt + 1 << ": " << L_out << std::endl);
        }
        __log("Turbo.decoded: " << L_out << std::endl);
        return L_out;
    }
}

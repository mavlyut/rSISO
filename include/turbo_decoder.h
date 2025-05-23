#ifndef TURBO_DECODER_H
#define TURBO_DECODER_H

#include "decoder.h"
#include "recursive_decoder.h"

namespace long_domain {
    struct turbo_decoder : public soft_decoder {
        turbo_decoder(short_domain::soft_decoder*, short_domain::soft_decoder*);

        binvector encode(binvector const& c) override;
        std::vector<double> decode_soft(std::vector<double> const&) override;

    private:
        static const unsigned ITER_CNT;

        unsigned n1, n2, k1, k2;
        short_domain::soft_decoder *dec_row, *dec_col;

        binvector ans;
        matrix enc1t;

        std::vector<double> colL;
        std::vector<std::vector<double>> dec2t;
        std::vector<double> L_out;
    };
}

#endif // TURBO_DECODER_H

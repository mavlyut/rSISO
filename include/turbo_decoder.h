#ifndef TURBO_DECODER_H
#define TURBO_DECODER_H

#include "decoder.h"

struct turbo_decoder : public soft_decoder {
    turbo_decoder(soft_decoder*, soft_decoder*);

    binvector encode(binvector const& c) const override;
    std::vector<double> decode_soft(std::vector<double> const&) override;

private:
    unsigned n1, n2, k1, k2;
    soft_decoder *dec_row, *dec_col;
    static const unsigned ITER_CNT = 10;
};

#endif // TURBO_DECODER_H

#ifndef TURBO_DECODER_H
#define TURBO_DECODER_H

#include "decoder.h"
#include "recursive_decoder.h"

struct turbo_decoder : public long_domain::soft_decoder {
    turbo_decoder(short_domain::soft_decoder*, short_domain::soft_decoder*);

    long_domain::binvector encode(long_domain::binvector const& c) override;
    std::vector<double> decode_soft(std::vector<double> const&) override;

private:
    static const unsigned ITER_CNT;

    unsigned n1, n2, k1, k2;
    short_domain::soft_decoder *dec_row, *dec_col;

    long_domain::binvector ans;
    short_domain::matrix enc1t;

    std::vector<double> colL;
    std::vector<std::vector<double>> dec2t;
    std::vector<double> L_out;
};

#endif // TURBO_DECODER_H

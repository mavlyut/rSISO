#ifndef TURBO_DECODER_H
#define TURBO_DECODER_H

#include "decoder.h"

class turbo_decoder : public soft_decoder {
public:
    turbo_decoder(soft_decoder const&, soft_decoder const&);

    binvector encode(binvector const&) const override;
    std::vector<double> soft_decode(std::vector<double> const&);

    unsigned length() const override;
    unsigned dim() const override;

private:
    static const unsigned ITER_CNT = 10;
    unsigned n1, n2, k1, k2;
    soft_decoder dec_row, dec_col;
};

#endif // TURBO_DECODER_H

#ifndef PLOTKIN_CONSTRUCTION_DECODER_H
#define PLOTKIN_CONSTRUCTION_DECODER_H

#include "decoder.h"

struct plotkin_construction_decoder : public soft_decoder {
    plotkin_construction_decoder(soft_decoder* dec1, soft_decoder* dec2);

    binvector encode(binvector const&) override;
    std::vector<double> decode_soft(std::vector<double> const&) override;

private:
    unsigned m;
    soft_decoder *dec1, *dec2;
    std::vector<double> Lq0, Lq1, Lq2, L_out;
    std::vector<std::array<double, 3>> Lr;

    static const unsigned MAX_ITER_COUNT = 5;
    static inline double phi(double);
};

#endif // PLOTKIN_CONSTRUCTION_DECODER_H

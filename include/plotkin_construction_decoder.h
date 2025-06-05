#ifndef PLOTKIN_CONSTRUCTION_DECODER_H
#define PLOTKIN_CONSTRUCTION_DECODER_H

#include "decoder.h"
#include "recursive_decoder.h"

namespace short_domain {
    struct plotkin_construction_decoder : soft_decoder {
        plotkin_construction_decoder(soft_decoder* dec1, soft_decoder* dec2);

        binvector encode(binvector const& x) override;
        std::vector<double> decode_soft(std::vector<double> const&) override;
		matrix generate_matrix() const override;

    private:
        static const unsigned MAX_ITER_COUNT;

        unsigned m;
        soft_decoder *dec1, *dec2;

        std::vector<double> Lq0, Lq1, Lq2, L_out;
        std::vector<double> Lr0, Lr1, Lr2;
        std::vector<double> L_ext_0, L_ext_2;

        static inline double phi(double);
    };
}

#endif // PLOTKIN_CONSTRUCTION_DECODER_H

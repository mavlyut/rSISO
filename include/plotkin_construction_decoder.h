#ifndef PLOTKIN_CONSTRUCTION_DECODER_H
#define PLOTKIN_CONSTRUCTION_DECODER_H

#include "decoder.h"
#include "recursive_decoder.h"

struct plotkin_construction_decoder {
    plotkin_construction_decoder(recursive_decoder* dec1, recursive_decoder* dec2);

    std::size_t encode(std::size_t const&) const;
    std::vector<double> decode_soft(std::vector<double> const&);
	std::pair<double, double> simulate(double snr, unsigned iter_cnt, unsigned max_error);

    unsigned length() const;
    unsigned dim() const;

    void printGraph(unsigned iter_num, unsigned state);

private:
    unsigned n, m, k;
    recursive_decoder *dec1, *dec2;
    std::vector<double> Lq0, Lq1, Lq2, L_out;
    std::vector<std::array<double, 3>> Lr;

    static const unsigned MAX_ITER_COUNT;
    static inline double phi(double);
};

#endif // PLOTKIN_CONSTRUCTION_DECODER_H

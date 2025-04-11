#ifndef TRUST_PROPAGATION_DECODER_H
#define TRUST_PROPAGATION_DECODER_H

#include "decoder.h"

struct trust_propagation_decoder : public linear_soft_decoder {
    trust_propagation_decoder(matrix const& H);

    std::vector<double> decode_soft(std::vector<double> const&) override;
    void print_graph(std::ofstream, binvector const&) const;

private:
    unsigned m;
    static const unsigned MAX_ITER_COUNT = 50;
    std::vector<std::vector<unsigned>> C, V;

    static inline double phi(double);
    std::ostream& print(std::ostream& out, std::vector<std::map<int, double>> const& a) const;
};

#endif // TRUST_PROPAGATION_DECODER_H

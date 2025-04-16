#ifndef DECODER_H
#define DECODER_H

#include "matrix.h"

struct decoder {
    decoder(unsigned n, unsigned k);
    virtual ~decoder();

    virtual binvector encode(binvector const&);
    virtual binvector decode(std::vector<double> const&);
    virtual std::pair<double, double> simulate(_Float64, unsigned, unsigned);

    unsigned length() const;
    unsigned dim() const;

protected:
    unsigned n, k;
};

struct soft_decoder : public decoder {
    soft_decoder(unsigned n, unsigned k);
    virtual ~soft_decoder();

    virtual std::vector<double> decode_soft(std::vector<double> const&);
    binvector decode(std::vector<double> const&) override;

protected:
    static const double MAX_TRUNC;
    static double truncate(double const&);
};

struct linear_soft_decoder : public soft_decoder {
    linear_soft_decoder(matrix const&);
    virtual ~linear_soft_decoder();

    binvector encode(binvector const&) override;

protected:
    matrix const G;
};

#endif // DECODER_H

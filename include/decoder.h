#ifndef DECODER_H
#define DECODER_H

#include "matrix.h"

class decoder {
public:
    virtual ~decoder() = default;

    virtual binvector encode(binvector const&) const;
    virtual binvector decode(std::vector<double> const&);
    std::pair<unsigned, unsigned> simulate(_Float64, unsigned, unsigned);

    virtual unsigned length() const;
    virtual unsigned dim() const;
};

class linear_decoder : public virtual decoder {
public:
    virtual ~linear_decoder() = default;

    linear_decoder(matrix const&);

    binvector encode(binvector const&) const override;
    unsigned length() const override;
    unsigned dim() const override;

protected:
    unsigned n, k;
    matrix const G;
};

class soft_decoder : public virtual decoder {
public:
    virtual ~soft_decoder() = default;
 
    virtual std::vector<double> decode_soft(std::vector<double> const&);
    binvector decode(std::vector<double> const&) override;
};

class linear_soft_decoder : public linear_decoder, public soft_decoder {
public:
    linear_soft_decoder(matrix const&);
};

#endif // DECODER_H

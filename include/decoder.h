#include "matrix.h"

class decoder {
public:
    decoder(matrix const&);

    binvector encode(binvector const&) const;
    virtual binvector decode(std::vector<double> const&) = 0;
    std::pair<unsigned, unsigned> simulate(_Float64, unsigned, unsigned);

    unsigned length() const;
    unsigned dim() const;

protected:
    unsigned n, k;
    matrix const G;
};

class soft_decoder : public decoder {
public:
    soft_decoder(matrix const&);
    
    virtual std::vector<double> decode_soft(std::vector<double> const&) = 0;
    binvector decode(std::vector<double> const&) override;
};

#ifndef DECODER_H
#define DECODER_H

#include "matrix.h"

#define MAX_TRUNC 16

static double truncate(double const& a) {
    if (a > MAX_TRUNC) {
        return MAX_TRUNC;
    }
    if (a < -MAX_TRUNC) {
        return -MAX_TRUNC;
    }
    return a;
}

namespace short_domain {
    struct soft_decoder {
        soft_decoder(unsigned n, unsigned k);
        virtual ~soft_decoder();

        virtual binvector encode(binvector const&);

        virtual binvector decode(std::vector<double> const&);
        virtual std::vector<double> decode_soft(std::vector<double> const&);

        std::pair<double, double> simulate(_Float64 snr, unsigned iter_cnt, unsigned max_error);

        virtual matrix generate_matrix() const;

        unsigned length() const;
        unsigned dim() const;

    protected:
        unsigned n, k;
    };
}

namespace long_domain {
    struct soft_decoder {
        soft_decoder(unsigned n, unsigned k);
        virtual ~soft_decoder();

        virtual binvector encode(binvector const&);
        virtual binvector decode(std::vector<double> const&);
        virtual std::pair<double, double> simulate(_Float64 snr, unsigned iter_cnt, unsigned max_error);
        virtual std::vector<double> decode_soft(std::vector<double> const&);

        unsigned length() const;
        unsigned dim() const;

    private:
        unsigned n, k;
    };
}

#endif // DECODER_H

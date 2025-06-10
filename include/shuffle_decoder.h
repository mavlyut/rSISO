#ifndef SHUFFLE_DECODER
#define SHUFFLE_DECODER

#include "decoder.h"

namespace short_domain {
	struct shuffle_decoder : soft_decoder {
		shuffle_decoder(std::vector<unsigned> const& perm, soft_decoder*);

		binvector encode(binvector const& c) override;
		std::vector<double> decode_soft(std::vector<double> const&) override;

		matrix generate_matrix() const override;

	private:
		std::vector<unsigned> perm;
		soft_decoder* inner_decoder;

		binvector encoded_perm;
		std::vector<double> input_perm, output_perm;
	};
}

#endif // SHUFFLE_DECODER

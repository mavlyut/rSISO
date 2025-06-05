#include "include/recursive_decoder.h"
#include "include/plotkin_construction_decoder.h"
#include "include/turbo_decoder.h"

unsigned iter_cnt = 10000, max_error = 100;

namespace short_domain {
    matrix generate_RM(unsigned r, unsigned m) {
        if (r == 0) {
            matrix G(1, get_empty(1ull << m));
            for (unsigned i = 0; i < (1ull << m); i++) {
                setbit(G[0], i, true);
            }
            return G;
        }
        if (r == m) {
            matrix G(1ull << m, get_empty(1ull << m));
            for (unsigned i = 0; i < (1ull << m); i++) {
                setbit(G[i], i, true);
            }
            return G;
        }
        matrix G1 = generate_RM(r, m - 1), G2 = generate_RM(r - 1, m - 1);
        unsigned sh = (1ull << (m - 1));
        matrix G;
        for (unsigned i = 0; i < G1.size(); i++) {
            binvector row = get_empty(1ull << m);
            for (unsigned j = 0; j < sh; j++) {
                setbit(row, j, getbit(G1[i], j));
                setbit(row, j + sh, getbit(G1[i], j));
            }
            G.push_back(row);
        }
        for (unsigned i = 0; i < G2.size(); i++) {
            binvector row = get_empty(1ull << m);
            for (unsigned j = 0; j < sh; j++) {
                setbit(row, j + sh, getbit(G2[i], j));
            }
            G.push_back(row);
        }
        return G;
    }
}

int main() {
	srand(time(NULL));

#if defined(LOCAL)
	std::ostream& fout = std::cout;
#else
	std::ofstream fout("output.txt");
#endif
	// fout << std::scientific;

    short_domain::matrix G = short_domain::generate_RM(2, 4);
    short_domain::printmatrix(std::cerr, 16, G);

	unsigned n = 16, k = G.size();
	auto* rsiso_dec = new short_domain::recursive_decoder(n, G);

    auto* c1 = new short_domain::recursive_decoder(n / 2, short_domain::generate_RM(2, 3));
    auto* c2 = new short_domain::recursive_decoder(n / 2, short_domain::generate_RM(1, 3));
    auto* pt_dec = new short_domain::plotkin_construction_decoder(c1, c2);

    turbo_decoder turbo_dec(rsiso_dec, pt_dec);

    for (double snr = 0.0; snr <= 5; snr += 0.5) {
		fail(iter_cnt > 0, "simulate: iter_cnt must be positive");
		fail(max_error >= 0, "simulate: max_err_cnt must be non-negative");
		_Float64 sigma = sqrt(0.5 * pow(10.0, -snr / 10.0) * n / k);
        std::normal_distribution<_Float64> norm(0.0, sigma);

		auto start = std::chrono::system_clock::now();
        auto [fer, ber] = turbo_dec.simulate(snr, iter_cnt, max_error);
		auto end = std::chrono::system_clock::now();
		auto time_in_ms = std::chrono::duration_cast<ms>(end - start).count();

		#ifdef TSV_FORMAT
			fout << std::to_string(fer).replace(1, 1, ",");
            std::cout << "snr = " << snr << "; time = " << time_in_ms / 1000.0 << "s" << std::endl;
		#else
			fout << "Snr = " << snr << "\n";
			fout << "\tPt   : FER = " << fer << "; BER = " << ber << "; time = " << time_in_ms / 1000.0 << "s\n";
		#endif
		fout << std::endl;
	}

	std::cout << "FINISH\n";
}

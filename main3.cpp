#include <string>
#include "include/recursive_decoder.h"
#include "include/turbo_decoder.h"
#include "include/plotkin_construction_decoder.h"

#define printbv(msg, n, v) 					\
	__log(msg);								\
	for (unsigned j = 0; j < n; ++j) { 		\
		if (j != 0) __log(" ");				\
		__log((((v >> j) & 1) ? 1 : 0));	\
	}
#define printmatrix(msg, n, M)					\
	__log(msg << "\n");							\
	for (std::size_t const& i : M) { 			\
		printbv("", n, i);						\
		__log(std::endl);						\
	}											\
	__log(std::endl);

unsigned iter_cnt = 10000, max_error = 100;

using namespace short_domain;

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

int main() {
	srand(time(NULL));

	std::ofstream fout("output.txt");

    unsigned R = 2, M = 5;
    matrix G1 = generate_RM(R, M), G2 = generate_RM(R - 1, M);
	unsigned m = (1ull << M), k1 = G1.size(), k2 = G2.size();
    std::cout << "Generated: (" << m << "," << k1 << ") and (" << m << "," << k2 << ")" << std::endl;
	unsigned n = 2 * m, k = k1 + k2;
	matrix G(k, 0);
	for (unsigned i = 0; i < k1; i++) {
		G[i] = (G1[i] | (G1[i] << m));
	}
	for (unsigned i = 0; i < k2; i++) {
		G[i + k1] = (G2[i] << m);
	}
	recursive_decoder* rd1 = new recursive_decoder(m, G1);
	recursive_decoder* rd2 = new recursive_decoder(m, G2);

	plotkin_construction_decoder pt_coder(rd1, rd2);	
	recursive_decoder rSISO_coder(n, G);

	fout << "RM (" << n << "," << k << "), rSISO\tRM (" << n << "," << k << "), Plotkin" << std::endl;

	for (double snr = 0.0; snr <= 5; snr += 0.5) {
		fail(iter_cnt > 0, "simulate: iter_cnt must be positive");
		fail(max_error >= 0, "simulate: max_err_cnt must be non-negative");
		_Float64 sigma = sqrt(0.5 * pow(10.0, -snr / 10.0) * n / k);
        std::normal_distribution<_Float64> norm(0.0, sigma);

		auto start = std::chrono::system_clock::now();

        auto [fer_pt, ber_pt] = pt_coder.simulate(snr, iter_cnt, max_error);

		auto end = std::chrono::system_clock::now();
		auto time_in_ms_pt = std::chrono::duration_cast<ms>(end - start).count();


		start = std::chrono::system_clock::now();

        auto [fer_rec, ber_rec] = std::pair{0.0, 0.0};//rSISO_coder.simulate(snr, iter_cnt, max_error);

		end = std::chrono::system_clock::now();
		auto time_in_ms_rec = std::chrono::duration_cast<ms>(end - start).count();

		#ifdef TSV_FORMAT
			fout
                << "0," << std::to_string(fer_rec).substr(2) << "\t"
				<< "0," << std::to_string(fer_pt).substr(2);
            std::cout << "snr = " << snr << "; time(plotkin) = " << time_in_ms_pt / 1000.0 << "s; time(rSISO) = " << time_in_ms_rec / 1000.0 << "s" << std::endl;
		#else
			fout << "Snr = " << snr << "\n";
			fout << "\trSISO: FER = " << fer_rec << "; BER = " << ber_rec << "; time = " << time_in_ms_rec / 1000.0 << "s\n";
			fout << "\tPt   : FER = " << fer_pt << "; BER = " << ber_pt << "; time = " << time_in_ms_pt / 1000.0 << "s\n";
		#endif
		fout << std::endl;
	}

	std::cout << "FINISH\n";
}

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

unsigned iter_cnt = 100000, max_error = 100;

int main() {
	srand(time(NULL));

	std::ifstream fin("input.txt");
	std::ofstream fout("output.txt");

	unsigned m, k1, k2;
	fin >> m >> k1 >> k2;
	unsigned n = 2 * m, k = k1 + k2;
	std::vector<std::size_t> G1(k1, 0), G2(k2, 0), G(k, 0);
	for (unsigned i = 0; i < k1; i++) {
		for (unsigned j = 0; j < m; j++) {
			bool b;
			fin >> b;
			if (b) {
				G1[i] ^= (1ull << j);
			}
		}
		G[i] = (G1[i] | (G1[i] << m));
	}
	for (unsigned i = 0; i < k2; i++) {
		for (unsigned j = 0; j < m; j++) {
			bool b;
			fin >> b;
			if (b) {
				G2[i] ^= (1ull << j);
			}
		}
		G[i + k1] = (G2[i] << m);
	}
	recursive_decoder* rd1 = new recursive_decoder(m, G1);
	recursive_decoder* rd2 = new recursive_decoder(m, G2);

	plotkin_construction_decoder pt_coder(rd1, rd2);	
	recursive_decoder rSISO_coder(n, G);

	fout << "RM (" << n << "," << k << ", rSISO\tRM (" << n << "," << k << ", Plotkin\n";

	for (double snr = 0.0; snr <= 5; snr += 0.5) {
		fail(iter_cnt > 0, "simulate: iter_cnt must be positive");
		fail(max_error >= 0, "simulate: max_err_cnt must be non-negative");

		auto start = std::chrono::system_clock::now();

		_Float64 sigma = sqrt(0.5 * pow(10.0, -snr / 10.0) * n / k);
			std::normal_distribution<_Float64> norm(0.0, sigma);
			unsigned errr_pt = 0, berr_pt = 0;
			unsigned errr_rec = 0, berr_rec = 0;
			unsigned sim_cnt = 0;
			while ((errr_pt < max_error || errr_rec < max_error) && sim_cnt < iter_cnt) {
				std::size_t x = 0;
				for (unsigned t = 0; t < k; ++t) {
					if (rand() % 2) {
						x |= (1ull << t);
					}
				}
				std::size_t enc = pt_coder.encode(x);
				fail(enc == rSISO_coder.encode(x), "Main2: encoded vectors aren't equals");
				_Float64 coef = 2.0f / (sigma * sigma);
				std::vector<_Float64> L_in(n);
				for (unsigned t = 0; t < n; ++t) {
					L_in[t] = coef * ((((enc >> t) & 1) ? -1 : 1) + norm(gen));
				}

				{
					std::vector<double> soft_dec = pt_coder.decode_soft(L_in);
					std::size_t dec = 0;
					for (unsigned t = 0; t < n; ++t) {
						if (soft_dec[t] < 0) {
							dec |= (1ull << t);
						}
						if ((soft_dec[t] < 0) != ((enc >> t) & 1)) {
							++berr_pt;
						}
					}
					errr_pt += (dec != enc);
				}

				{
					std::vector<double> soft_dec = rSISO_coder.decode_soft(L_in);
					std::size_t dec = 0;
					for (unsigned t = 0; t < n; ++t) {
						if (soft_dec[t] < 0) {
							dec |= (1ull << t);
						}
						if ((soft_dec[t] < 0) != ((enc >> t) & 1)) {
							++berr_rec;
						}
					}
					errr_rec += (dec != enc);
				}

				++sim_cnt;
			}

		auto end = std::chrono::system_clock::now();
		auto time_in_ms = std::chrono::duration_cast<ms>(end - start).count();

		#ifdef TSV_FORMAT
			fout << "0," << std::to_string((errr_rec + 0.0) / sim_cnt).substr(2) << "\t"
				 << "0," << std::to_string((errr_pt + 0.0) / sim_cnt).substr(2); 
		#else
			fout << "Snr = " << snr << "\n";
			fout << "\trSISO: FER = " << (errr_rec + 0.0) / sim_cnt << "; BER = " << (berr_rec + 0.0) / (sim_cnt * n) << "\n";
			fout << "\tPt   : FER = " << (errr_pt + 0.0) / sim_cnt << "; BER = " << (berr_pt + 0.0) / (sim_cnt * n) << "\n";
			fout << "\tTime = " << time_in_ms / 1000.0 << "s";
		#endif
		fout << std::endl;
	}

	std::cout << "FINISH\n";
}

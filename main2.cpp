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

int main() {
	srand(time(NULL));

	std::ifstream fin("input.txt");
	std::ofstream fout("output.txt");

	unsigned m, k1, k2;
	fin >> m >> k1 >> k2;
	unsigned n = 2 * m;
	std::vector<std::size_t> G1(k1, 0), G2(k2, 0), G(k1 + k2, 0);
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
		G[i + k1] = G2[i];
	}
	recursive_decoder* rd1 = new recursive_decoder(m, G1);
	recursive_decoder* rd2 = new recursive_decoder(m, G2);
	plotkin_construction_decoder coder(rd1, rd2);
	
	recursive_decoder coder_correct(n, G);

	printmatrix("G32: ", n, G);

	std::string command;
	bool skip = false;
	while (fin >> command) {
		if (command[0] == '-' && command[1] == '-') {
			std::getline(fin, command);
			continue;
		} else if (command == "{-") {
			skip = true;
			std::getline(fin, command);
			continue;
		} else if (command == "-}") {
			skip = false;
			std::getline(fin, command);
			continue;
		} else if (skip) {
			std::getline(fin, command);
			continue;
		} else if (command == "Encode") {
			binvector x(coder.dim());
			fin >> x;
			fout << coder.encode(x);
		} else if (command == "Decode") {
			std::vector<_Float64> y(coder.length());
			fin >> y;
			fout << coder.decode_soft(y);
		} else if (command == "Simulate") {
			_Float64 snr;
			int iter_cnt, max_error;
			fin >> snr >> iter_cnt >> max_error;
			fail(iter_cnt > 0, "simulate: iter_cnt must be positive");
			fail(max_error >= 0, "simulate: max_err_cnt must be non-negative");

			auto start = std::chrono::system_clock::now();

			_Float64 sigma = sqrt(0.5 * pow(10.0, -snr / 10.0) * coder.length() / coder.dim());
				std::normal_distribution<_Float64> norm(0.0, sigma);
				unsigned errr = 0, berr = 0, sim_cnt = 0;
				while (errr < max_error && sim_cnt < iter_cnt) {
					std::size_t x = 0;
					for (unsigned t = 0; t < coder.dim(); ++t) {
						if (rand() % 2) {
							x |= (1ull << t);
						}
					}
					std::size_t enc = coder.encode(x);
					fail(coder_correct.encode(x) == enc, "encode is incorrect");
					_Float64 coef = 2.0f / (sigma * sigma);
					std::vector<_Float64> L_in(n);
					for (unsigned t = 0; t < n; ++t) {
						L_in[t] = coef * ((((enc >> t) & 1) ? -1 : 1) + norm(gen));
					}
					std::vector<double> soft_dec = coder.decode_soft(L_in);
					std::vector<double> soft_dec_correct = coder_correct.decode_soft(L_in);
					std::size_t dec = 0;
					for (unsigned t = 0; t < n; ++t) {
						if (soft_dec[t] < 0) {
							dec |= (1ull << t);
						}
						if ((soft_dec[t] < 0) != ((enc >> t) & 1)) {
							++berr;
						}
						if ((soft_dec[t] < 0) != (soft_dec_correct[t] < 0)) {
							__log("Main: plotkin decoding error!" << std::endl);
							break;
						}
					}
					errr += (dec != enc);
					++sim_cnt;
				}

			auto end = std::chrono::system_clock::now();
			auto time_in_ms = std::chrono::duration_cast<ms>(end - start).count();

			fout << "Snr = " << snr << "; FER = " << (errr + 0.0) / sim_cnt << "; BER = " << (berr + 0.0) / (sim_cnt * coder.length());
			fout << "; time = " << time_in_ms / 1000.0 << "s";
		} else {
			fail(false, "main: incorrect command");
			continue;
		}
		fout << std::endl;
	}

	std::cout << "FINISH\n";
}

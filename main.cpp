#include "include/recursive_decoder.h"
#include "include/turbo_decoder.h"
#include "include/plotkin_construction_decoder.h"

int main() {
	srand(time(NULL));

	std::ifstream fin("input.txt");
#if defined(LOCAL)
	std::ostream& fout = std::cout;
#else
	std::ofstream fout("output.txt");
#endif
	// fout << std::scientific;

	// for (unsigned N = 1; N <= 32; N++) {
	// 	small_gray_code code1(N);
	// 	gray_code code2(N);
	// 	auto it1 = code1.begin();
	// 	auto it2 = code2.begin();
	// 	bool fl = false;
	// 	for (; it1 != code1.end() && it2 != code2.end(); ++it1, ++it2) {
	// 		if (*it1 != *it2) {
	// 			std::cerr << "Bad: " << N << " " << *it1 << " " << *it2 << std::endl;
	// 			fl = true;
	// 			break;
	// 		}
	// 	}
	// 	if (fl) {
	// 		continue;
	// 	}
	// 	if (!(it1 == code1.end() && it2 == code2.end())) {
	// 		std::cerr << "Bad: " << N << std::endl;
	// 		continue;
	// 	}
	// 	std::cout << "Good: " << N << std::endl;
	// }
	// return 0;

	// unsigned n, k1, k2;
	// fin >> n >> k1 >> k2;
	// std::vector<std::size_t> G1(k1, 0), G2(k2, 0);
	// for (unsigned i = 0; i < k1; i++) {
	// 	for (unsigned j = 0; j < n; j++) {
	// 		bool b;
	// 		fin >> b;
	// 		if (b) {
	// 			G1[i] ^= (1ull << j);
	// 		}
	// 	}
	// }
	// for (unsigned i = 0; i < k2; i++) {
	// 	for (unsigned j = 0; j < n; j++) {
	// 		bool b;
	// 		fin >> b;
	// 		if (b) {
	// 			G2[i] ^= (1ull << j);
	// 		}
	// 	}
	// }

	unsigned n, k;
	fin >> n >> k;
	std::vector<std::size_t> G(k, 0);
	for (unsigned i = 0; i < k; i++) {
		for (unsigned j = 0; j < n; j++) {
			bool b;
			fin >> b;
			if (b) {
				G[i] ^= (1ull << j);
			}
		}
	}

	recursive_decoder coder(n, G);
	// recursive_decoder* rd1 = new recursive_decoder(n, G1);
	// recursive_decoder* rd2 = new recursive_decoder(n, G2);
	// plotkin_construction_decoder coder(rd1, rd2);

	// unsigned n, k;
	// fin >> n >> k;
	// std::vector<std::size_t> G(k, 0);
	// for (unsigned i = 0; i < k; i++) {
	// 	for (unsigned j = 0; j < n; j++) {
	// 		bool b;
	// 		fin >> b;
	// 		if (b) {
	// 			G[i] ^= (1ull << j);
	// 		}
	// 	}
	// }
	// recursive_decoder coder(n, G);

	// #ifdef PRINT_SECTION_TREE
	// coder.printTree("Sectionalization.dot");
	// #endif

	#ifdef INIT_COUNTS
	std::cout << "Init: " << CNT_BIN << "\n";
	CNT_BIN = 0;
	#endif

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

			#ifdef TEST
			auto start = std::chrono::system_clock::now();
			#ifdef CNTLOG
			clear_cnt();
			clear_cnt_bin();
			#endif
			#endif

			auto [fer, ber] = coder.simulate(snr, iter_cnt, max_error);

			#ifdef TEST
			auto end = std::chrono::system_clock::now();
			auto time_in_ms = std::chrono::duration_cast<ms>(end - start).count();
			#endif
			fout << "Snr = " << snr << "; FER = " << fer << "; BER = " << ber;
			#if defined(TEST) && defined(TIMELOG)
			fout << "; time = " << time_in_ms / 1000.0 << "s";
			#endif
			#if defined(CNTLOG)
			fout << SUM_CNT << " " << MUL_CNT << " " << CNT_BIN;
			#endif
		} else {
			fail(false, "main: incorrect command");
			continue;
		}
		fout << std::endl;
	}

	std::cout << "FINISH\n";
}

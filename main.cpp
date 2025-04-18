// #include "include/trust_propagation_decoder.h"
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

	// int n, k1, k2;
	// fin >> n >> k1 >> k2;
	// matrix G1(k1, binvector(n));
	// matrix G2(k2, binvector(n));
	// fin >> G1 >> G2;
	// recursive_decoder* rd1 = new recursive_decoder(G1);
	// recursive_decoder* rd2 = new recursive_decoder(G2);

	int n, k;
	fin >> n >> k;
	matrix G(k, binvector(n));
	fin >> G;
	recursive_decoder* rd = new recursive_decoder(G);

	// plotkin_construction_decoder coder(rd, rd);
	turbo_decoder coder(rd, rd);

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
			fout << "Snr = " << snr << "; FER = " << fer << "; BER = " << ber;
			#endif
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

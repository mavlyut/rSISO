#include "RecursiveDecoder.cpp"

int main() {
	srand(time(NULL));

	std::ifstream fin("input.txt");
#if defined(LOCAL)
	std::ostream& fout = std::cout;
#else
	std::ofstream fout("output.txt");
#endif
	// fout << std::scientific;

	int n, k;
	fin >> n >> k;
	std::vector<binvector> G(k, binvector(n));
	fin >> G;

	RecursiveDecoder coder(G);

	#ifdef PRINT_SECTION_TREE
	coder.printTree("Sectionalization.dot");
	#endif

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
		} else if (command == "Decode" || command == "DecodeSISO") {
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

			auto [errr, sim_cnt] = simulate(coder, snr, iter_cnt, max_error);

			#ifdef TEST
			auto end = std::chrono::system_clock::now();
			auto time_in_ms = std::chrono::duration_cast<ms>(end - start).count();
			#endif
			#ifdef TEST
			fout << (errr + 0.0) / sim_cnt;
			#endif
			#if defined(TEST) && defined(TIMELOG)
			fout << " " << time_in_ms / 1000.0;
			#endif
			#if defined(CNTLOG)
			fout << SUM_CNT << " " << MUL_CNT << " " << CNT_BIN;
			#endif
		} else {
			fail(false, "main: incorrect command");
		}
		fout << std::endl;
	}

	#ifdef TEST
	std::cout << "FINISH\n";
	#endif
}

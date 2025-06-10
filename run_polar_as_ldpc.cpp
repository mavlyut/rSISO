#include <string>
#include "include/plotkin_construction_decoder.h"
#include "include/recursive_decoder.h"
#include "include/shuffle_decoder.h"
#include "include/turbo_decoder.h"

unsigned iter_cnt = 100000, max_error = 100;

using namespace short_domain;

static inline matrix Arikan_kernel{binvector(1), binvector(3)};

std::pair<matrix, std::vector<double>> __generate_Polar(unsigned m, double p, matrix const& Kernel) {
	if (m == 0) {
		return {Kernel, {p}};
	}
	auto [A, Z0] = __generate_Polar(m - 1, p, Kernel);
	unsigned sh = (1ull << (m - 1)), n = (1ull << m);
	matrix B(n, get_empty(n));
	std::vector<double> Z;
	for (unsigned i = 0; i < sh; i++) {
		for (unsigned j = 0; j < sh; j++) {
			setbit(B[i], j, getbit(A[i], j));
			setbit(B[i + sh], j, getbit(A[i], j));
			setbit(B[i + sh], j + sh, getbit(A[i], j));
		}
		Z.push_back(Z0[i] * (2 - Z0[i]));
		Z.push_back(Z0[i] * Z0[i]);
	}
	return {B, Z};
}

matrix permutate_columns(matrix const& A, std::vector<unsigned> const& pi) {
	matrix B = A;
	for (unsigned i = 0; i < A.size(); ++i) {
		for (unsigned j = 0; j < pi.size(); ++j) {
			setbit(B[i], j, getbit(A[i], pi[j]));
		}
	}
	return B;
}

std::vector<unsigned> get_bit_reversal_perm(unsigned n) {
	std::vector<unsigned> ans;
	for (unsigned i = 0; i < (1ull << n); ++i) {
		unsigned tmp = 0;
		for (unsigned j = 0; j < n; ++j) {
			if ((i >> (n - 1 - j)) & 1) {
				tmp |= (1ull << j);
			}
		}
		ans.push_back(tmp);
	}
	return ans;
}

std::pair<matrix, matrix> generate_polar_component(unsigned m, unsigned k, double p, matrix const& Kernel) {
	unsigned n = (1ull << m), n0 = n / 2;
	auto [A, Z] = __generate_Polar(m, p, Kernel);

	std::vector<std::pair<double, unsigned>> Z_paired;
	for (unsigned i = 0; i < n; i++) {
		Z_paired.emplace_back(Z[i], i);
	}
	std::sort(Z_paired.begin(), Z_paired.end());
	double threshold = Z_paired[k - 1].first;

	matrix G1, G2;
	for (unsigned j = 0; j < n && G1.size() + G2.size() < k; j++) {
		if (Z[j] <= threshold) {
			if (j < n0) {
				G1.push_back(subvector(A[j], 0, n0));
			} else {
				G2.push_back(subvector(A[j], 0, n0));
			}
		}
	}

	auto b = get_bit_reversal_perm(m - 1);
	return {permutate_columns(G1, b), permutate_columns(G2, b)};
}

matrix generate_polar(unsigned m, unsigned k, double p, matrix const& Kernel) {
	unsigned n = (1ull << m);
	auto [A, Z] = __generate_Polar(m, p, Kernel);

	std::vector<std::pair<double, unsigned>> Z_paired;
	for (unsigned i = 0; i < n; i++) {
		Z_paired.emplace_back(Z[i], i);
	}
	std::sort(Z_paired.begin(), Z_paired.end());
	double threshold = Z_paired[k - 1].first;

	matrix G;
	for (unsigned j = 0; j < n && G.size() < k; j++) {
		if (Z[j] <= threshold) {
			G.push_back(A[j]);
		}
	}
	return permutate_columns(G, get_bit_reversal_perm(m));
}

int main() {
	srand(time(NULL));

	std::ofstream fout("output.txt");

	fout << "FER,rSISO\tFER,Plotkin\n";
	unsigned M = 5;
	for (unsigned K = 1; K <= 32; ++K) {
		fout << "\nK = " << K << std::endl;
		unsigned n = (1ull << M), m = n / 2, k = K;

		auto [G1, G2] = generate_polar_component(M, K, 0.5, Arikan_kernel);

		std::vector<unsigned> bit_revs00(n), bit_revs0 = get_bit_reversal_perm(M - 1), bit_revs = get_bit_reversal_perm(M);
		std::vector<unsigned> swap_parts(n);
		for (unsigned i = 0; i < m; ++i) {
			swap_parts[i] = i + m;
			swap_parts[i + m] = i;
			bit_revs00[i] = bit_revs0[i];
			bit_revs00[i + m] = bit_revs0[i] + m;
		}

		std::vector<unsigned> perm(n);
		for (unsigned i = 0; i < n; ++i) {
			perm[i] = bit_revs00[swap_parts[bit_revs[i]]];
		}
		// std::cout << perm << std::endl;

		shuffle_decoder perm_pt_coder(
			perm,
			new plotkin_construction_decoder(
				new recursive_decoder(m, G2),
				new recursive_decoder(m, G1)
			)
		);
		recursive_decoder rSISO_coder(n, generate_polar(M, K, 0.5, Arikan_kernel));

		// printmatrix(std::cerr, "shuffle gen:", n, perm_pt_coder.generate_matrix()); std::cerr << "\n";
		// printmatrix(std::cerr, "rSISO gen:", n, rSISO_coder.generate_matrix()); std::cerr << "\n";
		// std::cerr << "real pt gen:\n";
		// for (unsigned i = 0; i < k; ++i) {
		//	 printbv(std::cerr, n, perm_pt_coder.encode(1ull << i)); std::cerr << "\n";
		// }

		#ifdef DEBUG
			double start_snr = 8;
			double end_snr = 8;
		#else
			double start_snr = 0;
			double end_snr = 5;
		#endif

		double metric = 0, denom = 0;
		for (double snr = start_snr; snr <= end_snr; snr += 0.5) {
			auto start = std::chrono::system_clock::now();

			auto [fer_pt, ber_pt] = perm_pt_coder.simulate(snr, iter_cnt, max_error);

			auto end = std::chrono::system_clock::now();
			auto time_in_ms_pt = std::chrono::duration_cast<ms>(end - start).count();


			start = std::chrono::system_clock::now();

			auto [fer_rec, ber_rec] = rSISO_coder.simulate(snr, iter_cnt, max_error);

			end = std::chrono::system_clock::now();
			auto time_in_ms_rec = std::chrono::duration_cast<ms>(end - start).count();

			#ifdef TSV_FORMAT
				fout
					<< std::to_string(fer_rec).replace(1, 1, ",") << "\t"
					<< std::to_string(fer_pt).replace(1, 1, ",");
				// std::cout << "snr = " << snr << "; time(plotkin) = " << time_in_ms_pt / 1000.0 << "s; time(rSISO) = " << time_in_ms_rec / 1000.0 << "s" << std::endl;
			#else
				fout << "Snr = " << snr << "\n";
				fout << "\trSISO: FER = " << fer_rec << "; BER = " << ber_rec << "; time = " << time_in_ms_rec / 1000.0 << "s\n";
				fout << "\tPt   : FER = " << fer_pt << "; BER = " << ber_pt << "; time = " << time_in_ms_pt / 1000.0 << "s\n";
			#endif
			fout << std::endl;

			metric += fer_pt / fer_rec;
			++denom;
		}

		std::cout << "K=" << K << "; metric=" << metric / denom << std::endl;
		fout << "metric=" << metric / denom << std::endl;
	}

	std::cout << "FINISH\n";
}

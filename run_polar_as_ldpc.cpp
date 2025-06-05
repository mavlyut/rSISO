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

std::pair<std::pair<matrix, std::pair<matrix, matrix>>, std::vector<unsigned>> 
generate_Polar_as_Pt(unsigned m, unsigned k, double p, matrix const& Kernel) {
    auto [A, Z] = __generate_Polar(m, p, Kernel);
    std::vector<std::pair<double, unsigned>> Z_paired;
    for (unsigned i = 0; i < Z.size(); i++) {
        Z_paired.emplace_back(Z[i], i);
    }
    std::sort(Z_paired.begin(), Z_paired.end());
    double threshold = Z_paired[k - 1].first;
    matrix G1, G2, G, G0;
    std::vector<unsigned> revs(1 << m);
    for (unsigned i = 0; i < revs.size(); i++) {
        unsigned rev_i = 0;
        for (unsigned j = 0; j < m; j++) {
            if ((i >> (m - 1 - j)) & 1) {
                rev_i |= (1 << j);
            }
        }
        revs[i] = rev_i;
    }
    for (unsigned j = 0; j < Z.size() && G1.size() + G2.size() < k; j++) {
        if (Z[j] <= threshold) {
            if (j < (1ull << (m - 1))) {
                G1.push_back(subvector(A[j], 0, (1ull << (m - 1))));
            } else {
                G2.push_back(subvector(A[j], 0, (1ull << (m - 1))));
            }
            binvector row = get_empty(Z.size());
            for (unsigned c = 0; c < Z.size(); c++) {
                setbit(row, c, getbit(A[j], revs[c]));
            }
            G.push_back(row);
            G0.push_back(A[j]);
        }
    }
    // printmatrix(std::cerr, (1ull << m), G0);
    return {{G, {G1, G2}}, revs};
}

int main() {
	srand(time(NULL));

	std::ofstream fout("output.txt");

    unsigned K = 32, M = 6;
    unsigned n = (1ull << M), m = n / 2, k = K;
    auto [Gs, bit_reversal] = generate_Polar_as_Pt(M, K, 0.5, Arikan_kernel);
    auto [G, G12] = Gs;
    auto [G1, G2] = G12;
    // printmatrix(std::cout, m, G1);
    // printmatrix(std::cout, m, G2);

    std::vector<unsigned> swap_parts(n);
    for (unsigned i = 0; i < m; i++) {
        swap_parts[i] = i + m;
        swap_parts[i + m] = i;
    }

    std::vector<unsigned> perm(n);
    for (unsigned i = 0; i < n; ++i) {
        perm[i] = swap_parts[bit_reversal[i]];
    }

    std::cout << perm << std::endl;
    shuffle_decoder perm_pt_coder(
        perm,
        new plotkin_construction_decoder(
            new recursive_decoder(m, G2),
            new recursive_decoder(m, G1)
        )
    );
	recursive_decoder rSISO_coder(n, G);

    // printmatrix(std::cerr, "shuffle gen:", n, perm_pt_coder.generate_matrix());

	for (double snr = 0.0; snr <= 5; snr += 0.5) {
		fail(iter_cnt > 0, "simulate: iter_cnt must be positive");
		fail(max_error >= 0, "simulate: max_err_cnt must be non-negative");
		_Float64 sigma = sqrt(0.5 * pow(10.0, -snr / 10.0) * n / k);
        std::normal_distribution<_Float64> norm(0.0, sigma);

		auto start = std::chrono::system_clock::now();

        auto [fer_pt, ber_pt] = perm_pt_coder.simulate(snr, iter_cnt, max_error);

		auto end = std::chrono::system_clock::now();
		auto time_in_ms_pt = std::chrono::duration_cast<ms>(end - start).count();


		start = std::chrono::system_clock::now();

        auto [fer_rec, ber_rec] = std::pair{0,0};//rSISO_coder.simulate(snr, iter_cnt, max_error);

		end = std::chrono::system_clock::now();
		auto time_in_ms_rec = std::chrono::duration_cast<ms>(end - start).count();

		#ifdef TSV_FORMAT
			fout
                << std::to_string(fer_rec).replace(1, 1, ",") << "\t"
				<< std::to_string(fer_pt).replace(1, 1, ",");
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

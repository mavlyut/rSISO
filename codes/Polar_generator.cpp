#include <iostream>

#include "../include/matrix.h"

static inline matrix Arikan_kernel{binvector(2, 1), binvector(2, 3)};

matrix __generate_Polar(int m, matrix const& Kernel) {
    if (m == 0) {
        return Kernel;
    }
    matrix A = __generate_Polar(m - 1, Kernel);
    std::size_t sh = (1ull << (m - 1)), n = (1ull << m);
    matrix B(n, binvector(n));
    for (std::size_t i = 0; i < sh; i++) {
        for (std::size_t j = 0; j < sh; j++) {
            B[i].set(j, A[i][j]);
            B[i + sh].set(j, A[i][j]);
            B[i + sh].set(j + sh, A[i][j]);
        }
    }
    return B;
}

matrix generate_Polar(int m, std::size_t k, double p, matrix const& Kernel) {
    matrix A = __generate_Polar(m, Kernel);
    std::vector<double> Z_prev(1, p), Z;
    for (int i = 0; i < m; i++) {
        Z.clear();
        for (double j : Z_prev) {
            Z.push_back(2 * j - j * j);
            Z.push_back(j * j);
        }
        Z_prev = Z;
    }
    std::vector<std::pair<double, std::size_t>> Z_paired;
    for (std::size_t i = 0; i < Z.size(); i++) {
        Z_paired.emplace_back(Z[i], i);
    }
    std::sort(Z_paired.begin(), Z_paired.end());
    double threshold = Z_paired[k - 1].first;
    matrix G;
    for (std::size_t j = 0; j < Z.size() && G.size() < k; j++) {
        if (Z[j] <= threshold) {
            G.push_back(A[j]);
        }
    }
    return G;
}

int main() {
    int m = 6, k = 32;
    std::ofstream fout("../Polar_" + std::to_string(1ull << m) + "_" + std::to_string(k) + ".txt");
    fout << (1ull << m) << " " << k << "\n"
         << generate_Polar(m, k, 0.5, Arikan_kernel);
}

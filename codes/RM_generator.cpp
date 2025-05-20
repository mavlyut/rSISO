#include <iostream>
#include <vector>

#include "../include/matrix.h"

matrix generate_RM(int r, int m) {
    if (r == 0) {
        matrix G(1, binvector(1ull << m));
        for (unsigned i = 0; i < (1ull << m); i++) {
            G[0].set(i, true);
        }
        return G;
    }
    if (r == m) {
        matrix G(1ull << m, binvector(1ull << m));
        for (unsigned i = 0, deg = (1ull << m); i < (1ull << m); i++, deg /= 2) {
            for (unsigned j = 0; j < (1ull << m); j++) {
            }
            G[i].set(i, true);
        }
        return G;
    }
    matrix G1 = generate_RM(r, m - 1), G2 = generate_RM(r - 1, m - 1);
    unsigned sh = (1ull << (m - 1));
    matrix G;
    for (unsigned i = 0; i < G1.size(); i++) {
        binvector row(1ull << m);
        for (unsigned j = 0; j < sh; j++) {
            row.set(j, G1[i][j]);
            row.set(j + sh, G1[i][j]);
        }
        G.push_back(row);
    }
    for (unsigned i = 0; i < G2.size(); i++) {
        binvector row(1ull << m);
        for (unsigned j = 0; j < sh; j++) {
            row.set(j + sh, G2[i][j]);
        }
        G.push_back(row);
    }
    return G;
}

int main() {
    matrix G = generate_RM(2, 2);
    unsigned k = G.size(), n = G.front().size();
    std::ofstream fout("codes/RM_" + std::to_string(n) + "_" + std::to_string(k) + ".txt");
    fout << n << " " << k << "\n" << G;
    std::cout << "RM (" << n << "," << k << ") generated\n";
}

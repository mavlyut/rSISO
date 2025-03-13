#include <iostream>
#include <fstream>
#include <vector>

#include "../include/matrix.h"

// p in P; 0 < r < l
matrix generateLDPMatrix(int p, int r, int l) {
    matrix ans(r * p, binvector(l * p));
    for (int i = 0; i < r; i++) {
        for (int j = 0; j < l; j++) {
            int shift = 0;
            if (i > 0 && j > 0) {
                shift = ((rand() % p) + p) % p;
            }
            for (int c = 0; c < p; c++) {
                ans[i * p + ((c + shift) % p)].set(j * p + c, 1);
            }
        }
    }
    return ans;
}

int main() {
    // int r = 10, l = 23, p = 37;
    int r = 3, l = 6, p = 17;

    std::vector<binvector> H = generateLDPMatrix(p, r, l);
    int n = H.front().size();
    int k = n - H.size();
    std::ofstream fout("LDPC_" + std::to_string(n) + "_" + std::to_string(k) + ".txt");
    fout << l * p << " " << r * p << "\n";
    for (binvector const& row : H) {
        fout << row << "\n";
    }
}

#include <iostream>
#include <vector>

#include "../include/matrix.h"

using namespace short_domain;

matrix generate_RM(int r, int m) {
	if (r == 0) {
		matrix G(1, get_empty(1ull << m));
		for (unsigned i = 0; i < (1ull << m); ++i) {
			setbit(G[0], i, true);
		}
		return G;
	}
	if (r == m) {
		matrix G(1ull << m, get_empty(1ull << m));
		for (unsigned i = 0; i < (1ull << m); ++i) {
			setbit(G[i], i, true);
		}
		return G;
	}
	matrix G1 = generate_RM(r, m - 1), G2 = generate_RM(r - 1, m - 1);
	unsigned sh = (1ull << (m - 1));
	matrix G;
	for (unsigned i = 0; i < G1.size(); ++i) {
		binvector row = get_empty(1ull << m);
		for (unsigned j = 0; j < sh; j++) {
			setbit(row, j, getbit(G1[i], j));
			setbit(row, j + sh, getbit(G1[i], j));
		}
		G.push_back(row);
	}
	for (unsigned i = 0; i < G2.size(); ++i) {
		binvector row = get_empty(1ull << m);
		for (unsigned j = 0; j < sh; j++) {
			setbit(row, j + sh, getbit(G2[i], j));
		}
		G.push_back(row);
	}
	return G;
}

int main() {
	unsigned R = 2, M = 4;
	matrix G = generate_RM(R, M);
	unsigned k = G.size(), n = (1ull << M);
	std::ofstream fout("codes/RM_" + std::to_string(n) + "_" + std::to_string(k) + ".txt");
	printmatrix(fout, std::to_string(n) + " " + std::to_string(k), n, G);
	std::cout << "RM (" << n << "," << k << ") generated\n";
}

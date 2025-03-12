#include "../include/matrix.h"

std::vector<int> gauss(matrix& a) {
	unsigned k = a.size();
	unsigned n = a[0].size();
	unsigned last_c = 0;
	std::vector<int> gamma;
	for (unsigned i = 0; i < k; ++i) {
		while (last_c < n) {
			bool fl = false;
			for (int j = i; j < k; ++j) {
				if (a[j][last_c]) {
					fl = true;
					std::swap(a[i], a[j]);
					break;
				}
			}
			if (fl) {
				break;
			}
			last_c++;
		}
		if (last_c == n) {
			for (unsigned ii = i; ii < k; ii++) {
				a.pop_back();
			}
			break;
		}
		gamma.push_back(last_c);
		for (unsigned j = i + 1; j < k; ++j) {
			if (a[j][last_c]) {
				for (unsigned t = 0; t < n; ++t) {
					a[j].set(t, a[j][t] ^ a[i][t]);
				}
			}
		}
		last_c++;
	}
	return gamma;
}

void full_gauss(matrix& a) {
	std::vector<int> gamma = gauss(a);
	for (unsigned i = gamma.size(); i-- > 0; ) {
		for (unsigned j = 0; j < i; j++) {
			if (a[j][gamma[i]]) {
				a[j] ^= a[i];
			}
		}
	}
}

std::vector<int> minimal_span_form(matrix& G) {
	std::vector<int> gamma = gauss(G);
	unsigned n = G[0].size(), k = G.size();
	for (unsigned i = k; i-- > 0; ) {
		unsigned j = n;
		while (j-- > 0 && G[i][j] == 0) {}
		for (unsigned r = 0; r < i; r++) {
			if (G[r][j] == 1) {
				G[r] ^= G[i];
			}
		}
	}
	return gamma;
}

unsigned rg(matrix M) {
	gauss(M);
	return M.size();
}

bool lin_indep(matrix const& M, binvector const& a) {
	matrix M_ext(M);
	M_ext.push_back(a);
	return rg(M_ext) == M_ext.size();
}

std::ostream& operator<<(std::ostream& out, matrix const& a) {
	for (binvector const& row : a) {
			for (unsigned i = 0; i < row.size(); i++) {
			if (i != 0) {
				out << ' ';
			}
			out << row[i];
		}
		out << '\n';
	}
	return out;
}

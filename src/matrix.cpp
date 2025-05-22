#include "../include/matrix.h"

namespace short_domain {
    binvector multiply(binvector const& bv, matrix const& M) {
		unsigned i = 0;
		binvector ans = 0, tmp = bv;
		while (tmp > 0) {
			if (tmp & 1) {
				ans ^= M[i];
			}
			++i;
			tmp >>= 1;
		}
		return ans;
	}

	std::vector<unsigned> gauss(unsigned n, matrix& a) {
		unsigned k = a.size();
		unsigned last_c = 0;
		std::vector<unsigned> gamma;
		for (unsigned i = 0; i < k; ++i) {
			while (last_c < n) {
				bool fl = false;
				for (unsigned j = i; j < k; ++j) {
					if (getbit(a[j], last_c)) {
						fl = true;
						std::swap(a[i], a[j]);
						break;
					}
				}
				if (fl) {
					break;
				}
				++last_c;
			}
			if (last_c == n) {
				for (unsigned ii = i; ii < k; ++ii) {
					a.pop_back();
				}
				break;
			}
			gamma.push_back(last_c);
			for (unsigned j = i + 1; j < k; ++j) {
				if (getbit(a[j], last_c)) {
					for (unsigned t = 0; t < n; ++t) {
						setbit(a[j], t, getbit(a[j], t) ^ getbit(a[i], t));
					}
				}
			}
			++last_c;
		}
		return gamma;
	}

	void full_gauss(unsigned n, matrix& a) {
		std::vector<unsigned> gamma = gauss(n, a);
		for (unsigned i = gamma.size(); i-- > 0; ) {
			for (unsigned j = 0; j < i; ++j) {
				if (getbit(a[j], gamma[i])) {
					a[j] ^= a[i];
				}
			}
		}
	}

	void minimal_span_form(unsigned n, matrix& G) {
		gauss(n, G);
		unsigned k = G.size();
		for (unsigned i = k; i-- > 0; ) {
			unsigned j = n;
			while (j-- > 0 && getbit(G[i], j) == 0) {}
			for (unsigned r = 0; r < i; ++r) {
				if (getbit(G[r], j)) {
					G[r] ^= G[i];
				}
			}
		}
	}

	unsigned rg(unsigned n, matrix M) {
		gauss(n, M);
		return M.size();
	}

	bool lin_indep(unsigned n, matrix const& M, binvector const& a) {
		matrix M_ext(M);
		M_ext.push_back(a);
		return rg(n, M_ext) == M_ext.size();
	}

    std::ostream& printmatrix(std::ostream& out, std::string const& msg, unsigned n, matrix const& M) {
		out << msg;
		for (binvector const& bv : M) {
			printbv(out, n, bv);
			out << "\n";
		}
		return out;
	}

    std::ostream& printmatrix(std::ostream& out, unsigned n, matrix const& M) {
		return printmatrix(out, "", n, M);
	}

    void _log_matrix(std::string const& msg, unsigned n, matrix const& M) {
		__log(msg << "\n");
		for (binvector const& bv : M) {
			_log_bv("", n, bv);
			__log(std::endl);
		}
		__log(std::endl);
	}
}

namespace long_domain {
    binvector multiply(binvector const& bv, matrix const& M) {
		binvector ans(bv.size());
		for (unsigned i = 0; i < bv.size(); ++i) {
			if (bv[i]) {
				ans ^= M[i];
			}
		}
		return ans;
	}

	std::vector<unsigned> gauss(matrix& a) {
		unsigned k = a.size();
		unsigned n = a[0].size();
		unsigned last_c = 0;
		std::vector<unsigned> gamma;
		for (unsigned i = 0; i < k; ++i) {
			while (last_c < n) {
				bool fl = false;
				for (unsigned j = i; j < k; ++j) {
					if (a[j][last_c]) {
						fl = true;
						std::swap(a[i], a[j]);
						break;
					}
				}
				if (fl) {
					break;
				}
				++last_c;
			}
			if (last_c == n) {
				for (unsigned ii = i; ii < k; ++ii) {
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
			++last_c;
		}
		return gamma;
	}

	void full_gauss(matrix& a) {
		std::vector<unsigned> gamma = gauss(a);
		for (unsigned i = gamma.size(); i-- > 0; ) {
			for (unsigned j = 0; j < i; ++j) {
				if (a[j][gamma[i]]) {
					a[j] ^= a[i];
				}
			}
		}
	}

	void minimal_span_form(matrix& G) {
		std::vector<unsigned> gamma = gauss(G);
		unsigned n = G[0].size(), k = G.size();
		for (unsigned i = k; i-- > 0; ) {
			unsigned j = n;
			while (j-- > 0 && G[i][j] == 0) {}
			for (unsigned r = 0; r < i; ++r) {
				if (G[r][j] == 1) {
					G[r] ^= G[i];
				}
			}
		}
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

    std::ostream& printmatrix(std::ostream& out, std::string const& msg, unsigned n, matrix const& M) {
		out << msg;
		for (binvector const& bv : M) {
			printbv(out, n, bv);
			out << "\n";
		}
		return out;
	}

    std::ostream& printmatrix(std::ostream& out, unsigned n, matrix const& M) {
		return printmatrix(out, "", n, M);
	}

    void _log_matrix(std::string const& msg, unsigned n, matrix const& M) {
		__log(msg);
		for (binvector const& bv : M) {
			_log_bv("", n, bv);
			__log(std::endl);
		}
	}
}

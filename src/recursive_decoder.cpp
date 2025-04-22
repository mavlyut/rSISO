#include "../include/gray_code.h"
#include "../include/recursive_decoder.h"

#define subvector(v, x, y) (((v) >> (x)) & ((1ull << ((y) - (x))) - 1))
#define getbit(v, i) (((v) >> (i)) & 1)
#define setbit(v, i, b) if (b) {v |= (1ull << (i));} else {v &= ~(1ull << (i));}
#define changebit(v, i) v ^= (1ull << (i))
#define printbv(msg, n, v) 					\
	__log(msg);								\
	for (unsigned j = 0; j < n; ++j) { 		\
		if (j != 0) __log(" ");				\
		__log((((v >> j) & 1) ? 1 : 0));	\
	}
#define printmatrix(msg, n, M)					\
	__log(msg << "\n");							\
	for (std::size_t const& i : M) { 			\
		printbv("", n, i);						\
		__log(std::endl);						\
	}											\
	__log(std::endl);

std::vector<unsigned> gauss(unsigned n, std::vector<std::size_t>& a) {
	// printmatrix("Gauss", n, a);
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
	// printmatrix("After Gauss: ", n, a);
	return gamma;
}

void full_gauss(unsigned n, std::vector<std::size_t>& a) {
	std::vector<unsigned> gamma = gauss(n, a);
	for (unsigned i = gamma.size(); i-- > 0; ) {
		for (unsigned j = 0; j < i; ++j) {
			if (getbit(a[j], gamma[i])) {
				a[j] ^= a[i];
			}
		}
	}
}

void minimal_span_form(unsigned n, std::vector<std::size_t>& G) {
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

unsigned rg(unsigned n, std::vector<std::size_t> M) {
	gauss(n, M);
	return M.size();
}

bool lin_indep(unsigned n, std::vector<std::size_t> const& M, std::size_t const& a) {
	std::vector<std::size_t> M_ext(M);
	M_ext.push_back(a);
	return rg(n, M_ext) == M_ext.size();
}

recursive_decoder::recursive_decoder(unsigned n, matrix const& _G) : n(n), k(_G.size()), G(_G), G_enc(_G), L_out(n) {
	time_measure(minimal_span_form(n, G));

	printmatrix("MSF:", n, G);

	p0.resize(n);
	p1.resize(n);
	_split = std::vector(n + 1, std::vector<unsigned>(n + 1, UNINIT));
	_ctor = std::vector(n + 1, std::vector<unsigned>(n + 1, 0));
	time_measure(sectionalization());
	time_measure(root = rec_build_section_tree(0, n));
}

unsigned recursive_decoder::length() const {
	return n;
}

unsigned recursive_decoder::dim() const {
	return k;
}

std::size_t recursive_decoder::encode(std::size_t const& c) {
	printbv("Encode: ", k, c); __log(std::endl);
    std::size_t ans = 0;
    for (unsigned i = 0; i < k; ++i) {
        if (getbit(c, i)) {
            ans ^= G_enc[i];
        }
    }
    return ans;
}

const double recursive_decoder::MAX_TRUNC = 16;
double recursive_decoder::truncate(double const &a) {
    if (a > MAX_TRUNC) {
        return MAX_TRUNC;
    }
    if (a < -MAX_TRUNC) {
        return -MAX_TRUNC;
    }
    return a;
}

std::pair<double, double> recursive_decoder::simulate(double snr, unsigned iter_cnt, unsigned max_error) {
    _Float64 sigma = sqrt(0.5 * pow(10.0, -snr / 10.0) * n / k);
    std::normal_distribution<_Float64> norm(0.0, sigma);
    unsigned errr = 0, berr = 0, sim_cnt = 0;
    while (errr < max_error && sim_cnt < iter_cnt) {
        std::size_t x = 0;
        for (unsigned t = 0; t < k; ++t) {
			if (rand() % 2) {
				x |= (1ull << t);
			}
        }
        std::size_t enc = encode(x);
        _Float64 coef = 2.0f / (sigma * sigma);
        std::vector<_Float64> L_in(n);
        for (unsigned t = 0; t < n; ++t) {
            L_in[t] = coef * ((((enc >> t) & 1) ? -1 : 1) + norm(gen));
        }
        std::vector<double> soft_dec = decode_soft(L_in);
		std::size_t dec = 0;
		for (unsigned t = 0; t < n; ++t) {
			if (soft_dec[t] < 0) {
				dec |= (1ull << t);
			}
			if ((soft_dec[t] < 0) != ((enc >> t) & 1)) {
				++berr;
			}
		}
		printbv("Encoded: ", n, enc); __log("\n");
		printbv("Decoded: ", n, dec); __log(std::endl);
        errr += (dec != enc);
        ++sim_cnt;
    }
    return {(errr + 0.0) / sim_cnt, (berr + 0.0) / (sim_cnt * n)};
}

void recursive_decoder::sectionalization() {
	std::vector<std::vector<unsigned>> _Gp_size(n + 1, std::vector<unsigned>(n + 1, 0));
	std::vector<std::vector<unsigned>> _Gs_size(n + 1, std::vector<unsigned>(n + 1, 0));
	std::vector<std::vector<std::vector<std::size_t>>> phi_u_c(n + 1, std::vector<std::vector<std::size_t>>(n + 1, std::vector<std::size_t>(n + 1, UNINIT)));
	std::vector<std::vector<std::vector<std::size_t>>> phi_d_c(n + 1, std::vector<std::vector<std::size_t>>(n + 1, std::vector<std::size_t>(n + 1, UNINIT)));
	std::vector<std::vector<std::size_t>> phi_u_l(n + 1, std::vector<std::size_t>(n + 1, UNINIT));
	std::vector<std::vector<std::size_t>> phi_d_l(n + 1, std::vector<std::size_t>(n + 1, UNINIT));
	std::vector<std::vector<unsigned>> phi(n + 1, std::vector<unsigned>(n + 1, UNINIT));

	auto update_phi = [&](std::vector<std::vector<std::size_t>>& __phi, unsigned x, unsigned y, std::size_t const& val, unsigned ctorNo) {
		if (val < __phi[x][y]) {
			__phi[x][y] = val;
			_ctor[x][y] = ctorNo;
		}
	};

	auto update_uninit = [](std::size_t& x, std::size_t const& y) {
		if (x == UNINIT) {
			x = y;
		}
	};

	for (unsigned x = 0; x <= n; ++x) {
		for (unsigned y = x + 1; y <= n; ++y) {
			matrix Gp;
			unsigned k3 = 0;
			for (std::size_t const& row : G) {
				Gp.push_back(subvector(row, x, y));
				if (subvector(row, 0, x) == 0ull && subvector(row, y, n) == 0ull) {
					++k3;
				}
			}
			// __log("segment " << x << " " << y << std::endl);
			gauss(y - x, Gp);
			_Gp_size[x][y] = Gp.size();
			_Gs_size[x][y] = k3;

			#ifdef ENABLE_OPT_1
			if (Gp.size() == 1 && Gp[0] == ((1ull << (y - x)) - 1) && k3 == 1) {
				update_phi(phi_u_l, x, y, 2 * (y - x), 1);
				update_phi(phi_d_l, x, y, 0, 1);
			}
			#endif

			#ifdef ENABLE_OPT_2
			if (Gp.size() == 1 && Gp[0] == ((1ull << (y - x)) - 1) && k3 == 0) {
				update_phi(phi_u_l, x, y, 2 * (y - x), 2);
				update_phi(phi_d_l, x, y, 3, 2);
			}
			#endif

			#ifdef ENABLE_OPT_3
			if (Gp.size() == 2 && y - x == 2) {
				update_phi(phi_u_l, x, y, 6, 3);
				update_phi(phi_d_l, x, y, 8, 3);
			}
			#endif

			#ifdef ENABLE_OPT_5
			if (k3 == 0 && Gp.size() == y - x && y - x >= 2) {
				update_phi(phi_u_l, x, y, (1ull << (y - x)) + 2 * (y - x) - 1, 5);
				update_phi(phi_d_l, x, y, (1ull << (y - x + 1)) * (y - x + 1) + 2 * (y - x), 5);
			}
			#endif
		}
	}

	for (unsigned x = 0; x <= n; ++x) {
		for (unsigned y = x + 1; y <= n; ++y) {
			update_phi(phi_u_l, x, y, (1ull << _Gp_size[x][y]) * (2 * (y - x) - 1), 0);
			update_phi(phi_d_l, x, y, ((1ull << (_Gp_size[x][y] - _Gs_size[x][y] + 1))) * (y - x + 1) + 2 * (y - x), 0);
			for (unsigned z = x + 1; z < y - 1; ++z) {
				unsigned k1_k2 = _Gp_size[x][y] - _Gs_size[x][z] - _Gs_size[z][y];
				update_uninit(phi_u_c[x][y][z], (1ull << (k1_k2 + 1)));
				update_uninit(phi_d_c[x][y][z], (1ull << (k1_k2 + 2)));
			}
		}
	}

	std::function<unsigned(unsigned, unsigned)> get_phi = [&](unsigned x, unsigned y) -> unsigned {
		if (phi[x][y] == UNINIT) {
			unsigned ind_min = UNINIT;
			std::size_t min_score = phi_u_l[x][y] + phi_d_l[x][y];
			for (unsigned z = x + 1; z < y - 1; ++z) {
				std::size_t score = phi_u_c[x][y][z] + phi_d_c[x][y][z] + get_phi(x, z) + get_phi(z, y);
				if (score < min_score) {
					min_score = score;
					ind_min = z;
				}
			}
			_split[x][y] = ind_min;
			phi[x][y] = min_score;
		}
		return phi[x][y];
	};
	get_phi(0, n);

	#ifdef INIT_COUNTS
	std::cout << "Expected operations: " << get_phi(0, n) << std::endl;
	#endif

	__log("Gp_sz:  " << _Gp_size << "\nGs_sz:  " << _Gs_size << std::endl);
	__log("split: " << _split << std::endl);
}

recursive_decoder::node* recursive_decoder::rec_build_section_tree(unsigned x, unsigned y) {
	// __log("build, handle " << x << " " << y << " | " << _split[x][y] << std::endl);
	fail(x < y && y <= n, "rec_build: incorrect bounds");
	unsigned z = _split[x][y];
	if (z == UNINIT) {
		matrix Gp, Gs, Gl;
		for (unsigned i = 0; i < k; ++i) {
			std::size_t row = G[i];
			std::size_t append = subvector(row, x, y);
			if (subvector(row, 0, x) == 0ull && subvector(row, y, n) == 0ull) {
				Gs.push_back(append);
				Gp.push_back(append);
			} else {
				Gl.push_back(append);
			}
		}
		unsigned k1 = 0;
		full_gauss(y - x, Gl);
		for (std::size_t const& row : Gl) {
			if (lin_indep(y - x, Gp, row)) {
				Gp.push_back(row);
				++k1;
			}
		}

		unsigned optimizationNo = _ctor[x][y];
		#ifdef ENABLE_OPT_1
		if (optimizationNo == 1) {
			fail(k1 == 0 && Gp.size() == 1 && Gp[0] == ((1ull << (y - x)) - 1), "incorrect pre-conditions (1)");
			__log("Simplify 1: " << x << " " << y << std::endl);
			return new leaf_simplify_1(x, y, Gp);
		}
		#endif

		#ifdef ENABLE_OPT_2
		if (optimizationNo == 2) {
			fail(k1 == 1 && Gp.size() == 1 && Gp[0] == ((1ull << (y - x)) - 1), "incorrect pre-conditions (2)");
			__log("Simplify 2: " << x << " " << y << std::endl);
			return new leaf_simplify_2(x, y, Gp);
		}
		#endif

		#ifdef ENABLE_OPT_3
		if (optimizationNo == 3) {
			fail(y - x == 2 && Gp.size() == 2, "incorrect pre-conditions (3)");
			__log("Simplify 3: " << x << " " << y << std::endl);
			return new leaf_simplify_3(x, y, k1, Gp);
		}
		#endif

		#ifdef ENABLE_OPT_5
		if (optimizationNo == 5) {
			fail(y - x == k1, "incorrect pre-conditions (5)");
			__log("Simplify 5: " << x << " " << y << std::endl);
			return new leaf_simplify_5(x, y, Gp);
		}
		#endif

		__log("No simplify: " << x << " " << y << std::endl);
		return new leaf_no_simplify(x, y, k1, Gp);
	}

	matrix Gp, _G0, _G1; unsigned k3 = 0;
	for (unsigned i = 0; i < k; ++i) {
		std::size_t row = G[i];
		std::size_t append = subvector(row, x, y);
		if (subvector(row, 0, x) == 0ull && subvector(row, y, n) == 0ull) {
			if (subvector(row, x, z) == 0ull || subvector(row, z, y) == 0ull) {
				Gp.push_back(append);
			} else {
				_G0.push_back(append);
			}
			++k3;
		} else {
			_G1.push_back(append);
		}
	}

	// __log("build, inner, splitting matrix\n");
	// printmatrix("Gp", y - x, Gp);
	// printmatrix("G0", y - x, _G0);
	// printmatrix("G1", y - x, _G1);

	matrix G_0, G_1;
	unsigned k1 = 0, k2 = 0;
	for (std::size_t const& row : _G0) {
		if (lin_indep(y - x, Gp, row)) {
			Gp.push_back(row);
			G_0.push_back(subvector(row, 0, z - x));
			G_1.push_back(subvector(row, z - x, y - x));
			++k2;
		}
	}
	for (std::size_t const& row : _G1) {
		if (lin_indep(y - x, Gp, row)) {
			Gp.push_back(row);
			G_0.push_back(subvector(row, 0, z - x));
			G_1.push_back(subvector(row, z - x, y - x));
			++k1;
		}
	}

	__log("inner: " << x << " " << y << "\n");
	printmatrix("Gp:", y - x, Gp);
	printmatrix("G_0:", z - x, G_0);
	printmatrix("G_1:", y - z, G_1);

	node* left = rec_build_section_tree(x, z);
	node* right = rec_build_section_tree(z, y);

	matrix G_hat(k1 + k2, 0);
	matrix G_tilda(k1 + k2, 0);
	{
		matrix G_hat_ext(z - x, 0);
		for (unsigned i = 0; i < z - x; ++i) {
			unsigned j = 0;
			for (unsigned t = 0; t < left->Gp.size(); ++t, ++j) {
				setbit(G_hat_ext[i], j, getbit(left->Gp[t], i));
			}
			for (unsigned t = 0; t < G_0.size(); ++t, ++j) {
				setbit(G_hat_ext[i], j, getbit(G_0[t], i));
			}
		}
		full_gauss(left->k1 + left->k3 + k1 + k2, G_hat_ext);

		printmatrix("G_hat_ext:", left->k1 + left->k3 + k1 + k2, G_hat_ext);

		for (unsigned i = 0; i < k1 + k2; ++i) {
			for (unsigned j = 0; j < left->k1; ++j) {
				setbit(G_hat[i], left->k1 - j - 1, getbit(G_hat_ext[G_hat_ext.size() - j - 1], i + left->k1 + left->k3));
			}
		}
	}
	{
		matrix G_tilda_ext(y - z, 0);
		for (unsigned i = 0; i < y - z; ++i) {
			unsigned j = 0;
			for (unsigned t = 0; t < right->Gp.size(); ++t, ++j) {
				setbit(G_tilda_ext[i], j, getbit(right->Gp[t], i));
			}
			for (unsigned t = 0; t < G_1.size(); ++t, ++j) {
				setbit(G_tilda_ext[i], j, getbit(G_1[t], i));
			}
		}
		full_gauss(right->k1 + right->k3 + k1 + k2, G_tilda_ext);

		printmatrix("G_tilda_ext:", right->k1 + right->k3 + k1 + k2, G_tilda_ext);

		for (unsigned i = 0; i < k1 + k2; ++i) {
			for (unsigned j = 0; j < right->k1; ++j) {
				setbit(G_tilda[i], right->k1 - j - 1, getbit(G_tilda_ext[G_tilda_ext.size() - j - 1], i + right->k1 + right->k3));
			}
		}
	}

	printmatrix("G_hat:", left->k1, G_hat);
	printmatrix("G_tilda:", right->k1, G_tilda);

	return new inner(x, y, z, left, right, k1, k2, k3, G_hat, G_tilda, Gp);
}

std::vector<double> recursive_decoder::decode_soft(std::vector<_Float64> const& L_in) {
	__log("Decode: " << L_in << std::endl);
	for (unsigned i = 0; i < n; ++i) {
		_Float64 L_exp = exp(truncate(L_in[i]));
		_Float64 z = 1.0 + L_exp;
		p0[i] = L_exp / z;
		p1[i] = 1.0 / z;
	}

	time_measure(root->upward_pass(p0, p1));
	time_measure(root->downward_pass(L_out));

	__log("Ans   : " << L_out << std::endl);
	return L_out;
}

void recursive_decoder::printTree(std::string const& fileName) const {
	std::ofstream fout(fileName);
	fout << "digraph G {\n";
	root->print(fout);
	fout << "}";
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

recursive_decoder::node::node(unsigned x, unsigned y, unsigned k1, unsigned k2, unsigned k3, matrix const& Gp)
	: x(x), y(y), k1(k1), k2(k2), k3(k3), A(1ull << k1, I), B(1ull << k1, I), Gp(Gp), traverse_order(k1 + k2) {
	__log("Node ctor: x=" << x << ", y=" << y << "; k=(" << k1 << ", " << k2 << ", " << k3 << ")\n");
}

void recursive_decoder::node::print(std::ostream& out) const {
	out << "\t" << getName() << " [label=\"[" << x << ";" << y << ")\"]\n";
	__print(out);
}

std::string recursive_decoder::node::getName() const {
	return "node_" + std::to_string(x) + "_" + std::to_string(y);
}

recursive_decoder::leaf::leaf(unsigned x, unsigned y, unsigned k1, matrix const& Gp)
	: node(x, y, k1, Gp.size() - k1, Gp.size() - k1, Gp) {}

double recursive_decoder::leaf::F(std::size_t const& c, std::vector<double> const& p0, std::vector<double> const& p1) const {
	double ans = 1.0;
	for (unsigned i = x, j = 0; i < y; ++i, ++j) {
		ans *= (getbit(c, j) ? p1[i] : p0[i]);
	}
	return ans;
}

double recursive_decoder::leaf::LLR(double M0, double M1) const {
	// if (M0 < EPS && M1 < EPS) {
	// 	// __log("Leaf, LLR: nan" << std::endl);
	// 	return NAN;
	// } else if (M0 < EPS) {
	// 	// __log("Leaf, LLR: -inf" << std::endl);
	// 	return -INF;
	// } else if (M1 < EPS) {
	// 	// __log("Leaf, LLR: inf" << std::endl);
	// 	return INF;
	// }
	return truncate(log(M0 / M1));
}

void recursive_decoder::leaf::__print(std::ostream&) const {}

recursive_decoder::leaf_no_simplify::leaf_no_simplify(unsigned x, unsigned y, unsigned k1, matrix const& Gp)
	: leaf(x, y, k1, Gp), P0(y - x), P1(y - x)
	, A0(y - x, std::vector<double>(1ull << k1)), A1(y - x, std::vector<double>(1ull << k1)) {}

void recursive_decoder::leaf_no_simplify::upward_pass(std::vector<double> const& p0, std::vector<double> const& p1) {
	std::size_t v = 0, c = 0;
	double T;
	for (unsigned ind : traverse_order) {
		if (ind != UNINIT) {
			if (ind >= k2) {
				changebit(v, ind - k2);
				A[v] = B[v] = I;
				for (unsigned i = 0; i < y - x; ++i) {
					A0[i][v] = A1[i][v] = I;
				}
			}
			c ^= Gp[ind];
		} else {
			A[v] = B[v] = I;
			for (unsigned i = 0; i < y - x; ++i) {
				A0[i][v] = A1[i][v] = I;
			}
		}
		T = F(c, p0, p1);
		A[v] += T;
		for (unsigned i = 0; i < y - x; ++i) {
			if (getbit(c, i)) {
				A1[i][v] += T;
			} else {
				A0[i][v] += T;
			}
		}
	}
}

void recursive_decoder::leaf_no_simplify::downward_pass(std::vector<double>& L) {
	for (unsigned i = 0; i < y - x; ++i) {
		P0[i] = P1[i] = 0.0;
	}
	for (std::size_t v = 0; v < (1ull << k1); ++v) {
		for (unsigned i = 0; i < y - x; ++i) {
			P0[i] += A0[i][v] * B[v];
			P1[i] += A1[i][v] * B[v];
		}
	}
	for (unsigned i = x, j = 0; i < y; ++i, ++j) {
		L[i] = LLR(P0[j], P1[j]);
	}
}

#ifdef ENABLE_OPT_1

recursive_decoder::leaf_simplify_1::leaf_simplify_1(unsigned x, unsigned y, matrix const& Gp) : leaf(x, y, 0, Gp) {
	fail(Gp.size() == 1 && Gp[0] == ((1ull << (y - x)) - 1), "leaf_s1: inappropriate optimization");
}

void recursive_decoder::leaf_simplify_1::upward_pass(std::vector<double> const& p0, std::vector<double> const& p1) {
	ext_l = I;
	for (unsigned i = x; i < y; ++i) {
		ext_l += LLR(p0[i], p1[i]);
	}
}

void recursive_decoder::leaf_simplify_1::downward_pass(std::vector<double>& L) {
	for (unsigned i = x; i < y; ++i) {
		L[i] = ext_l;
	}
}

#endif // ENABLE_OPT_1

#ifdef ENABLE_OPT_2

recursive_decoder::leaf_simplify_2::leaf_simplify_2(unsigned x, unsigned y, matrix const& Gp) : leaf(x, y, 1, Gp) {
	fail(Gp.size() == 1 && Gp[0] == ((1ull << (y - x)) - 1), "leaf_s2: inappropriate optimization");
}

void recursive_decoder::leaf_simplify_2::upward_pass(std::vector<double> const& p0, std::vector<double> const& p1) {
	A[0] = A[1] = 1.0;
	for (unsigned i = x; i < y; ++i) {
		A[0] *= p0[i];
		A[1] *= p1[i];
	}
}

void recursive_decoder::leaf_simplify_2::downward_pass(std::vector<double>& L) {
	ext_l = LLR(B[0] * A[0], B[1] * A[1]);
	for (unsigned i = x; i < y; ++i) {
		L[i] = ext_l;
	}
}

#endif // ENABLE_OPT_2

#ifdef ENABLE_OPT_3

recursive_decoder::leaf_simplify_3::leaf_simplify_3(unsigned x, unsigned y, unsigned k1, matrix const& Gp) : leaf(x, y, k1, Gp) {
	fail(y == x + 2, "leaf_s3: inappropriate optimization");
	fail(k1 + k2 == 2 && k1 != 0, "leaf_s3: inappropriate optimization");
}

void recursive_decoder::leaf_simplify_3::upward_pass(std::vector<double> const& p0, std::vector<double> const& p1) {
	Phi00 = p0[x] * p0[x + 1];
	Phi10 = p0[x + 1] - Phi00;
	Phi01 = p0[x] - Phi00;
	Phi11 = p1[x + 1] - Phi01;
	if (k1 == 1) {
		A[0] = Phi00 + Phi11;
		A[1] = Phi01 + Phi10;
		B[0] = B[1] = I;
	} else {
		A[0] = Phi00;
		A[1] = Phi10;
		A[2] = Phi01;
		A[3] = Phi11;
		B[0] = B[1] = B[2] = B[3] = I;
	}
}

void recursive_decoder::leaf_simplify_3::downward_pass(std::vector<double>& L) {
	if (k1 == 1) {
		double Phi00B0 = Phi00 * B[0];
		double Phi01B1 = Phi01 * B[1];
		double Phi10B1 = Phi10 * B[1];
		double Phi11B0 = Phi11 * B[0];
		L[x]     = LLR(Phi00B0 + Phi01B1, Phi11B0 + Phi10B1);
		L[x + 1] = LLR(Phi00B0 + Phi10B1, Phi11B0 + Phi01B1);
	} else {
		double Phi00B0 = Phi00 * B[0];
		double Phi10B1 = Phi10 * B[1];
		double Phi01B2 = Phi01 * B[2];
		double Phi11B3 = Phi11 * B[3];
		L[x]     = LLR(Phi00B0 + Phi01B2, Phi10B1 + Phi11B3);
		L[x + 1] = LLR(Phi00B0 + Phi10B1, Phi01B2 + Phi11B3);
	}
}

#endif // ENABLE_OPT_3

#ifdef ENABLE_OPT_5

recursive_decoder::leaf_simplify_5::leaf_simplify_5(unsigned x, unsigned y, matrix const& Gp) : leaf(x, y, y - x, Gp) {
	fail(Gp.size() == y - x && y - x >= 2, "leaf_s5: inappropriate optimization");
}

void recursive_decoder::leaf_simplify_5::upward_pass(std::vector<double> const& p0, std::vector<double> const& p1) {
	std::vector<double> diffs(k1);
	double T = 1;
	for (unsigned i = x, j = 0; i < y; ++i, ++j) {
		T *= p0[i];
		diffs[j] = p0[i] / p1[i];
	}
	std::size_t v = 0;
	for (unsigned ind : traverse_order) {
		if (ind != UNINIT) {
			if (getbit(v, ind)) {
				T *= diffs[ind];
			} else {
				T /= diffs[ind];
			}
			changebit(v, ind);
		}
		A[v] = T;
	}
}

void recursive_decoder::leaf_simplify_5::downward_pass(std::vector<double>& L_out) {
	for (unsigned i = 0; i < k1; ++i) {
		P0[i] = P1[i] = I;
	}
	for (std::size_t v = 1; v < (1ull << k1); ++v) {
		for (unsigned i = 0; i < k1; ++i) {
			if ((v >> i) & 1) {
				P1[i] += A[v] * B[v];
			} else {
				P0[i] += A[v] * B[v];
			}
		}
	}
	for (unsigned i = x, j = 0; i < y; ++i, ++j) {
		L_out[i] = LLR(P0[j], P1[j]);
	}
}

#endif // ENABLE_OPT_5

recursive_decoder::inner::inner(unsigned x, unsigned y, unsigned z, node* left, node* right, unsigned k1, unsigned k2, unsigned k3, matrix const& G_hat, matrix const& G_tilda, matrix const& Gp)
	: node(x, y, k1, k2, k3, Gp), z(z), left(left), right(right), G_hat(G_hat), G_tilda(G_tilda) {
	fail(k3 == left->k3 + right->k3 + k2, "inner-ctor: incorrect dims");
	fail(x < z && z < y, "inner-ctor: z must be in (x,y)");
	fail(left != nullptr && right != nullptr, "innert-ctor: empty left or right");
}

void recursive_decoder::inner::upward_pass(std::vector<double> const& p0, std::vector<double> const& p1) {
	left->upward_pass(p0, p1);
	right->upward_pass(p0, p1);
	if (k1 != 0) {
		std::size_t a = 0, b = 0, v = 0;
		for (unsigned ind : traverse_order) {
			if (ind != UNINIT) {
				if (ind >= k2) {
					changebit(v, ind - k2);
					A[v] = B[v] = I;
				}
				a ^= G_hat[ind];
				b ^= G_tilda[ind];
			} else {
				A[v] = B[v] = I;
			}
			A[v] += left->A[a] * right->A[b];
		}
	}
}

void recursive_decoder::inner::downward_pass(std::vector<double>& L) {
	std::size_t a = 0, b = 0;
	if (k1 == 0) {
		for (unsigned ind : traverse_order) {
			if (ind != UNINIT) {
				a ^= G_hat[ind];
				b ^= G_tilda[ind];
			}
			left->B[a] = right->A[b];
			right->B[b] = left->A[a];
		}
	} else {
		std::size_t v = 0;
		for (unsigned ind : traverse_order) {
			if (ind != UNINIT) {
				if (ind >= k2) {
					changebit(v, ind - k2);
				}
				a ^= G_hat[ind];
				b ^= G_tilda[ind];
			}
			left->B[a] += B[v] * right->A[b];
			right->B[b] += B[v] * left->A[a];
		}
	}
	left->downward_pass(L);
	right->downward_pass(L);
}

void recursive_decoder::inner::__print(std::ostream& out) const {
	out << "\t" << getName() << " -> " << left->getName() << "\n";
	out << "\t" << getName() << " -> " << right->getName() << "\n";
	left->print(out);
	right->print(out);
}

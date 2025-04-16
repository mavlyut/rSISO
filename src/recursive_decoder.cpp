#include "../include/gray_code.h"
#include "../include/recursive_decoder.h"

recursive_decoder::recursive_decoder(matrix const& _G) : linear_soft_decoder(_G), G(_G), L_out(n) {
	time_measure(minimal_span_form(G));

	__log("MSF:\n" << G << "\n")

	p0.resize(n);
	p1.resize(n);
	_split = std::vector(n + 1, std::vector<unsigned>(n + 1, UNINIT));
	_ctor = std::vector(n + 1, std::vector<unsigned>(n + 1, 0));
	if (k > 1) {
		time_measure(sectionalization());
		time_measure(root = rec_build_section_tree(0, n));
	} else {
		root = new leaf_simplify_2(0, n, G);
	}
}

void recursive_decoder::sectionalization() {
	std::vector<std::vector<unsigned>> _Gp_size(n + 1, std::vector<unsigned>(n + 1, 0));
	std::vector<std::vector<unsigned>> _Gs_size(n + 1, std::vector<unsigned>(n + 1, 0));
	std::vector<std::vector<std::vector<std::size_t>>> phi_u_c(n + 1, std::vector<std::vector<std::size_t>>(n + 1, std::vector<std::size_t>(n + 1, UNINIT)));
	std::vector<std::vector<std::vector<std::size_t>>> phi_d_c(n + 1, std::vector<std::vector<std::size_t>>(n + 1, std::vector<std::size_t>(n + 1, UNINIT)));
	std::vector<std::vector<std::size_t>> phi_u_l(n + 1, std::vector<std::size_t>(n + 1, UNINIT));
	std::vector<std::vector<std::size_t>> phi_d_l(n + 1, std::vector<std::size_t>(n + 1, UNINIT));
	std::vector<std::vector<unsigned>> phi(n + 1, std::vector<unsigned>(n + 1, UNINIT));

	std::function<void(std::vector<std::vector<std::size_t>>&, unsigned, unsigned, std::size_t const&, unsigned)> update_phi
	= [&](std::vector<std::vector<std::size_t>>& __phi, unsigned x, unsigned y, std::size_t const& val, unsigned ctorNo) {
		if (val < __phi[x][y]) {
			__phi[x][y] = val;
			_ctor[x][y] = ctorNo;
		}
	};

	std::function<void(std::size_t&, std::size_t const&)> update_uninit = [](std::size_t& x, std::size_t const& y) {
		if (x == UNINIT) {
			x = y;
		}
	};

	for (unsigned x = 0; x <= n; x++) {
		for (unsigned y = x + 1; y <= n; y++) {
			matrix Gp;
			int k3 = 0;
			for (binvector const& row : G) {
				Gp.push_back(row.subvector(x, y));
				if (row.subvector(0, x).isZero() && row.subvector(y, n).isZero()) {
					k3++;
				}
			}
			gauss(Gp);
			_Gp_size[x][y] = Gp.size();
			_Gs_size[x][y] = k3;

			#ifdef ENABLE_OPT_1
			if (Gp.size() == 1 && G[0].isOnes() && k3 == 1) {
				update_phi(phi_u_l, x, y, 2 * (y - x), 1);
				update_phi(phi_d_l, x, y, 0, 1);
			}
			#endif

			#ifdef ENABLE_OPT_2
			if (Gp.size() == 1 && G[0].isOnes() && k3 == 0) {
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

	for (unsigned x = 0; x <= n; x++) {
		for (unsigned y = x + 1; y <= n; y++) {
			update_phi(phi_u_l, x, y, (1ull << _Gp_size[x][y]) * (2 * (y - x) - 1), 0);
			update_phi(phi_d_l, x, y, ((1ull << (_Gp_size[x][y] - _Gs_size[x][y] + 1))) * (y - x + 1) + 2 * (y - x), 0);
			for (unsigned z = x + 1; z < y - 1; z++) {
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
			for (unsigned z = x + 1; z < y - 1; z++) {
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

	__log("Gp_sz:  " << _Gp_size << "\nGs_sz:  " << _Gs_size << "\n");
	__log("split: " << _split << "\n");
}

recursive_decoder::node* recursive_decoder::rec_build_section_tree(unsigned x, unsigned y) {
	fail(x < y && y <= n, "rec_build: incorrect bounds");
	unsigned z = _split[x][y];
	if (z == UNINIT) {
		matrix Gp, Gs, Gl;
		for (unsigned i = 0; i < k; i++) {
			binvector row = G[i];
			binvector append = row.subvector(x, y);
			if (row.subvector(0, x).isZero() && row.subvector(y, n).isZero()) {
				Gs.push_back(append);
				Gp.push_back(append);
			} else {
				Gl.push_back(append);
			}
		}
		unsigned k1 = 0;
		full_gauss(Gl);
		for (binvector const& row : Gl) {
			if (lin_indep(Gp, row)) {
				Gp.push_back(row);
				k1++;
			}
		}

		unsigned optimizationNo = _ctor[x][y];
		#ifdef ENABLE_OPT_1
		if (optimizationNo == 1) {
			fail(k1 == 0 && Gp.size() == 1 && Gp[0].isOnes(), "incorrect pre-conditions (1)");
			__log("Simplify 1: " << x << " " << y << std::endl);
			return new leaf_simplify_1(x, y, Gp);
		}
		#endif

		#ifdef ENABLE_OPT_2
		if (optimizationNo == 2) {
			fail(k1 == 1 && Gp.size() == 1 && Gp[0].isOnes(), "incorrect pre-conditions (2)");
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
	for (unsigned i = 0; i < k; i++) {
		binvector row = G[i];
		binvector append = row.subvector(x, y);
		if (row.subvector(0, x).isZero() && row.subvector(y, n).isZero()) {
			if (row.subvector(x, z).isZero() || row.subvector(z, y).isZero()) {
				Gp.push_back(append);
			} else {
				_G0.push_back(append);
			}
			k3++;
		} else {
			_G1.push_back(append);
		}
	}

	matrix G_0, G_1;
	unsigned k1 = 0, k2 = 0;
	for (binvector const& row : _G0) {
		if (lin_indep(Gp, row)) {
			Gp.push_back(row);
			G_0.push_back(row.subvector(0, z - x));
			G_1.push_back(row.subvector(z - x, y - x));
			k2++;
		}
	}
	for (binvector const& row : _G1) {
		if (lin_indep(Gp, row)) {
			Gp.push_back(row);
			G_0.push_back(row.subvector(0, z - x));
			G_1.push_back(row.subvector(z - x, y - x));
			k1++;
		}
	}

	__log("inner: " << x << " " << y << "\nGp:\n" << Gp << "\nG_0:\n" << G_0 << "\nG_1:\n" << G_1 << "\n");

	node* left = rec_build_section_tree(x, z);
	node* right = rec_build_section_tree(z, y);

	matrix G_hat(k1 + k2, binvector(left->k1));
	matrix G_tilda(k1 + k2, binvector(right->k1));
	{
		matrix G_hat_ext(z - x, binvector(left->k1 + left->k3 + k1 + k2));
		for (unsigned i = 0; i < z - x; i++) {
			unsigned j = 0;
			for (unsigned t = 0; t < left->Gp.size(); t++, j++) {
				G_hat_ext[i].set(j, left->Gp[t][i]);
			}
			for (unsigned t = 0; t < G_0.size(); t++, j++) {
				G_hat_ext[i].set(j, G_0[t][i]);
			}
		}
		full_gauss(G_hat_ext);

		__log("G_hat_ext:\n" << G_hat_ext << "\n");

		for (unsigned i = 0; i < k1 + k2; i++) {
			for (unsigned j = 0; j < left->k1; j++) {
				G_hat[i].set(left->k1 - j - 1, G_hat_ext[G_hat_ext.size() - j - 1][i + left->k1 + left->k3]);
			}
		}
	}
	{
		matrix G_tilda_ext(y - z, binvector(right->k1 + right->k3 + k1 + k2));
		for (unsigned i = 0; i < y - z; i++) {
			unsigned j = 0;
			for (unsigned t = 0; t < right->Gp.size(); t++, j++) {
				G_tilda_ext[i].set(j, right->Gp[t][i]);
			}
			for (unsigned t = 0; t < G_1.size(); t++, j++) {
				G_tilda_ext[i].set(j, G_1[t][i]);
			}
		}

		__log("G_tilda_ext (before Gauss):\n" << G_tilda_ext << "\n");

		full_gauss(G_tilda_ext);

		__log("G_tilda_ext:\n" << G_tilda_ext << "\n");

		for (unsigned i = 0; i < k1 + k2; i++) {
			for (unsigned j = 0; j < right->k1; j++) {
				G_tilda[i].set(right->k1 - j - 1, G_tilda_ext[G_tilda_ext.size() - j - 1][i + right->k1 + right->k3]);
			}
		}
	}

	__log("G_hat:\n" << G_hat << "\nG_tilda:\n" << G_tilda << "\n");

	return new inner(x, y, z, left, right, k1, k2, k3, G_hat, G_tilda, Gp);
}

std::vector<double> recursive_decoder::decode_soft(std::vector<_Float64> const& L_in) {
	for (unsigned i = 0; i < n; i++) {
		_Float64 L_exp = exp(truncate(L_in[i]));
		_Float64 z = 1.0 + L_exp;
		p0[i] = L_exp / z;
		p1[i] = 1.0 / z;
	}

	root->clear();
	time_measure(root->upward_pass(p0, p1));
	time_measure(root->downward_pass(L_out));

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
	: x(x), y(y), k1(k1), k2(k2), k3(k3), A(1ull << k1, I), B(1ull << k1, I), Gp(Gp) {
	__log("Node ctor: x=" << x << ", y=" << y << "; k=(" << k1 << ", " << k2 << ", " << k3 << ")\n");
}

void recursive_decoder::node::clear() {
	__clear();
	for (double& i : A) {
		i = I;
	}
	for (double& i : B) {
		i = I;
	}
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

double recursive_decoder::leaf::F(binvector const& c, std::vector<double> const& p0, std::vector<double> const& p1) const {
	double ans = 1.0;
	for (unsigned i = x, j = 0; i < y; i++, j++) {
		ans *= (c[j] ? p1[i] : p0[i]);
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
	: leaf(x, y, k1, Gp), A1(y - x, std::vector<double>(1ull << k1, I)) {}

void recursive_decoder::leaf_no_simplify::__clear() {
	for (auto& i : A1) {
		for (double& j : i) {
			j = I;
		}
	}
}

void recursive_decoder::leaf_no_simplify::upward_pass(std::vector<double> const& p0, std::vector<double> const& p1) {
	binvector v(k1), c(y - x);
	double T;
	for (unsigned ind : gray_code(k1 + k2)) {
		if (ind != UNINIT) {
			if (ind >= k2) {
				v.change(ind - k2);
			}
			c ^= Gp[ind];
		}
		T = F(c, p0, p1);
		A[v] += T;
		for (unsigned i = 0; i < y - x; i++) {
			if (c[i]) {
				A1[i][v] += T;
			}
		}
	}
}

void recursive_decoder::leaf_no_simplify::downward_pass(std::vector<double>& L) const {
	std::vector<double> P1(y - x, 0.0);
	double P_summary = 0.0;
	for (std::size_t v = 0; v < (1ull << k1); v++) {
		P_summary += A[v] * B[v];
		for (unsigned i = 0; i < y - x; i++) {
			P1[i] += A1[i][v] * B[v];
		}
	}
	for (unsigned i = x, j = 0; i < y; i++, j++) {
		L[i] = LLR((P_summary - P1[j]), P1[j]);
	}
}

#ifdef ENABLE_OPT_1

recursive_decoder::leaf_simplify_1::leaf_simplify_1(unsigned x, unsigned y, matrix const& Gp) : leaf(x, y, 0, Gp) {
	fail(Gp.size() == 1 && Gp[0].isOnes(), "leaf_s1: inappropriate optimization");
}

void recursive_decoder::leaf_simplify_1::upward_pass(std::vector<double> const& p0, std::vector<double> const& p1) {
	for (unsigned i = x; i < y; i++) {
		ext_l += LLR(p0[i], p1[i]);
	}
}

void recursive_decoder::leaf_simplify_1::downward_pass(std::vector<double>& L) const {
	for (unsigned i = x; i < y; i++) {
		L[i] = ext_l;
	}
}

void recursive_decoder::leaf_simplify_1::__clear() {
	ext_l = I;
}

#endif // ENABLE_OPT_1

#ifdef ENABLE_OPT_2

recursive_decoder::leaf_simplify_2::leaf_simplify_2(unsigned x, unsigned y, matrix const& Gp) : leaf(x, y, 1, Gp) {
	fail(Gp.size() == 1 && Gp[0].isOnes(), "leaf_s2: inappropriate optimization");
}

void recursive_decoder::leaf_simplify_2::upward_pass(std::vector<double> const& p0, std::vector<double> const& p1) {
	for (unsigned i = x; i < y; i++) {
		A[0] *= p0[i];
		A[1] *= p1[i];
	}
}

void recursive_decoder::leaf_simplify_2::downward_pass(std::vector<double>& L) const {
	double ext_l = LLR(B[0] * A[0], B[1] * A[1]);
	for (unsigned i = x; i < y; i++) {
		L[i] = ext_l;
	}
}

void recursive_decoder::leaf_simplify_2::__clear() {
	A[0] = 1.0;
	A[1] = 1.0;
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
	} else {
		A[0] = Phi00;
		A[1] = Phi10;
		A[2] = Phi01;
		A[3] = Phi11;
	}
}

void recursive_decoder::leaf_simplify_3::downward_pass(std::vector<double>& L) const {
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

void recursive_decoder::leaf_simplify_3::__clear() {}

#endif // ENABLE_OPT_3

#ifdef ENABLE_OPT_5

recursive_decoder::leaf_simplify_5::leaf_simplify_5(unsigned x, unsigned y, matrix const& Gp) : leaf(x, y, y - x, Gp) {
	fail(Gp.size() == y - x && y - x >= 2, "leaf_s5: inappropriate optimization");
}

void recursive_decoder::leaf_simplify_5::upward_pass(std::vector<double> const& p0, std::vector<double> const& p1) {
	std::vector<double> diffs(k1);
	double T = 1;
	for (unsigned i = x, j = 0; i < y; i++, j++) {
		T *= p0[i];
		diffs[j] = p0[i] / p1[i];
	}
	binvector v(k1);
	for (unsigned ind : gray_code(k1)) {
		if (ind != UNINIT) {
			if (v[ind]) {
				T *= diffs[ind];
			} else {
				T /= diffs[ind];
			}
			v.change(ind);
		}
		A[v] = T;
	}
}

void recursive_decoder::leaf_simplify_5::downward_pass(std::vector<double>& L_out) const {
	std::vector<double> P1(k1, 0.0);
	double P_summary = 0;
	for (std::size_t v = 0; v < (1ull << k1); v++) {
		P_summary += A[v] * B[v];
		for (unsigned i = 0; i < k1; i++) {
			if ((v >> i) & 1) {
				P1[i] += A[v] * B[v];
			}
		}
	}
	for (unsigned i = x, j = 0; i < y; i++, j++) {
		L_out[i] = LLR((P_summary - P1[j]), P1[j]);
	}
}

void recursive_decoder::leaf_simplify_5::__clear() {}

#endif // ENABLE_OPT_5

recursive_decoder::inner::inner(unsigned x, unsigned y, unsigned z, node* left, node* right, unsigned k1, unsigned k2, unsigned k3, matrix const& G_hat, matrix const& G_tilda, matrix const& Gp)
	: node(x, y, k1, k2, k3, Gp), z(z), left(left), right(right), G_hat(G_hat), G_tilda(G_tilda) {
	fail(k3 == left->k3 + right->k3 + k2, "inner-ctor: incorrect dims");
	fail(x < z && z < y, "inner-ctor: z must be in (x,y)");
	fail(left != nullptr && right != nullptr, "innert-ctor: empty left or right");
	fail(G_hat.size() == k1 + k2 && G_hat.front().size() == left->k1, "inner-ctor: incorrecy size of \\hat{G}");
	fail(G_tilda.size() == k1 + k2 && G_tilda.front().size() == right->k1, "inner-ctor: incorrecy size of \\tilda{G}");
}

void recursive_decoder::inner::upward_pass(std::vector<double> const& p0, std::vector<double> const& p1) {
	left->upward_pass(p0, p1);
	right->upward_pass(p0, p1);
	if (k1 != 0) {
		binvector a(left->k1), b(right->k1), v(k1);
		for (unsigned ind : gray_code(k1 + k2)) {
			if (ind != UNINIT) {
				if (ind >= k2) {
					v.change(ind - k2);
				}
				a ^= G_hat[ind];
				b ^= G_tilda[ind];
			}
			A[v] += left->A[a] * right->A[b];
		}
	}
}

void recursive_decoder::inner::downward_pass(std::vector<double>& L) const {
	binvector a(left->k1), b(right->k1);
	if (k1 == 0) {
		for (unsigned ind : gray_code(k1 + k2)) {
			if (ind != UNINIT) {
				a ^= G_hat[ind];
				b ^= G_tilda[ind];
			}
			left->B[a] = right->A[b];
			right->B[b] = left->A[a];
		}
	} else {
		binvector v(k1);
		for (unsigned ind : gray_code(k1 + k2)) {
			if (ind != UNINIT) {
				if (ind >= k2) {
					v.change(ind - k2);
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

void recursive_decoder::inner::__clear() {
	left->clear();
	right->clear();
}

void recursive_decoder::inner::__print(std::ostream& out) const {
	out << "\t" << getName() << " -> " << left->getName() << "\n";
	out << "\t" << getName() << " -> " << right->getName() << "\n";
	left->print(out);
	right->print(out);
}

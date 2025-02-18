#include <algorithm>
#include <array>
#include "binpoly.cpp"
#include <chrono>
#include <cmath>
#include <ctime>
#include <fstream>
#include <functional>
#include <iostream>
#include <iomanip>
#include <map>
#include <queue>
#include <random>
#include <set>
#include <unordered_map>
#include <vector>

#define fail(b, msg)						\
	if (!(b)) { 							\
		std::ofstream fout("output.txt");	\
		fout << "FAILURE: " << msg << "\n";	\
		fout.close();						\
		throw 2;							\
	}

#define LOG
#define TIMELOG
#define TEST
#define EPS 1e-10
static const double INF = INFINITY;
static constexpr double COEF = M_2_SQRTPI * M_SQRT1_2 / 2;
static std::mt19937 gen(time(0));
std::ofstream _log_time("_log_time");
std::ofstream _log_out("_log_out");

typedef std::chrono::milliseconds ms;
#define __time_measure_with_msg(msg, line)									\
	auto start = std::chrono::system_clock::now();							\
	line;																	\
	auto end = std::chrono::system_clock::now();							\
	auto time_in_ms = std::chrono::duration_cast<ms>(end - start).count();	\
	_log_time << msg << ": "												\
			  << (time_in_ms > 1000 ? (time_in_ms / 1000.0) : time_in_ms)	\
			  << (time_in_ms > 1000 ? "s" : "ms") << std::endl;
#define __time_measure(line) \
	__time_measure_with_msg("Instruction: \"" #line "\"", line)

#if defined(TIMELOG) && !defined(TEST)
#define time_measure(line) { __time_measure(line) }
#else
#define time_measure(line) line
#endif

#if defined(LOG) && !defined(TEST)
#define __log(msg) _log_out << msg;
#else
#define __log(msg)
#endif

using matrix = std::vector<binvector>;

template <typename T>
std::ostream& operator<<(std::ostream& out, std::vector<T> const& a) {
	for (int i = 0; i < a.size(); i++) {
		if (i != 0) {
			out << ' ';
		}
		out << a[i];
	}
	return out;
}

template <typename T>
std::ostream& operator<<(std::ostream& out, std::vector<std::vector<T>> const& a) {
	for (int i = 0; i < a.size(); i++) {
		if (i != 0) {
			out << "\n\t";
		}
		out << a[i];
	}
	return out;
}

std::ostream& operator<<(std::ostream& out, matrix const& a) {
	for (binvector const& row : a) {
			for (int i = 0; i < row.size(); i++) {
			if (i != 0) {
				out << ' ';
			}
			out << row[i];
		}
		out << '\n';
	}
	return out;
}

template <typename T>
std::istream& operator>>(std::istream& in, std::vector<T>& a) {
	for (T& i : a) {
		in >> i;
	}
	return in;
}

binvector operator*(binvector const& x, matrix const& A) {
	fail(x.size() == A.size(), "mul: incompatible sizes");
	binvector ans(A.front().size());
	for (int i = 0; i < ans.size(); i++) {
		bool val = false;
		for (int j = 0; j < x.size(); j++) {
			val ^= (x[j] & A[j][i]);
		}
		ans.set(i, val);
	}
	return ans;
}

std::vector<int> gauss(matrix& a) {
	int k = a.size();
	int n = a[0].size();
	int last_c = 0;
	std::vector<int> gamma;
	for (int i = 0; i < k; ++i) {
		while (last_c < n) {
			bool fl = false;
			for (int j = i; j < k; ++j) {
				if (a[j][last_c] == 1) {
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
			for (int ii = i; ii < k; ii++) {
				a.pop_back();
			}
			break;
		}
		gamma.push_back(last_c);
		for (int j = i + 1; j < k; ++j) {
			if (a[j][last_c] == 1) {
				for (int t = 0; t < n; ++t) {
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
	for (int i = gamma.size() - 1; i >= 0; i--) {
		for (int j = 0; j < i; j++) {
			if (a[j][gamma[i]]) {
				a[j] = a[j] + a[i];
			}
		}
	}
}

std::vector<int> minimal_span_form(matrix& G) {
	std::vector<int> gamma = gauss(G);
	int n = G[0].size(), k = G.size();
	for (int i = k - 1; i >= 0; i--) {
		int j = n - 1;
		while (j >= 0 && G[i][j] == 0) {
			j--;
		}
		for (int r = 0; r < i; r++) {
			if (G[r][j] == 1) {
				for (int c = 0; c < n; c++) {
					G[r].set(c, (G[r][c] + G[i][c]) % 2);
				}
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

matrix operator*(matrix const& A, matrix const& B) {
	std::cout << "mul:\n" << A << "*\n" << B << "\n";
	fail(A.front().size() == B.size(), "mulMatrix: incompatible sizes");
	matrix C(A.size(), binvector(B.front().size()));
	for (int i = 0; i < C.size(); i++) {
		for (int j = 0; j < C.front().size(); j++) {
			for (int k = 0; k < B.size(); k++) {
				C[i].set(j, C[i][j] ^ (A[i][k] & B[k][j]));
			}
		}
	}
	return C;
}

class RecursiveDecoder {
public:
	RecursiveDecoder(matrix const& _G) : G(_G) {
		time_measure(minimal_span_form(G));

		__log("MSF:\n" << G << "\n")
		
		k = G.size();
		n = G.front().size();
		_split = std::vector(n + 1, std::vector(n + 1, -1));
		time_measure(root = sectionalization());
	}

private: // Sectionalization
	class node {
	protected:
		node(int x, int y, int k1, int k2, int k3, matrix const& Gp)
			: x(x), y(y), k1(k1), k2(k2), k3(k3), A(1 << k1, I), B(1 << k1, I), Gp(Gp) {}

	public:
		virtual bool is_leaf() const = 0;

		void clear() {
			_clear();
			for (double& i : A) {
				i = I;
			}
			for (double& i : B) {
				i = I;
			}
		}

	public:
		const int x, y;
		const int k1, k2, k3;
		std::vector<double> A, B;
		matrix Gp;

	protected:
		static constexpr double I = 0;

		// (u v) * M
		binvector combineAndMul(binvector const& u, binvector const& v, matrix const& M) const {
			fail(u.size() + v.size() == M.size(), "combineAndMul: incompatible sizes");
			binvector ans(M.front().size());
			for (int i = 0; i < ans.size(); i++) {
				bool val = false;
				int t = 0;
				for (int j = 0; j < u.size(); j++, t++) {
					val ^= (u[j] & M[t][i]);
				}
				for (int j = 0; j < v.size(); j++, t++) {
					val ^= (v[j] & M[t][i]);
				}
				ans.set(i, val);
			}
			return ans;
		}

		virtual void _clear() = 0;
	};

	class leaf : public node {
	public:
		leaf(int x, int y, int k1, matrix const& Gp)
			: node(x, y, k1, Gp.size() - k1, Gp.size() - k1, Gp)
			, A0(y - x, std::vector<double>(1 << k1, I))
			, A1(y - x, std::vector<double>(1 << k1, I)) {}

		bool is_leaf() const override {
			return true;
		}

		binvector dot(binvector const& u, binvector const& v) const {
			return combineAndMul(u, v, Gp);
		}

	public:
		std::vector<std::vector<double>> A0, A1;
	
	protected:
		void _clear() {
			for (auto& i : A0) {
				for (double& j : i) {
					j = I;
				}
			}
			for (auto& i : A1) {
				for (double& j : i) {
					j = I;
				}
			}
		}
	};

	class inner : public node {
	public:
		inner(int x, int y, int z, node* left, node* right, int k1, int k2, int k3, matrix const& G_hat, matrix const& G_tilda, matrix const& Gp)
			: node(x, y, k1, k2, k3, Gp), z(z), left(left), right(right), G_hat(G_hat), G_tilda(G_tilda) {
			fail(k3 == left->k3 + right->k3 + k2, "inner-ctor: incorrect dims");
			fail(x < z && z < y, "inner-ctor: z must be in (x,y)");
			fail(left != nullptr && right != nullptr, "innert-ctor: empty left or right");
			fail(G_hat.size() == k1 + k2 && G_hat.front().size() == left->k1, "inner-ctor: incorrecy size of \\hat{G}");
			fail(G_tilda.size() == k1 + k2 && G_tilda.front().size() == right->k1, "inner-ctor: incorrecy size of \\tilda{G}");
		}

		bool is_leaf() const override {
			return false;
		}

		std::pair<binvector, binvector> mapping(binvector const& u, binvector const& v) const {
			return {combineAndMul(u, v, G_hat), combineAndMul(u, v, G_tilda)};
		}

	public:
		int z;
		node *left, *right;

	private:
		matrix G_hat, G_tilda;
	
	protected:
		void _clear() {
			left->clear();
			right->clear();
		}
	};

	node* sectionalization() {
		std::vector<std::vector<int>> _Gp_size(n + 1, std::vector<int>(n + 1, 0));
		std::vector<std::vector<int>> _Gs_size(n + 1, std::vector<int>(n + 1, 0));
		std::vector<std::vector<std::vector<std::size_t>>> phi_u_c(n + 1, std::vector<std::vector<std::size_t>>(n + 1, std::vector<std::size_t>(n + 1, 0)));
		std::vector<std::vector<std::vector<std::size_t>>> phi_d_c(n + 1, std::vector<std::vector<std::size_t>>(n + 1, std::vector<std::size_t>(n + 1, 0)));
		std::vector<std::vector<std::size_t>> phi_u_l(n + 1, std::vector<std::size_t>(n + 1, 0));
		std::vector<std::vector<std::size_t>> phi_d_l(n + 1, std::vector<std::size_t>(n + 1, 0));
		std::vector<std::vector<std::size_t>> phi(n + 1, std::vector<std::size_t>(n + 1, -1));

		for (int x = 0; x <= n; x++) {
			for (int y = x + 1; y <= n; y++) {
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
			}
		}

		for (int x = 0; x <= n; x++) {
			for (int y = x + 1; y <= n; y++) {
				phi_u_l[x][y] = (1ull << _Gp_size[x][y]) * (y - x - 1);
				phi_d_l[x][y] = ((1ull << (_Gp_size[x][y] - _Gs_size[x][y] + 1)) + 1) * (y - x);
				for (int z = x + 1; z < y - 1; z++) {
					int k1_k2 = _Gp_size[x][y] - _Gs_size[x][z] - _Gs_size[z][y];
					phi_u_c[x][y][z] = (1ull << k1_k2);
					phi_d_c[x][y][z] = (1ull << (k1_k2 + 1));
				}
			}
		}

		std::function<int(int, int)> get_phi = [&](int x, int y) -> int {
			if (phi[x][y] == -1) {
				std::size_t min_score = phi_u_l[x][y] + phi_d_l[x][y], ind_min = -1;
				for (int z = x + 1; z < y - 1; z++) {
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
		for (int x = 0; x <= n; x++) {
			for (int y = x + 1; y <= n; y++) {
				get_phi(x, y);
			}
		}

		__log("Gp_sz:  " << _Gp_size << "\nGs_sz:  " << _Gs_size << "\n");
		__log("split: " << _split << "\n");

		std::vector<int> sections = {0, n};
		std::function<void(int, int)> rec_log = [&](int x, int y) {
			int z = _split[x][y];

			__log("[" << x << ";" << y << ") at " << z << "\n");

			if (z == -1) {
				return;
			}
			sections.push_back(z);
			rec_log(x, z);
			rec_log(z, y);
		};
		rec_log(0, n);
		std::sort(sections.begin(), sections.end());
		std::cout << sections << "\n\n";

		return rec_build_section_tree(0, n);
	}

	node* rec_build_section_tree(int x, int y) {
		fail(0 <= x && x < y && y <= n, "rec_build: incorrect bounds");
		int z = _split[x][y];
		if (z == -1) {
			matrix Gp, Gs, Gl;
			for (int i = 0; i < k; i++) {
				binvector row = G[i];
				binvector append = row.subvector(x, y);
				if (row.subvector(0, x).isZero() && row.subvector(y, n).isZero()) {
					Gs.push_back(append);
					Gp.push_back(append);
				} else {
					Gl.push_back(append);
				}
			}
			int k1 = 0;
			for (binvector const& row : Gl) {
				if (lin_indep(Gp, row)) {
					Gp.push_back(row);
					k1++;
				}
			}

			__log("leaf: " << x << " " << y << "\nGp:\n" << Gp << "\n");

			return new leaf(x, y, k1, Gp);
		}

		matrix Gp, _G0, _G1; int k3 = 0;
		for (int i = 0; i < k; i++) {
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
		int k1 = 0, k2 = 0;
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
			for (int i = 0; i < z - x; i++) {
				int j = 0;
				for (int t = 0; t < left->Gp.size(); t++, j++) {
					G_hat_ext[i].set(j, left->Gp[t][i]);
				}
				for (int t = 0; t < G_0.size(); t++, j++) {
					G_hat_ext[i].set(j, G_0[t][i]);
				}
			}
			full_gauss(G_hat_ext);

			__log("G_hat_ext:\n" << G_hat_ext << "\n");

			for (int i = 0; i < k1 + k2; i++) {
				for (int j = 0; j < left->k1; j++) {
					G_hat[i].set(left->k1 - j - 1, G_hat_ext[G_hat_ext.size() - j - 1][i + left->k1 + left->k3]);
				}
			}
		}
		{
			matrix G_tilda_ext(y - z, binvector(right->k1 + right->k3 + k1 + k2));
			for (int i = 0; i < y - z; i++) {
				int j = 0;
				for (int t = 0; t < right->Gp.size(); t++, j++) {
					G_tilda_ext[i].set(j, right->Gp[t][i]);
				}
				for (int t = 0; t < G_1.size(); t++, j++) {
					G_tilda_ext[i].set(j, G_1[t][i]);
				}
			}
			
			__log("G_tilda_ext (before Gauss):\n" << G_tilda_ext << "\n");

			full_gauss(G_tilda_ext);

			__log("G_tilda_ext:\n" << G_tilda_ext << "\n");

			for (int i = 0; i < k1 + k2; i++) {
				for (int j = 0; j < right->k1; j++) {
					G_tilda[i].set(right->k1 - j - 1, G_tilda_ext[G_tilda_ext.size() - j - 1][i + right->k1 + right->k3]);
				}
			}
		}

		__log("G_hat:\n" << G_hat << "\nG_tilda:\n" << G_tilda << "\n");

		return new inner(x, y, z, left, right, k1, k2, k3, G_hat, G_tilda, Gp);
	}

public:
	binvector encode(binvector const& c) const {
		fail(c.size() == k, "encode: size is incorrect");
		return c * G;
	}

public:
	std::vector<double> decode_soft(std::vector<double> const& y) const {
		root->clear();
		std::vector<double> L(length(), 0);
		time_measure(upward_pass(root, y));
		time_measure(downward_pass(root, L));

		printLog(root);

		return L;
	}

	binvector decode(std::vector<double> const& y) const {
		auto L = decode_soft(y);
		binvector ans(length());
		for (int i = 0; i < n; i++) {
			ans.set(i, L[i] < 0);
		}
		return ans;
	}

private:
	void printLog(node* rt) const {
		__log("[" << rt->x << ", " << rt->y << ")\n");
		__log("k = " << rt->k1 << " " << rt->k2 << " " << rt->k3 << "\n");
		__log("Gp:\n" << rt->Gp << "\n");
		if (rt->is_leaf()) {
			leaf* nd = reinterpret_cast<leaf*>(rt);
			__log("A/0: " << nd->A0 << "\nA/1: " << nd->A1
				<< "\nA: " << nd->A << "\nB: " << nd->B << "\n");
			return;
		}
		inner* nd = reinterpret_cast<inner*>(rt);
		__log("A: " << nd->A << "\nB: " << nd->B << "\n");
		printLog(nd->left);
		printLog(nd->right);
	}

	void upward_pass(node* rt, std::vector<double> const& R) const {
		if (rt->is_leaf()) {
			upward_pass_non_rec(reinterpret_cast<leaf*>(rt), R);
			return;
		}
		inner* nd = reinterpret_cast<inner*>(rt);
		upward_pass(nd->left, R);
		upward_pass(nd->right, R);
		if (nd->k1 != 0) {
			binvector w_zero(nd->k2);
			for (std::size_t v = 0; v < (1ull << nd->k1); v++) {
				binvector vv(nd->k1, v);
				for (std::size_t w = 0; w < (1ull << nd->k2); w++) {
					binvector ww(nd->k2, w);
					auto [a, b] = nd->mapping(ww, vv);
					nd->A[v] += nd->left->A[a] * nd->right->A[b];
				}
			}
		}
	}

	void downward_pass(node* rt, std::vector<double>& L) const {
		if (rt->is_leaf()) {
			downward_pass_non_rec(reinterpret_cast<leaf*>(rt), L);
			return;
		}
		inner* nd = reinterpret_cast<inner*>(rt);
		if (nd->k1 == 0) {
			binvector vv_zero(0);
			for (std::size_t w = 0; w < (1ull << nd->k2); w++) {
				binvector vv(nd->k2, w);
				auto [a, b] = nd->mapping(binvector(nd->k2, w), vv_zero);
				nd->left->B[a] = nd->right->A[b];
				nd->right->B[b] = nd->left->A[a];
			}
		} else {
			binvector w_zero(nd->k2);
			for (std::size_t v = 0; v < (1ull << nd->k1); v++) {
				binvector vv(nd->k1, v);
				for (std::size_t w = 0; w < (1ull << nd->k2); w++) {
					binvector ww(nd->k2, w);
					auto [a, b] = nd->mapping(ww, vv);
					nd->left->B[a] += nd->B[v] * nd->right->A[b];
					nd->right->B[b] += nd->B[v] * nd->left->A[a];
				}
			}
		}
		downward_pass(nd->left, L);
		downward_pass(nd->right, L);
	}

	void upward_pass_non_rec(leaf* nd, std::vector<double> const& R) const {
		for (std::size_t v = 0; v < (1ull << nd->k1); v++) {
			binvector vv(nd->k1, v);
			for (std::size_t w = 0; w < (1 << nd->k2); w++) {
				binvector ww(nd->k2, w);
				binvector c = nd->dot(ww, vv);
				double T = F(c, R, nd->x, nd->y);
				nd->A[v] += T;
				for (int i = 0; i < nd->y - nd->x; i++) {
					if (c[i]) {
						nd->A1[i][v] += T;
					} else {
						nd->A0[i][v] += T;
					}
				}
			}
		}
	}

	void downward_pass_non_rec(leaf* nd, std::vector<double>& L) const {
		int sz = nd->y - nd->x;
		std::vector<double> P0(sz, 0.0), P1(sz, 0.0);
		for (std::size_t v = 0; v < (1ull << nd->k1); v++) {
			for (int i = 0; i < sz; i++) {
				P0[i] += nd->A0[i][v] * nd->B[v];
				P1[i] += nd->A1[i][v] * nd->B[v];
			}
		}

		__log("P0 " << P0 << "\nP1 " << P1 << "\n");

		for (int i = nd->x, j = 0; j < sz; i++, j++) {
			if (P0[j] < EPS && P1[j] < EPS) {
				L[i] = NAN;
			} else if (P0[j] < EPS) {
				L[i] = -INF;
			} else if (P1[j] < EPS) {
				L[i] = INF;
			} else {
				L[i] = log(P0[j] / P1[j]);
			}
		}
	}

public:
	friend std::pair<int, int> simulate(RecursiveDecoder const& coder, double snr, int iter_cnt, int max_error) {
		double sigma = sqrt(0.5 * pow(10.0, -snr / 10.0) * coder.length() / coder.dim());
		std::normal_distribution norm(0.0, sigma);
		int errr = 0, sim_cnt = 0;
		while (errr < max_error && sim_cnt < iter_cnt) {
			binvector x(coder.dim());
			for (int t = 0; t < coder.dim(); t++) {
				x.set(t, rand() % 2);
			}
			binvector enc = coder.encode(x);

			std::vector<double> y(coder.length());
			double coef = 2.0f / (sigma * sigma);
			for (int t = 0; t < coder.length(); t++) {
				y[t] = ((enc[t] ? -1 : 1) + norm(gen));
			}
			binvector dec = coder.decode(y);
			errr += (dec != enc);
			sim_cnt++;
		}
		return {errr, sim_cnt};
	}

public:
	int length() const {
		return n;
	}

	int dim() const {
		return k;
	}

private:
	int n, k, d = -1;
	std::vector<std::vector<int>> _split;
	std::vector<binvector> G;
	node* root;

	double F(binvector const& c, std::vector<double> const& R, int x, int y) const {
		double ans = 0;
		for (int i = x, j = 0; i < y; i++, j++) {
			ans += (R[i] - (c[j] ? -1 : 1)) * (R[i] - (c[j] ? -1 : 1));
		}
		return exp(- ans / 2);
	}
};

int main() {
	srand(time(NULL));

	std::ifstream fin("input.txt");
#if defined(LOCAL)
	std::ostream& fout = std::cout;
#else
	std::ofstream fout("output.txt");
#endif
	// fout << std::scientific;

	int n, k;
	fin >> n >> k;
	std::vector<binvector> G(k, binvector(n));
	fin >> G;

	RecursiveDecoder coder(G);
	fout << CNT << std::endl;

	std::string command;
	bool skip = false;
	while (fin >> command) {
		if (command[0] == '-' && command[1] == '-') {
			std::getline(fin, command);
			continue;
		} else if (command == "{-") {
			skip = true;
			std::getline(fin, command);
			continue;
		} else if (command == "-}") {
			skip = false;
			std::getline(fin, command);
			continue;
		} else if (skip) {
			std::getline(fin, command);
			continue;
		} else if (command == "Encode") {
			binvector x(coder.dim());
			fin >> x;
			fout << coder.encode(x);
		} else if (command == "Decode") {
			std::vector<double> y(coder.length());
			fin >> y;
			fout << coder.decode(y);
		} else if (command == "DecodeSO") {
			std::vector<double> y(coder.length());
			bool c;
			for (int i = 0; i < y.size(); i++) {
				fin >> c;
				y[i] = (c ? -INF : INF);
			}
			fout << coder.decode_soft(y);
		} else if (command == "DecodeSISO") {
			std::vector<double> y(coder.length());
			fin >> y;
			fout << coder.decode_soft(y);
		} else if (command == "Simulate") {
			double snr;
			int iter_cnt, max_error;
			fin >> snr >> iter_cnt >> max_error;
			fail(iter_cnt > 0, "simulate: iter_cnt must be positive");
			fail(max_error >= 0, "simulate: max_err_cnt must be non-negative");

			#ifdef TEST
			auto start = std::chrono::system_clock::now();
			clear_cnt();
			#endif

			auto [errr, sim_cnt] = simulate(coder, snr, iter_cnt, max_error);

			#ifdef TEST
			auto end = std::chrono::system_clock::now();
			auto time_in_ms = std::chrono::duration_cast<ms>(end - start).count();
			#endif

			fout << (errr + 0.0) / sim_cnt;

			#if defined(TEST) && defined(TIMELOG)
			fout << " " << time_in_ms / 1000.0;
			#endif
			#if defined(TEST) && defined(LOG)
			fout << " " << CNT;
			#endif
		} else {
			fail(false, "main: incorrect command");
		}
		fout << std::endl;
	}

	std::cout << "FINISH\n";
}

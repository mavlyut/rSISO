#include <algorithm>
#include <array>
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

#include "defines.h"
#include "binvector.h"
#include "gray_code.cpp"
#ifdef CNTLOG
#include "count_ops.cpp"
#endif

#define fail(b, msg)						\
	if (!(b)) { 							\
		std::ofstream fout("output.txt");	\
		fout << "FAILURE: " << msg << "\n";	\
		fout.close();						\
		throw 2;							\
	}

static const _Float64 INF = INFINITY;
static constexpr _Float64 COEF = M_2_SQRTPI * M_SQRT1_2 / 2;
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
	for (unsigned i = 0; i < a.size(); i++) {
		if (i != 0) {
			out << ' ';
		}
		out << a[i];
	}
	return out;
}

template <typename T>
std::ostream& operator<<(std::ostream& out, std::vector<std::vector<T>> const& a) {
	for (unsigned i = 0; i < a.size(); i++) {
		if (i != 0) {
			out << "\n\t";
		}
		out << a[i];
	}
	return out;
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

template <typename T>
std::istream& operator>>(std::istream& in, std::vector<T>& a) {
	for (T& i : a) {
		in >> i;
	}
	return in;
}

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

class RecursiveDecoder {
public:
	RecursiveDecoder(matrix const& _G) : G(_G), G_enc(G) {
		time_measure(minimal_span_form(G));

		__log("MSF:\n" << G << "\n")

		k = G.size();
		n = G.front().size();
		p0.resize(n);
		p1.resize(n);
		_split = std::vector(n + 1, std::vector<unsigned>(n + 1, UNINIT));
		_ctor = std::vector(n + 1, std::vector<unsigned>(n + 1, 0));
		time_measure(root = sectionalization());
	}

private: // Sectionalization
	class node {
	protected:
		node(unsigned x, unsigned y, unsigned k1, unsigned k2, unsigned k3, matrix const& Gp)
			: x(x), y(y), k1(k1), k2(k2), k3(k3), A(1ull << k1, I), B(1ull << k1, I), Gp(Gp) {
			__log("Node ctor: x=" << x << ", y=" << y << "; k=(" << k1 << ", " << k2 << ", " << k3 << ")\n");
		}

	public:
		void clear() {
			_clear();
			for (double& i : A) {
				i = I;
			}
			for (double& i : B) {
				i = I;
			}
		}

		virtual void upward_pass(std::vector<double> const& p0, std::vector<double> const& p1) = 0;

		virtual void downward_pass(std::vector<double>& L) const = 0;

	public:
		const unsigned x, y;
		const unsigned k1, k2, k3;
		std::vector<double> A, B;
		matrix Gp;

	protected:
		inline static const double I = 0;
		virtual void _clear() = 0;
	};

	class leaf : public node {
	public:
		leaf(unsigned x, unsigned y, unsigned k1, matrix const& Gp)
			: node(x, y, k1, Gp.size() - k1, Gp.size() - k1, Gp) {}

	protected:
		virtual void _clear() = 0;

		double F(binvector const& c, std::vector<double> const& p0, std::vector<double> const& p1) const {
			double ans = 1.0;
			for (unsigned i = x, j = 0; i < y; i++, j++) {
				ans *= (c[j] ? p1[i] : p0[i]);
			}
			return ans;
		}

		double LLR(double M0, double M1) const {
			if (M0 < EPS && M1 < EPS) {
				return NAN;
			} else if (M0 < EPS) {
				return -INF;
			} else if (M1 < EPS) {
				return INF;
			} else {
				return log(M0 / M1);
			}
		}
	};

	class leaf_no_simplify : public leaf {
	public:
		leaf_no_simplify(unsigned x, unsigned y, unsigned k1, matrix const& Gp)
			: leaf(x, y, k1, Gp)
			, A0(y - x, std::vector<double>(1ull << k1, I))
			, A1(y - x, std::vector<double>(1ull << k1, I)) {}

	protected:
		std::vector<std::vector<double>> A0, A1;

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

		void upward_pass(std::vector<double> const& p0, std::vector<double> const& p1) override {
			binvector v(k1), c(y - x);
			for (unsigned ind : GrayCode(k1 + k2)) {
				if (ind != UNINIT) {
					if (ind >= k2) {
						v.change(ind - k2);
					}
					c ^= Gp[ind];
				}
				double T = F(c, p0, p1);
				A[v] += T;
				for (unsigned i = 0; i < y - x; i++) {
					if (c[i]) {
						A1[i][v] += T;
					} else {
						A0[i][v] += T;
					}
				}
			}
		}

		void downward_pass(std::vector<double>& L) const override {
			unsigned sz = y - x;
			std::vector<double> P0(sz, 0.0), P1(sz, 0.0);
			for (std::size_t v = 0; v < (1ull << k1); v++) {
				for (unsigned i = 0; i < sz; i++) {
					P0[i] += A0[i][v] * B[v];
					P1[i] += A1[i][v] * B[v];
				}
			}

			__log("P0 " << P0 << "\nP1 " << P1 << "\n");

			for (unsigned i = x, j = 0; j < sz; i++, j++) {
				L[i] = LLR(P0[j], P1[j]);
			}
		}
	};

#ifdef ENABLE_OPT_1
	class leaf_simplify_1 : public leaf {
	public:
		leaf_simplify_1(unsigned x, unsigned y, matrix const& Gp) : leaf(x, y, 0, Gp) {
			fail(Gp.size() == 1 && Gp[0].isOnes(), "leaf_s1: inappropriate optimization");
		}

		void upward_pass(std::vector<double> const& p0, std::vector<double> const& p1) override {
			for (unsigned i = x; i < y; i++) {
				ext_l += LLR(p0[i], p1[i]);
			}
		}

		void downward_pass(std::vector<double>& L) const override {
			for (unsigned i = x; i < y; i++) {
				L[i] = ext_l;
			}
		}

	private:
		double ext_l;

		void _clear() override {
			ext_l = I;
		}
	};
#endif // ENABLE_OPT_1

#ifdef ENABLE_OPT_2
	class leaf_simplify_2 : public leaf {
	public:
		leaf_simplify_2(unsigned x, unsigned y, matrix const& Gp) : leaf(x, y, 1, Gp) {
			fail(Gp.size() == 1 && Gp[0].isOnes(), "leaf_s2: inappropriate optimization");
		}

		void upward_pass(std::vector<double> const& p0, std::vector<double> const& p1) override {
			for (unsigned i = x; i < y; i++) {
				A[0] *= p0[i];
				A[1] *= p1[i];
			}
		}

		void downward_pass(std::vector<double>& L) const override {
			double ext_l = LLR(B[0] * A[0], B[1] * A[1]);
			for (unsigned i = x; i < y; i++) {
				L[i] = ext_l;
			}
		}

	private:
		void _clear() override {
			A[0] = 1.0;
			A[1] = 1.0;
		}
	};
#endif // ENABLE_OPT_2

#ifdef ENABLE_OPT_3
	class leaf_simplify_3 : public leaf {
	public:
		leaf_simplify_3(unsigned x, unsigned y, unsigned k1, matrix const& Gp) : leaf(x, y, k1, Gp) {
			fail(y == x + 2, "leaf_s3: inappropriate optimization");
			fail(k1 + k2 == 2 && k1 != 0, "leaf_s3: inappropriate optimization");
		}

		void upward_pass(std::vector<double> const& p0, std::vector<double> const& p1) override {
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

		void downward_pass(std::vector<double>& L) const override {
			std::vector<double> P0(2, I), P1(2, I);
			if (k1 == 1) {
				P0[0] = Phi00 * B[0] + Phi01 * B[1];
				P1[0] = Phi11 * B[0] + Phi10 * B[1];
				P0[1] = Phi00 * B[0] + Phi10 * B[1];
				P1[1] = Phi11 * B[0] + Phi01 * B[1];
			} else {
				P0[0] = Phi00 * B[0] + Phi01 * B[2];
				P0[1] = Phi00 * B[0] + Phi10 * B[1];
				P1[0] = Phi10 * B[1] + Phi11 * B[3];
				P1[1] = Phi01 * B[2] + Phi11 * B[3];
			}
			for (unsigned i = x, j = 0; i < y; i++, j++) {
				L[i] = LLR(P0[j], P1[j]);
			}
		}

	private:
		double Phi00, Phi10, Phi01, Phi11;

		void _clear() override {}
	};
#endif // ENABLE_OPT_3

#ifdef ENABLE_OPT_5
	class leaf_simplify_5 : public leaf {
	public:
		leaf_simplify_5(unsigned x, unsigned y, matrix const& Gp) : leaf(x, y, y - x, Gp) {
			fail(Gp.size() == y - x && y - x >= 2, "leaf_s5: inappropriate optimization");
		}

		void upward_pass(std::vector<double> const& p0, std::vector<double> const& p1) override {
			std::vector<double> diffs(k1);
			double T = 1;
			for (unsigned i = x, j = 0; i < y; i++, j++) {
				T *= p0[i];
				diffs[j] = p0[i] / p1[i];
			}
			binvector v(k1);
			for (unsigned ind : GrayCode(k1)) {
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

		void downward_pass(std::vector<double>& L_out) const override {
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

		private:
			void _clear() override {}
	};
#endif // ENABLE_OPT_5

	class inner : public node {
	public:
		inner(unsigned x, unsigned y, unsigned z, node* left, node* right, unsigned k1, unsigned k2, unsigned k3, matrix const& G_hat, matrix const& G_tilda, matrix const& Gp)
			: node(x, y, k1, k2, k3, Gp), z(z), left(left), right(right), G_hat(G_hat), G_tilda(G_tilda) {
			fail(k3 == left->k3 + right->k3 + k2, "inner-ctor: incorrect dims");
			fail(x < z && z < y, "inner-ctor: z must be in (x,y)");
			fail(left != nullptr && right != nullptr, "innert-ctor: empty left or right");
			fail(G_hat.size() == k1 + k2 && G_hat.front().size() == left->k1, "inner-ctor: incorrecy size of \\hat{G}");
			fail(G_tilda.size() == k1 + k2 && G_tilda.front().size() == right->k1, "inner-ctor: incorrecy size of \\tilda{G}");
		}

		void upward_pass(std::vector<double> const& p0, std::vector<double> const& p1) override {
			left->upward_pass(p0, p1);
			right->upward_pass(p0, p1);
			if (k1 != 0) {
				binvector a(left->k1), b(right->k1), v(k1);
				for (unsigned ind : GrayCode(k1 + k2)) {
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

		void downward_pass(std::vector<double>& L) const override {
			binvector a(left->k1), b(right->k1);
			if (k1 == 0) {
				for (unsigned ind : GrayCode(k1 + k2)) {
					if (ind != UNINIT) {
						a ^= G_hat[ind];
						b ^= G_tilda[ind];
					}
					left->B[a] = right->A[b];
					right->B[b] = left->A[a];
				}
			} else {
				binvector v(k1);
				for (unsigned ind : GrayCode(k1 + k2)) {
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

	public:
		unsigned z;
		node *left, *right;
		matrix G_hat, G_tilda;

	protected:
		void _clear() {
			left->clear();
			right->clear();
		}
	};

#ifdef ENABLE_OPT_7

#endif // ENABLE_OPT_7

	node* sectionalization() {
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
					update_phi(phi_d_l, x, y, 14, 3);
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
				update_phi(phi_d_l, x, y, ((1ull << (_Gp_size[x][y] - _Gs_size[x][y]))) * (y - x) * 4, 0);
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
		std::cout << "Expected operations: " << get_phi(0, n) << std::endl;

		__log("Gp_sz:  " << _Gp_size << "\nGs_sz:  " << _Gs_size << "\n");
		__log("split: " << _split << "\n");

		std::vector<unsigned> sections = {0, n};
		std::function<void(unsigned, unsigned)> rec_log = [&](unsigned x, unsigned y) {
			unsigned z = _split[x][y];

			__log("[" << x << ";" << y << ") at " << z << "\n");

			if (z == UNINIT) {
				return;
			}
			sections.push_back(z);
			rec_log(x, z);
			rec_log(z, y);
		};
		rec_log(0, n);
		std::sort(sections.begin(), sections.end());
		#ifdef LOG_SECTIONS
		std::cout << sections << std::endl;
		#endif

		return rec_build_section_tree(0, n);
	}

	node* rec_build_section_tree(unsigned x, unsigned y) {
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
				_log_out << "Simplify 1: " << x << " " << y << std::endl;
				return new leaf_simplify_1(x, y, Gp);
			}
			#endif

			#ifdef ENABLE_OPT_2
			if (optimizationNo == 2) {
				fail(k1 == 1 && Gp.size() == 1 && Gp[0].isOnes(), "incorrect pre-conditions (2)");
				_log_out << "Simplify 2: " << x << " " << y << std::endl;
				return new leaf_simplify_2(x, y, Gp);
			}
			#endif

			#ifdef ENABLE_OPT_3
			if (optimizationNo == 3) {
				fail(y - x == 2 && Gp.size() == 2, "incorrect pre-conditions (3)");
				_log_out << "Simplify 3: " << x << " " << y << std::endl;
				return new leaf_simplify_3(x, y, k1, Gp);
			}
			#endif

			#ifdef ENABLE_OPT_5
			if (optimizationNo == 5) {
				fail(y - x == k1, "incorrect pre-conditions (5)");
				_log_out << "Simplify 5: " << x << " " << y << std::endl;
				return new leaf_simplify_5(x, y, Gp);
			}
			#endif

			_log_out << "No simplify: " << x << " " << y << std::endl;
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

public:
	binvector encode(binvector const& c) const {
		fail(c.size() == k, "encode: size is incorrect");
		binvector ans(n);
		for (unsigned i = 0; i < k; i++) {
			if (c[i]) {
				ans ^= G_enc[i];
			}
		}
		return ans;
	}

public:
	std::vector<double> decode_soft(std::vector<_Float64> const& L_in) {
		for (unsigned i = 0; i < n; i++) {
			_Float64 L_exp = exp(L_in[i]);
			_Float64 z = 1.0 + L_exp;
			p0[i] = L_exp / z;
			p1[i] = 1.0 / z;
		}

		root->clear();
		std::vector<double> L_out(length(), 0.0);
		time_measure(root->upward_pass(p0, p1));
		time_measure(root->downward_pass(L_out));

		return L_out;
	}

public:
	friend std::pair<int, int> simulate(RecursiveDecoder coder, _Float64 snr, int iter_cnt, int max_error) {
		_Float64 sigma = sqrt(0.5 * pow(10.0, -snr / 10.0) * coder.length() / coder.dim());
		std::normal_distribution<_Float64> norm(0.0, sigma);
		int errr = 0, sim_cnt = 0;
		while (errr < max_error && sim_cnt < iter_cnt) {
			binvector x(coder.dim());
			for (unsigned t = 0; t < coder.dim(); t++) {
				x.set(t, rand() % 2);
			}
			binvector enc = coder.encode(x);
			_Float64 coef = 2.0f / (sigma * sigma);
			std::vector<_Float64> L_in(coder.length());
			for (unsigned t = 0; t < coder.length(); t++) {
				L_in[t] = coef * ((enc[t] ? -1 : 1) + norm(gen));
			}
			std::vector<double> L_out = coder.decode_soft(L_in);
			binvector dec(coder.length());
			for (unsigned t = 0; t < coder.length(); t++) {
				dec.set(t, static_cast<_Float64>(L_out[t]) < 0);
			}
			errr += (dec != enc);
			sim_cnt++;
		}
		return {errr, sim_cnt};
	}

public:
	unsigned length() const {
		return n;
	}

	unsigned dim() const {
		return k;
	}

private:
	unsigned n, k;
	std::vector<double> p0, p1;
	std::vector<std::vector<unsigned>> _split;
	std::vector<std::vector<unsigned>> _ctor;
	std::vector<binvector> G, G_enc;
	node* root;
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
	#if defined(CNTLOG)
	std::cout << "Init: " << CNT_BIN << "\n";
	CNT_BIN = 0;
	#endif

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
		} else if (command == "Decode" || command == "DecodeSISO") {
			std::vector<_Float64> y(coder.length());
			fin >> y;
			fout << coder.decode_soft(y);
		} else if (command == "Simulate") {
			_Float64 snr;
			int iter_cnt, max_error;
			fin >> snr >> iter_cnt >> max_error;
			fail(iter_cnt > 0, "simulate: iter_cnt must be positive");
			fail(max_error >= 0, "simulate: max_err_cnt must be non-negative");

			#ifdef TEST
			auto start = std::chrono::system_clock::now();
			#ifdef CNTLOG
			clear_cnt();
			clear_cnt_bin();
			#endif
			#endif

			auto [errr, sim_cnt] = simulate(coder, snr, iter_cnt, max_error);

			#ifdef TEST
			auto end = std::chrono::system_clock::now();
			auto time_in_ms = std::chrono::duration_cast<ms>(end - start).count();
			#endif
			#ifdef TEST
			fout << (errr + 0.0) / sim_cnt;
			#endif
			#if defined(TEST) && defined(TIMELOG)
			fout << " " << time_in_ms / 1000.0;
			#endif
			#if defined(CNTLOG)
			fout << SUM_CNT << " " << MUL_CNT << " " << CNT_BIN;
			#endif
		} else {
			fail(false, "main: incorrect command");
		}
		fout << std::endl;
	}

	#ifdef TEST
	std::cout << "FINISH\n";
	#endif
}

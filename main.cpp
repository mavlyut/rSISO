#include <algorithm>
#include <array>
#include "binpoly.cpp"
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

#define fail(b, msg) 						\
	if (!(b)) { 							\
		std::ofstream fout("output.txt"); 	\
		fout << "FAILURE: " << msg << "\n";	\
		fout.close();						\
		exit(135);							\
	}

using matrix = std::vector<binvector>;

template <typename T>
std::ostream& operator<<(std::ostream& out, std::vector<T> const& a) {
	for (int i = 0; i < a.size(); i++) {
		if (i != 0) {
			out << ' ';
		}
		out << a[i];
	}
	return out << '\n';
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
	std::vector<int> gamma(k);
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
		gamma[i] = last_c;
        for (int j = i + 1; j < k; ++j) {
            if (a[j][last_c] == 1) {
                for (int t = 0; t < n; ++t) {
                    a[j].set(t, a[j][t] ^ a[i][t]);
                }
            }
        }
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
		minimal_span_form(G);
		std::cout << "MSF:\n" << G << "\n";
		k = G.size();
		n = G.front().size();
		ZERO = binvector(n);
		_split = std::vector(n + 1, std::vector(n + 1, -1));
		root = sectionalization();
	}

private: // Sectionalization
	class node {
	protected:
		node(int x, int y, int k1, int k2, int k3, matrix const& Gp)
			: x(x), y(y), k1(k1), k2(k2), k3(k3), A(1 << k1), B(1 << k1), Gp(Gp) {}

	public:
		virtual bool is_leaf() const = 0;

	public:
		const int x, y;
		const int k1, k2, k3;
		std::vector<double> A, B;
		matrix Gp;

	protected:
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
	};

	class leaf : public node {
	public:
		leaf(int x, int y, int k1, matrix const& Gp)
			: node(x, y, k1, Gp.size() - k1, Gp.size() - k1, Gp)
			, A0(y - x, std::vector<double>(1 << k1, 0))
			, A1(y - x, std::vector<double>(1 << k1, 0)) {}

		bool is_leaf() const override {
			return true;
		}

		binvector dot(binvector const& u, binvector const& v) const {
			return combineAndMul(u, v, Gp);
		}

	public:
		std::vector<std::vector<double>> A0, A1;
	};

	class inner : public node {
	public:
		inner(int x, int y, int z, node* left, node* right, int k1, int k2, int k3, matrix const& G_hat, matrix const& G_tilda, matrix const& Gp)
			: node(x, y, k1, k2, k3, Gp), z(z), left(left), right(right), G_hat(G_hat), G_tilda(G_tilda) {}

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
	};

	node* sectionalization() {
		// std::vector<std::vector<int>> phi_u_c(n, std::vector<int>(n, 0));
		// std::vector<std::vector<int>> phi_d_c(n, std::vector<int>(n, 0));
		// std::vector<std::vector<int>> phi_u_l(n, std::vector<int>(n, 0));
		// std::vector<std::vector<int>> phi_d_l(n, std::vector<int>(n, 0));

		// for (int x = 0; x <= n; x++) {
		// 	for (int y = x + 1; y <= n; y++) {
		// 		_split[x][y] = (x + y) / 2;
		// 	}
		// }

		_split[0][8] = 4;
		_split[0][4] = 2;
		_split[4][8] = 6;

		return rec_build_section_tree(0, n);
	}

	node* rec_build_section_tree(int x, int y) {
		if (x >= y) {
			return nullptr;
		}
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
			std::cout << "leaf: " << x << " " << y << "\n" << Gp << "\n";
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

		node* left = rec_build_section_tree(x, z);
		node* right = rec_build_section_tree(z, y);

		std::cout << "inner: " << x << " " << y << "\n" << Gp << "\n" << G_0 << "\n" << G_1 << "\n";
		
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
			std::cout << "G_hat_ext\n" << G_hat_ext << "\n";
			for (int i = 0; i < k1 + k2; i++) {
				for (int j = 0; j < left->k1; j++) {
					G_hat[i].set(left->k1 - j - 1, G_hat_ext[G_hat_ext.size() - j - 1][i + left->k1 + left->k3]);
				}
			}
		}
		{
			matrix G_tilda_ext(z - x, binvector(right->k1 + right->k3 + k1 + k2));
			for (int i = 0; i < z - x; i++) {
				int j = 0;
				for (int t = 0; t < right->Gp.size(); t++, j++) {
					G_tilda_ext[i].set(j, right->Gp[t][i]);
				}
				for (int t = 0; t < G_1.size(); t++, j++) {
					G_tilda_ext[i].set(j, G_1[t][i]);
				}
			}
			full_gauss(G_tilda_ext);
			std::cout << "G_tilda_ext\n" << G_tilda_ext << "\n";
			for (int i = 0; i < k1 + k2; i++) {
				for (int j = 0; j < right->k1; j++) {
					G_tilda[i].set(right->k1 - j - 1, G_tilda_ext[G_tilda_ext.size() - j - 1][i + right->k1 + right->k3]);
				}
			}
		}

		std::cout << G_hat << "\n" << G_tilda << "\n";

		return new inner(x, y, z, left, right, k1, k2, k3, G_hat, G_tilda, Gp);
	}

public:
	binvector encode(binvector const& c) const {
		fail(c.size() == k, "encode: size is incorrect");
		return c * G;
	}

public:
	std::vector<double> decode_soft(std::vector<double> const& y) const {
		std::vector<double> L(length(), 0);
		upward_pass(root, y);
		downward_pass(root, L);
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
	void upward_pass(node* rt, std::vector<double> const& R) const {
		if (rt->is_leaf()) {
			upward_pass_non_rec(reinterpret_cast<leaf*>(rt), R);
		}
		inner* nd = reinterpret_cast<inner*>(rt);
		if (nd->k1 != 0) {
			upward_pass(nd->left, R);
			upward_pass(nd->right, R);
			binvector w_zero(nd->k2);
			for (std::size_t v = 0; v < (1ull << nd->k1); v++) {
				binvector vv(nd->k1, v);
				auto [a, b] = nd->mapping(w_zero, vv);
				nd->A[v] = nd->left->A[a] * nd->right->A[b];
				for (std::size_t w = 1; w < (1ull << nd->k2); w++) {
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
		binvector w_zero(nd->k1);
		if (nd->k1 == 0) {
			for (std::size_t w = 0; w < (1ull << nd->k2); w++) {
				binvector vv(nd->k2, w);
				auto [a, b] = nd->mapping(w_zero, vv);
				nd->left->B[a] = nd->right->A[b];
				nd->right->B[b] = nd->left->A[a];
			}
		} else {
			for (std::size_t v = 0; v < (1ull << nd->k1); v++) {
				binvector vv(nd->k1, v);
				auto [a, b] = nd->mapping(w_zero, vv);
				nd->left->B[a] = nd->B[v] * nd->right->A[b];
				nd->right->B[b] = nd->B[v] * nd->left->A[a];
				for (std::size_t w = 1; w < (1ull << nd->k2); w++) {
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
			binvector c = nd->dot(binvector(nd->k2), vv);
			nd->A[v] = F(c, R, nd->x, nd->y);
			for (std::size_t w = 1; w < (1 << nd->k2); w++) {
				binvector ww(nd->k2, w);
				c = nd->dot(ww, vv);
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
		std::vector<double> P0(sz, 0), P1(sz, 0);
		for (std::size_t v = 0; v < (1ull << nd->k1); v++) {
			for (int i = 0; i < sz; i++) {
				P0[i] += nd->A0[i][v] * nd->B[v];
				P1[i] += nd->A1[i][v] * nd->B[v];
			}
		}
		for (int i = nd->x, j = 0; j < sz; i++, j++) {
			L[i] = log(P0[j] / P1[j]);
		}
	}

public:
	int length() const {
		return n;
	}

	int dim() const {
		return k;
	}

	int dist() {
		if (d == -1) {
			d = n + 1;
			for (int i = 1; i < k; i++) {
				binvector bv(k, i);
				d = std::min(d, encode(bv).wt());
			}
		}
		return d;
	}

private:
	int n, k, d = -1;
	std::vector<std::vector<int>> _split;
	std::vector<binvector> G;
	node* root;
	binvector ZERO;

	double F(binvector c, std::vector<double> const& R, int x, int y) const {
		double ans = 0;
		for (int i = x, j = 0; i < y; i++, j++) {
			ans += (R[i] - (c[j] ? -1 : 1)) * (R[i] - (c[j] ? -1 : 1));
		}
		return ans;
	}
};

int main() {
	srand(time(NULL));
	std::mt19937 gen(time(0));

	std::ifstream fin("input.txt");
#if defined(LOCAL)
	std::ostream& fout = std::cout;
#else
	std::ofstream fout("output.txt");
#endif
	fout << std::scientific;

	int n, k;
	fin >> n >> k;
	std::vector<binvector> G(k, binvector(n));
	fin >> G;

	RecursiveDecoder coder(G);

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
		} else if (command == "Simulate") {
			double snr;
			int iter_cnt, max_error;
			fin >> snr >> iter_cnt >> max_error;
			fail(iter_cnt > 0, "simulate: iter_cnt must be positive");
			fail(max_error >= 0, "simulate: max_err_cnt must be non-negative");

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
				for (int t = 0; t < coder.length(); t++) {
					y[t] = (enc[t] ? -1 : 1) + norm(gen);
				}
				binvector dec = coder.decode(y);
				errr += (dec != enc);
				sim_cnt++;
			}
			fout << (errr + 0.0) / sim_cnt;
		} else {
			fail(false, "main: incorrect command");
		}
		fout << std::endl;
	}

	std::cout << "FINISH\n";
}

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

std::ostream& operator<<(std::ostream& out, std::vector<binvector> const& a) {
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

class RecursiveDecoder {
public:
	RecursiveDecoder(matrix const G) : _split(G[0].size(), std::vector(G[0].size(), -1)), G(G) {
		minimal_span_form(this->G);
		this->k = G.size();
		this->n = G.front().size();
		this->ZERO = binvector(n);

		sectionalization();
	}

private: // Sectionalization
	struct node {
		int x, y, z;
		int k1, k2, k3;
		node* left;
		node* right;
		matrix Gp, G_tilda, G_hat;
		std::vector<double> A, B;
	};

	void sectionalization() {
		// std::vector<std::vector<int>> phi_u_c(n, std::vector<int>(n, 0));
		// std::vector<std::vector<int>> phi_d_c(n, std::vector<int>(n, 0));
		// std::vector<std::vector<int>> phi_u_l(n, std::vector<int>(n, 0));
		// std::vector<std::vector<int>> phi_d_l(n, std::vector<int>(n, 0));

		for (int x = 0; x < n; x++) {
			for (int y = x + 1; y < n; y++) {
				_split[x][y] = (x + y) / 2;
			}
		}
	}

	int split(int x, int y) const {
		return _split[x][y];
	}

public:
	binvector encode(binvector const& c) const {
		fail(c.size() == k, "encode: size is incorrect");
		binvector ans(n);
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < k; j++) {
				ans.set(i, ans[i] ^ (c[j] & G[j][i]));
			}
		}
		return ans;
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
	void upward_pass(node* nd, std::vector<double> const& R) const {
		if (nd->z == -1) {
			upward_pass_non_rec(nd, R);
		}
		nd->A.resize(1LL << nd->k1);
		nd->B.resize(1LL << nd->k1);
		binvector a, b;
		if (nd->k1 != 0) {
			upward_pass(nd->left, R);
			upward_pass(nd->right, R);
			binvector vv(nd->k1 + nd->k2);
			for (std::size_t v = 0; v < (1ll << nd->k1); v++) {
				for (int t = 0; t < nd->k1; t++) {
					vv.set(vv.size() - t, (v >> (nd->k1 - t)) & 1);
				}
				a = vv * nd->G_hat;
				b = vv * nd->G_tilda;
				nd->A[v] = nd->left->A[a] * nd->right->A[b];
				for (std::size_t w = 1; w < (1ll << nd->k2); w++) {
					for (int t = 0; t < nd->k2; t++) {
						vv.set(t, (w >> t) & 1);
					}
					a = vv * nd->G_hat;
					b = vv * nd->G_tilda;
					nd->A[v] += nd->left->A[a] * nd->right->A[b];
				}
			}
		}
	}

	void downward_pass(node* nd, std::vector<double>& L) const {
		if (nd->z == -1) {
			downward_pass_non_rec(nd, L);
			return;
		}
		binvector a, b;
		if (nd->k1 == 0) {
			binvector vv(nd->k2);
			for (std::size_t w = 0; w < (1ll << nd->k2); w++) {
				binvector vv(nd->k2, w);
				a = vv * nd->G_hat;
				b = vv * nd->G_tilda;
				nd->left->B[a] = nd->right->A[b];
				nd->right->B[b] = nd->left->A[a];
			}
		} else {
			binvector vv(nd->k1 + nd->k2);
			for (std::size_t v = 0; v < (1ll << nd->k1); v++) {
				for (int t = 0; t < nd->k1; t++) {
					vv.set(vv.size() - t, (v >> (nd->k1 - t)) & 1);
				}
				a = vv * nd->G_hat;
				b = vv * nd->G_tilda;
				nd->left->B[a] = nd->B[v] * nd->right->A[b];
				nd->right->B[b] = nd->B[v] * nd->left->A[a];
				for (std::size_t w = 1; w < (1ll << nd->k2); w++) {
					for (int t = 0; t < nd->k2; t++) {
						vv.set(t, (w >> t) & 1);
					}
					a = vv * nd->G_hat;
					b = vv * nd->G_tilda;
					nd->left->B[a] += nd->B[v] * nd->right->A[b];
					nd->right->B[b] += nd->B[v] * nd->left->A[a];
				}
			}
		}
		downward_pass(nd->left, L);
		downward_pass(nd->right, L);
	}

	void upward_pass_non_rec(node* nd, std::vector<double> const& R) const {

	}

	void downward_pass_non_rec(node* nd, std::vector<double>& L) const {

	}

public:
	matrix Gp(int x, int y) const {
		matrix ans(dim(), binvector(y - x));
		for (int i = 0; i < dim(); i++) {
			for (int j = x, jj = 0; j < y; j++, jj++) {
				ans[i].set(jj, G[i][j]);
			}
		}
		int i = 0;
		while (i < ans.size()) {
			if (ans[i].isZero()) {
				ans.erase(ans.begin() + i);
			} else {
				i++;
			}
		}
		minimal_span_form(ans);
		return ans;
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

private:
	binvector ZERO;
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
	std::cout << coder.Gp(4, 8);

	return 0;

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
			std::cout << sigma << '\n';
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
}
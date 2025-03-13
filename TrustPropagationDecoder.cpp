#include <fstream>
#include <iostream>
#include <queue>

#include "include/matrix.h"

std::vector<binvector> getParityMatrix(std::vector<binvector> a) {
    int k0 = a.size();
    int n = a[0].size();
    int last_c = 0;
    std::vector<int> gamma;
    binvector in_gamma(n);
    for (int i = 0; i < k0; ++i) {
        while (last_c < n) {
            bool fl = false;
            for (int j = i; j < k0; ++j) {
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
            break;
        }
        in_gamma.set(last_c, true);
        gamma.push_back(last_c);
        for (int j = i + 1; j < k0; ++j) {
            if (a[j][last_c] == 1) {
                for (int t = 0; t < n; ++t) {
                    a[j].set(t, a[j][t] ^ a[i][t]);
                }
            }
        }
    }
    int k = gamma.size();
    for (int c = k - 1; c >= 0; c--) {
        for (int i = 0; i < c; i++) {
            if (a[i][gamma[c]] == 1) {
                for (int j = 0; j < n; j++) {
                    a[i].set(j, a[i][j] ^ a[c][j]);
                }
            }
        }
    }
    std::vector<binvector> H(n - k, binvector(n));
    for (int j = 0, cnt = 0; j < n; j++) {
        if (!in_gamma[j]) {
            for (int i = 0; i < k; i++) {
                H[cnt].set(gamma[i], a[i][j]);
            }
            H[cnt].set(j, true);
            cnt++;
        }
    }
    return H;
}

class TrustPropagationDecoder {
public:
    TrustPropagationDecoder(std::vector<binvector> const& H) : n(H[0].size()), G(getParityMatrix(H)) {
        C.resize(n), V.resize(m);
        for (int j = 0; j < m; j++) {
            for (int i = 0; i < n; i++) {
                if (H[j][i]) {
                    C[i].push_back(j);
                    V[j].push_back(i);
                }
            }
        }
    }

public:
    binvector encode(binvector const& x) const {
        fail(x.size() == dim(), "encode: invalid size");
        binvector c(n);
        for (int j = 0; j < dim(); j++) {
            if (x[j]) {
                c ^= G[j];
            }
        }
        return c;
    }

public:
    std::vector<double> decode(std::vector<double> const& L) const {
        fail(L.size() == length(), "decode: invalid size");
        std::vector<double> Lc(n, 0);
        std::vector<std::map<int, double>> Lq(n);
        std::vector<std::map<int, double>> Lr(m);
        for (int i = 0; i < n; i++) {
            for (int j : C[i]) {
                Lq[i][j] = L[i];
            }
        }
        binvector c(n);
        for (int _t = 0; _t < MAX_ITER_COUNT; _t++) {
            for (int j = 0; j < m; j++) {
                double sum_b = 0.0, sum_a = 1.0;
                for (int i : V[j]) {
                    auto l = Lq[i][j];
                    sum_b += phi(std::abs(l));
                    sum_a *= (l > 0 ? 1.0 : -1.0);
                }
                for (int i : V[j]) {
                    auto l = Lq[i][j];
                    Lr[j][i] = (l > 0 ? 1.0 : -1.0) * sum_a * phi(sum_b - phi(std::abs(l)));
                }
            }

            for (int i = 0; i < n; i++) {
                double sum_L = 0;
                for (int j : C[i]) {
                    sum_L += Lr[j][i];
                }
                for (int j : C[i]) {
                    Lq[i][j] = L[i] + sum_L - Lr[j][i];
                }
            }

            for (int i = 0; i < n; i++) {
                Lc[i] = L[i];
                for (int j : C[i]) {
                    Lc[i] += Lr[j][i];
                }
                c.set(i, Lc[i] < 0);
            }

            bool is_codeword = true;
            for (int j = 0; j < m; j++) {
                bool bit = 0;
                for (int i : V[j]) {
                    bit ^= c[i];
                }
                if (bit != 0) {
                    is_codeword = false;
                    break;
                }
            }
            if (is_codeword) {
                break;
            }
        }
        return Lc;
    }

private:
    static inline double phi(double x) {
        return -log(tanh(x / 2.0));
    }

    std::ostream& print(std::ostream& out, std::vector<std::map<int, double>> const& a) const {
        for (int i = 0; i < a.size(); i++) {
            for (int j = 0; j < n + m - a.size(); j++) {
                auto p = a[i].find(j);
                auto pp = ((p == a[i].end()) ? 0.0 : p->second);
                out << pp << " ";
            }
            out << "\n";
        }
        return out << "\n";
    }

    std::ostream& print(std::ostream& out, std::vector<double>& a) const {
        for (double i : a) {
            out << i << " ";
        }
        return out << "\n";
    }

private:
    void printGraph(std::ofstream fgraph, binvector const& y) const {
        fgraph << "digraph G {\n";
        fgraph << "\tranksep=1\n\n";
        fgraph << "\t subgraph checks {\n";
        fgraph << "\t\trankdir=LR\n";
        fgraph << "\t\trank=same\n";
        fgraph << "\t\trankC [style=invisible]\n";
        for (int j = 0; j < m; j++) {
            bool bit = 0;
            for (int i : V[j]) {
                bit ^= y[i];
            }
            fgraph << "\t\tc" << j << " [shape=box, label=\"\", color=" << (bit ? "red" : "green") << "]" << "\n";
            fgraph << "\t\t";
            if (j > 0) {
                fgraph << "c" << j - 1;
            } else {
                fgraph << "rankC";
            }
            fgraph << " -> c" << j << " [style=invis]\n";
        }
        fgraph << "\t}\n\n";
        fgraph << "\t subgraph symbols {\n";
        fgraph << "\t\trankdir=LR\n";
        fgraph << "\t\trank=same\n";
        fgraph << "\t\trankV [style=invisible]\n";
        for (int i = 0; i < n; i++) {
            fgraph << "\t\tv" << i << " [shape=circle, label=\"" << (y[i] ? '1' : '0') << "\"]" << "\n";
            fgraph << "\t\t";
            if (i > 0) {
                fgraph << "v" << i - 1;
            } else {
                fgraph << "rankV";
            }
            fgraph << " -> v" << i << " [style=invis]\n";
        }
        fgraph << "\t}\n\n";
        for (int j = 0; j < m; j++) {
            for (int i : V[j]) {
                fgraph << "\tc" << j << " -> v" << i << " [dir=none]\n";
            }
        }
        fgraph << "}\n";
        fgraph.close();
    }

public:
    int length() const {
        return n;
    }

    int dim() const {
        return G.size();
    }

private:
    int n, m;
    const int MAX_ITER_COUNT = 50;
    std::vector<std::vector<int>> C, V;
    std::vector<binvector> G;
};

int main() {
    std::mt19937 gen(time(0));
    std::ifstream fin("input.txt");
    std::ofstream fout("output.txt");

    // int n, m, z;
    // fin >> n >> m >> z >> z;
    // std::vector<binvector> H(m, binvector(n));
    // for (int i = 0; i < n; i++) {
    //     int w;
    //     fin >> w;
    //     for (int t = 0; t < w; t++) {
    //         int j;
    //         fin >> j;
    //         H[j].set(i, true);
    //     }
    // }
    // for (int j = 0; j < m; j++) {
    //     int w;
    //     fin >> w;
    //     for (int t = 0; t < w; t++) {
    //         int i;
    //         fin >> i;
    //         fail(H[j][i], "read: " + std::to_string(j) + " " + std::to_string(i));
    //     }
    // }

    int n, m;
    fin >> n >> m;
    std::vector<binvector> H(m, binvector(n));
    fin >> H;

    TrustPropagationDecoder coder(H);
    fout << coder.length() << " " << coder.dim() << std::endl;

    std::string command;
    while (fin >> command) {
        if (command == "Encode") {
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
            fail(iter_cnt > 0, "simulate: iter_cnt must be positive integer");
            fail(max_error >= 0, "simulate: max_error must be non-negative integer");

            double sigma = sqrt(0.5 * pow(10.0, -snr / 10.0) * coder.length() / coder.dim());
            std::cout << sigma << std::endl;
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
                std::vector<double> L = coder.decode(y);
                binvector dec(n);
                for (int i = 0; i < n; i++) {
                    dec.set(i, L[i] < 0);
                }
                errr += (dec != enc);
                sim_cnt++;
            }
            fout << errr << " " << sim_cnt;
        } else {
            fail(false, "main: unknown command");
        }
        fout << std::endl;
    }
}

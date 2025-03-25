#include "../include/trust_propagation_decoder.h"

matrix getParityMatrix(matrix a) {
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

trust_propagation_decoder::trust_propagation_decoder(matrix const& H) : linear_soft_decoder(getParityMatrix(H)), m(H.size()) {
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

std::vector<double> trust_propagation_decoder::decode_soft(std::vector<double> const& L) {
    fail(L.size() == length(), "ldpc, decode: invalid dim");
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

double trust_propagation_decoder::phi(double x) {
    if (x >= 7) {
        return 2.0 / (exp(x) + 1);
    }
    if (x <= 5e-2) {
        return log(2.0 / x);
    }
    return -log(tanh(x / 2.0));
}

std::ostream& trust_propagation_decoder::print(std::ostream& out, std::vector<std::map<int, double>> const& a) const {
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

void trust_propagation_decoder::print_graph(std::ofstream fgraph, binvector const& y) const {
    bool empty = (y.size() == 0);
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
        fgraph << "\t\tc" << j << " [shape=box, label=\"\", color="
                << (empty ? "black" : (bit ? "red" : "green")) << "]" << "\n";
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
        fgraph << "\t\tv" << i << " [shape=circle, label=\""
                << (empty ? "" : (y[i] ? "1" : "0")) << "\"]" << "\n";
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

#include <cmath>

#include "../include/plotkin_construction_decoder.h"

plotkin_construction_decoder::plotkin_construction_decoder(recursive_decoder* dec1, recursive_decoder* dec2)
        : n(dec1->length() * 2), m(dec1->length()), k(dec1->dim() + dec2->dim())
        , Lq0(m, 0), Lq1(m), Lq2(m), L_out(n), Lr(m)
        , dec1(dec1), dec2(dec2) {
    fail(dec1->length() == dec2->length(), "Pt, ctor: incompatible sizes");
    __log("Pt, created" << std::endl);
}

#define getbit(v, i) (((v) >> (i)) & 1)
#define subvector(v, x, y) (((v) >> (x)) & ((1ull << ((y) - (x))) - 1))
#define setbit(v, i, b) if (b) {v |= (1ull << (i));} else {v &= ~(1ull << (i));}
#define printbv(msg, n, v) 					\
	__log(msg);								\
	for (unsigned j = 0; j < n; ++j) { 		\
		if (j != 0) __log(" ");				\
		__log((((v >> j) & 1) ? 1 : 0));	\
	}                                       \
    __log(std::endl);

std::size_t plotkin_construction_decoder::encode(std::size_t const& c) const {
    printbv("Pt, encode: ", k, c);
    auto u1 = dec1->encode(subvector(c, 0, dec1->dim()));
    auto u2 = dec2->encode(subvector(c, dec1->dim(), dec1->dim() + dec2->dim()));
    printbv("Pt, encode, half1-encoded: ", m, u1);
    printbv("Pt, encode, half2-encoded: ", m, u2);
    printbv("Pt, encode, half1^half2  : ", m, (u1 ^ u2));
    std::size_t ans = (((u1 ^ u2) << m) | u1);
    printbv("Pt, encoded: ", n, ans);
    return ans;
}

double sign(double x) {
    if (x < 0) {
        return -1;
    }
    if (x > 0) {
        return 1;
    }
    return 0;
}

const double MAX_TRUNC = 32;
double truncate(double const &a) {
    if (a > MAX_TRUNC) {
        return MAX_TRUNC;
    }
    if (a < -MAX_TRUNC) {
        return -MAX_TRUNC;
    }
    return a;
}

const unsigned plotkin_construction_decoder::MAX_ITER_COUNT = 20;
std::vector<double> plotkin_construction_decoder::decode_soft(std::vector<double> const& L0) {
    fail(L0.size() == length(), "Pt, decode: incorrect size");
    __log("Pt, decode: " << L0 << std::endl);
    bool is_codeword;
    for (unsigned i = 0; i < m; ++i) {
        Lq0[i] = truncate(L0[i]);
        Lq1[i] = truncate(L0[i + m]);
        Lq2[i] = 0.0;
        Lr[i][0] = Lr[i][1] = Lr[i][2] = 0.0;
    }
    // 0: Lq0
    // 1: Lq1
    // 2: Lq2
    // 3: Lr_0
    // 4: Lr_1
    // 5: Lr_2
    // 6: checks
    printGraph(0, 0);
    for (unsigned _t = 1; _t <= MAX_ITER_COUNT; ++_t) {
        Lq0 = dec1->decode_soft(Lq0);
        printGraph(_t, 1);
        Lq2 = dec2->decode_soft(Lq2);
        printGraph(_t, 4);
    
        for (unsigned j = 0; j < m; ++j) {
            Lr[j][0] = truncate((sign(Lq1[j]) * sign(Lq2[j])) * phi(phi(std::abs(Lq1[j])) + phi(std::abs(Lq2[j]))));
        }
        printGraph(_t, 8);

        for (unsigned j = 0; j < m; ++j) {
            Lr[j][1] = truncate((sign(Lq2[j]) * sign(Lq0[j])) * phi(phi(std::abs(Lq2[j])) + phi(std::abs(Lq0[j]))));
        }
        printGraph(_t, 16);

        for (unsigned j = 0; j < m; ++j) {
            Lr[j][2] = truncate((sign(Lq0[j]) * sign(Lq1[j])) * phi(phi(std::abs(Lq0[j])) + phi(std::abs(Lq1[j]))));
        }
        printGraph(_t, 32);

        for (unsigned i = 0; i < m; ++i) {
            Lq0[i] = truncate(Lq0[i] + Lr[i][0]);
        }
        printGraph(_t, 128+1);

        for (unsigned i = 0; i < m; ++i) {
            Lq1[i] = truncate(Lq1[i] + Lr[i][1]);
        }
        printGraph(_t, 128+2);

        for (unsigned i = 0; i < m; ++i) {
            Lq2[i] = truncate(Lq2[i] + Lr[i][2]);
        }
        printGraph(_t, 128+4);

        __log("After " << _t + 1 << " iteration: " << Lq1 << " " << Lq2 << std::endl);

        // is_codeword = true;
        // for (unsigned j = 0; j < m; ++j) {
        //     if ((Lq0[j] < 0) ^ (Lq1[j] < 0) ^ (Lq2[j] < 0)) {
        //         is_codeword = false;
        //         break;
        //     }
        // }
        printGraph(_t, 128+64);
        // if (is_codeword) {
        //     break;
        // }
    }
    for (unsigned i = 0; i < m; ++i) {
        L_out[i] = Lq0[i] + Lr[i][0];
        L_out[i + m] = Lq1[i] + Lr[i][1];
    }
    return L_out;
}

std::string fill_color(bool b) {
    if (!b) {
        return "";
    }
    return ",style=\"filled\",fillcolor=\"lightgreen\"";
}

void plotkin_construction_decoder::printGraph(unsigned iter_num, unsigned state) {
    #ifdef VISUALIZE
    std::ofstream out("img/DecodeState_" + std::to_string(iter_num) + "_" + std::to_string(state) + ".dot");
    bool __lq0 = getbit(state, 0);
    bool __lq1 = getbit(state, 1);
    bool __lq2 = getbit(state, 2);
    bool __lr0 = getbit(state, 3);
    bool __lr1 = getbit(state, 4);
    bool __lr2 = getbit(state, 5);
    bool __chs = getbit(state, 6);
    bool __lr[]{__lr0, __lr1, __lr2};
    bool __lr_any = (__lr0 || __lr1 || __lr2);

    out << std::fixed;
    out << std::setprecision(1);
    out << "digraph G {\n";

    out << "\tsubgraph value_symbols {\n";
    out << "\t\trankdir=LR\n\t\trank=same\n";
    for (unsigned i = 0; i < m; i++) {
        out << "\t\tv" << i << " [label=\"" << Lq0[i] << "\"]\n";
    }
    for (unsigned i = 0; i < m; i++) {
        out << "\t\tv" << i + m << " [label=\"" << Lq1[i] << "\"" << fill_color(__lq1) << "]\n";
    }
    for (unsigned i = 0; i < n; i++) {
        if (i != 0) {
            out << " -> ";
        } else {
            out << "\t\t";
        }
        out << "v" << i;
    }
    out << " [style=invis]\n";
    out << "\t}\n\n";

    out << "\tsubgraph check_symbols {\n";
    out << "\t\trankdir=LR\n\t\trank=same\n";
    for (unsigned i = 0; i < m; i++) {
        out << "\t\tch" << i << " [shape=square,label=\"+\"" << fill_color(__lr_any);
        if (__chs) {
            out << ",color=" << (((Lq0[i] < 0) ^ (Lq1[i] < 0) ^ (Lq2[i] < 0)) ? "red" : "green");
        }
        out << "]\n";
    }
    for (unsigned i = 0; i < m; i++) {
        if (i != 0) {
            out << " -> ";
        } else {
            out << "\t\t";
        }
        out << "ch" << i;
    }
    out << " [style=invis]\n";
    out << "\t}\n\n";

    for (unsigned i = 0; i < n; i++) {
        out << "\tv" << i << " -> ch" << (i % m) << " [dir=none, headlabel=<" << (__lr[i / m] ? "<b>" : "") << Lr[i % m][i / m] << (__lr[i / m] ? "</b>" : "") << ">]\n";
    }
    out << "\n";

    out << "\tsubgraph subcodes {\n";
    out << "\t\trankdir=LR\n\t\trank=same\n";
    for (unsigned i = 0; i < m; i++) {
        out << "\t\th" << i << " [color=blue,label=\"" << Lq0[i] << "\"" << fill_color(__lq0) << "]\n";
    }
    for (unsigned i = 0; i < m; i++) {
        out << "\t\th" << i + m << " [color=purple,label=\"" << Lq2[i] << "\"" << fill_color(__lq2) << "]\n";
    }
    for (unsigned i = 0; i < n; i++) {
        if (i != 0) {
            out << " -> ";
        } else {
            out << "\t\t";
        }
        out << "h" << i;
    }
    out << " [style=invis]\n";
    out << "\t}\n\n";

    for (unsigned i = 0; i < m; i++) {
        out << "\th" << i << " -> v" << i << " [dir=none]\n";
    }
    for (unsigned i = 0; i < m; i++) {
        out << "\th" << i + m << " -> ch" << i << " [dir=none, headlabel=<" << (__lr[2] ? "<b>" : "") << Lr[i][2] << (__lr[2] ? "</b>" : "") << ">]\n";
    }
    out << "\n";

    out << "\t{\n";
    out << "\t\tnode [style=invis]\n";
    out << "\t\tA;B;C;D\n";
    out << "\t\tedge [weight=10, style=invis]\n";
    out << "\t\tA -> B [minlen=3]\n";
    out << "\t\tB -> C [minlen=3]\n";
    out << "\t\tC -> D\n";
    out << "\t}\n\n";

    out << "\t{ rank=same; v" << n - 1 << " -> A [weight=10,style=invis];}\n";
    out << "\t{ rank=same; ch" << m - 1 << " -> B [weight=10,style=invis];}\n";
    out << "\t{ rank=same; h" << n - 1 << " -> C [weight=10,style=invis];}\n";

    out << "}\n";

    __log("Pt, printGraph: finish" << std::endl);
    #endif
}

double plotkin_construction_decoder::phi(double x) {
    if (x >= 7) {
        return 2.0 / (exp(x) + 1);
    }
    if (x <= 5e-2) {
        return log(2.0 / x);
    }
    return -log(tanh(x / 2.0));  
}

std::pair<double, double> plotkin_construction_decoder::simulate(double snr, unsigned iter_cnt, unsigned max_error) {
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
        errr += (dec != enc);
        ++sim_cnt;
    }
    return {(errr + 0.0) / sim_cnt, (berr + 0.0) / (sim_cnt * n)};
}

unsigned plotkin_construction_decoder::length() const {
    return n;
}

unsigned plotkin_construction_decoder::dim() const {
    return k;
}

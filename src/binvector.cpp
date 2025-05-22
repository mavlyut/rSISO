#include "../include/binvector.h"

namespace short_domain {
    bool getbit(binvector const& v, unsigned i) {
        return ((v >> i) & 1);
    }

    binvector& setbit(binvector& v, unsigned i, bool b) {
        return (b ? (v |= (1ull << i)) : (v &= (~(1ull << i))));
    }

	binvector& changebit(binvector& v, unsigned i) {
        return (v ^= (1ull << i));
    }

	binvector subvector(binvector const& v, unsigned x, unsigned y) {
        return ((v >> x) & ((1ull << (y - x)) - 1));
    }

    std::ostream& printbv(std::ostream& out, std::string const& msg, unsigned n, binvector const& bv) {
        out << msg;
        for (unsigned i = 0; i < n; ++i) {
            if (i != 0) {
                out << " ";
            }
            out << getbit(bv, i);
        }
        return out;
    }

    std::ostream& printbv(std::ostream& out, unsigned n, binvector const& bv) {
        return printbv(out, "", n, bv);
    }

	void _log_bv(std::string const& msg, unsigned n, binvector const& bv) {
        __log(msg);
        for (unsigned i = 0; i < n; ++i) {
            if (i != 0) {
                __log(" ");
            }
            __log(getbit(bv, i));
        }
    }

	unsigned wt(binvector const& v) {
        binvector tmp = v;
        unsigned ans = 0;
        while (tmp > 0) {
            if (tmp & 1) {
                ++ans;
            }
            tmp >>= 1;
        }
        return ans;
    }
}

namespace long_domain {
    bool getbit(binvector const& v, unsigned i) {
        return v[i];
    }

	binvector& setbit(binvector& v, unsigned i, bool b) {
        v.set(i, b);
        return v;
    }

	binvector& changebit(binvector& v, unsigned i) {
        v.change(i);
        return v;
    }

	binvector subvector(binvector const& v, unsigned x, unsigned y) {
        return v.subvector(x, y);
    }

    std::ostream& printbv(std::ostream& out, std::string const& msg, unsigned n, binvector const& bv) {
        out << msg;
        for (unsigned i = 0; i < n; ++i) {
            if (i != 0) {
                out << " ";
            }
            out << bv[i];
        }
        return out;
    }

    std::ostream& printbv(std::ostream& out, unsigned n, binvector const& bv) {
        return printbv(out, "", n, bv);
    }

	void _log_bv(std::string const& msg, unsigned n, binvector const& bv) {
        __log(msg);
        for (unsigned i = 0; i < n; ++i) {
            if (i != 0) {
                __log(" ");
            }
            __log(bv[i]);
        }
    }

	unsigned wt(binvector const& v) {
        return v.wt();
    }
}

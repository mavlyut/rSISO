#ifndef BINVECTOR_H
#define BINVECTOR_H

#include <algorithm>
#include <array>
#include <chrono>
#include <climits>
#include <cmath>
#include <ctime>
#include <fstream>
#include <functional>
#include <iostream>
#include <iomanip>
#include <map>
#include <ratio>
#include <vector>

#include "defines.h"

namespace short_domain {
	using binvector = std::uint64_t;

	binvector get_empty(unsigned sz);
	bool getbit(binvector const& v, unsigned i);
	binvector& setbit(binvector& v, unsigned i, bool b);
	binvector& changebit(binvector& v, unsigned i);
	binvector subvector(binvector const& v, unsigned x, unsigned y);
	std::ostream& printbv(std::ostream& out, unsigned n, binvector const& bv);
	std::ostream& printbv(std::ostream& out, std::string const& msg, unsigned n, binvector const& bv);
	void _log_bv(std::string const& msg, unsigned n, binvector const& bv);
	unsigned wt(binvector const& v);
}

namespace long_domain {
	template <typename _int_t, std::size_t chunk_size>
	class __binvector {
		static inline constexpr _int_t __min(_int_t const& a, _int_t const& b) {
			return (a < b ? a : b);
		}

		static constexpr const _int_t MAX_DIGIT
			= (chunk_size < 64)
			? __min(std::numeric_limits<_int_t>::max(), (1ull << chunk_size))
			: std::numeric_limits<_int_t>::max();

	public:
		__binvector() = default;

		__binvector(unsigned m) : __binvector(m, 0) {}

		explicit __binvector(unsigned m, _int_t x) : sz(m), coefs((m + chunk_size - 1) / chunk_size, 0) {
			for (unsigned i = 0; i < m; ++i, x >>= 1) {
				set(i, x & 1);
			}
		}

		static __binvector getEmpty(int m) {
			return __binvector(m);
		}

		__binvector(__binvector const& other) : __binvector(other.sz, other.coefs) {}

		__binvector& operator=(__binvector const& other) {
			sz = other.sz;
			coefs = other.coefs;
			return *this;
		}

		bool operator[](unsigned i) const {
			if (i < sz) {
				return ((coefs[i / chunk_size] >> (i % chunk_size)) & 1);
			}
			return false;
		}

		void set(unsigned i, bool val) {
			++CNT_BIN;
			if (val) {
				coefs[i / chunk_size] |= (1ull << (i % chunk_size));
			} else {
				coefs[i / chunk_size] &= ~(1ull << (i % chunk_size));
			}
		}

		void change(unsigned i) {
			++CNT_BIN;
			coefs[i / chunk_size] ^= (1ull << (i % chunk_size));
		}

		__binvector subvector(unsigned x, unsigned y) const {
			if (x >= y) {
				return __binvector(0);
			}
			__binvector ans(y - x);
			for (unsigned i = 0; i < ans.size(); ++i) {
				ans.set(i, operator[](i + x));
			}
			return ans;
		}

		unsigned size() const {
			return sz;
		}

		_int_t to_integer() const {
			return coefs.size() > 0 ? coefs[0] : 0;
		}

		operator _int_t() const {
			return to_integer();
		}

		int wt() const {
			int ans = 0;
			for (unsigned i = 0; i < size(); ++i) {
				ans += operator[](i);
			}
			return ans;
		}

	public:
		__binvector& operator^=(__binvector const& r) {
			fail(r.size() <= size(), "^=: invalid size");
			CNT_BIN += r.coefs.size();
			for (unsigned i = 0; i < r.coefs.size(); ++i) {
				coefs[i] ^= r.coefs[i];
			}
			return *this;
		}

		friend __binvector operator^(__binvector a, __binvector const& b) {
			return a ^= b;
		}

	public:
		friend std::istream& operator>>(std::istream& in, __binvector& a) {
			int b;
			for (unsigned i = 0; i < a.size(); ++i) {
				in >> b;
				a.set(i, b);
			}
			return in;
		}

	public:
		friend std::ostream& operator<<(std::ostream& out, __binvector const& a) {
			for (unsigned i = 0; i < a.size(); ++i) {
				if (i != 0) {
					out << ' ';
				}
				out << a[i];
			}
			return out;
		}

	public:
		friend bool operator<(__binvector const& a, __binvector const& b) {
			if (a.size() != b.size()) {
				return a.size() < b.size();
			}
			for (unsigned i = a.size(); i-- > 0; ) {
				if (a[i] < b[i]) {
					return true;
				}
				if (a[i] > b[i]) {
					return false;
				}
			}
			return false;
		}

		friend bool operator>(__binvector const& a, __binvector const& b) {
			return b < a;
		}

		friend bool operator<=(__binvector const& a, __binvector const& b) {
			return !(a > b);
		}

		friend bool operator>=(__binvector const& a, __binvector const& b) {
			return !(a < b);
		}

		friend bool operator==(__binvector const& a, __binvector const& b) {
			return a.sz == b.sz && a.coefs == b.coefs;
		}

		friend bool operator!=(__binvector const& a, __binvector const& b) {
			return !(a == b);
		}

		bool isZero() const {
			for (auto const& i : coefs) {
				if (i != 0) {
					return false;
				}
			}
			return true;
		}

		bool isOnes() const {
			for (auto const& i : coefs) {
				if (i != MAX_DIGIT) {
					return false;
				}
			}
			return true;
		}

	private:
		unsigned sz;
		std::vector<_int_t> coefs;

		explicit __binvector(unsigned m, std::vector<_int_t> const& coefs) : sz(m), coefs(coefs) {}

		explicit __binvector(unsigned m, std::vector<_int_t> const& coefs, bool) : sz(m), coefs(coefs) {
			while (!operator[](sz - 1)) {
				sz--;
			}
		}
	};

	using binvector = __binvector<std::size_t, 32>;

	binvector get_empty(unsigned sz);
	bool getbit(binvector const& v, unsigned i);
	binvector& setbit(binvector& v, unsigned i, bool b);
	binvector& changebit(binvector& v, unsigned i);
	binvector subvector(binvector const& v, unsigned x, unsigned y);
	std::ostream& printbv(std::ostream& out, unsigned n, binvector const& bv);
	std::ostream& printbv(std::ostream& out, std::string const& msg, unsigned n, binvector const& bv);
	void _log_bv(std::string const& msg, unsigned n, binvector const& bv);
	unsigned wt(binvector const& v);
}

#endif // BINVECTOR_H

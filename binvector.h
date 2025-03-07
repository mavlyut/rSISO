#pragma once
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

static std::size_t CNT_BIN = 0;
void clear_cnt_bin() {
	CNT_BIN = 0;
}

template <typename _int_t, std::size_t chunk_size>
class binvector_ {
	static inline constexpr _int_t __min(_int_t const& a, _int_t const& b) {
		return (a < b ? a : b);
	}

	static constexpr const _int_t MAX_DIGIT
		= (chunk_size < 64)
		? __min(std::numeric_limits<_int_t>::max(), (1ull << chunk_size))
		: std::numeric_limits<_int_t>::max();

public:
	binvector_() = default;

	binvector_(unsigned m) : binvector_(m, 0) {}

	explicit binvector_(unsigned m, _int_t x) : sz(m), coefs((m + chunk_size - 1) / chunk_size, 0) {
		for (unsigned i = 0; i < m; i++, x >>= 1) {
			set(i, x & 1);
		}
	}

	static binvector_ getEmpty(int m) {
		return binvector_(m);
	}

	binvector_(binvector_ const& other) : binvector_(other.sz, other.coefs) {}

	binvector_& operator=(binvector_ const& other) {
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
		CNT_BIN++;
		if (val) {
			coefs[i / chunk_size] |= (1ull << (i % chunk_size));
		} else {
			coefs[i / chunk_size] &= ~(1ull << (i % chunk_size));
		}
	}

	void change(unsigned i) {
		CNT_BIN++;
		coefs[i / chunk_size] ^= (1ull << (i % chunk_size));
	}

	binvector_ subvector(unsigned x, unsigned y) const {
		if (x >= y) {
			return binvector_(0);
		}
		binvector_ ans(y - x);
		for (unsigned i = 0; i < ans.size(); i++) {
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
		for (unsigned i = 0; i < size(); i++) {
			ans += operator[](i);
		}
		return ans;
	}

public:
	binvector_& operator^=(binvector_ const& r) {
		// fail(r.size() <= size(), "^=: invalid size");
		CNT_BIN += r.coefs.size();
		for (unsigned i = 0; i < r.coefs.size(); i++) {
			coefs[i] ^= r.coefs[i];
		}
		return *this;
	}

	friend binvector_ operator^(binvector_ a, binvector_ const& b) {
		return a ^= b;
	}

public:
	friend std::istream& operator>>(std::istream& in, binvector_& a) {
		int b;
		for (unsigned i = 0; i < a.size(); i++) {
			in >> b;
			a.set(i, b);
		}
		return in;
	}

public:
	friend std::ostream& operator<<(std::ostream& out, binvector_ const& a) {
		for (unsigned i = 0; i < a.size(); i++) {
			if (i != 0) {
				out << ' ';
			}
			out << a[i];
		}
		return out;
	}

public:
	friend bool operator<(binvector_ const& a, binvector_ const& b) {
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

	friend bool operator>(binvector_ const& a, binvector_ const& b) {
		return b < a;
	}

	friend bool operator<=(binvector_ const& a, binvector_ const& b) {
		return !(a > b);
	}

	friend bool operator>=(binvector_ const& a, binvector_ const& b) {
		return !(a < b);
	}

	friend bool operator==(binvector_ const& a, binvector_ const& b) {
		return a.sz == b.sz && a.coefs == b.coefs;
	}

	friend bool operator!=(binvector_ const& a, binvector_ const& b) {
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

	explicit binvector_(unsigned m, std::vector<_int_t> const& coefs) : sz(m), coefs(coefs) {}

	explicit binvector_(unsigned m, std::vector<_int_t> const& coefs, bool) : sz(m), coefs(coefs) {
		while (!operator[](sz - 1)) {
			sz--;
		}
	}

	friend struct std::hash<binvector_>;
};

template <typename T, std::size_t N>
struct std::hash<binvector_<T, N>> {
	std::size_t operator()(binvector_<T, N> const& p) const noexcept {
		auto int_hasher = std::hash<T>{};
		std::size_t ans = N * int_hasher(p.sz);
		for (T const& x : p.coefs) {
			ans = ((ans << 1) | int_hasher(x));
		}
		return ans;
	}
};

using binvector = binvector_<std::size_t, 32>;

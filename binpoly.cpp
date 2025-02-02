#pragma once
#include <algorithm>
#include <array>
#include <bitset>
#include <chrono>
#include <climits>
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
#include <unordered_set>
#include <vector>

template <typename _int_t>
class binvector_ {
	static const int chunk_size = sizeof(_int_t) * CHAR_BIT;

public:
	binvector_() = default;

	explicit binvector_(int m, int x) : sz(m), coefs((m + 1 + chunk_size) / chunk_size, 0) {
		for (int i = 0; i <= m; i++) {
			set(i, (x >> i) & 1);
		}
	}

	binvector_(int m) : binvector_(m, 0) {}

	explicit binvector_(int m, _int_t x) : sz(m), coefs((m + 1 + chunk_size) / chunk_size, 0) {
		coefs[0] = x;
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

	bool operator[](int i) const {
		// fail(i >= 0, "index must be non-negative");
		if (i < sz) {
			return ((coefs[i / chunk_size] >> (i % chunk_size)) & 1);
		}
		return 0;
	}

	void set(int i, bool val) {
		// fail(i >= 0 && i < sz, "index out of bounds");
		if (val) {
			coefs[i / chunk_size] |= (1ull << (i % chunk_size));
		} else {
			coefs[i / chunk_size] &= ~(1ull << (i % chunk_size));
		}
	}

	binvector_ subvector(int x, int y) const {
		binvector_ ans(std::max(0, y - x));
		for (int i = 0; i < ans.size(); i++) {
			ans.set(i, operator[](i + x));
		}
		return ans;
	}

	std::size_t size() const {
		return sz;
	}

	_int_t to_integer() const {
		return coefs[0];
	}

	operator _int_t() const {
		return to_integer();
	}

	int wt() const {
		int ans = 0;
		for (int i = 0; i < size(); i++) {
			ans += operator[](i);
		}
		return ans;
	}

public:
	friend binvector_ operator+(binvector_ const& a, binvector_ const& b) {
		binvector_ ans = getEmpty(std::max(a.size(), b.size()));
		for (int i = 0; i < ans.size(); i++) {
			ans.set(i, a[i] ^ b[i]);
		}
		return ans;
	}

	friend binvector_ mulAsPoly(binvector_ const& a, binvector_ const& b) {
		binvector_ ans = getEmpty(a.size() + b.size());
		for (int i = 0; i < ans.size(); i++) {
			bool val = false;
			for (int j = 0; j <= i; j++) {
				val ^= (a[j] & b[i - j]);
			}
			ans.set(i, val);
		}
		return ans;
	}

	friend bool operator*(binvector_ const& a, binvector_ const& b) {
		bool ans = false;
		for (int i = 0; i < std::min(a.size(), b.size()); i++) {
			ans ^= (a[i] & b[i]);
		}
		return ans;
	}

	friend binvector_ operator%(binvector_ a, binvector_ const& b) {
		while (true) {
			int diff = a.size() - b.size();
			if (diff < 0 || a.isZero()) {
				return a;
			}
			for (int i = diff; i < a.size(); i++) {
				a.set(i, a[i] ^ b[i - diff]);
			}
		}
	}

public:
	friend std::istream& operator>>(std::istream& in, binvector_& a) {
		int b;
		for (int i = 0; i < a.size(); i++) {
			in >> b;
			a.set(i, b);
		}
		return in;
	}

public:
	friend std::ostream& operator<<(std::ostream& out, binvector_ const& a) {
		for (int i = 0; i < a.size(); i++) {
			if (i != 0) {
				out << ' ';
			}
			out << a[i];
		}
		return out;
	}

	std::ostream& printAsPoly(std::ostream& out) const {
		for (int i = size() - 1, start = 1; i >= 0; i--) {
			if (operator[](i)) {
				if (start) {
					start = 0;
				} else {
					out << " + ";
				}
				if (i == 0) {
					out << "1";
				} else {
					out << "x";
					if (i > 1) {
						out << "^" << i;
					}
				}
			}
		}
		return out;
	}

public:
	friend bool operator<(binvector_ const& a, binvector_ const& b) {
		if (a.size() != b.size()) {
				return a.size() < b.size();
		}
		for (int i = a.size() - 1; i >= 0; i--) {
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
		return (a + b).isZero();
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

private:
	int sz;
	std::vector<_int_t> coefs;

	explicit binvector_(int m, std::vector<_int_t> const& coefs) : sz(m), coefs(coefs) {}

	friend struct std::hash<binvector_>;
};

template <typename T>
struct std::hash<binvector_<T>> {
	std::size_t operator()(binvector_<T> const& p) const noexcept {
		auto int_hasher = std::hash<T>{};
		std::size_t ans = int_hasher(p.sz);
		for (T x : p.coefs) {
			ans = ((ans << 1) | int_hasher(x));
		}
		return ans;
	}
};

using binvector = binvector_<std::size_t>;

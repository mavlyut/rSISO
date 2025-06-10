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

#include "../include/matrix.h"

class binpoly {
	using _int_t = std::size_t;
	static const int chunk_size = sizeof(_int_t) * CHAR_BIT;

public:
	binpoly() = default;

	binpoly(int m, int x) : binpoly(m, x, true) {}

private:
	binpoly(int m, int x, bool need_canonize) : deg(m), coefs((m * 2 + 1 + chunk_size) / chunk_size, 0) {
		for (int i = 0; i <= m; i++) {
			set(i, (x >> i) & 1);
		}
		if (need_canonize) {
			canonize();
		}
	}

public:
	static binpoly getEmpty(int m) {
		return binpoly(m, 0, false);
	}

	static binpoly getEmptyBySize(int m) {
		return getEmpty(m - 1);
	}

public:
	binpoly(binpoly const& other) : binpoly(other.deg, other.coefs) {}

	binpoly& operator=(binpoly const& other) {
		deg = other.deg;
		coefs = other.coefs;
		return *this;
	}

	bool operator[](int i) const {
		fail(i >= 0, "index must be non-negative");
		if (i <= deg) {
			return ((coefs[i / chunk_size] >> (i % chunk_size)) & 1);
		}
		return 0;
	}

	void set(int i, bool val) {
		fail(i >= 0 && i <= deg, "index out of bounds");
		if (val) {
			coefs[i / chunk_size] |= (1ull << (i % chunk_size));
		} else {
			coefs[i / chunk_size] &= ~(1ull << (i % chunk_size));
		}
	}

	int getDeg() const {
		return deg;
	}

	int size() const {
		return deg + 1;
	}

	_int_t to_integer() const {
		return coefs[0];
	}

public:
	friend binpoly operator+(binpoly const& a, binpoly const& b) {
		binpoly ans = getEmpty(std::max(a.getDeg(), b.getDeg()));
		for (int i = 0; i <= ans.getDeg(); i++) {
			ans.set(i, a[i] ^ b[i]);
		}
		return ans.canonize();
	}

	friend binpoly operator*(binpoly const& a, binpoly const& b) {
		binpoly ans = getEmpty(a.getDeg() + b.getDeg());
		for (int i = 0; i <= ans.getDeg(); i++) {
			bool val = false;
			for (int j = 0; j <= i; j++) {
				val ^= (a[j] & b[i - j]);
			}
			ans.set(i, val);
		}
		return ans.canonize();
	}

	friend binpoly operator%(binpoly a, binpoly const& b) {
		while (true) {
			int diff = a.getDeg() - b.getDeg();
			if (diff < 0 || a.isZero()) {
				return a;
			}
			for (int i = diff; i <= a.getDeg(); i++) {
				a.set(i, a[i] ^ b[i - diff]);
			}
			a.canonize();
		}
	}

public:
	friend std::istream& operator>>(std::istream& in, binpoly& a) {
		int b;
		for (int i = 0; i <= a.getDeg(); i++) {
			in >> b;
			a.set(i, b);
		}
		return in;
	}

public:
	friend std::ostream& operator<<(std::ostream& out, binpoly const& a) {
		for (int i = 0; i <= a.getDeg(); i++) {
			if (i != 0) {
				out << ' ';
			}
			out << a[i];
		}
		return out;
	}

	std::ostream& printPoly(std::ostream& out) const {
		if (getDeg() == 0) {
			return out << operator[](0);
		}
		for (int i = getDeg(), start = 1; i >= 0; i--) {
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
	friend bool operator<(binpoly const& a, binpoly const& b) {
		if (a.getDeg() != b.getDeg()) {
			return a.getDeg() < b.getDeg();
		}
		for (int i = a.getDeg(); i >= 0; i--) {
			if (a[i] < b[i]) {
				return true;
			}
			if (a[i] > b[i]) {
				return false;
			}
		}
		return false;
	}

	friend bool operator>(binpoly const& a, binpoly const& b) {
		return b < a;
	}

	friend bool operator<=(binpoly const& a, binpoly const& b) {
		return !(a > b);
	}

	friend bool operator>=(binpoly const& a, binpoly const& b) {
		return !(a < b);
	}

	friend bool operator==(binpoly const& a, binpoly const& b) {
		return a <= b && a >= b;
	}

	friend bool operator!=(binpoly const& a, binpoly const& b) {
		return !(a == b);
	}

	bool isZero() const {
		return deg == 0 && coefs[0] == 0;
	}

private:
	int deg;
	std::vector<_int_t> coefs;

	explicit binpoly(int m, std::vector<_int_t> const& coefs) : deg(m), coefs(coefs) {}

	binpoly& canonize() {
		while (deg > 0 && !operator[](deg)) {
			deg--;
		}
		return *this;
	}

	friend struct std::hash<binpoly>;
};

template <>
struct std::hash<binpoly> {
	std::size_t operator()(binpoly const& p) const noexcept {
		auto int_hasher = std::hash<binpoly::_int_t>{};
		std::size_t ans = int_hasher(p.deg);
		for (binpoly::_int_t x : p.coefs) {
			ans = ((ans << 1) | int_hasher(x));
		}
		return ans;
	}
};

using field_elem = std::uint16_t;

class GF {
private:
	field_elem __mulPoly(field_elem a, field_elem b, field_elem const& c) {
		field_elem ans = 0ull;
		while (b > 0) {
			if (b % 2) {
				ans ^= a;
			}
			b >>= 1;
			a <<= 1;
			if (a > n1) {
				a ^= c;
			}
		}
		return ans;
	}

public:
	GF(int m, int prim)
		: m(m), n1((1ull << m) - 1)
		, prim(prim)
		, ZERO(0ull), ONE(1ull)
		, mapDegToPoly(1ull << m, ZERO)
		, mapPolyToDeg(1ull << m, ZERO)
		, CayleyTable(1ull << m, std::vector<field_elem>(1ull << m))
		{

		for (field_elem a = 0; a < (1ull << m); a++) {
			for (field_elem b = a; b < (1ull << m); b++) {
				CayleyTable[a][b] = __mulPoly(a, b, prim);
				if (a != b) {
					CayleyTable[b][a] = CayleyTable[a][b];
				}
			}
		}

		mapPolyToDeg[ZERO] = 0;
		mapDegToPoly[0] = ZERO;
		for (field_elem a = 2; a <= n1; a++) {
			field_elem tmp = ONE;
			bool restart = false;
			for (field_elem i = 0; i < n1; i++) {
				mapPolyToDeg[tmp] = i;
				mapDegToPoly[i] = tmp;
				tmp = mulPoly(tmp, a);
				if (tmp == ONE && i < n1 - 1) {
					restart = true;
					break;
				}
			}
			if (!restart) {
				break;
			}
			fail(a != n1, "primitive element not found");
		}
	}

public:
	field_elem invPoly(field_elem const& x) const {
		return mapDegToPoly[(n1 - mapPolyToDeg[x]) % n1];
	}

	field_elem plusPoly(field_elem const& a, field_elem const& b) const {
		return a ^ b;
	}

	field_elem mulPoly(field_elem const& a, field_elem const& b) const {
		return CayleyTable[a][b];
	}

public:
	binpoly genMinimal(int deg) const {
		std::vector<field_elem> M(1, ONE);
		field_elem beta = mapDegToPoly[deg];
		field_elem tmp = beta;
		while (true) {
			M.push_back(ZERO);
			for (int i = M.size() - 1; i >= 0; i--) {
				M[i] = mulPoly(M[i], tmp);
				if (i > 0) {
					M[i] = plusPoly(M[i], M[i - 1]);
				}
			}
			tmp = mulPoly(tmp, tmp);
			if (tmp == beta) {
				break;
			}
		}
		binpoly ans = binpoly::getEmptyBySize(M.size());
		for (int i = 0; i < M.size(); i++) {
			ans.set(i, M[i] != ZERO);
		}
		return ans;
	}

public:
	int getMulGroupSize() const {
		return n1;
	}

	field_elem getPolyByDeg(int deg) const {
		return mapDegToPoly[deg];
	}

private:
	int m, n1;
	field_elem prim;
	field_elem const ZERO, ONE;
	std::vector<field_elem> mapDegToPoly;
	std::vector<field_elem> mapPolyToDeg;
	std::vector<std::vector<field_elem>> CayleyTable;
};

int main() {
	int N = 64, delta = 9;
	int prim = 64 + 2 + 1, m = -1;

	int n = N;
	while (n > 0) {
		m++;
		n /= 2;
	}

	GF field(m, prim);
	binpoly g(0, 1);
	std::unordered_set<binpoly> used_minimal;
	for (int j = 1; j <= delta - 1; j++) {
		binpoly M = field.genMinimal(j);
		if (used_minimal.find(M) == used_minimal.end()) {
			g = g * M;
			used_minimal.insert(M);
		}
	}

	std::cout << g << "\n";

	int K = N - 1 - g.getDeg();

	std::ofstream fout("eBCH_" + std::to_string(N) + "_" + std::to_string(K) + ".txt");

	fout << N << " " << K << "\n";

	for (int i = 0; i < K; i++) {
		long_domain::binvector row(N);
		bool parity = false;
		for (int j = 0; i + j < N - 1 && j <= g.getDeg(); j++) {
			row.set(i + j, g[j]);
			parity ^= g[j];
		}
		row.set(N - 1, parity);
		fout << row << "\n";
	}
}

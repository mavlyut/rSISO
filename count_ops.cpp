#include <istream>
#include <ostream>

#define CMP_AS_OP

static std::size_t SUM_CNT = 0ull;
static std::size_t MUL_CNT = 0ull;

void clear_cnt() {
	SUM_CNT = 0ull;
	MUL_CNT = 0ull;
}

struct __double_cnt {
	__double_cnt() : x(0) {}

	__double_cnt(_Float64 x) : x(x) {}

	__double_cnt(__double_cnt const& b) : x(b.x) {}

	__double_cnt& operator=(__double_cnt const& b) {
		x = b.x;
		return *this;
	}

	friend __double_cnt& operator+=(__double_cnt& a, __double_cnt const& b) {
		SUM_CNT++;
		a.x += b.x;
		return a;
	}

	template<typename T>
	friend __double_cnt& operator+=(__double_cnt& a, T const& b) {
		SUM_CNT++;
		a.x += ((_Float64)b);
		return a;
	}

	friend __double_cnt& operator-=(__double_cnt& a, __double_cnt const& b) {
		SUM_CNT++;
		a.x -= b.x;
		return a;
	}

	template <typename T>
	friend __double_cnt& operator-=(__double_cnt& a, T const& b) {
//		SUM_CNT++;
		a.x -= ((_Float64)b);
		return a;
	}

	friend __double_cnt& operator*=(__double_cnt& a, __double_cnt const& b) {
		MUL_CNT++;
		a.x *= b.x;
		return a;
	}

	template <typename T>
	friend __double_cnt& operator*=(__double_cnt& a, T const& b) {
//		MUL_CNT++;
		a.x *= ((_Float64)b);
		return a;
	}

	friend __double_cnt& operator/=(__double_cnt& a, __double_cnt const& b) {
		MUL_CNT++;
		a.x /= b.x;
		return a;
	}

	template <typename T>
	friend __double_cnt& operator/=(__double_cnt& a, T const& b) {
//		MUL_CNT++;
		a.x /= ((_Float64)b);
		return a;
	}

	friend __double_cnt operator+(__double_cnt const& a, __double_cnt const& b) {
		SUM_CNT++;
		return __double_cnt(a.x + b.x);
	}

	template <typename T>
	friend __double_cnt operator+(T const& a, __double_cnt const& b) {
//		SUM_CNT++;
		return ((_Float64)a) + b.x;
	}

	template <typename T>
	friend __double_cnt operator+(__double_cnt const& a, T const& b) {
//		SUM_CNT++;
		return a.x + ((_Float64)b);
	}

	friend __double_cnt operator-(__double_cnt const& a, __double_cnt const& b) {
		SUM_CNT++;
		return __double_cnt(a.x - b.x);
	}

	template <typename T>
	friend __double_cnt operator-(T const& a, __double_cnt const& b) {
//		SUM_CNT++;
		return ((_Float64)a) - b.x;
	}

	template <typename T>
	friend __double_cnt operator-(__double_cnt const& a, T const& b) {
//		SUM_CNT++;
		return a.x - ((_Float64)b);
	}

	friend __double_cnt operator*(__double_cnt const& a, __double_cnt const& b) {
		MUL_CNT++;
		return __double_cnt(a.x * b.x);
	}

	template <typename T>
	friend __double_cnt operator*(T const& a, __double_cnt const& b) {
//		MUL_CNT++;
		return ((_Float64)a) * b.x;
	}

	template <typename T>
	friend __double_cnt operator*(__double_cnt const& a, T const& b) {
//		MUL_CNT++;
		return a.x * ((_Float64)b);
	}

	friend __double_cnt operator/(__double_cnt const& a, __double_cnt const& b) {
		MUL_CNT++;
		return __double_cnt(a.x / b.x);
	}

	template <typename T>
	friend __double_cnt operator/(T const& a, __double_cnt const& b) {
//		MUL_CNT++;
		return ((_Float64)a) / b.x;
	}

	template <typename T>
	friend __double_cnt operator/(__double_cnt const& a, T const& b) {
//		MUL_CNT++;
		return a.x / ((_Float64)b);
	}

	friend bool operator==(__double_cnt const& a, __double_cnt const& b) {
		return a.x == b.x;
	}

	template <typename T>
	friend bool operator==(__double_cnt const& a, T const& b) {
		return a.x == ((_Float64)b);
	}

	template <typename T>
	friend bool operator==(T const& a, __double_cnt const& b) {
		return ((_Float64)a) == b.x;
	}

	friend bool operator!=(__double_cnt const& a, __double_cnt const& b) {
		return a.x != b.x;
	}

	template <typename T>
	friend bool operator!=(__double_cnt const& a, T const& b) {
		return a.x != ((_Float64)b);
	}

	template <typename T>
	friend bool operator!=(T const& a, __double_cnt const& b) {
		return ((_Float64)a) != b.x;
	}

	friend bool operator<(__double_cnt const& a, __double_cnt const& b) {
		return a.x < b.x;
	}

	template <typename T>
	friend bool operator<(__double_cnt const& a, T const& b) {
		return a.x < ((_Float64)b);
	}

	template <typename T>
	friend bool operator<(T const& a, __double_cnt const& b) {
		return ((_Float64)a) < b.x;
	}

	friend bool operator>(__double_cnt const& a, __double_cnt const& b) {
		return a.x > b.x;
	}

	template <typename T>
	friend bool operator>(__double_cnt const& a, T const& b) {
		return a.x > ((_Float64)b);
	}

	template <typename T>
	friend bool operator>(T const& a, __double_cnt const& b) {
		return ((_Float64)a) > b.x;
	}

	friend bool operator<=(__double_cnt const& a, __double_cnt const& b) {
		return a.x <= b.x;
	}

	template <typename T>
	friend bool operator<=(__double_cnt const& a, T const& b) {
		return a.x <= ((_Float64)b);
	}

	template <typename T>
	friend bool operator<=(T const& a, __double_cnt const& b) {
		return ((_Float64)a) <= b.x;
	}

	friend bool operator>=(__double_cnt const& a, __double_cnt const& b) {
		return a.x >= b.x;
	}

	template <typename T>
	friend bool operator>=(__double_cnt const& a, T const& b) {
		return a.x >= ((_Float64)b);
	}

	template <typename T>
	friend bool operator>=(T const& a, __double_cnt const& b) {
		return ((_Float64)a) >= b.x;
	}

	operator _Float64() const {
		return x;
	}

	__double_cnt operator+() const {
		return *this;
	}

	__double_cnt operator-() const {
		return __double_cnt(-x);
	}

	friend std::istream& operator>>(std::istream& in, __double_cnt& a) {
		return in >> a.x;
	}

	friend std::ostream& operator<<(std::ostream& out, __double_cnt const& a) {
		return out << a.x;
	}

private:
	_Float64 x;
};

#define double __double_cnt

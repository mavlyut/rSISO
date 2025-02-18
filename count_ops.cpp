#include <istream>
#include <ostream>

#define CMP_AS_OP

static std::size_t CNT = 0ull;

void clear_cnt() {
	CNT = 0ull;
}

struct __int_cnt {
	__int_cnt() : x(0) {}

	__int_cnt(signed x) : x(x) {}

	__int_cnt(__int_cnt const& b) : x(b.x) {}

	__int_cnt& operator=(__int_cnt const& b) {
		x = b.x;
		return *this;
	}

	friend __int_cnt& operator+=(__int_cnt& a, __int_cnt const& b) {
		CNT++;
		a.x += b.x;
		return a;
	}

	template<typename T>
	friend __int_cnt& operator+=(__int_cnt& a, T const& b) {
		CNT++;
		a.x += ((signed)b);
		return a;
	}

	friend __int_cnt& operator-=(__int_cnt& a, __int_cnt const& b) {
		CNT++;
		a.x -= b.x;
		return a;
	}

	template <typename T>
	friend __int_cnt& operator-=(__int_cnt& a, T const& b) {
		CNT++;
		a.x -= ((signed)b);
		return a;
	}

	friend __int_cnt& operator*=(__int_cnt& a, __int_cnt const& b) {
		CNT++;
		a.x *= b.x;
		return a;
	}

	template <typename T>
	friend __int_cnt& operator*=(__int_cnt& a, T const& b) {
		CNT++;
		a.x *= ((signed)b);
		return a;
	}

	friend __int_cnt& operator/=(__int_cnt& a, __int_cnt const& b) {
		CNT++;
		a.x /= b.x;
		return a;
	}

	template <typename T>
	friend __int_cnt& operator/=(__int_cnt& a, T const& b) {
		CNT++;
		a.x /= ((signed)b);
		return a;
	}

	friend __int_cnt operator+(__int_cnt const& a, __int_cnt const& b) {
		CNT++;
		return a.x + b.x;
	}

	template <typename T>
	friend __int_cnt operator+(T const& a, __int_cnt const& b) {
		CNT++;
		return ((signed)a) + b.x;
	}

	template <typename T>
	friend __int_cnt operator+(__int_cnt const& a, T const& b) {
		CNT++;
		return a.x + ((signed)b);
	}

	friend __int_cnt operator-(__int_cnt const& a, __int_cnt const& b) {
		CNT++;
		return a.x - b.x;
	}

	template <typename T>
	friend __int_cnt operator-(T const& a, __int_cnt const& b) {
		CNT++;
		return ((signed)a) - b.x;
	}

	template <typename T>
	friend __int_cnt operator-(__int_cnt const& a, T const& b) {
		CNT++;
		return a.x - ((signed)b);
	}

	friend __int_cnt operator*(__int_cnt const& a, __int_cnt const& b) {
		CNT++;
		return a.x * b.x;
	}

	template <typename T>
	friend __int_cnt operator*(T const& a, __int_cnt const& b) {
		CNT++;
		return ((signed)a) * b.x;
	}

	template <typename T>
	friend __int_cnt operator*(__int_cnt const& a, T const& b) {
		CNT++;
		return a.x * ((signed)b);
	}

	friend __int_cnt operator/(__int_cnt const& a, __int_cnt const& b) {
		CNT++;
		return a.x / b.x;
	}

	template <typename T>
	friend __int_cnt operator/(T const& a, __int_cnt const& b) {
		CNT++;
		return ((signed)a) / b.x;
	}

	template <typename T>
	friend __int_cnt operator/(__int_cnt const& a, T const& b) {
		CNT++;
		return a.x / ((signed)b);
	}

	friend bool operator==(__int_cnt const& a, __int_cnt const& b) {
		#ifdef CMP_AS_OP
		CNT++;
		#endif
		return a.x == b.x;
	}

	template <typename T>
	friend bool operator==(__int_cnt const& a, T const& b) {
		#ifdef CMP_AS_OP
		CNT++;
		#endif
		return a.x == ((signed)b);
	}

	template <typename T>
	friend bool operator==(T const& a, __int_cnt const& b) {
		#ifdef CMP_AS_OP
		CNT++;
		#endif
		return ((signed)a) == b.x;
	}

	friend bool operator!=(__int_cnt const& a, __int_cnt const& b) {
		#ifdef CMP_AS_OP
		CNT++;
		#endif
		return a.x != b.x;
	}

	template <typename T>
	friend bool operator!=(__int_cnt const& a, T const& b) {
		#ifdef CMP_AS_OP
		CNT++;
		#endif
		return a.x != ((signed)b);
	}

	template <typename T>
	friend bool operator!=(T const& a, __int_cnt const& b) {
		#ifdef CMP_AS_OP
		CNT++;
		#endif
		return ((signed)a) != b.x;
	}

	friend bool operator<(__int_cnt const& a, __int_cnt const& b) {
		#ifdef CMP_AS_OP
		CNT++;
		#endif
		return a.x < b.x;
	}

	template <typename T>
	friend bool operator<(__int_cnt const& a, T const& b) {
		#ifdef CMP_AS_OP
		CNT++;
		#endif
		return a.x < ((signed)b);
	}

	template <typename T>
	friend bool operator<(T const& a, __int_cnt const& b) {
		#ifdef CMP_AS_OP
		CNT++;
		#endif
		return ((signed)a) < b.x;
	}

	friend bool operator>(__int_cnt const& a, __int_cnt const& b) {
		#ifdef CMP_AS_OP
		CNT++;
		#endif
		return a.x > b.x;
	}

	template <typename T>
	friend bool operator>(__int_cnt const& a, T const& b) {
		#ifdef CMP_AS_OP
		CNT++;
		#endif
		return a.x > ((signed)b);
	}

	template <typename T>
	friend bool operator>(T const& a, __int_cnt const& b) {
		#ifdef CMP_AS_OP
		CNT++;
		#endif
		return ((signed)a) > b.x;
	}

	friend bool operator<=(__int_cnt const& a, __int_cnt const& b) {
		#ifdef CMP_AS_OP
		CNT++;
		#endif
		return a.x <= b.x;
	}

	template <typename T>
	friend bool operator<=(__int_cnt const& a, T const& b) {
		#ifdef CMP_AS_OP
		CNT++;
		#endif
		return a.x <= ((signed)b);
	}

	template <typename T>
	friend bool operator<=(T const& a, __int_cnt const& b) {
		#ifdef CMP_AS_OP
		CNT++;
		#endif
		return ((signed)a) <= b.x;
	}

	friend bool operator>=(__int_cnt const& a, __int_cnt const& b) {
		#ifdef CMP_AS_OP
		CNT++;
		#endif
		return a.x >= b.x;
	}

	template <typename T>
	friend bool operator>=(__int_cnt const& a, T const& b) {
		#ifdef CMP_AS_OP
		CNT++;
		#endif
		return a.x >= ((signed)b);
	}

	template <typename T>
	friend bool operator>=(T const& a, __int_cnt const& b) {
		#ifdef CMP_AS_OP
		CNT++;
		#endif
		return ((signed)a) >= b.x;
	}

	__int_cnt& operator++() {
		CNT++;
		this->x++;
		return *this;
	}

	__int_cnt operator++(signed) {
		__int_cnt tmp(*this);
		CNT++;
		this->x++;
		return tmp;
	}

	__int_cnt& operator--() {
		CNT++;
		this->x--;
		return *this;
	}

	__int_cnt operator--(signed) {
		__int_cnt tmp(*this);
		CNT++;
		this->x--;
		return tmp;
	}

	operator signed() const {
		return x;
	}

	friend std::istream& operator>>(std::istream& in, __int_cnt& a) {
		return in >> a.x;
	}

	friend std::ostream& operator<<(std::ostream& out, __int_cnt const& a) {
		return out << a.x;
	}

private:
	signed x;
};

struct __double_cnt {
	__double_cnt() : x(0) {}

	__double_cnt(_Float64 x) : x(x) {}

	__double_cnt(__double_cnt const& b) : x(b.x) {}

	__double_cnt& operator=(__double_cnt const& b) {
		x = b.x;
		return *this;
	}

	friend __double_cnt& operator+=(__double_cnt& a, __double_cnt const& b) {
		CNT++;
		a.x += b.x;
		return a;
	}

	template<typename T>
	friend __double_cnt& operator+=(__double_cnt& a, T const& b) {
		CNT++;
		a.x += ((_Float64)b);
		return a;
	}

	friend __double_cnt& operator-=(__double_cnt& a, __double_cnt const& b) {
		CNT++;
		a.x -= b.x;
		return a;
	}

	template <typename T>
	friend __double_cnt& operator-=(__double_cnt& a, T const& b) {
		CNT++;
		a.x -= ((_Float64)b);
		return a;
	}

	friend __double_cnt& operator*=(__double_cnt& a, __double_cnt const& b) {
		CNT++;
		a.x *= b.x;
		return a;
	}

	template <typename T>
	friend __double_cnt& operator*=(__double_cnt& a, T const& b) {
		CNT++;
		a.x *= ((_Float64)b);
		return a;
	}

	friend __double_cnt& operator/=(__double_cnt& a, __double_cnt const& b) {
		CNT++;
		a.x /= b.x;
		return a;
	}

	template <typename T>
	friend __double_cnt& operator/=(__double_cnt& a, T const& b) {
		CNT++;
		a.x /= ((_Float64)b);
		return a;
	}

	friend __double_cnt operator+(__double_cnt const& a, __double_cnt const& b) {
		CNT++;
		return a.x + b.x;
	}

	template <typename T>
	friend __double_cnt operator+(T const& a, __double_cnt const& b) {
		CNT++;
		return ((_Float64)a) + b.x;
	}

	template <typename T>
	friend __double_cnt operator+(__double_cnt const& a, T const& b) {
		CNT++;
		return a.x + ((_Float64)b);
	}

	friend __double_cnt operator-(__double_cnt const& a, __double_cnt const& b) {
		CNT++;
		return a.x - b.x;
	}

	template <typename T>
	friend __double_cnt operator-(T const& a, __double_cnt const& b) {
		CNT++;
		return ((_Float64)a) - b.x;
	}

	template <typename T>
	friend __double_cnt operator-(__double_cnt const& a, T const& b) {
		CNT++;
		return a.x - ((_Float64)b);
	}

	friend __double_cnt operator*(__double_cnt const& a, __double_cnt const& b) {
		CNT++;
		return a.x * b.x;
	}

	template <typename T>
	friend __double_cnt operator*(T const& a, __double_cnt const& b) {
		CNT++;
		return ((_Float64)a) * b.x;
	}

	template <typename T>
	friend __double_cnt operator*(__double_cnt const& a, T const& b) {
		CNT++;
		return a.x * ((_Float64)b);
	}

	friend __double_cnt operator/(__double_cnt const& a, __double_cnt const& b) {
		CNT++;
		return a.x / b.x;
	}

	template <typename T>
	friend __double_cnt operator/(T const& a, __double_cnt const& b) {
		CNT++;
		return ((_Float64)a) / b.x;
	}

	template <typename T>
	friend __double_cnt operator/(__double_cnt const& a, T const& b) {
		CNT++;
		return a.x / ((_Float64)b);
	}

	friend bool operator==(__double_cnt const& a, __double_cnt const& b) {
		#ifdef CMP_AS_OP
		CNT++;
		#endif
		return a.x == b.x;
	}

	template <typename T>
	friend bool operator==(__double_cnt const& a, T const& b) {
		#ifdef CMP_AS_OP
		CNT++;
		#endif
		return a.x == ((_Float64)b);
	}

	template <typename T>
	friend bool operator==(T const& a, __double_cnt const& b) {
		#ifdef CMP_AS_OP
		CNT++;
		#endif
		return ((_Float64)a) == b.x;
	}

	friend bool operator!=(__double_cnt const& a, __double_cnt const& b) {
		#ifdef CMP_AS_OP
		CNT++;
		#endif
		return a.x != b.x;
	}

	template <typename T>
	friend bool operator!=(__double_cnt const& a, T const& b) {
		#ifdef CMP_AS_OP
		CNT++;
		#endif
		return a.x != ((_Float64)b);
	}

	template <typename T>
	friend bool operator!=(T const& a, __double_cnt const& b) {
		#ifdef CMP_AS_OP
		CNT++;
		#endif
		return ((_Float64)a) != b.x;
	}

	friend bool operator<(__double_cnt const& a, __double_cnt const& b) {
		#ifdef CMP_AS_OP
		CNT++;
		#endif
		return a.x < b.x;
	}

	template <typename T>
	friend bool operator<(__double_cnt const& a, T const& b) {
		#ifdef CMP_AS_OP
		CNT++;
		#endif
		return a.x < ((_Float64)b);
	}

	template <typename T>
	friend bool operator<(T const& a, __double_cnt const& b) {
		#ifdef CMP_AS_OP
		CNT++;
		#endif
		return ((_Float64)a) < b.x;
	}

	friend bool operator>(__double_cnt const& a, __double_cnt const& b) {
		#ifdef CMP_AS_OP
		CNT++;
		#endif
		return a.x > b.x;
	}

	template <typename T>
	friend bool operator>(__double_cnt const& a, T const& b) {
		#ifdef CMP_AS_OP
		CNT++;
		#endif
		return a.x > ((_Float64)b);
	}

	template <typename T>
	friend bool operator>(T const& a, __double_cnt const& b) {
		#ifdef CMP_AS_OP
		CNT++;
		#endif
		return ((_Float64)a) > b.x;
	}

	friend bool operator<=(__double_cnt const& a, __double_cnt const& b) {
		#ifdef CMP_AS_OP
		CNT++;
		#endif
		return a.x <= b.x;
	}

	template <typename T>
	friend bool operator<=(__double_cnt const& a, T const& b) {
		#ifdef CMP_AS_OP
		CNT++;
		#endif
		return a.x <= ((_Float64)b);
	}

	template <typename T>
	friend bool operator<=(T const& a, __double_cnt const& b) {
		#ifdef CMP_AS_OP
		CNT++;
		#endif
		return ((_Float64)a) <= b.x;
	}

	friend bool operator>=(__double_cnt const& a, __double_cnt const& b) {
		#ifdef CMP_AS_OP
		CNT++;
		#endif
		return a.x >= b.x;
	}

	template <typename T>
	friend bool operator>=(__double_cnt const& a, T const& b) {
		#ifdef CMP_AS_OP
		CNT++;
		#endif
		return a.x >= ((_Float64)b);
	}

	template <typename T>
	friend bool operator>=(T const& a, __double_cnt const& b) {
		#ifdef CMP_AS_OP
		CNT++;
		#endif
		return ((_Float64)a) >= b.x;
	}

	operator _Float64() const {
		return x;
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

// #define int __int_cnt
// #define double __double_cnt

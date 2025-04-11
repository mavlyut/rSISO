#ifndef COUNT_OPS_H
#define COUNT_OPS_H

#include <istream>
#include <ostream>

#include "defines.h"

template <typename __type>
struct __count_type {
	__count_type() : x(0) {}

	__count_type(__type x) : x(x) {}

	__count_type(__count_type const& b) : x(b.x) {}

	__count_type& operator=(__count_type const& b) {
		x = b.x;
		return *this;
	}

	friend __count_type& operator+=(__count_type& a, __count_type const& b) {
		SUM_CNT++;
		a.x += b.x;
		return a;
	}

	template<typename T>
	friend __count_type& operator+=(__count_type& a, T const& b) {
		SUM_CNT++;
		a.x += ((__type)b);
		return a;
	}

	friend __count_type& operator-=(__count_type& a, __count_type const& b) {
		SUM_CNT++;
		a.x -= b.x;
		return a;
	}

	template <typename T>
	friend __count_type& operator-=(__count_type& a, T const& b) {
		a.x -= ((__type)b);
		return a;
	}

	friend __count_type& operator*=(__count_type& a, __count_type const& b) {
		MUL_CNT++;
		a.x *= b.x;
		return a;
	}

	template <typename T>
	friend __count_type& operator*=(__count_type& a, T const& b) {
		a.x *= ((__type)b);
		return a;
	}

	friend __count_type& operator/=(__count_type& a, __count_type const& b) {
		MUL_CNT++;
		a.x /= b.x;
		return a;
	}

	template <typename T>
	friend __count_type& operator/=(__count_type& a, T const& b) {
		a.x /= ((__type)b);
		return a;
	}

	friend __count_type operator+(__count_type const& a, __count_type const& b) {
		SUM_CNT++;
		return __count_type(a.x + b.x);
	}

	template <typename T>
	friend __count_type operator+(T const& a, __count_type const& b) {
		return ((__type)a) + b.x;
	}

	template <typename T>
	friend __count_type operator+(__count_type const& a, T const& b) {
		return a.x + ((__type)b);
	}

	friend __count_type operator-(__count_type const& a, __count_type const& b) {
		SUM_CNT++;
		return __count_type(a.x - b.x);
	}

	template <typename T>
	friend __count_type operator-(T const& a, __count_type const& b) {
		return ((__type)a) - b.x;
	}

	template <typename T>
	friend __count_type operator-(__count_type const& a, T const& b) {
		return a.x - ((__type)b);
	}

	friend __count_type operator*(__count_type const& a, __count_type const& b) {
		MUL_CNT++;
		return __count_type(a.x * b.x);
	}

	template <typename T>
	friend __count_type operator*(T const& a, __count_type const& b) {
		return ((__type)a) * b.x;
	}

	template <typename T>
	friend __count_type operator*(__count_type const& a, T const& b) {
		return a.x * ((__type)b);
	}

	friend __count_type operator/(__count_type const& a, __count_type const& b) {
		MUL_CNT++;
		return __count_type(a.x / b.x);
	}

	template <typename T>
	friend __count_type operator/(T const& a, __count_type const& b) {
		return ((__type)a) / b.x;
	}

	template <typename T>
	friend __count_type operator/(__count_type const& a, T const& b) {
		return a.x / ((__type)b);
	}

	friend bool operator==(__count_type const& a, __count_type const& b) {
		return a.x == b.x;
	}

	template <typename T>
	friend bool operator==(__count_type const& a, T const& b) {
		return a.x == ((__type)b);
	}

	template <typename T>
	friend bool operator==(T const& a, __count_type const& b) {
		return ((__type)a) == b.x;
	}

	friend bool operator!=(__count_type const& a, __count_type const& b) {
		return a.x != b.x;
	}

	template <typename T>
	friend bool operator!=(__count_type const& a, T const& b) {
		return a.x != ((__type)b);
	}

	template <typename T>
	friend bool operator!=(T const& a, __count_type const& b) {
		return ((__type)a) != b.x;
	}

	friend bool operator<(__count_type const& a, __count_type const& b) {
		return a.x < b.x;
	}

	template <typename T>
	friend bool operator<(__count_type const& a, T const& b) {
		return a.x < ((__type)b);
	}

	template <typename T>
	friend bool operator<(T const& a, __count_type const& b) {
		return ((__type)a) < b.x;
	}

	friend bool operator>(__count_type const& a, __count_type const& b) {
		return a.x > b.x;
	}

	template <typename T>
	friend bool operator>(__count_type const& a, T const& b) {
		return a.x > ((__type)b);
	}

	template <typename T>
	friend bool operator>(T const& a, __count_type const& b) {
		return ((__type)a) > b.x;
	}

	friend bool operator<=(__count_type const& a, __count_type const& b) {
		return a.x <= b.x;
	}

	template <typename T>
	friend bool operator<=(__count_type const& a, T const& b) {
		return a.x <= ((__type)b);
	}

	template <typename T>
	friend bool operator<=(T const& a, __count_type const& b) {
		return ((__type)a) <= b.x;
	}

	friend bool operator>=(__count_type const& a, __count_type const& b) {
		return a.x >= b.x;
	}

	template <typename T>
	friend bool operator>=(__count_type const& a, T const& b) {
		return a.x >= ((__type)b);
	}

	template <typename T>
	friend bool operator>=(T const& a, __count_type const& b) {
		return ((__type)a) >= b.x;
	}

	operator __type() const {
		return x;
	}

	__count_type operator+() const {
		return *this;
	}

	__count_type operator-() const {
		return __count_type(-x);
	}

	friend std::istream& operator>>(std::istream& in, __count_type& a) {
		return in >> a.x;
	}

	friend std::ostream& operator<<(std::ostream& out, __count_type const& a) {
		return out << a.x;
	}

private:
	__type x;
};

#define double __count_type<double>

#endif // COUNT_OPS_H

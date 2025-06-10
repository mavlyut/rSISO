#ifndef MATRIX_H
#define MATRIX_H

#include <vector>

#include "binvector.h"

namespace short_domain {
	using matrix = std::vector<binvector>;

	binvector multiply(binvector const& bv, matrix const& M);
	std::vector<unsigned> gauss(unsigned n, matrix& a);
	void full_gauss(unsigned n, matrix& a);
	void minimal_span_form(unsigned n, matrix& G);
	unsigned rg(unsigned n, matrix M);
	bool lin_indep(unsigned n, matrix const& M, binvector const& a);
	std::ostream& printmatrix(std::ostream& out, unsigned n, matrix const& M);
	std::ostream& printmatrix(std::ostream& out, std::string const& msg, unsigned n, matrix const& M);
	void _log_matrix(std::string const& msg, unsigned n, matrix const& M);
}

namespace long_domain {
	using matrix = std::vector<binvector>;

	binvector multiply(binvector const& bv, matrix const& M);
	std::vector<unsigned> gauss(matrix& a);
	void full_gauss(matrix& a);
	void minimal_span_form(matrix& G);
	unsigned rg(matrix M);
	bool lin_indep(matrix const& M, binvector const& a);
	std::ostream& printmatrix(std::ostream& out, unsigned n, matrix const& M);
	std::ostream& printmatrix(std::ostream& out, std::string const& msg, unsigned n, matrix const& M);
	void _log_matrix(std::string const& msg, unsigned n, matrix const& M);
}

#endif // MATRIX_H

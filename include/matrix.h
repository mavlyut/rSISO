#ifndef MATRIX_H
#define MATRIX_H

#include <vector>

#include "binvector.h"

using matrix = std::vector<binvector>;

std::vector<int> gauss(matrix& a);
void full_gauss(matrix& a);
std::vector<int> minimal_span_form(matrix& G);
unsigned rg(matrix M);
bool lin_indep(matrix const& M, binvector const& a);
std::ostream& operator<<(std::ostream& out, matrix const& a);

#endif // MATRIX_H

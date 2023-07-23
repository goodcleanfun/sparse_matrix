#ifndef DOUBLE_CSR_MATRIX_H
#define DOUBLE_CSR_MATRIX_H

#include "sparse_matrix.h"
#include "num_arrays/double_array.h"
#include "num_arrays/uint32_array.h"

COMPRESSED_SPARSE_MATRIX_INIT(double_csr_matrix, double, double_array, row)

#endif
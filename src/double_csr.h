#ifndef DOUBLE_CSR_MATRIX_H
#define DOUBLE_CSR_MATRIX_H

#include "csr.h"
#include "num_arrays/double_array.h"
#include "hashes/uint32_double.h"
#include "heaps/uint32_minheap.h"
#include "num_arrays/uint32_array.h"

COMPRESSED_SPARSE_MATRIX_INIT(double_csr_matrix, double, double_array, row)
CSR_DOT_SPARSE(double_csr_matrix, double, uint32_t, uint32_array, row, uint32_double, uint32_minheap)

#endif
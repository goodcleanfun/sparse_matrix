#ifndef FLOAT_CSR_MATRIX_H
#define FLOAT_CSR_MATRIX_H

#include "csr.h"
#include "num_arrays/float_array.h"
#include "hashes/uint32_float.h"
#include "heaps/uint32_minheap.h"
#include "num_arrays/uint32_array.h"

COMPRESSED_SPARSE_MATRIX_INIT(float_csr_matrix, float, float_array, row)
CSR_DOT_SPARSE(float_csr_matrix, float, uint32_t, uint32_array, row, uint32_float, uint32_minheap)

#endif
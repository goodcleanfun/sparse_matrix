#ifndef DOUBLE_CSR_MATRIX_H
#define DOUBLE_CSR_MATRIX_H

#include "num_array/double_array.h"
#include "hashes/uint32_double.h"
#include "heap/uint32_minheap.h"
#include "num_array/uint32_array.h"

#define SPARSE_TYPE_NAME double_csr_matrix
#define SPARSE_INDEX_TYPE uint32_t
#define SPARSE_INDEX_TYPE_NAME uint32
#define SPARSE_ORIENTATION row
#define SPARSE_INDEX_ARRAY_TYPE uint32_array
#define SPARSE_DATA_TYPE double
#define SPARSE_DATA_ARRAY_TYPE double_array
#include "compressed.h"
#undef SPARSE_TYPE_NAME
#undef SPARSE_INDEX_TYPE
#undef SPARSE_ORIENTATION
#undef SPARSE_INDEX_ARRAY_TYPE
#undef SPARSE_DATA_TYPE
#undef SPARSE_DATA_ARRAY_TYPE

#endif
#ifndef FLOAT_CSR_MATRIX_H
#define FLOAT_CSR_MATRIX_H

#include "num_arrays/float_array.h"
#include "hashes/uint32_float.h"
#include "heaps/uint32_minheap.h"
#include "num_arrays/uint32_array.h"

#define SPARSE_TYPE_NAME float_csr
#define SPARSE_INDEX_TYPE uint32_t
#define SPARSE_INDEX_TYPE_NAME uint32
#define SPARSE_ORIENTATION row
#define SPARSE_INDEX_ARRAY_TYPE uint32_array
#define SPARSE_DATA_TYPE float
#define SPARSE_DATA_ARRAY_TYPE float_array
#include "compressed.h"
#undef SPARSE_TYPE_NAME
#undef SPARSE_INDEX_TYPE
#undef SPARSE_ORIENTATION
#undef SPARSE_INDEX_ARRAY_TYPE
#undef SPARSE_DATA_TYPE
#undef SPARSE_DATA_ARRAY_TYPE

#endif
/*
sparse/compressed.h
-------------------

Dynamic compressed sparse row (CSR) or compressed sparse column (CSC)
sparse matrix.

These types of matrices arise often when representing graphs,
term-document matrices in text collections, etc.

The compressed sparse row format stores the following 3x7
dense matrix:

{   1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 1.0  }

with the following 3 arrays:

indptr = { 0, 2, 4, 6 }
indices = { 0, 4, 1, 5, 3, 6 }
data = { 1.0, 1.0, 1.0, 1.0, 2.0, 1.0 }

For a given row i, the indices indptr[i] through indptr[i+1]
denotes the number of nonzero columns in row i. The column
indices can be found at that contiguous location in indices
and the data values can be found at the same location in the
data array.

This implementation is generic, abstracting both the data type
and the differences between row and column (they are the same,
just switch out the meaning of indptr and indices).

Sparse matrices can be constructed dynamically by appending
entire rows/columns at a time or by appending individual values
to the last row/column and then finalizing the row/column, though
this requires the indices to arrive in sorted order.

Sparse matrix row/col iteration, row/col indexing, scalar arithmetic
and dot products with vectors dense matrices, and other sparse
matrices are efficient.
*/

#ifndef COMPRESSED_SPARSE_H
#define COMPRESSED_SPARSE_H

#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>

#include "file_utils/file_utils.h"

typedef enum {
    ROW,
    COLUMN
} sparse_matrix_orientation_t;


#define compressed_sparse_matrix_foreach(sp, indptr_var, start_var, length_var, code) do {      \
    __typeof__(indptr_var) ind_start = 0;                                                       \
    __typeof__(indptr_var) ind_end = 0;                                                         \
    __typeof__(indptr_var) *indptr = sp->indptr->a;                                             \
    __typeof__(indptr_var) m = sp->m;                                                           \
                                                                                                \
    for (__typeof__(indptr_var) i = 0; i < m; i++) {                                            \
        (indptr_var) = i;                                                                       \
        ind_start = indptr[i];                                                                  \
        ind_end = indptr[i + 1];                                                                \
        (start_var) = ind_start;                                                                \
        (length_var) = ind_end - ind_start;                                                     \
        code;                                                                                   \
    }                                                                                           \
} while(0);

#define compressed_sparse_matrix_foreach_value(sp, x_var, y_var, data_var, code) do {           \
    __typeof__(x_var) index, length;                                                            \
                                                                                                \
    compressed_sparse_matrix_foreach(sp, x_var, index, length, {                                \
        for (__typeof__(x_var) j = index; j < index + length; j++) {                            \
            (y_var) = sp->indices->a[j];                                                        \
            (data_var) = sp->data->a[j];                                                        \
            code;                                                                               \
        }                                                                                       \
    })                                                                                          \
} while(0);

#endif // COMPRESSED_SPARSE_H

#ifndef SPARSE_TYPE_NAME
#error "SPARSE_TYPE_NAME must be defined"
#endif

#ifndef SPARSE_INDEX_TYPE
#error "SPARSE_INDEX_TYPE must be defined"
#endif

#ifndef SPARSE_INDEX_TYPE_NAME
#define SPARSE_INDEX_TYPE_NAME SPARSE_INDEX_TYPE
#define SPARSE_INDEX_TYPE_NAME_DEFINED
#endif

#ifndef SPARSE_ORIENTATION
#error "SPARSE_ORIENTATION must be defined"
#endif

#if (SPARSE_ORIENTATION == row)
#define SPARSE_INVERSE_ORIENTATION col
#elif (SPARSE_ORIENTATION == col)
#define SPARSE_INVERSE_ORIENTATION row
#else
#error "SPARSE_ORIENTATION must be row or col"
#endif

#ifndef SPARSE_DATA_TYPE
#error "SPARSE_DATA_TYPE must be defined"
#endif

#ifndef SPARSE_DATA_TYPE_NAME
#define SPARSE_DATA_TYPE_NAME SPARSE_DATA_TYPE
#define SPARSE_DATA_TYPE_NAME_DEFINED
#endif

#define SPARSE_COMPRESSED_CONCAT_(a, b) a ## b
#define SPARSE_COMPRESSED_CONCAT(a, b) SPARSE_COMPRESSED_CONCAT_(a, b)
#define SPARSE_COMPRESSED_CONCAT3_(a, b, c) a ## b ## c
#define SPARSE_COMPRESSED_CONCAT3(a, b, c) SPARSE_COMPRESSED_CONCAT3_(a, b, c)
#define SPARSE_FUNC(name) SPARSE_COMPRESSED_CONCAT(SPARSE_TYPE_NAME, _##name)
#define SPARSE_TYPE_NAME_ORIENTATION SPARSE_COMPRESSED_CONCAT3(SPARSE_TYPE_NAME, _, SPARSE_ORIENTATION)
#define SPARSE_TYPE_NAME_INVERSE_ORIENTATION SPARSE_COMPRESSED_CONCAT3(SPARSE_TYPE_NAME, _, SPARSE_INVERSE_ORIENTATION)
#define SPARSE_INDICES_NAME SPARSE_COMPRESSED_CONCAT(SPARSE_TYPE_NAME_ORIENTATION, s)
#define SPARSE_INDICES_INVERSE_NAME SPARSE_COMPRESSED_CONCAT(SPARSE_TYPE_NAME_INVERSE_ORIENTATION, s)
#define SPARSE_INDICES_FUNC(name) SPARSE_COMPRESSED_CONCAT(SPARSE_INDICES_NAME, _##name)
#define SPARSE_ORIENTATION_FUNC(name) SPARSE_COMPRESSED_CONCAT3(SPARSE_TYPE_NAME, _##name##_, SPARSE_ORIENTATION)

#ifndef SPARSE_HASH_TYPE
#define SPARSE_HASH_TYPE SPARSE_COMPRESSED_CONCAT3(SPARSE_INDEX_TYPE_NAME, _, SPARSE_DATA_TYPE_NAME)
#define SPARSE_HASH_TYPE_DEFINED
#endif

#ifndef SPARSE_HEAP_TYPE
#define SPARSE_HEAP_TYPE SPARSE_COMPRESSED_CONCAT(SPARSE_INDEX_TYPE_NAME, _minheap)
#define SPARSE_HEAP_TYPE_DEFINED
#endif

#define HASH_NAME SPARSE_COMPRESSED_CONCAT(SPARSE_HASH_TYPE, _hash)
#define HASH_FUNC(name) SPARSE_COMPRESSED_CONCAT(SPARSE_HASH_TYPE, _hash_##name)
#define HEAP_FUNC(name) SPARSE_COMPRESSED_CONCAT(SPARSE_HEAP_TYPE, _##name)

#ifndef SPARSE_INDEX_ARRAY_TYPE
#define SPARSE_INDEX_ARRAY_TYPE SPARSE_COMPRESSED_CONCAT(SPARSE_INDEX_TYPE_NAME, _array)
#define SPARSE_INDEX_ARRAY_TYPE_DEFINED
#endif

#define INDEX_VECTOR_FUNC(name) SPARSE_COMPRESSED_CONCAT3(SPARSE_INDEX_TYPE_NAME, _vector, _##name)

#ifndef SPARSE_DATA_ARRAY_TYPE
#define SPARSE_DATA_ARRAY_TYPE SPARSE_COMPRESSED_CONCAT(SPARSE_DATA_TYPE_NAME, _array)
#define SPARSE_DATA_ARRAY_TYPE_DEFINED
#endif

#define DATA_VECTOR_FUNC(name) SPARSE_COMPRESSED_CONCAT3(SPARSE_DATA_TYPE_NAME, _vector, _##name)
#define INDEX_ARRAY_FUNC(name) SPARSE_COMPRESSED_CONCAT(SPARSE_INDEX_ARRAY_TYPE, _##name)
#define DATA_ARRAY_FUNC(name) SPARSE_COMPRESSED_CONCAT(SPARSE_DATA_ARRAY_TYPE, _##name)
#define FILE_INDEX_ARRAY_FUNC(name) SPARSE_COMPRESSED_CONCAT3(name##_, SPARSE_INDEX_TYPE_NAME, _array)
#define FILE_DATA_ARRAY_FUNC(name) SPARSE_COMPRESSED_CONCAT3(name##_, SPARSE_DATA_TYPE_NAME, _array)


typedef struct {
    size_t m;
    size_t n;
    SPARSE_INDEX_ARRAY_TYPE *indptr;
    SPARSE_INDEX_ARRAY_TYPE *indices;
    SPARSE_DATA_ARRAY_TYPE *data;
    size_t max_nnz;
} SPARSE_TYPE_NAME;


static inline void SPARSE_FUNC(destroy)(SPARSE_TYPE_NAME *self) {
    if (self == NULL) return;

    if (self->indptr != NULL) {
        INDEX_ARRAY_FUNC(destroy)(self->indptr);
    }

    if (self->indices != NULL) {
        INDEX_ARRAY_FUNC(destroy)(self->indices);
    }

    if (self->data != NULL) {
        DATA_ARRAY_FUNC(destroy)(self->data);
    }

    free(self);
}

static inline SPARSE_TYPE_NAME *SPARSE_FUNC(new_shape)(SPARSE_INDEX_TYPE m, SPARSE_INDEX_TYPE n) {
    SPARSE_TYPE_NAME *sp = calloc(1, sizeof(SPARSE_TYPE_NAME));
    if (sp == NULL) return NULL;
    sp->m = m;
    sp->n = n;
    sp->max_nnz = 0;
    sp->indptr = INDEX_ARRAY_FUNC(new_size)(m + 1);
    if (sp->indptr == NULL) {
        goto exit_sparse_allocated;
    }
    INDEX_ARRAY_FUNC(push)(sp->indptr, 0);

    sp->indices = INDEX_ARRAY_FUNC(new)();
    if (sp->indices == NULL) {
        goto exit_sparse_allocated;
    }

    sp->data = DATA_ARRAY_FUNC(new)();
    if (sp->data == NULL) {
        goto exit_sparse_allocated;
    }

    return sp;

exit_sparse_allocated:
    SPARSE_FUNC(destroy)(sp);
    return NULL;
}

static inline SPARSE_TYPE_NAME *SPARSE_FUNC(new)(void) {
    return SPARSE_FUNC(new_shape)(0, 0);
}


static inline void SPARSE_FUNC(clear)(SPARSE_TYPE_NAME *self) {
    INDEX_ARRAY_FUNC(clear)(self->indptr);
    INDEX_ARRAY_FUNC(push)(self->indptr, 0);

    INDEX_ARRAY_FUNC(clear)(self->indices);
    DATA_ARRAY_FUNC(clear)(self->data);
}

static inline void SPARSE_ORIENTATION_FUNC(finalize)(SPARSE_TYPE_NAME *self) {
    INDEX_ARRAY_FUNC(push)(self->indptr, (SPARSE_INDEX_TYPE)self->indices->n);
    if (self->indptr->n > self->m + 1) {
        self->m++;
    }
    SPARSE_INDEX_TYPE ind_len = self->indptr->a[self->m] - self->indptr->a[self->m - 1];
    if (ind_len > self->max_nnz) {
        self->max_nnz = ind_len;
    }
}

// rows for CSR, cols for CSC
static inline SPARSE_INDEX_TYPE SPARSE_INDICES_NAME(SPARSE_TYPE_NAME *self) {
    return self->m;
}


// cols for CSR, rows for CSC
static inline SPARSE_INDEX_TYPE SPARSE_INDICES_INVERSE_NAME(SPARSE_TYPE_NAME *self) {
    return self->n;
}

// nonzero elements overall (this assumes no deletions)
static inline SPARSE_INDEX_TYPE SPARSE_FUNC(nonzeros)(SPARSE_TYPE_NAME *self) {
    return self->data->n;
}

// row nonzero elements for CSR, col nonzero elements for CSC
static inline SPARSE_INDEX_TYPE SPARSE_ORIENTATION_FUNC(nonzeros)(SPARSE_TYPE_NAME *self, SPARSE_INDEX_TYPE ind) {
    if (ind > self->m) return 0;
    return self->indptr->a[ind + 1] - self->indptr->a[ind];
}


static inline SPARSE_INDEX_TYPE SPARSE_ORIENTATION_FUNC(len)(SPARSE_TYPE_NAME *self, SPARSE_INDEX_TYPE ind) {
    if (ind > self->m) return 0;
    return self->indptr->a[ind + 1] - self->indptr->a[ind];
}

static inline void SPARSE_FUNC(append)(SPARSE_TYPE_NAME *self, SPARSE_INDEX_TYPE index, SPARSE_DATA_TYPE val) {
    INDEX_ARRAY_FUNC(push)(self->indices, index);
    DATA_ARRAY_FUNC(push)(self->data, val);
    if (index >= self->n) self->n = index + 1;
}

static inline void SPARSE_ORIENTATION_FUNC(append)(SPARSE_TYPE_NAME *self, SPARSE_INDEX_TYPE *indices, SPARSE_DATA_TYPE *values, SPARSE_INDEX_TYPE n) {
    INDEX_ARRAY_FUNC(extend)(self->indices, indices, n);
    DATA_ARRAY_FUNC(extend)(self->data, values, n);
    SPARSE_INDEX_TYPE max_index = INDEX_VECTOR_FUNC(max)(indices, n);
    if (max_index >= self->n) self->n = max_index + 1;
    SPARSE_ORIENTATION_FUNC(finalize)(self);
}

static inline void SPARSE_FUNC(concat)(SPARSE_TYPE_NAME *self, SPARSE_TYPE_NAME *other) {
    SPARSE_ORIENTATION_FUNC(append)(self, other->indices->a, other->data->a, other->data->n);
    if (other->n > self->n) self->n = other->n;
    size_t ind = self->indptr->a[self->m - 1];
    for (SPARSE_INDEX_TYPE i = 0; i < other->m; i++) {
        INDEX_ARRAY_FUNC(push)(self->indptr, ind);
        size_t ind_len = SPARSE_ORIENTATION_FUNC(len)(other, i);
        ind += ind_len;
    }
}

#define SPARSE_INDEX_VALUE_TYPE SPARSE_COMPRESSED_CONCAT(SPARSE_TYPE_NAME, _index_value)
#define SPARSE_INDEX_VALUE_SORT_FUNC(name) SPARSE_COMPRESSED_CONCAT(name##_, SPARSE_INDEX_VALUE_TYPE)

typedef struct {
    SPARSE_INDEX_TYPE index;
    SPARSE_DATA_TYPE value;
} SPARSE_INDEX_VALUE_TYPE;

#define INTROSORT_TYPE SPARSE_INDEX_VALUE_TYPE
#define INTROSORT_LT(a, b) ((a).index < (b).index)
#include "sorting/introsort.h"
#undef INTROSORT_TYPE
#undef INTROSORT_LT

static inline void SPARSE_FUNC(sort_indices)(SPARSE_TYPE_NAME *self) {
    SPARSE_INDEX_TYPE x, ind_start, ind_len;

    size_t max_nnz = self->max_nnz;
    SPARSE_INDEX_VALUE_TYPE *tmp_index_values = malloc(sizeof(SPARSE_INDEX_VALUE_TYPE) * max_nnz);

    compressed_sparse_matrix_foreach(self, x, ind_start, ind_len, {
        for (size_t j = 0; j < ind_len; j++) {
            SPARSE_INDEX_VALUE_TYPE tmp;
            tmp.index = self->indices->a[ind_start + j];
            tmp.value = self->data->a[ind_start + j];
            tmp_index_values[j] = tmp;
        }
        SPARSE_INDEX_VALUE_SORT_FUNC(introsort)(ind_len, tmp_index_values);
        for (size_t j = 0; j < ind_len; j++) {
            self->indices->a[ind_start + j] = tmp_index_values[j].index;
            self->data->a[ind_start + j] = tmp_index_values[j].value;
        }
    })
    free(tmp_index_values);

}

static inline bool SPARSE_FUNC(dot_vector)(SPARSE_TYPE_NAME *self, SPARSE_DATA_TYPE *vec, size_t n, SPARSE_DATA_TYPE *result, size_t result_size) {
    if (n != self->n || result_size != self->m) return false;

    SPARSE_INDEX_TYPE x, ind_start, ind_len;
    SPARSE_DATA_TYPE val;
    SPARSE_INDEX_TYPE *indices = self->indices->a;
    SPARSE_DATA_TYPE *data = self->data->a;

    #pragma omp parallel for
    compressed_sparse_matrix_foreach(self, x, ind_start, ind_len, {
        SPARSE_DATA_TYPE sum = result[x];
        for (SPARSE_INDEX_TYPE y = ind_start; y < ind_start + ind_len; y++) {
            SPARSE_INDEX_TYPE col = indices[y];
            sum += data[y] * vec[col];
        }
        result[x] = sum;
    })
    return true;
}

/* rows_dot_vector or cols_dot_vector */
static inline bool SPARSE_INDICES_FUNC(dot_vector)(SPARSE_TYPE_NAME *self, SPARSE_INDEX_TYPE *xs, size_t m, SPARSE_DATA_TYPE *vec, size_t n, SPARSE_DATA_TYPE *result, size_t result_size) {
    if (n != self->n || m != result_size) return false;

    SPARSE_INDEX_TYPE *indptr = self->indptr->a;
    SPARSE_INDEX_TYPE *indices = self->indices->a;
    SPARSE_DATA_TYPE *data = self->data->a;

    #pragma omp parallel for
    for (SPARSE_INDEX_TYPE i = 0; i < m; i++) {
        SPARSE_INDEX_TYPE index = xs[i];

        SPARSE_DATA_TYPE sum = result[i];
        if (index >= self->m) return false;

        for (SPARSE_INDEX_TYPE j = indptr[index]; j < indptr[index + 1]; j++) {
            SPARSE_INDEX_TYPE col = indices[j];
            sum += data[j] * vec[col];
        }

        result[i] = sum;
    }
    return true;
}


static inline bool SPARSE_FUNC(dot_sparse)(
    SPARSE_TYPE_NAME *a,
    SPARSE_TYPE_NAME *b,
    SPARSE_TYPE_NAME *c,
    bool sort_indices
) {
    SPARSE_INDEX_TYPE a_rows = a->m;
    SPARSE_INDEX_TYPE a_cols = a->n;
    SPARSE_INDEX_TYPE b_rows = b->m;
    SPARSE_INDEX_TYPE b_cols = b->n;

    if (a_cols != b_rows) {
        return false;
    }

    SPARSE_INDEX_TYPE *a_indptr = a->indptr->a;
    SPARSE_INDEX_TYPE *a_indices = a->indices->a;
    SPARSE_DATA_TYPE *a_values = a->data->a;
    SPARSE_INDEX_TYPE *b_indptr = b->indptr->a;
    SPARSE_INDEX_TYPE *b_indices = b->indices->a;
    SPARSE_DATA_TYPE *b_values = b->data->a;

    if (c == NULL) return false;

    if (c->indices->n > 0) {
        SPARSE_FUNC(clear)(c);
    }
    if (c->indptr->n > 0) {
        INDEX_ARRAY_FUNC(clear)(c->indptr);
        INDEX_ARRAY_FUNC(push)(c->indptr, 0);
    }
    if (c->data->n > 0) {
        DATA_ARRAY_FUNC(clear)(c->data);
    }
    c->m = a_rows;
    c->n = b_cols;
    c->max_nnz = 0;

    SPARSE_INDEX_TYPE nnz = 0;

    HASH_NAME *sums = HASH_FUNC(new)();
    SPARSE_INDEX_ARRAY_TYPE *sorted_indices;
    if (sort_indices) {
        sorted_indices = INDEX_ARRAY_FUNC(new)();
    }

    for (SPARSE_INDEX_TYPE i = 0; i < a_rows; i++) {
        SPARSE_INDEX_TYPE head = 0;
        SPARSE_INDEX_TYPE length = 0;

        SPARSE_INDEX_TYPE a_ind_start = a_indptr[i];
        SPARSE_INDEX_TYPE a_ind_end = a_indptr[i+1];

        for (SPARSE_INDEX_TYPE j = a_ind_start; j < a_ind_end; j++) {
            SPARSE_INDEX_TYPE a_index = a_indices[j];
            SPARSE_DATA_TYPE a_value = a_values[j];

            SPARSE_INDEX_TYPE b_ind_start = b_indptr[a_index];
            SPARSE_INDEX_TYPE b_ind_end = b_indptr[a_index+1];

            for (SPARSE_INDEX_TYPE k = b_ind_start; k < b_ind_end; k++) {
                SPARSE_INDEX_TYPE b_index = b_indices[k];
                SPARSE_DATA_TYPE b_value = b_values[k];
                bool index_seen = false;
                HASH_FUNC(incr_by_exists)(sums, b_index, a_value * b_value, &index_seen);

                if (sort_indices && !index_seen) {
                    HEAP_FUNC(push)(sorted_indices, b_index);
                }
            }
        }

        if (sort_indices) {
            SPARSE_INDEX_TYPE nonzero_index;
            while (HEAP_FUNC(pop)(sorted_indices, &nonzero_index)) {
                SPARSE_DATA_TYPE sum_value;
                if (!(HASH_FUNC(get)(sums, nonzero_index, &sum_value))) {
                    break;
                }
                SPARSE_FUNC(append)(c, nonzero_index, sum_value);
                nnz++;
            }
        } else {
            size_t i = 0;
            SPARSE_INDEX_TYPE key;
            SPARSE_DATA_TYPE value;
            kh_foreach(sums, key, value, {
                SPARSE_FUNC(append)(c, key, value);
                nnz++;
            });
        }
        HASH_FUNC(clear)(sums);
        SPARSE_ORIENTATION_FUNC(finalize)(c);
    }
    HASH_FUNC(destroy)(sums);
    if (sort_indices) {
        INDEX_ARRAY_FUNC(destroy)(sorted_indices);
    }

    return true;
}




static inline SPARSE_TYPE_NAME *SPARSE_FUNC(read)(FILE *f) {
    SPARSE_TYPE_NAME *sp = malloc(sizeof(SPARSE_TYPE_NAME));
    if (sp == NULL) return NULL;

    sp->indptr = NULL;
    sp->indices = NULL;
    sp->data = NULL;

    uint64_t m = 0;
    uint64_t n = 0;

    if (!file_read_uint64(f, &m) ||
        !file_read_uint64(f, &n)) {
        goto exit_read_sparse_allocated;
    }

    sp->m = (size_t)m;
    sp->n = (size_t)n;

    uint64_t len_indptr;

    if (!file_read_uint64(f, &len_indptr)) {
        goto exit_read_sparse_allocated;
    }

    SPARSE_INDEX_ARRAY_TYPE *indptr = INDEX_ARRAY_FUNC(new_size)((size_t)len_indptr);
    if (indptr == NULL) {
        goto exit_read_sparse_allocated;
    }

    if (!FILE_INDEX_ARRAY_FUNC(file_read)(f, indptr->a, len_indptr)) {
        INDEX_ARRAY_FUNC(destroy)(indptr);
        goto exit_read_sparse_allocated;
    }

    indptr->n = (size_t)len_indptr;
    sp->indptr = indptr;

    uint64_t len_indices;

    if (!file_read_uint64(f, &len_indices)) {
        goto exit_read_sparse_allocated;
    }

    SPARSE_INDEX_ARRAY_TYPE *indices = INDEX_ARRAY_FUNC(new_size)((size_t)len_indices);
    if (indices == NULL) {
        goto exit_read_sparse_allocated;
    }

    if (!FILE_INDEX_ARRAY_FUNC(file_read)(f, indices->a, len_indices)) {
        INDEX_ARRAY_FUNC(destroy)(indices);
        goto exit_read_sparse_allocated;
    }

    indices->n = (size_t)len_indices;
    sp->indices = indices;

    uint64_t len_data;

    if (!file_read_uint64(f, &len_data)) {
        goto exit_read_sparse_allocated;
    }

    SPARSE_DATA_ARRAY_TYPE *data = DATA_ARRAY_FUNC(new_size)((size_t)len_data);
    if (data == NULL) {
        goto exit_read_sparse_allocated;
    }

    if (!FILE_DATA_ARRAY_FUNC(file_read)(f, data->a, len_data)) {
        DATA_ARRAY_FUNC(destroy)(data);
        goto exit_read_sparse_allocated;
    }

    data->n = (size_t)len_data;
    sp->data = data;

    return sp;

exit_read_sparse_allocated:
    SPARSE_FUNC(destroy)(sp);
    return NULL;
}

static inline bool SPARSE_FUNC(write)(SPARSE_TYPE_NAME *self, FILE *f) {
    if (self == NULL || self->indptr == NULL || self->indices == NULL || self->data == NULL) {
        return false;
    }

    if (!file_write_uint64(f, self->m) ||
        !file_write_uint64(f, self->n)) {
        return false;
    }

    uint64_t len_indptr = self->indptr->n;

    if (!file_write_uint64(f, len_indptr)) {
        return false;
    }

    if (!FILE_INDEX_ARRAY_FUNC(file_write)(f, self->indptr->a, len_indptr)) {
        return false;
    }

    uint64_t len_indices = (uint64_t)self->indices->n;

    if (!file_write_uint64(f, len_indices)) {
        return false;
    }

    if (!FILE_INDEX_ARRAY_FUNC(file_write)(f, self->indices->a, len_indices)) {
        return false;
    }

    uint64_t len_data = (uint64_t)self->data->n;

    if (!file_write_uint64(f, len_data)) {
        return false;
    }

    if (!FILE_DATA_ARRAY_FUNC(file_write)(f, self->data->a, len_data)) {
        return false;
    }

    return true;
}



#ifndef SPARSE_TYPE_NAME
#error "SPARSE_TYPE_NAME must be defined"
#endif

#ifndef SPARSE_INDEX_TYPE
#error "SPARSE_INDEX_TYPE must be defined"
#endif

#ifndef SPARSE_INDEX_TYPE_NAME
#define SPARSE_INDEX_TYPE_NAME SPARSE_INDEX_TYPE
#define SPARSE_INDEX_TYPE_NAME_DEFINED
#endif

#ifndef SPARSE_ORIENTATION
#error "SPARSE_ORIENTATION must be defined"
#endif

#ifndef SPARSE_DATA_TYPE
#error "SPARSE_DATA_TYPE must be defined"
#endif

#ifndef SPARSE_DATA_TYPE_NAME
#define SPARSE_DATA_TYPE_NAME SPARSE_DATA_TYPE
#define SPARSE_DATA_TYPE_NAME_DEFINED
#endif

#ifndef SPARSE_HASH_TYPE
#define SPARSE_HASH_TYPE SPARSE_COMPRESSED_CONCAT3(SPARSE_INDEX_TYPE_NAME, _, SPARSE_DATA_TYPE_NAME)
#define SPARSE_HASH_TYPE_DEFINED
#endif

#ifndef SPARSE_HEAP_TYPE
#define SPARSE_HEAP_TYPE SPARSE_COMPRESSED_CONCAT(SPARSE_INDEX_TYPE_NAME, _minheap)
#define SPARSE_HEAP_TYPE_DEFINED
#endif

#undef SPARSE_INVERSE_ORIENTATION

#undef SPARSE_COMPRESSED_CONCAT_
#undef SPARSE_COMPRESSED_CONCAT
#undef SPARSE_COMPRESSED_CONCAT3_
#undef SPARSE_COMPRESSED_CONCAT3
#undef SPARSE_FUNC
#undef SPARSE_TYPE_NAME_ORIENTATION
#undef SPARSE_TYPE_NAME_INVERSE_ORIENTATION
#undef SPARSE_INDICES_NAME
#undef SPARSE_INDICES_INVERSE_NAME
#undef SPARSE_INDICES_FUNC
#undef SPARSE_ORIENTATION_FUNC

#ifdef SPARSE_INDEX_TYPE_NAME_DEFINED
#undef SPARSE_INDEX_TYPE_NAME
#undef SPARSE_INDEX_TYPE_NAME_DEFINED
#endif

#ifdef SPARSE_DATA_TYPE_NAME_DEFINED
#undef SPARSE_DATA_TYPE_NAME
#undef SPARSE_DATA_TYPE_NAME_DEFINED
#endif

#ifdef SPARSE_HASH_TYPE_DEFINED
#undef SPARSE_HASH_TYPE
#undef SPARSE_HASH_TYPE_DEFINED
#endif

#ifdef SPARSE_HEAP_TYPE_DEFINED
#undef SPARSE_HEAP_TYPE
#undef SPARSE_HEAP_TYPE_DEFINED
#endif

#undef SPARSE_INDEX_VALUE_TYPE
#undef HASH_NAME
#undef HASH_FUNC
#undef INDEX_ARRAY_FUNC
#undef DATA_ARRAY_FUNC
#undef FILE_INDEX_ARRAY_FUNC
#undef FILE_DATA_ARRAY_FUNC

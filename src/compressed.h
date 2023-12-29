/*
sparse_matrices/compressed.h
----------------------------

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

#ifndef COMPRESSED_SPARSE_MATRIX_H
#define COMPRESSED_SPARSE_MATRIX_H

#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>

#include "vector/vector.h"
#include "file_utils/file_utils.h"
#include "sorting/introsort.h"

#define ks_lt_index_value(a, b) ((a).index < (b).index)

typedef enum {
    ROW,
    COLUMN
} sparse_matrix_orientation_t;


#define compressed_sparse_matrix_foreach(sp, indptr_var, index_var, length_var, code) {                                 \
    __typeof__(indptr_var) _ind_start = 0;                                                                              \
    __typeof__(indptr_var) _ind_end = 0;                                                                                \
    __typeof__(indptr_var) *_indptr = sp->indptr->a;                                                                    \
    __typeof__(indptr_var) _m = sp->m;                                                                                  \
                                                                                                                        \
    for (__typeof__(indptr_var) _i = 0; _i < _m; _i++) {                                                                \
        (indptr_var) = _i;                                                                                              \
        _ind_start = _indptr[_i];                                                                                       \
        _ind_end = _indptr[_i + 1];                                                                                     \
        (index_var) = _ind_start;                                                                                       \
        (length_var) = _ind_end - _ind_start;                                                                           \
        code;                                                                                                           \
    }                                                                                                                   \
}

#define compressed_sparse_matrix_foreach_value(sp, x_var, y_var, data_var, code) {                                      \
    __typeof__(y_var) *_indices = sp->indices->a;                                                                       \
    __typeof__(data_var) *_data = sp->data->a;                                                                          \
    __typeof__(x_var) _index, _length;                                                                                  \
    compressed_sparse_matrix_foreach(sp, x_var, _index, _length, {                                                      \
        for (__typeof__(x_var) _j = _index; _j < _index + _length; _j++) {                                              \
            (y_var) = _indices[_j];                                                                                     \
            (data_var) = _data[_j];                                                                                     \
            code;                                                                                                       \
        }                                                                                                               \
    })                                                                                                                  \
}                                                                                                                       \


#define COMPRESSED_SPARSE_MATRIX_INIT_SIZE(name, data_type, data_array_type, index_size, index_type, index_array_type, index_name)      \
    typedef struct {                                                                                                        \
        index_type m;                                                                                                       \
        index_type n;                                                                                                       \
        index_array_type *indptr;                                                                                           \
        index_array_type *indices;                                                                                          \
        data_array_type *data;                                                                                              \
    } name;                                                                                                                 \
                                                                                                                            \
    static inline void name##_destroy(name *self) {                                                                         \
        if (self == NULL) return;                                                                                           \
                                                                                                                            \
        if (self->indptr != NULL) {                                                                                         \
            index_array_type##_destroy(self->indptr);                                                                       \
        }                                                                                                                   \
                                                                                                                            \
        if (self->indices != NULL) {                                                                                        \
            index_array_type##_destroy(self->indices);                                                                      \
        }                                                                                                                   \
                                                                                                                            \
        if (self->data != NULL) {                                                                                           \
            data_array_type##_destroy(self->data);                                                                          \
        }                                                                                                                   \
                                                                                                                            \
        free(self);                                                                                                         \
    }                                                                                                                       \
                                                                                                                            \
    static inline name *name##_new_shape(index_type m, index_type n) {                                                      \
        name *matrix = calloc(1, sizeof(name));                                                                             \
        if (matrix == NULL) return NULL;                                                                                    \
        matrix->m = m;                                                                                                      \
        matrix->n = n;                                                                                                      \
        matrix->indptr = index_array_type##_new_size(m + 1);                                                                \
        if (matrix->indptr == NULL) {                                                                                       \
            goto exit_##name##_allocated;                                                                                   \
        }                                                                                                                   \
        index_array_type##_push(matrix->indptr, 0);                                                                         \
                                                                                                                            \
        matrix->indices = index_array_type##_new();                                                                         \
        if (matrix->indices == NULL) {                                                                                      \
            goto exit_##name##_allocated;                                                                                   \
        }                                                                                                                   \
                                                                                                                            \
        matrix->data = data_array_type##_new();                                                                             \
        if (matrix->data == NULL) {                                                                                         \
            goto exit_##name##_allocated;                                                                                   \
        }                                                                                                                   \
                                                                                                                            \
        return matrix;                                                                                                      \
                                                                                                                            \
    exit_##name##_allocated:                                                                                                \
        name##_destroy(matrix);                                                                                             \
        return NULL;                                                                                                        \
    }                                                                                                                       \
                                                                                                                            \
    name *name##_new(void) {                                                                                                \
        return name##_new_shape(0, 0);                                                                                      \
    }                                                                                                                       \
                                                                                                                            \
                                                                                                                            \
    static inline void name##_clear(name *self) {                                                                           \
        index_array_type##_clear(self->indptr);                                                                             \
        index_array_type##_push(self->indptr, 0);                                                                           \
                                                                                                                            \
        index_array_type##_clear(self->indices);                                                                            \
        data_array_type##_clear(self->data);                                                                                \
    }                                                                                                                       \
                                                                                                                            \
    static inline void name##_finalize_##index_name(name *self) {                                                           \
        index_array_type##_push(self->indptr, (index_type)self->indices->n);                                                \
        if (self->indptr->n > self->m + 1) {                                                                                \
            self->m++;                                                                                                      \
        }                                                                                                                   \
    }                                                                                                                       \
                                                                                                                            \
    static inline void name##_append(name *self, index_type index, data_type val) {                                         \
        index_array_type##_push(self->indices, index);                                                                      \
        data_array_type##_push(self->data, val);                                                                            \
        if (index >= self->n) self->n = index + 1;                                                                          \
    }                                                                                                                       \
                                                                                                                            \
    static inline void name##_append_##index_name(name *self, index_type *indices, data_type *values, index_type n) {       \
        for (index_type i = 0; i < n; i++) {                                                                                \
            name##_append(self, indices[i], values[i]);                                                                     \
        }                                                                                                                   \
        name##_finalize_##index_name(self);                                                                                 \
    }                                                                                                                       \
                                                                                                                            \
    typedef struct name##_index_value {                                                                                     \
        index_type index;                                                                                                   \
        data_type val;                                                                                                      \
    } name##_index_value_t;                                                                                                 \
                                                                                                                            \
    VECTOR_INIT(name##_index_value_array, name##_index_value_t)                                                             \
                                                                                                                            \
    INTROSORT_INIT(name##_index_value_array, name##_index_value_t, ks_lt_index_value)                                       \
                                                                                                                            \
    static inline void name##_sort_indices(name *self) {                                                                    \
        index_type x, ind_start, ind_len, j;                                                                                \
                                                                                                                            \
        name##_index_value_array *index_vals = name##_index_value_array_new();                                              \
                                                                                                                            \
        compressed_sparse_matrix_foreach(self, x, ind_start, ind_len, {                                                     \
            for (j = ind_start; j < ind_start + ind_len; j++) {                                                             \
                name##_index_value_array_push(index_vals, (name##_index_value_t){self->indices->a[j], self->data->a[j]});   \
            }                                                                                                               \
            ks_introsort(name##_index_value_array, index_vals->n, index_vals->a);                                           \
                                                                                                                            \
            for (j = 0; j < index_vals->n; j++) {                                                                           \
                name##_index_value_t ind_val = index_vals->a[j];                                                            \
                self->indices->a[ind_start + j] = ind_val.index;                                                            \
                self->data->a[ind_start + j] = ind_val.val;                                                                 \
            }                                                                                                               \
        })                                                                                                                  \
                                                                                                                            \
    }                                                                                                                       \
                                                                                                                            \
                                                                                                                            \
    static inline bool name##_dot_vector(name *self, data_type *vec, index_type n, data_type *result) {                     \
        if (n != self->n) return false;                                                                                     \
                                                                                                                            \
        index_type x, ind_start, ind_len;                                                                                   \
        data_type val;                                                                                                      \
        data_type *data = self->data->a;                                                                                    \
                                                                                                                            \
        compressed_sparse_matrix_foreach(self, x, ind_start, ind_len, {                                                     \
            data_type sum = result[x];                                                                                      \
            for (index_type y = ind_start; y < ind_start + ind_len; y++) {                                                  \
                sum += data[y] * vec[y];                                                                                    \
            }                                                                                                               \
            result[x] = sum;                                                                                                \
        })                                                                                                                  \
        return true;                                                                                                        \
    }                                                                                                                       \
                                                                                                                            \
    static inline bool name##_##index_name##s_dot_vector(name *self, index_type *xs, index_type m, data_type *vec, index_type n, data_type *result) { \
        if (n != self->n) return false;                                                                                     \
                                                                                                                            \
        index_type *indptr = self->indptr->a;                                                                               \
        index_type *indices = self->indices->a;                                                                             \
        data_type *data = self->data->a;                                                                                    \
                                                                                                                            \
        for (index_type i = 0; i < m; i++) {                                                                                \
            index_type index = indices[i];                                                                                  \
                                                                                                                            \
            data_type sum = result[i];                                                                                      \
            if (index >= self->m) return false;                                                                             \
                                                                                                                            \
            for (index_type j = indptr[index]; j < indptr[index+1]; j++) {                                                  \
                sum += data[j] * vec[indices[j]];                                                                           \
            }                                                                                                               \
                                                                                                                            \
            result[i] = sum;                                                                                                \
                                                                                                                            \
        }                                                                                                                   \
        return 0;                                                                                                           \
    }                                                                                                                       \
                                                                                                                            \
    static inline name *name##_read(FILE *f) {                                                                              \
        name *sp = malloc(sizeof(name));                                                                                    \
        if (sp == NULL) return NULL;                                                                                        \
                                                                                                                            \
        sp->indptr = NULL;                                                                                                  \
        sp->indices = NULL;                                                                                                 \
        sp->data = NULL;                                                                                                    \
                                                                                                                            \
        if (!file_read_uint##index_size(f, &sp->m) ||                                                                       \
            !file_read_uint##index_size(f, &sp->n)) {                                                                       \
            goto exit_##name##_allocated;                                                                                   \
        }                                                                                                                   \
                                                                                                                            \
        uint64_t len_indptr;                                                                                                \
                                                                                                                            \
        if (!file_read_uint64(f, &len_indptr)) {                                                                            \
            goto exit_##name##_allocated;                                                                                   \
        }                                                                                                                   \
                                                                                                                            \
        index_array_type *indptr = index_array_type##_new_size((size_t)len_indptr);                                         \
        if (indptr == NULL) {                                                                                               \
            goto exit_##name##_allocated;                                                                                   \
        }                                                                                                                   \
                                                                                                                            \
        if (!file_read_##index_array_type(f, indptr->a, len_indptr)) {                                                      \
            index_array_type##_destroy(indptr);                                                                             \
            goto exit_##name##_allocated;                                                                                   \
        }                                                                                                                   \
                                                                                                                            \
        indptr->n = (size_t)len_indptr;                                                                                     \
        sp->indptr = indptr;                                                                                                \
                                                                                                                            \
        uint64_t len_indices;                                                                                               \
                                                                                                                            \
        if (!file_read_uint64(f, &len_indices)) {                                                                           \
            goto exit_##name##_allocated;                                                                                   \
        }                                                                                                                   \
                                                                                                                            \
        index_array_type *indices = index_array_type##_new_size((size_t)len_indices);                                       \
        if (indices == NULL) {                                                                                              \
            goto exit_##name##_allocated;                                                                                   \
        }                                                                                                                   \
                                                                                                                            \
        if (!file_read_##index_array_type(f, indices->a, len_indices)) {                                                    \
            index_array_type##_destroy(indices);                                                                            \
            goto exit_##name##_allocated;                                                                                   \
        }                                                                                                                   \
                                                                                                                            \
        indices->n = (size_t)len_indices;                                                                                   \
        sp->indices = indices;                                                                                              \
                                                                                                                            \
        uint64_t len_data;                                                                                                  \
                                                                                                                            \
        if (!file_read_uint64(f, &len_data)) {                                                                              \
            goto exit_##name##_allocated;                                                                                   \
        }                                                                                                                   \
                                                                                                                            \
        data_array_type *data = data_array_type##_new_size((size_t)len_data);                                               \
        if (data == NULL) {                                                                                                 \
            goto exit_##name##_allocated;                                                                                   \
        }                                                                                                                   \
                                                                                                                            \
        if (!file_read_##data_type##_array(f, data->a, len_data)) {                                                         \
            data_array_type##_destroy(data);                                                                                \
            goto exit_##name##_allocated;                                                                                   \
        }                                                                                                                   \
                                                                                                                            \
        data->n = (size_t)len_data;                                                                                         \
        sp->data = data;                                                                                                    \
                                                                                                                            \
        return sp;                                                                                                          \
                                                                                                                            \
    exit_##name##_allocated:                                                                                                \
        name##_destroy(sp);                                                                                                 \
        return NULL;                                                                                                        \
    }                                                                                                                       \
                                                                                                                            \
    static inline bool name##_write(name *self, FILE *f) {                                                                  \
        if (self == NULL || self->indptr == NULL || self->indices == NULL || self->data == NULL) {                          \
            return false;                                                                                                   \
        }                                                                                                                   \
                                                                                                                            \
        if (!file_write_uint##index_size(f, self->m) ||                                                                     \
            !file_write_uint##index_size(f, self->n)) {                                                                     \
            return false;                                                                                                   \
        }                                                                                                                   \
                                                                                                                            \
        uint64_t len_indptr = self->indptr->n;                                                                              \
                                                                                                                            \
        if (!file_write_uint64(f, len_indptr)) {                                                                            \
            return false;                                                                                                   \
        }                                                                                                                   \
                                                                                                                            \
        for (int i = 0; i < len_indptr; i++) {                                                                              \
            if (!file_write_uint##index_size(f, self->indptr->a[i])) {                                                      \
                return false;                                                                                               \
            }                                                                                                               \
        }                                                                                                                   \
                                                                                                                            \
        uint64_t len_indices = (uint64_t)self->indices->n;                                                                  \
                                                                                                                            \
        if (!file_write_uint64(f, len_indices)) {                                                                           \
            return false;                                                                                                   \
        }                                                                                                                   \
                                                                                                                            \
        for (int i = 0; i < len_indices; i++) {                                                                             \
            if (!file_write_uint##index_size(f, self->indices->a[i])) {                                                     \
                return false;                                                                                               \
            }                                                                                                               \
        }                                                                                                                   \
                                                                                                                            \
        uint64_t len_data = (uint64_t)self->data->n;                                                                        \
                                                                                                                            \
        if (!file_write_uint64(f, len_data)) {                                                                              \
            return false;                                                                                                   \
        }                                                                                                                   \
                                                                                                                            \
        for (int i = 0; i < len_data; i++) {                                                                                \
            if (!file_write_##data_type(f, self->data->a[i])) {                                                             \
                return false;                                                                                               \
            }                                                                                                               \
        }                                                                                                                   \
                                                                                                                            \
        return true;                                                                                                        \
    }                                                                                                                       \



#define COMPRESSED_SPARSE_MATRIX_INIT(name, data_type, data_array_type, index_name) \
COMPRESSED_SPARSE_MATRIX_INIT_SIZE(name, data_type, data_array_type, 32, uint32_t, uint32_array, index_name)


#endif

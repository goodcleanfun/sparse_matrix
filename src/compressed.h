/*
sparse_matrix.h
---------------

Dynamic compressed sparse row (CSR) sparse matrix.

A sparse matrix is a matrix where an overwhelming number of the
entries are 0. Thus we get significant space/time advantages
from only storing the non-zero values.

These data_types of matrices arise often when representing graphs
and term-document matrices in text collections.

The compressed sparse row format stores the following 3x4
dense matrix:

{   1.0, 0.0, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 2.0  }

with the following 3 arrays:

indptr = { 0, 1, 2, 3 }
indices = { 0, 1, 3 }
data = { 1.0, 1.0, 2.0 }

For a given row i, the indices indptr[i] through indptr[i+1]
denotes the number of nonzero columns in row i. The column
indices can be found at that contiguous location in indices
and the data values can be found at the same location in the
data array.

Sparse matrix row iteration, row indexing, scalar arithmetic
and dot products with vectors and dense matrices are efficient.
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
    size_t _m = sp->m;                                                                                                  \
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


#define COMPRESSED_SPARSE_MATRIX_INIT_SIZE(name, data_type, data_array_type, index_size, index_name)                        \
    typedef struct {                                                                                                        \
        uint##index_size##_t m;                                                                                            \
        uint##index_size##_t n;                                                                                             \
        uint##index_size##_array *indptr;                                                                                  \
        uint##index_size##_array *indices;                                                                                  \
        data_array_type *data;                                                                                              \
    } name##_t;                                                                                                             \
                                                                                                                            \
    static inline void name##_destroy(name##_t *self) {                                                                     \
        if (self == NULL) return;                                                                                           \
                                                                                                                            \
        if (self->indptr != NULL) {                                                                                         \
            uint##index_size##_array_destroy(self->indptr);                                                                \
        }                                                                                                                   \
                                                                                                                            \
        if (self->indices != NULL) {                                                                                        \
            uint##index_size##_array_destroy(self->indices);                                                                \
        }                                                                                                                   \
                                                                                                                            \
        if (self->data != NULL) {                                                                                           \
            data_array_type##_destroy(self->data);                                                                          \
        }                                                                                                                   \
                                                                                                                            \
        free(self);                                                                                                         \
    }                                                                                                                       \
                                                                                                                            \
    static inline name##_t *name##_new_shape(size_t m, size_t n) {                                                          \
        name##_t *matrix = calloc(1, sizeof(name##_t));                                                                     \
        if (matrix == NULL) return NULL;                                                                                    \
        matrix->m = m;                                                                                                      \
        matrix->n = n;                                                                                                      \
        matrix->indptr = uint##index_size##_array_new_size(m + 1);                                                         \
        if (matrix->indptr == NULL) {                                                                                       \
            goto exit_##name##_allocated;                                                                                   \
        }                                                                                                                   \
        uint##index_size##_array_push(matrix->indptr, 0);                                                                  \
                                                                                                                            \
        matrix->indices = uint##index_size##_array_new();                                                                   \
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
    name##_t *name##_new(void) {                                                                                            \
        return name##_new_shape(0, 0);                                                                                      \
    }                                                                                                                       \
                                                                                                                            \
                                                                                                                            \
    static inline void name##_clear(name##_t *self) {                                                                       \
        uint##index_size##_array_clear(self->indptr);                                                                      \
        uint##index_size##_array_push(self->indptr, 0);                                                                    \
                                                                                                                            \
        uint##index_size##_array_clear(self->indices);                                                                      \
        data_array_type##_clear(self->data);                                                                                \
    }                                                                                                                       \
                                                                                                                            \
    static inline void name##_finalize_##index_name(name##_t *self) {                                                       \
        uint##index_size##_array_push(self->indptr, (uint##index_size##_t)self->indices->n);                               \
        if (self->indptr->n > self->m + 1) {                                                                                \
            self->m++;                                                                                                      \
        }                                                                                                                   \
    }                                                                                                                       \
                                                                                                                            \
    static inline void name##_append(name##_t *self, uint##index_size##_t index, data_type val) {                           \
        uint##index_size##_array_push(self->indices, index);                                                                \
        data_array_type##_push(self->data, val);                                                                            \
        if (index >= self->n) self->n = index + 1;                                                                          \
    }                                                                                                                       \
                                                                                                                            \
    static inline void name##_append_##index_name(name##_t *self, uint##index_size##_t *indices, data_type *values, size_t n) {    \
        for (size_t i = 0; i < n; i++) {                                                                                    \
            name##_append(self, indices[i], values[i]);                                                                     \
        }                                                                                                                   \
        name##_finalize_##index_name(self);                                                                                 \
    }                                                                                                                       \
                                                                                                                            \
    typedef struct name##_index_value {                                                                                     \
        uint##index_size##_t index;                                                                                         \
        data_type val;                                                                                                      \
    } name##_index_value_t;                                                                                                 \
                                                                                                                            \
    VECTOR_INIT(name##_index_value_array, name##_index_value_t)                                                             \
                                                                                                                            \
    INTROSORT_INIT(name##_index_value_array, name##_index_value_t, ks_lt_index_value)                                       \
                                                                                                                            \
    static inline void name##_sort_indices(name##_t *self) {                                                                \
        uint##index_size##_t x, ind_start, ind_len, j;                                                                     \
                                                                                                                            \
        name##_index_value_array *index_vals = name##_index_value_array_new();                                              \
                                                                                                                            \
        compressed_sparse_matrix_foreach(self, x, ind_start, ind_len, {                                               \
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
    static inline int name##_dot_vector(name##_t *self, data_type *vec, size_t n, data_type *result) {                      \
        if (n != self->n) return -1;                                                                                        \
                                                                                                                            \
        uint##index_size##_t x, ind_start, ind_len;                                                                        \
        data_type val;                                                                                                      \
        data_type *data = self->data->a;                                                                                    \
                                                                                                                            \
        compressed_sparse_matrix_foreach(self, x, ind_start, ind_len, {                                               \
            data_type sum = result[x];                                                                                      \
            for (uint##index_size##_t y = ind_start; y < ind_start + ind_len; y++) {                                       \
                sum += data[y] * vec[y];                                                                                    \
            }                                                                                                               \
            result[x] = sum;                                                                                                \
        })                                                                                                                  \
        return 0;                                                                                                           \
    }                                                                                                                       \
                                                                                                                            \
    static inline int name##_##index_name##s_dot_vector(name##_t *self, uint##index_size##_t *xs, size_t m, data_type *vec, size_t n, data_type *result) { \
        if (n != self->n) return -1;                                                                                            \
                                                                                                                                \
        uint##index_size##_t *indptr = self->indptr->a;                                                                        \
        uint##index_size##_t *indices = self->indices->a;                                                                       \
        data_type *data = self->data->a;                                                                                        \
                                                                                                                                \
        for (uint##index_size##_t i = 0; i < m; i++) {                                                                         \
            uint##index_size##_t index = indices[i];                                                                            \
                                                                                                                                \
            data_type sum = result[i];                                                                                          \
            if (index >= self->m) return -1;                                                                                    \
                                                                                                                                \
            for (uint##index_size##_t j = indptr[index]; j < indptr[index+1]; j++) {                                            \
                sum += data[j] * vec[indices[j]];                                                                               \
            }                                                                                                                   \
                                                                                                                                \
            result[i] = sum;                                                                                                    \
                                                                                                                                \
        }                                                                                                                       \
        return 0;                                                                                                               \
    }                                                                                                                           \
                                                                                                                                \
    static inline name##_t *name##_read(FILE *f) {                                                                              \
        name##_t *sp = malloc(sizeof(name##_t));                                                                                \
        if (sp == NULL) return NULL;                                                                                            \
                                                                                                                                \
        sp->indptr = NULL;                                                                                                      \
        sp->indices = NULL;                                                                                                     \
        sp->data = NULL;                                                                                                        \
                                                                                                                                \
        if (!file_read_uint##index_size(f, &sp->m) ||                                                                           \
            !file_read_uint##index_size(f, &sp->n)) {                                                                           \
            goto exit_##name##_allocated;                                                                                       \
        }                                                                                                                       \
                                                                                                                                \
        uint64_t len_indptr;                                                                                                    \
                                                                                                                                \
        if (!file_read_uint64(f, &len_indptr)) {                                                                                \
            goto exit_##name##_allocated;                                                                                       \
        }                                                                                                                       \
                                                                                                                                \
        uint##index_size##_array *indptr = uint##index_size##_array_new_size((size_t)len_indptr);                               \
        if (indptr == NULL) {                                                                                                   \
            goto exit_##name##_allocated;                                                                                       \
        }                                                                                                                       \
                                                                                                                                \
        if (!file_read_uint##index_size##_array(f, indptr->a, len_indptr)) {                                                    \
            uint##index_size##_array_destroy(indptr);                                                                           \
            goto exit_##name##_allocated;                                                                                       \
        }                                                                                                                       \
                                                                                                                                \
        indptr->n = (size_t)len_indptr;                                                                                         \
        sp->indptr = indptr;                                                                                                    \
                                                                                                                                \
        uint64_t len_indices;                                                                                                   \
                                                                                                                                \
        if (!file_read_uint64(f, &len_indices)) {                                                                               \
            goto exit_##name##_allocated;                                                                                       \
        }                                                                                                                       \
                                                                                                                                \
        uint##index_size##_array *indices = uint##index_size##_array_new_size(len_indices);                                     \
        if (indices == NULL) {                                                                                                  \
            goto exit_##name##_allocated;                                                                                       \
        }                                                                                                                       \
                                                                                                                                \
        if (!file_read_uint##index_size##_array(f, indices->a, len_indices)) {                                                  \
            uint##index_size##_array_destroy(indices);                                                                          \
            goto exit_##name##_allocated;                                                                                       \
        }                                                                                                                       \
                                                                                                                                \
        indices->n = (size_t)len_indices;                                                                                       \
        sp->indices = indices;                                                                                                  \
                                                                                                                                \
        uint64_t len_data;                                                                                                      \
                                                                                                                                \
        if (!file_read_uint64(f, &len_data)) {                                                                                  \
            goto exit_##name##_allocated;                                                                                       \
        }                                                                                                                       \
                                                                                                                                \
        data_array_type *data = data_array_type##_new_size(len_data);                                                           \
        if (data == NULL) {                                                                                                     \
            goto exit_##name##_allocated;                                                                                       \
        }                                                                                                                       \
                                                                                                                                \
        if (!file_read_##data_type##_array(f, data->a, len_data)) {                                                             \
            data_array_type##_destroy(data);                                                                                    \
            goto exit_##name##_allocated;                                                                                       \
        }                                                                                                                       \
                                                                                                                                \
        data->n = (size_t)len_data;                                                                                             \
        sp->data = data;                                                                                                        \
                                                                                                                                \
        return sp;                                                                                                              \
                                                                                                                                \
    exit_##name##_allocated:                                                                                                    \
        name##_destroy(sp);                                                                                                     \
        return NULL;                                                                                                            \
    }                                                                                                                           \
                                                                                                                                \
    static inline bool name##_write(name##_t *self, FILE *f) {                                                                  \
        if (self == NULL || self->indptr == NULL || self->indices == NULL || self->data == NULL) {                              \
            return false;                                                                                                       \
        }                                                                                                                       \
                                                                                                                                \
        if (!file_write_uint##index_size(f, self->m) ||                                                                         \
            !file_write_uint##index_size(f, self->n)) {                                                                         \
            return false;                                                                                                       \
        }                                                                                                                       \
                                                                                                                                \
        uint64_t len_indptr = self->indptr->n;                                                                                  \
                                                                                                                                \
        if (!file_write_uint64(f, len_indptr)) {                                                                                \
            return false;                                                                                                       \
        }                                                                                                                       \
                                                                                                                                \
        for (int i = 0; i < len_indptr; i++) {                                                                                  \
            if (!file_write_uint##index_size(f, self->indptr->a[i])) {                                                          \
                return false;                                                                                                   \
            }                                                                                                                   \
        }                                                                                                                       \
                                                                                                                                \
        uint64_t len_indices = (uint64_t)self->indices->n;                                                                      \
                                                                                                                                \
        if (!file_write_uint64(f, len_indices)) {                                                                               \
            return false;                                                                                                       \
        }                                                                                                                       \
                                                                                                                                \
        for (int i = 0; i < len_indices; i++) {                                                                                 \
            if (!file_write_uint##index_size(f, self->indices->a[i])) {                                                         \
                return false;                                                                                                   \
            }                                                                                                                   \
        }                                                                                                                       \
                                                                                                                                \
        uint64_t len_data = (uint64_t)self->data->n;                                                                            \
                                                                                                                                \
        if (!file_write_uint64(f, len_data)) {                                                                                  \
            return false;                                                                                                       \
        }                                                                                                                       \
                                                                                                                                \
        for (int i = 0; i < len_data; i++) {                                                                                    \
            if (!file_write_##data_type(f, self->data->a[i])) {                                                                 \
                return false;                                                                                                   \
            }                                                                                                                   \
        }                                                                                                                       \
                                                                                                                                \
        return true;                                                                                                            \
    }                                                                                                                           \



#define COMPRESSED_SPARSE_MATRIX_INIT(name, data_type, data_array_type, index_name)                \
COMPRESSED_SPARSE_MATRIX_INIT_SIZE(name, data_type, data_array_type, 32, index_name)


#endif

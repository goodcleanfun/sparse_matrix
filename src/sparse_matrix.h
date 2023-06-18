/*
sparse_matrix.h
---------------

Dynamic compressed sparse row (CSR) sparse matrix.

A sparse matrix is a matrix where an overwhelming number of the
entries are 0. Thus we get significant space/time advantages
from only storing the non-zero values.

These types of matrices arise often when representing graphs
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

#ifndef SPARSE_MATRIX_H
#define SPARSE_MATRIX_H

#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>

#include "vector/vector.h"
#include "num_arrays/uint32_array.h"

#define SPARSE_MATRIX_INIT(name, type, array_type)                                                              \
    typedef struct {                                                                                            \
        uint32_t m;                                                                                             \
        uint32_t n;                                                                                             \
        uint32_array *indptr;                                                                                   \
        uint32_array *indices;                                                                                  \
        array_type *data;                                                                                       \
    } name##_t;                                                                                                 \
                                                                                                                \
    name##_t *name##_new(void);                                                                                 \
    name##_t *name##_new_shape(size_t m, size_t n);                                                             \
    void name##_destroy(name##_t *self);                                                                        \
                                                                                                                \
    void name##_clear(name##_t *self);                                                                          \
                                                                                                                \
    void name##_append(name##_t *self, uint32_t col, type val);                                                 \
    void name##_append_row(name##_t *self, uint32_t *col, type *val, size_t n);                                 \
                                                                                                                \
    void name##_finalize_row(name##_t *self);                                                                   \
                                                                                                                \
    void name##_sort_indices(name##_t *self);                                                                   \
                                                                                                                \
    int name##_dot_vector(name##_t *self, type *vec, size_t n, type *result);                                   \
    int name##_rows_dot_vector(name##_t *self, uint32_t *rows, size_t m, type *vec, size_t n, type *result);    \
                                                                                                                \
    // No need to allocate vector. Equivalent to doing a dot product with a vector of all ones                  \
    int name##_sum_cols(name##_t *self, type *result, size_t n);                                                \
    int name##_rows_sum_cols(name##_t *self, uint32_t *rows, size_t m, type *result, size_t n);                 \
                                                                                                                \
    int name##_sum_all_rows(name##_t *self, type *result, size_t n);                                            \
    int name##_sum_rows(name##_t *self, uint32_t *rows, size_t m, type *result, size_t n);                      \
                                                                                                                \
    bool name##_write(name##_t *self, FILE *f);                                                                 \
    name##_t *name##_read(FILE *f);                                                                             \
                                                                                                                \
    #define name##_foreach_row(sp, row_var, index_var, length_var, code) {                                      \
        uint32_t _row_start = 0, _row_end = 0;                                                                  \
        uint32_t *_indptr = sp->indptr->a;                                                                      \
        size_t _m = sp->m;                                                                                      \
                                                                                                                \
        for (uint32_t _i = 0; _i < _m; _i++) {                                                                  \
            (row_var) = _i;                                                                                     \
            _row_start = _indptr[_i];                                                                           \
            _row_end = _indptr[_i + 1];                                                                         \
            (index_var) = _row_start;                                                                           \
            (length_var) = _row_end - _row_start;                                                               \
            code;                                                                                               \
        }                                                                                                       \
    }                                                                                                           \
    #define name##_foreach(sp, row_var, col_var, data_var, code) {                                              \
        uint32_t *_indices = sp->indices->a;                                                                    \
        type *_data = sp->data->a;                                                                              \
        uint32_t _index, _length;                                                                               \
        name##_foreach_row(sp, row_var, _index, _length, {                                                      \
            for (uint32_t _j = _index; _j < _index + _length; _j++) {                                           \
                (col_var) = _indices[_j];                                                                       \
                (data_var) = _data[_j];                                                                         \
                code;                                                                                           \
            }                                                                                                   \
        })                                                                                                      \
    }                                                                                                           \
                                                                                                                \
    name##_t *name##_new_shape(size_t m, size_t n) {                                                            \
        name##_t *matrix = calloc(1, sizeof(name##_t));                                                         \
        if (matrix == NULL) return NULL;                                                                        \
        matrix->m = m;                                                                                          \
        matrix->n = n;                                                                                          \
        matrix->indptr = uint32_array_new_size(m + 1);                                                          \
        if (matrix->indptr == NULL) {                                                                           \
            goto exit_sparse_matrix_created;                                                                    \
        }                                                                                                       \
        uint32_array_push(matrix->indptr, 0);                                                                   \
                                                                                                                \
        matrix->indices = uint32_array_new();                                                                   \
        if (matrix->indices == NULL) {                                                                          \
            goto exit_sparse_matrix_created;                                                                    \
        }                                                                                                       \
                                                                                                                \
        matrix->data = array_type##_new();                                                                      \
        if (matrix->data == NULL) {                                                                             \
            goto exit_sparse_matrix_created;                                                                    \
        }                                                                                                       \
                                                                                                                \
        return matrix;                                                                                          \
                                                                                                                \
    exit_sparse_matrix_created:                                                                                 \
        name##_destroy(matrix);                                                                                 \
        return NULL;                                                                                            \
    }                                                                                                           \
                                                                                                                \
    name##_t *name##_new(void) {                                                                                \
        return name##_new_shape(0, 0);                                                                          \
    }                                                                                                           \
                                                                                                                \
                                                                                                                \
    void name##_destroy(name##_t *self) {                                                                       \
        if (self == NULL) return;                                                                               \
                                                                                                                \
        if (self->indptr != NULL) {                                                                             \
            uint32_array_destroy(self->indptr);                                                                 \
        }                                                                                                       \
                                                                                                                \
        if (self->indices != NULL) {                                                                            \
            uint32_array_destroy(self->indices);                                                                \
        }                                                                                                       \
                                                                                                                \
        if (self->data != NULL) {                                                                               \
            array_type##_destroy(self->data);                                                                   \
        }                                                                                                       \
                                                                                                                \
        free(self);                                                                                             \
    }                                                                                                           \
                                                                                                                \
    inline void name##_clear(name##_t *self) {                                                                  \
        uint32_array_clear(self->indptr);                                                                       \
        uint32_array_push(self->indptr, 0);                                                                     \
                                                                                                                \
        uint32_array_clear(self->indices);                                                                      \
        array_type##_clear(self->data);                                                                         \
    }                                                                                                           \
                                                                                                                \
    inline void name##_finalize_row(name##_t *self) {                                                           \
        uint32_array_push(self->indptr, (uint32_t)self->indices->n);                                            \
        if (self->indptr->n > self->m + 1) {                                                                    \
            self->m++;                                                                                          \
        }                                                                                                       \
    }                                                                                                           \
                                                                                                                \
    inline void name##_append(name##_t *self, uint32_t col, type val) {                                         \
        uint32_array_push(self->indices, col);                                                                  \
        array_type##_push(self->data, val);                                                                     \
        if (col >= self->n) self->n = col + 1;                                                                  \
    }                                                                                                           \
                                                                                                                \
    inline void name##_append_row(name##_t *self, uint32_t *col, type *val, size_t n) {                         \
        for (int i = 0; i < n; i++) {                                                                           \
            name##_append(self, col[i], val[i]);                                                                \
        }                                                                                                       \
        name##_finalize_row(self);                                                                              \
    }                                                                                                           \
                                                                                                                \
    typedef struct column_value {                                                                               \
        uint32_t col;                                                                                           \
        type val;                                                                                               \
    } column_value_t;                                                                                           \
                                                                                                                \
    VECTOR_INIT(column_value_array, column_value_t)                                                             \
                                                                                                                \
    #define ks_lt_column_value(a, b) ((a).col < (b).col)                                                        \
                                                                                                                \
    KSORT_INIT(column_value_array, column_value_t, ks_lt_column_value)                                          \
                                                                                                                \
    void name##_sort_indices(name##_t *self) {                                                                  \
        uint32_t row, row_start, row_len, i;                                                                    \
                                                                                                                \
        column_value_array *col_vals = column_value_array_new();                                                \
                                                                                                                \
        name##_foreach_row(self, row, row_start, row_len, {                                                     \
            for (i = row_start; i < row_start + row_len; i++) {                                                 \
                column_value_array_push(col_vals, (column_value_t){self->indices->a[i], self->data->a[i]});     \
            }                                                                                                   \
            ks_introsort(column_value_array, col_vals->n, col_vals->a);                                         \
                                                                                                                \
            for (i = 0; i < col_vals->n; i++) {                                                                 \
                column_value_t col_val = col_vals->a[i];                                                        \
                self->indices->a[row_start + i] = col_val.col;                                                  \
                self->data->a[row_start + i] = col_val.val;                                                     \
            }                                                                                                   \
        })                                                                                                      \
                                                                                                                \
    }                                                                                                           \
                                                                                                                \
                                                                                                                \
    inline int name##_dot_vector(name##_t *self, type *vec, size_t n, type *result) {                           \
        if (n != self->n) return -1;                                                                            \
                                                                                                                \
        uint32_t row, row_start, row_len;                                                                       \
        type val;                                                                                               \
        type *data = self->data->a;                                                                             \
                                                                                                                \
        name##_foreach_row(self, row, row_start, row_len, {                                                     \
            type sum = result[row];                                                                             \
            for (uint32_t col = row_start; col < row_start + row_len; col++) {                                  \
                sum += data[col] * vec[col];                                                                    \
            }                                                                                                   \
            result[row] = sum;                                                                                  \
        })                                                                                                      \
        return 0;                                                                                               \
    }                                                                                                           \
                                                                                                                \
    int name##_rows_dot_vector(name##_t *self, uint32_t *rows, size_t m, type *vec, size_t n, type *result) {   \
        if (n != self->n) return -1;                                                                            \
                                                                                                                \
        uint32_t *indptr = self->indptr->a;                                                                     \
        uint32_t *indices = self->indices->a;                                                                   \
        type *data = self->data->a;                                                                             \
                                                                                                                \
        for (int i = 0; i < m; i++) {                                                                           \
            uint32_t row = rows[i];                                                                             \
                                                                                                                \
            type sum = result[i];                                                                               \
            if (row >= self->m) return -1;                                                                      \
                                                                                                                \
            for (int j = indptr[row]; j < indptr[row+1]; j++) {                                                 \
                sum += data[j] * vec[indices[j]];                                                               \
            }                                                                                                   \
                                                                                                                \
            result[i] = sum;                                                                                    \
                                                                                                                \
        }                                                                                                       \
        return 0;                                                                                               \
    }                                                                                                           \
                                                                                                                \
    int name##_sum_cols(name##_t *self, type *result, size_t n) {                                               \
        if (n != self->m) return -1;                                                                            \
                                                                                                                \
        uint32_t row, row_start, row_len;                                                                       \
        type val;                                                                                               \
        type *data = self->data->a;                                                                             \
                                                                                                                \
        name##_foreach_row(self, row, row_start, row_len, {                                                     \
            type sum = result[row];                                                                             \
            for (uint32_t col = row_start; col < row_start + row_len; col++) {                                  \
                sum += data[col];                                                                               \
            }                                                                                                   \
            result[row] = sum;                                                                                  \
        })                                                                                                      \
        return 0;                                                                                               \
                                                                                                                \
    }                                                                                                           \
                                                                                                                \
    // No need to allocate actual vector for values, sum rather than a dot product                              \
    int name##_rows_sum_cols(name##_t *self, uint32_t *rows, size_t m, type *result, size_t n) {                \
        if (m != n) return -1;                                                                                  \
                                                                                                                \
        uint32_t *indptr = self->indptr->a;                                                                     \
        uint32_t *indices = self->indices->a;                                                                   \
        type *data = self->data->a;                                                                             \
                                                                                                                \
        for (int i = 0; i < m; i++) {                                                                           \
            uint32_t row = rows[i];                                                                             \
                                                                                                                \
            type sum = result[i];                                                                               \
            if (row >= self->m) return -1;                                                                      \
                                                                                                                \
            for (int j = indptr[row]; j < indptr[row+1]; j++) {                                                 \
                sum += data[j];                                                                                 \
            }                                                                                                   \
                                                                                                                \
            result[i] = sum;                                                                                    \
        }                                                                                                       \
        return 0;                                                                                               \
    }                                                                                                           \
                                                                                                                \
                                                                                                                \
    int name##_sum_all_rows(name##_t *self, type *result, size_t n) {                                           \
        if (n != self->n) return -1;                                                                            \
                                                                                                                \
        uint32_t row, row_start, row_len;                                                                       \
        type val;                                                                                               \
        type *data = self->data->a;                                                                             \
                                                                                                                \
        name##_foreach_row(self, row, row_start, row_len, {                                                     \
            for (uint32_t col = row_start; col < row_start + row_len; col++) {                                  \
                result[col] += data[col];                                                                       \
            }                                                                                                   \
        })                                                                                                      \
        return 0;                                                                                               \
                                                                                                                \
    }                                                                                                           \
                                                                                                                \
    int name##_sum_rows(name##_t *self, uint32_t *rows, size_t m, type *result, size_t n) {                     \
        if (n != self->n) return -1;                                                                            \
                                                                                                                \
        uint32_t *indptr = self->indptr->a;                                                                     \
        uint32_t *indices = self->indices->a;                                                                   \
        type *data = self->data->a;                                                                             \
                                                                                                                \
        for (int i = 0; i < m; i++) {                                                                           \
            uint32_t row = rows[i];                                                                             \
                                                                                                                \
            if (row >= self->m) return -1;                                                                      \
                                                                                                                \
            for (int j = indptr[row]; j < indptr[row+1]; j++) {                                                 \
                result[j] += data[j];                                                                           \
            }                                                                                                   \
        }                                                                                                       \
        return 0;                                                                                               \
    }                                                                                                           \
                                                                                                                \
    name##_t *name##_read(FILE *f) {                                                                            \
        name##_t *sp = malloc(sizeof(name##_t));                                                                \
        if (sp == NULL) return NULL;                                                                            \
                                                                                                                \
        sp->indptr = NULL;                                                                                      \
        sp->indices = NULL;                                                                                     \
        sp->data = NULL;                                                                                        \
                                                                                                                \
        if (!file_read_uint32(f, &sp->m) ||                                                                     \
            !file_read_uint32(f, &sp->n)) {                                                                     \
            goto exit_sparse_matrix_allocated;                                                                  \
        }                                                                                                       \
                                                                                                                \
        uint64_t len_indptr;                                                                                    \
                                                                                                                \
        if (!file_read_uint64(f, &len_indptr)) {                                                                \
            goto exit_sparse_matrix_allocated;                                                                  \
        }                                                                                                       \
                                                                                                                \
        uint32_array *indptr = uint32_array_new_size((size_t)len_indptr);                                       \
        if (indptr == NULL) {                                                                                   \
            goto exit_sparse_matrix_allocated;                                                                  \
        }                                                                                                       \
                                                                                                                \
        if (!file_read_uint32_array(f, indptr->a, len_indptr)) {                                                \
            uint32_array_destroy(indptr);                                                                       \
            goto exit_sparse_matrix_allocated;                                                                  \
        }                                                                                                       \
                                                                                                                \
        indptr->n = (size_t)len_indptr;                                                                         \
        sp->indptr = indptr;                                                                                    \
                                                                                                                \
        uint64_t len_indices;                                                                                   \
                                                                                                                \
        if (!file_read_uint64(f, &len_indices)) {                                                               \
            goto exit_sparse_matrix_allocated;                                                                  \
        }                                                                                                       \
                                                                                                                \
        uint32_array *indices = uint32_array_new_size(len_indices);                                             \
        if (indices == NULL) {                                                                                  \
            goto exit_sparse_matrix_allocated;                                                                  \
        }                                                                                                       \
                                                                                                                \
        if (!file_read_uint32_array(f, indices->a, len_indices)) {                                              \
            uint32_array_destroy(indices);                                                                      \
            goto exit_sparse_matrix_allocated;                                                                  \
        }                                                                                                       \
                                                                                                                \
        indices->n = (size_t)len_indices;                                                                       \
        sp->indices = indices;                                                                                  \
                                                                                                                \
        uint64_t len_data;                                                                                      \
                                                                                                                \
        if (!file_read_uint64(f, &len_data)) {                                                                  \
            goto exit_sparse_matrix_allocated;                                                                  \
        }                                                                                                       \
                                                                                                                \
        array_type *data = array_type##_new_size(len_data);                                                     \
        if (data == NULL) {                                                                                     \
            goto exit_sparse_matrix_allocated;                                                                  \
        }                                                                                                       \
                                                                                                                \
        if (!file_read_##type##_array(f, data->a, len_data)) {                                                  \
            array_type##_destroy(data);                                                                         \
            goto exit_sparse_matrix_allocated;                                                                  \
        }                                                                                                       \
                                                                                                                \
        data->n = (size_t)len_data;                                                                             \
        sp->data = data;                                                                                        \
                                                                                                                \
        return sp;                                                                                              \
                                                                                                                \
    exit_sparse_matrix_allocated:                                                                               \
        name##_destroy(sp);                                                                                     \
        return NULL;                                                                                            \
    }                                                                                                           \
                                                                                                                \
    bool name##_write(name##_t *self, FILE *f) {                                                                \
        if (self == NULL || self->indptr == NULL || self->indices == NULL || self->data == NULL) {              \
            return false;                                                                                       \
        }                                                                                                       \
                                                                                                                \
        if (!file_write_uint32(f, self->m) ||                                                                   \
            !file_write_uint32(f, self->n)) {                                                                   \
            return false;                                                                                       \
        }                                                                                                       \
                                                                                                                \
        uint64_t len_indptr = self->indptr->n;                                                                  \
                                                                                                                \
        if (!file_write_uint64(f, len_indptr)) {                                                                \
            return false;                                                                                       \
        }                                                                                                       \
                                                                                                                \
        for (int i = 0; i < len_indptr; i++) {                                                                  \
            if (!file_write_uint32(f, self->indptr->a[i])) {                                                    \
                return false;                                                                                   \
            }                                                                                                   \
        }                                                                                                       \
                                                                                                                \
        uint64_t len_indices = (uint64_t)self->indices->n;                                                      \
                                                                                                                \
        if (!file_write_uint64(f, len_indices)) {                                                               \
            return false;                                                                                       \
        }                                                                                                       \
                                                                                                                \
        for (int i = 0; i < len_indices; i++) {                                                                 \
            if (!file_write_uint32(f, self->indices->a[i])) {                                                   \
                return false;                                                                                   \
            }                                                                                                   \
        }                                                                                                       \
                                                                                                                \
        uint64_t len_data = (uint64_t)self->data->n;                                                            \
                                                                                                                \
        if (!file_write_uint64(f, len_data)) {                                                                  \
            return false;                                                                                       \
        }                                                                                                       \
                                                                                                                \
        for (int i = 0; i < len_data; i++) {                                                                    \
            if (!file_write_type(f, self->data->a[i])) {                                                        \
                return false;                                                                                   \
            }                                                                                                   \
        }                                                                                                       \
                                                                                                                \
        return true;                                                                                            \
    }                                                                                                           \

#endif

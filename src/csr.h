
#ifndef CSR_MATRIX_H
#define CSR_MATRIX_H

#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>

#include "compressed.h"


#define CSR_DOT_SPARSE(name, data_type, index_type, index_array_type, index_name, hash_type, heap_type)    \
    static inline bool name##_dot_sparse( \
        name *a, \
        name *b, \
        name **c \
    ) { \
        index_type a_rows = a->m; \
        index_type a_cols = a->n; \
        index_type b_rows = b->m; \
        index_type b_cols = b->n; \
\
        if (a_cols != b_rows) { \
            return false; \
        } \
 \
        index_type *a_indptr = a->indptr->a; \
        index_type *a_indices = a->indices->a; \
        data_type *a_values = a->data->a; \
        index_type *b_indptr = b->indptr->a; \
        index_type *b_indices = b->indices->a; \
        data_type *b_values = b->data->a; \
 \
        name *result = name##_new_shape(a_rows, b_cols); \
         \
        index_type nnz = 0; \
 \
        hash_type##_hash *sums = hash_type##_hash_new(); \
        index_array_type *sorted_indices = index_array_type##_new(); \
 \
        for (index_type i = 0; i < a_rows; i++) { \
            index_type head = 0; \
            index_type length = 0; \
 \
            index_type a_ind_start = a_indptr[i]; \
            index_type a_ind_end = a_indptr[i+1]; \
 \
            for (index_type j = a_ind_start; j < a_ind_end; j++) { \
                index_type a_index = a_indices[j]; \
                data_type a_value = a_values[j]; \
 \
                index_type b_ind_start = b_indptr[a_index]; \
                index_type b_ind_end = b_indptr[a_index+1]; \
 \
                for (index_type k = b_ind_start; k < b_ind_end; k++) { \
                    index_type b_index = b_indices[k]; \
                    data_type b_value = b_values[k]; \
                    bool index_seen = false; \
                    hash_type##_hash_incr_by_exists(sums, b_index, a_value * b_value, &index_seen); \
 \
                    if (!index_seen) { \
                        heap_type##_push(sorted_indices, b_index); \
                    } \
                } \
            } \
 \
            index_type nonzero_index; \
            while (heap_type##_pop(sorted_indices, &nonzero_index)) { \
                data_type sum_value; \
                if (!hash_type##_hash_get(sums, nonzero_index, &sum_value)) { \
                    break; \
                } \
                hash_type##_hash_del(sums, nonzero_index); \
                name##_append(result, nonzero_index, sum_value); \
                nnz++; \
            } \
            name##_finalize_##index_name(result); \
        } \
        hash_type##_hash_destroy(sums); \
        index_array_type##_destroy(sorted_indices); \
        *c = result; \
 \
        return true; \
    }


#endif
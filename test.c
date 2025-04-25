#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>
#include <float.h>

#include "greatest/greatest.h"
#include "float_csr.h"

TEST test_float_csr_matrix(void) {
    float_csr_matrix *sp = float_csr_matrix_new();
    uint32_t indices_1[] = {0, 1, 2, 3};
    float values_1[] = {1.0f, 1.0f, 1.0f, 1.0f};
    float_csr_matrix_append_row(sp, (uint32_t *)indices_1, (float *)values_1, sizeof(indices_1) / sizeof(uint32_t));

    ASSERT_EQ(sp->m, 1);
    ASSERT_EQ(sp->n, 4);
    ASSERT_EQ(sp->indices->n, 4);
    for (uint32_t i = 0; i < sp->n; i++) {
        ASSERT_EQ(sp->indices->a[i], indices_1[i]);
        ASSERT_EQ(sp->data->a[i], values_1[i]);
    }

    uint32_t indices_2[] = {1, 4, 5};
    float values_2[] = {1.0f, 1.0f, 1.0f};

    float_csr_matrix_append_row(sp, indices_2, values_2, sizeof(indices_2) / sizeof(uint32_t));

    ASSERT_EQ(sp->m, 2);
    ASSERT_EQ(float_csr_matrix_rows(sp), 2);
    ASSERT_EQ(sp->n, 6);
    ASSERT_EQ(float_csr_matrix_cols(sp), 6);
    ASSERT_EQ(sp->max_nnz, 4);
    ASSERT_EQ(sp->indices->n, 7);
    ASSERT_EQ(sp->data->n, 7);
    for (uint32_t i = sizeof(indices_1) / sizeof(uint32_t); i < sp->n; i++) {
        ASSERT_EQ(sp->indices->a[i], indices_2[i - 4]);
        ASSERT_EQ(sp->data->a[i], values_2[i - 4]);
    }

    float_csr_matrix_append(sp, 4, 1.0f);
    float_csr_matrix_append(sp, 5, 1.0f);
    float_csr_matrix_append(sp, 6, 1.0f);
    float_csr_matrix_append(sp, 7, 1.0f);

    float_csr_matrix_finalize_row(sp);

    ASSERT_EQ(sp->m, 3);
    ASSERT_EQ(sp->n, 8);

    float_array *vec = float_array_new_size(sp->n);
    for (size_t i = 0; i < sp->n; i++) {
        float_array_push(vec, 1.0f);
    }

    float_array *res = float_array_new_size(sp->m);
    float_csr_matrix_dot_vector(sp, vec->a, vec->n, res->a, sp->m);
    uint32_t row = 0, row_start = 0, row_len = 0;
    compressed_sparse_matrix_foreach(sp, row, row_start, row_len, {
        ASSERT_EQ(res->a[row], row_len * 1.0f);
    })

    size_t num_rows = 2;
    uint32_t rows[] = {0, 2};
    float sparse_result[] = {0.0f, 0.0f};
    ASSERT(float_csr_matrix_rows_dot_vector(sp, rows, num_rows, vec->a, vec->n, sparse_result, num_rows));
    for (size_t i = 0; i < num_rows; i++) {
        uint32_t row = rows[i];
        uint32_t row_len = float_csr_matrix_len_row(sp, row);
        ASSERT_EQ(sparse_result[i], row_len * 1.0f);
    }

    float_csr_matrix *sp2 = float_csr_matrix_new();

    float_csr_matrix_append_row(sp2, (uint32_t[]){0}, (float[]){3.0f}, 1);
    float_csr_matrix_append_row(sp2, (uint32_t[]){1}, (float[]){5.0f}, 1);
    float_csr_matrix_append_row(sp2, (uint32_t[]){1}, (float[]){2.0f}, 1);
    float_csr_matrix_append_row(sp2, (uint32_t[]){0, 1},(float[]) {3.0f, 5.0f}, 2);
    float_csr_matrix_append_row(sp2, (uint32_t[]){1}, (float[]){8.0f}, 1);
    float_csr_matrix_append_row(sp2, (uint32_t[]){0}, (float[]){4.0f}, 1);
    float_csr_matrix_append_row(sp2, (uint32_t[]){0, 1}, (float[]){1.0f, 2.0f}, 2);
    float_csr_matrix_append_row(sp2, (uint32_t[]){1}, (float[]){7.0f}, 1);

    float_csr_matrix *sp_prod_sp2 = float_csr_matrix_new();

    ASSERT(float_csr_matrix_dot_sparse(sp, sp2, sp_prod_sp2, false));

    ASSERT_EQ(sp_prod_sp2->m, 3);
    ASSERT_EQ(sp_prod_sp2->n, 2);

    row = 0;
    uint32_t col = 0;
    float data = 0.0f;

    uint32_t expected_indptr[] = {0, 2, 4, 6};
    uint32_t expected_indices[] = {0, 1, 0, 1, 0, 1};
    float expected_data[] = {6.0f, 12.0f, 4.0f, 13.0f, 5.0f, 17.0f};

    size_t indptr_size = sizeof(expected_indptr) / sizeof(uint32_t);
    size_t indices_size = sizeof(expected_indices) / sizeof(uint32_t);
    size_t data_size = sizeof(expected_data) / sizeof(float);

    for (size_t i = 0; i < indptr_size; i++) {
        ASSERT_EQ(sp_prod_sp2->indptr->a[i], expected_indptr[i]);
    }
    for (size_t i = 0; i < indices_size; i++) {
        ASSERT_EQ(sp_prod_sp2->indices->a[i], expected_indices[i]);
    }
    for (size_t i = 0; i < data_size; i++) {
        ASSERT_IN_RANGE(sp_prod_sp2->data->a[i] - expected_data[i], 0.0f, FLT_EPSILON);
    }

    float_array_destroy(vec);
    float_array_destroy(res);
    float_csr_matrix_destroy(sp);
    float_csr_matrix_destroy(sp2);
    float_csr_matrix_destroy(sp_prod_sp2);

    PASS();
}

/* Add definitions that need to be in the test runner's main file. */
GREATEST_MAIN_DEFS();

int main(int argc, char **argv) {
    GREATEST_MAIN_BEGIN();      /* command-line options, initialization. */

    RUN_TEST(test_float_csr_matrix);

    GREATEST_MAIN_END();        /* display results */
}

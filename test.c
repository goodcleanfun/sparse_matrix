#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>

#include "greatest/greatest.h"
#include "float_csr.h"

TEST test_float_csr_matrix(void) {
    float_csr_matrix *sp = float_csr_matrix_new();
    uint32_t indices_1[] = {0, 1, 2, 3};
    float values_1[] = {1.0, 1.0, 1.0, 1.0};
    float_csr_matrix_append_row(sp, (uint32_t *)indices_1, (float *)values_1, sizeof(indices_1) / sizeof(uint32_t));

    ASSERT_EQ(sp->m, 1);
    ASSERT_EQ(sp->n, 4);
    ASSERT_EQ(sp->indices->n, 4);
    for (uint32_t i = 0; i < sp->n; i++) {
        ASSERT_EQ(sp->indices->a[i], indices_1[i]);
        ASSERT_EQ(sp->data->a[i], values_1[i]);
    }

    uint32_t indices_2[] = {1, 4, 5};
    float values_2[] = {1.0, 2.0, 1.0};

    float_csr_matrix_append_row(sp, indices_2, values_2, sizeof(indices_2) / sizeof(uint32_t));

    ASSERT_EQ(sp->m, 2);
    ASSERT_EQ(sp->n, 6);
    ASSERT_EQ(sp->indices->n, 7);
    ASSERT_EQ(sp->data->n, 7);
    for (uint32_t i = sizeof(indices_1) / sizeof(uint32_t); i < sp->n; i++) {
        ASSERT_EQ(sp->indices->a[i], indices_2[i - 4]);
        ASSERT_EQ(sp->data->a[i], values_2[i - 4]);
    }

    float_csr_matrix_append(sp, 4, 1.0);
    float_csr_matrix_append(sp, 5, 1.0);
    float_csr_matrix_append(sp, 6, 1.0);
    float_csr_matrix_append(sp, 7, 1.0);

    float_csr_matrix_finalize_row(sp);

    ASSERT_EQ(sp->m, 3);

    float_csr_matrix_destroy(sp);

    PASS();
}


/* Add definitions that need to be in the test runner's main file. */
GREATEST_MAIN_DEFS();

int main(int argc, char **argv) {
    GREATEST_MAIN_BEGIN();      /* command-line options, initialization. */

    RUN_TEST(test_float_csr_matrix);

    GREATEST_MAIN_END();        /* display results */
}

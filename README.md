# sparse
Generic dynamic sparse matrices in CSR (compressed sparse row) format.

Implements an efficient sparse matrix multiply which uses a hash table and a heap instead of an N-vector.

## General sparse matrix multiplication
Typically, a sparse dot product between two matrices is accomplished using some version of Gustavson's row-wise SpGEMM:

```
for a_i in matrix A do
    for a_ik in row a_i do
        for b_kj in row b_k do
            value ← a_ik * b_kj
            if c_ij < c_i
                then insert c_ij ← value
            else
                c_ij ← c_ij + value
```

One widely-used implementation is scipy.sparse (sparsetools), which allocates two intermediate vectors of size |N| (the number of columns in the output matrix C i.e. the number of rows in B if B is a CSC matrix). One array, called sums, is used to accumulate the summed values for all columns in the current row, and an array called next which keeps pointers to the next non-zero entry in the sums array. A counter is incremented each time a new column is encountered, the head is set to that new column, and at the end of the row, the algorithm starts at head and follows the pointers to the non-zero sums. At the end of each row, the values then have to be overwritten before the next step.

Instead of storing these two |N| vectors (where N, the number of columns of the matrix, could potentially be quite large, and depending on the structure each row may only need a small subset of the columns), this implementation makes the following observations:

1. the column-wise sums can be accumulated in a hashtable which is cleared (zeroed out without deallocating, so that the memory is reused) after each row and will not grow larger than the max number of non-zero columns in a given row, which in a high-sparsity scenario will tend to be closer to 1 than to |N|. If the indices do not need to be guaranteed sorted, we can just iterate through the hashtable and output the indices and values.
2. the output's indices can be sorted by using a heap/priority queue. Maintaining the heap invariant costs O(lg n) per insert where n is the number of non-zeros b_kj, which will again tend to be closer to 1 than |N| in the highly sparse case. Heap inserts/pops are only needed once per non-zero in b_k (insert when a new column is encountered, then pop from the heap at the end to retrieve the ordering).

This version should be lower memory and probably a bit faster since clearing the hashtable requires a number of operations proportional to the number of buckets used rather than the number of occupied slots (in the array case, clearing takes time proportional to the non-zeros for the current row).

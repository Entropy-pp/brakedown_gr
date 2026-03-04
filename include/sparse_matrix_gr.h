#ifndef SPARSE_MATRIX_GR_H
#define SPARSE_MATRIX_GR_H

#include <NTL/ZZ_pE.h>
#include <vector>
#include <brakedown_params.h>

using namespace NTL;

// One non-zero entry: (column_index, coefficient)
struct SparseEntry {
    long col;
    ZZ_pE val;
};

// One row of the sparse matrix
struct SparseRow {
    std::vector<SparseEntry> entries;
};

// Sparse matrix over GR(2^k, d)
struct SparseMatrixGR {
    long n_rows; // number of rows (= input dimension)
    long n_cols; // number of columns (= output dimension)
    long sparsity; // entries per row
    std::vector<SparseRow> rows;
};

// Create a random sparse matrix with given dimensions
SparseMatrixGR createRandomSparseMatrix(const LevelDim& dim);

// Compute: target[col] += sum_i input[i] * M[i][j]  (sparse matrix-vector product)
// input has length M.n_rows, target has length M.n_cols
void sparseMatVecMul(const SparseMatrixGR& M, const ZZ_pE* input, ZZ_pE* target);

// Same but add into target
void sparseMatVecMulAdd(const SparseMatrixGR& M, const ZZ_pE* input, ZZ_pE* target);

#endif
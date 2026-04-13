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

// ============================================================
// 小环专用：稀疏矩阵系数来自 base ring
// ============================================================

// 小环稀疏矩阵条目：系数是 ZZ_pX（base ring 多项式表示）
struct SparseEntrySmallRing {
    long col;
    ZZ_pX val;  // base ring 元素的多项式表示
};

struct SparseRowSmallRing {
    std::vector<SparseEntrySmallRing> entries;
};

struct SparseMatrixSmallRing {
    long n_rows;
    long n_cols;
    long sparsity;
    std::vector<SparseRowSmallRing> rows;
};

// 创建小环稀疏矩阵（系数来自 base ring）
// 调用前应在 base ring context
SparseMatrixSmallRing createRandomSparseMatrixSmallRing(const LevelDim& dim, long base_degree);

// 小环稀疏矩阵乘法：使用分量独立乘法
// input 是 packed 向量（每个元素是 ext ring 元素，但有分量结构）
// 结果也保持分量结构
void sparseMatVecMulSmallRing(const SparseMatrixSmallRing& M,
                               const ZZ_pE* input, ZZ_pE* target,
                               long pf, long base_r, const ZZ_pX& base_mod);

void sparseMatVecMulAddSmallRing(const SparseMatrixSmallRing& M,
                                  const ZZ_pE* input, ZZ_pE* target,
                                  long pf, long base_r, const ZZ_pX& base_mod);

#endif
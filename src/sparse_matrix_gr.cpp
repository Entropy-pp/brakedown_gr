#include <sparse_matrix_gr.h>
#include <brakedown_code_gr.h>
#include <gr.h>
#include <set>
#include <NTL/ZZ_pE.h>

using namespace NTL;

SparseMatrixGR createRandomSparseMatrix(const LevelDim& dim) {
    SparseMatrixGR M;
    M.n_rows = dim.n;
    M.n_cols = dim.m;
    M.sparsity = dim.d;
    M.rows.resize(dim.n);

    for (long i = 0; i < dim.n; i++) {
        std::set<long> chosen;
        // Pick 'd' distinct random column indices in [0, m)
        while ((long)chosen.size() < dim.d) {
            long col = RandomBnd(dim.m);
            chosen.insert(col);
        }
        M.rows[i].entries.reserve(dim.d);
        for (long col : chosen) {
            // Random invertible coefficient in GR*
            ZZ_pE val = random_ZZ_pE();
            M.rows[i].entries.push_back({col, val});
        }
    }
    return M;
}

void sparseMatVecMul(const SparseMatrixGR& M, const ZZ_pE* input, ZZ_pE* target) {
    // Zero out target
    for (long j = 0; j < M.n_cols; j++) clear(target[j]);
    sparseMatVecMulAdd(M, input, target);
}

void sparseMatVecMulAdd(const SparseMatrixGR& M, const ZZ_pE* input, ZZ_pE* target) {
    for (long i = 0; i < M.n_rows; i++) {
        for (const auto& entry : M.rows[i].entries) {
            // target[col] += input[i] * coeff
            target[entry.col] += input[i] * entry.val;
        }
    }
}

// ============================================================
// 小环专用函数
// ============================================================

SparseMatrixSmallRing createRandomSparseMatrixSmallRing(const LevelDim& dim, long base_degree) {
    SparseMatrixSmallRing M;
    M.n_rows = dim.n;
    M.n_cols = dim.m;
    M.sparsity = dim.d;
    M.rows.resize(dim.n);

    for (long i = 0; i < dim.n; i++) {
        std::set<long> chosen;
        while ((long)chosen.size() < dim.d) {
            long col = RandomBnd(dim.m);
            chosen.insert(col);
        }
        M.rows[i].entries.reserve(dim.d);
        for (long col : chosen) {
            // ★ 关键：系数来自 base ring（随机非零可逆元素）
            // 生成随机的 base ring 元素（在当前 ZZ_pE context 下）
            ZZ_pE val = randomNonZeroInExceptionalSet();
            ZZ_pX val_poly = rep(val);
            M.rows[i].entries.push_back({col, val_poly});
        }
    }
    return M;
}

// 小环稀疏矩阵乘法：使用分量独立乘法
void sparseMatVecMulSmallRing(const SparseMatrixSmallRing& M,
                               const ZZ_pE* input, ZZ_pE* target,
                               long pf, long base_r, const ZZ_pX& base_mod)
{
    // Zero out target
    for (long j = 0; j < M.n_cols; j++) clear(target[j]);
    sparseMatVecMulAddSmallRing(M, input, target, pf, base_r, base_mod);
}

void sparseMatVecMulAddSmallRing(const SparseMatrixSmallRing& M,
                                  const ZZ_pE* input, ZZ_pE* target,
                                  long pf, long base_r, const ZZ_pX& base_mod)
{
    for (long i = 0; i < M.n_rows; i++) {
        for (const auto& entry : M.rows[i].entries) {
            // ★ 使用分量独立乘法：entry.val 是 base ring 多项式
            ZZ_pE product = component_wise_scalar_mul(
                entry.val, input[i], pf, base_r, base_mod);
            target[entry.col] += product;
        }
    }
}
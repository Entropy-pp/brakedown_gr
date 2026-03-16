#include <sparse_matrix_gr.h>
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
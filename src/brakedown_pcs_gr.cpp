#include <brakedown_pcs_gr.h>
#include <gr.h>
#include <cassert>
#include <cstdlib>

using namespace NTL;

// ============================================================
// Inner product over GR
// ============================================================
static ZZ_pE inner_product_gr(const ZZ_pE* a, const ZZ_pE* b, long len) {
    ZZ_pE result;
    clear(result);
    for (long i = 0; i < len; i++) {
        result += a[i] * b[i];
    }
    return result;
}

static ZZ_pE inner_product_gr(const std::vector<ZZ_pE>& a,
                               const std::vector<ZZ_pE>& b) {
    assert(a.size() == b.size());
    return inner_product_gr(a.data(), b.data(), (long)a.size());
}

// ============================================================
// 辅助: 计算 Small Ring 下实际的 proximity test 次数
// 论文 Lemma 9: 重复 ceil(λ/r) 次
// ============================================================
static long effective_num_prox(const BrakedownCodeGR& code, long num_rows) {
    if (num_rows <= 1) return 0;
    long num_prox = code.num_prox_test;
    if (code.is_small_ring) {
        long repetitions = (code.lambda + code.base_degree - 1) / code.base_degree;
        num_prox = std::max(num_prox, repetitions);
    }
    return num_prox;
}

// ============================================================
// 辅助: Small Ring 下对单行做 pack → encode
// combined_row 是 row_len 个 GR(p^s,r) 元素
// 输出: codeword_len 个 GR(p^s,kr) 元素
// 调用前 NTL context 在 base ring; 调用后切到 ext ring
// ============================================================
static std::vector<ZZ_pE> encode_row_small_ring(
    const BrakedownCodeGR& code,
    const std::vector<ZZ_pE>& row)
{
    std::vector<ZZ_pE> result(code.codeword_len);
    brakedown_encode_small_ring(code, row.data(), (long)row.size(), result.data());
    return result;
}

// ============================================================
// 辅助: 将 NTL context 切换到 base ring GR(p^s, base_degree)
// ============================================================
static void switch_to_base_ring(long base_degree) {
    ZZ_pX base_mod = primitiveIrredPoly(base_degree);
    ZZ_pE::init(base_mod);
}

// ============================================================
// 辅助: 将 NTL context 切换到 ext ring GR(p^s, ext_degree)
// ============================================================
static void switch_to_ext_ring(long ext_degree) {
    ZZ_pX ext_mod = primitiveIrredPoly(ext_degree);
    ZZ_pE::init(ext_mod);
}

// ============================================================
// COMMIT
//
// Large Ring: 和之前一样，直接编码
// Small Ring:
//   - 系数矩阵 S 的每行有 row_len 个 GR(p^s,r) 元素
//   - 编码: pack 每行 → GR(p^s,kr) 上 packed_row_len 个元素 → Enc
//   - encoded_rows 存储编码后的 GR(p^s,kr) 元素
//   - Merkle tree 在 GR(p^s,kr) 元素上构建
// ============================================================
BrakedownCommitment brakedown_commit(
    const BrakedownCodeGR& code,
    const ZZ_pE* poly_coeffs,
    long n,
    long num_rows)
{
    BrakedownCommitment comm;
    comm.num_rows = num_rows;
    comm.codeword_len = code.codeword_len;
    long row_len = code.row_len;

    assert(n <= num_rows * row_len);

    if (!code.is_small_ring) {
        // ====== Large Ring: 原有逻辑 ======
        comm.encoded_rows.resize((long)num_rows * code.codeword_len);
        for (long i = 0; i < num_rows; i++) {
            for (long j = 0; j < row_len; j++) {
                long idx = i * row_len + j;
                if (idx < n)
                    comm.encoded_rows[i * code.codeword_len + j] = poly_coeffs[idx];
                else
                    clear(comm.encoded_rows[i * code.codeword_len + j]);
            }
            for (long j = row_len; j < code.codeword_len; j++) {
                clear(comm.encoded_rows[i * code.codeword_len + j]);
            }
            brakedown_encode(code, &comm.encoded_rows[i * code.codeword_len]);
        }
    } else {
        // ====== Small Ring: pack → encode in ext ring ======
        // 先在 base ring context 中提取每行系数，
        // 然后每行调 brakedown_encode_small_ring
        comm.encoded_rows.resize((long)num_rows * code.codeword_len);

        for (long i = 0; i < num_rows; i++) {
            // 确保在 base ring context 下
            switch_to_base_ring(code.base_degree);

            // 构建本行的 base ring 元素
            std::vector<ZZ_pE> row_base(row_len);
            for (long j = 0; j < row_len; j++) {
                long idx = i * row_len + j;
                if (idx < n)
                    row_base[j] = poly_coeffs[idx];
                else
                    clear(row_base[j]);
            }

            // 提取 ZZ_pX 系数 (在 base ring context 下)
            std::vector<ZZ_pX> row_polys(row_len);
            for (long j = 0; j < row_len; j++) {
                row_polys[j] = rep(row_base[j]);
            }

            // 切换到 ext ring 并执行 pack + encode
            switch_to_ext_ring(code.ext_degree);

            long pf = code.packing_factor;
            long base_r = code.base_degree;
            long packed_len = code.packed_row_len;
            long cw_len = code.codeword_len;

            // Pack
            ZZ_pE* cw_ptr = &comm.encoded_rows[i * cw_len];
            for (long j = 0; j < packed_len; j++) {
                std::vector<ZZ_pX> group(pf);
                for (long t = 0; t < pf; t++) {
                    long src = j * pf + t;
                    if (src < row_len) group[t] = row_polys[src];
                }
                cw_ptr[j] = pack_elements_pub(group, pf, base_r);
            }
            for (long j = packed_len; j < cw_len; j++) clear(cw_ptr[j]);

            // Encode in ext ring
            brakedown_encode(code, cw_ptr);
        }

        // 确保最终 context 在 ext ring (Merkle tree 需要)
        switch_to_ext_ring(code.ext_degree);
    }

    // Build Merkle tree (列哈希)
    std::vector<std::vector<unsigned char>> leaf_hashes(code.codeword_len);
    for (long col = 0; col < code.codeword_len; col++) {
        std::vector<ZZ_pE> column(num_rows);
        for (long row = 0; row < num_rows; row++) {
            column[row] = comm.encoded_rows[row * code.codeword_len + col];
        }
        leaf_hashes[col] = hash_column(column.data(), num_rows);
    }

    comm.tree = build_merkle_tree(leaf_hashes);
    comm.root = comm.tree.root;

    return comm;
}

// ============================================================
// PROVE
//
// Large Ring: 原有逻辑
// Small Ring:
//   - challenge 从 GR(p^s,r) 的 exceptional set 选取 (Lemma 4)
//   - 线性组合在 GR(p^s,r) 上做 → combined_row 仍是 row_len 个 GR(p^s,r) 元素
//   - 重复 ceil(λ/r) 次 (Lemma 9)
//   - column_items 是 ext ring 元素
// ============================================================
BrakedownEvalProof brakedown_prove(
    const BrakedownCodeGR& code,
    const BrakedownCommitment& comm,
    const ZZ_pE* poly_coeffs,
    long n,
    const std::vector<ZZ_pE>& q1,
    const std::vector<ZZ_pE>& q2)
{
    BrakedownEvalProof proof;
    long row_len  = code.row_len;
    long num_rows = comm.num_rows;

    // combine_rows 在 base ring 上做线性组合
    // 对 small ring: 确保在 base ring context
    if (code.is_small_ring) {
        switch_to_base_ring(code.base_degree);
    }

    auto combine_rows = [&](const std::vector<ZZ_pE>& coeffs) {
        std::vector<ZZ_pE> result(row_len);
        for (long j = 0; j < row_len; j++) clear(result[j]);
        for (long i = 0; i < num_rows; i++) {
            for (long j = 0; j < row_len; j++) {
                long idx = i * row_len + j;
                ZZ_pE val;
                if (idx < n) val = poly_coeffs[idx]; else clear(val);
                result[j] += coeffs[i] * val;
            }
        }
        return result;
    };

    long num_prox = effective_num_prox(code, num_rows);

    // Proximity tests
    proof.prox_coeffs.resize(num_prox);
    for (long t = 0; t < num_prox; t++) {
        if (code.is_small_ring) switch_to_base_ring(code.base_degree);

        proof.prox_coeffs[t].resize(num_rows);
        for (long i = 0; i < num_rows; i++) {
            if (code.is_small_ring) {
                // Small ring: challenge 来自 GR(p^s,r) 的 exceptional set
                proof.prox_coeffs[t][i] = randomNonZeroInExceptionalSet();
            } else {
                proof.prox_coeffs[t][i] = randomInvertible();
            }
        }
        proof.combined_rows.push_back(combine_rows(proof.prox_coeffs[t]));
    }

    // Evaluation consistency: combine with q1
    if (code.is_small_ring) switch_to_base_ring(code.base_degree);
    auto q1_combined = combine_rows(q1);
    proof.eval_value = inner_product_gr(q1_combined, q2);
    proof.combined_rows.push_back(q1_combined);

    // Column openings (在 ext ring context 下读 encoded_rows)
    if (code.is_small_ring) switch_to_ext_ring(code.ext_degree);

    long num_open = code.num_col_open;
    proof.column_indices.resize(num_open);
    proof.column_items.resize(num_open);
    proof.merkle_paths.resize(num_open);

    for (long t = 0; t < num_open; t++) {
        long col_idx = RandomBnd(code.codeword_len);
        proof.column_indices[t] = col_idx;

        std::vector<ZZ_pE> column(num_rows);
        for (long row = 0; row < num_rows; row++) {
            column[row] = comm.encoded_rows[row * comm.codeword_len + col_idx];
        }
        proof.column_items[t] = column;
        proof.merkle_paths[t] = get_merkle_path(comm.tree, col_idx);
    }

    return proof;
}

// ============================================================
// VERIFY
//
// Large Ring: 原有逻辑
// Small Ring:
//   - num_prox 同步使用 effective_num_prox (与 Prove 一致)
//   - 重新编码 combined rows: pack → encode in ext ring
//   - column consistency: 需要在 ext ring 下比较
//     · 对 column_items (ext ring 元素)，用 prox_coeffs (base ring)
//       乘以每行的打开列 → 在 ext ring 下做内积
//     · 论文保证: c ∈ GR(p^s,r), b ∈ GR(p^s,kr)
//       c*b 等价于对 b 的 k 个分量分别乘 c (Lemma 4)
// ============================================================
bool brakedown_verify(
    const BrakedownCodeGR& code,
    const BrakedownCommitment& comm,
    const BrakedownEvalProof& proof,
    const std::vector<ZZ_pE>& q1,
    const std::vector<ZZ_pE>& q2)
{
    long row_len      = code.row_len;
    long codeword_len = code.codeword_len;
    long num_rows     = comm.num_rows;

    // ====== 关键修复: 与 Prove 使用相同的 num_prox ======
    long num_prox = effective_num_prox(code, num_rows);

    // Check proof structure
    if ((long)proof.combined_rows.size() != num_prox + 1) {
        std::cout << "VERIFY FAIL: wrong number of combined_rows ("
                  << proof.combined_rows.size() << " vs expected "
                  << num_prox + 1 << ")\n";
        return false;
    }

    // 1. Evaluation consistency (在 base ring 下)
    if (code.is_small_ring) switch_to_base_ring(code.base_degree);
    const auto& q1_row = proof.combined_rows.back();
    ZZ_pE computed_eval = inner_product_gr(q1_row, q2);
    if (computed_eval != proof.eval_value) {
        std::cout << "VERIFY FAIL: evaluation consistency\n";
        return false;
    }

    // 2. Build coefficient lists
    std::vector<const std::vector<ZZ_pE>*> all_coeffs;
    for (long t = 0; t < num_prox; t++) {
        all_coeffs.push_back(&proof.prox_coeffs[t]);
    }
    all_coeffs.push_back(&q1);

    // 3. Re-encode all combined rows
    std::vector<std::vector<ZZ_pE>> encoded_combined(proof.combined_rows.size());
    for (size_t t = 0; t < proof.combined_rows.size(); t++) {
        if (!code.is_small_ring) {
            // Large Ring: 直接编码
            encoded_combined[t].resize(codeword_len);
            for (long j = 0; j < row_len; j++) {
                encoded_combined[t][j] = proof.combined_rows[t][j];
            }
            for (long j = row_len; j < codeword_len; j++) {
                clear(encoded_combined[t][j]);
            }
            brakedown_encode(code, encoded_combined[t].data());
        } else {
            // Small Ring: pack → encode in ext ring
            // combined_rows[t] 有 row_len 个 base ring 元素
            switch_to_base_ring(code.base_degree);

            // 提取 ZZ_pX 系数
            std::vector<ZZ_pX> row_polys(row_len);
            for (long j = 0; j < row_len; j++) {
                row_polys[j] = rep(proof.combined_rows[t][j]);
            }

            // 切换到 ext ring
            switch_to_ext_ring(code.ext_degree);

            long pf = code.packing_factor;
            long base_r = code.base_degree;
            long packed_len = code.packed_row_len;

            encoded_combined[t].resize(codeword_len);
            for (long j = 0; j < packed_len; j++) {
                std::vector<ZZ_pX> group(pf);
                for (long tt = 0; tt < pf; tt++) {
                    long src = j * pf + tt;
                    if (src < row_len) group[tt] = row_polys[src];
                }
                encoded_combined[t][j] = pack_elements_pub(group, pf, base_r);
            }
            for (long j = packed_len; j < codeword_len; j++) {
                clear(encoded_combined[t][j]);
            }
            brakedown_encode(code, encoded_combined[t].data());
        }
    }

    // 确保在 ext ring context 下做列一致性检查
    if (code.is_small_ring) switch_to_ext_ring(code.ext_degree);

    // 4. Verify Merkle openings + column consistency
    for (long t = 0; t < code.num_col_open; t++) {
        long col = proof.column_indices[t];
        const auto& items = proof.column_items[t];
        const auto& path  = proof.merkle_paths[t];

        if ((long)items.size() != num_rows) {
            std::cout << "VERIFY FAIL: wrong column_items length\n";
            return false;
        }

        // Merkle path
        auto leaf_hash = hash_column(items.data(), num_rows);
        if (!verify_merkle_path(leaf_hash, col, path, comm.root)) {
            std::cout << "VERIFY FAIL: Merkle path at column " << col << "\n";
            return false;
        }

        // Column consistency
        for (size_t s = 0; s < all_coeffs.size(); s++) {
            ZZ_pE expected;
            if (num_rows > 1) {
                if (!code.is_small_ring) {
                    // Large Ring: 直接内积
                    const auto& coeffs = *all_coeffs[s];
                    expected = inner_product_gr(coeffs.data(), items.data(), num_rows);
                } else {
                    // Small Ring: coeffs 在 base ring,  items 在 ext ring
                    // 但由于 GR(p^s,r) ⊂ GR(p^s,kr), 且 NTL 中
                    // base ring 元素被视为 ext ring 中低次元素,
                    // 需要将 coeffs 的 ZZ_pX 提升到 ext ring 后做内积
                    //
                    // 论文 Lemma 4: c ∈ GR(p^s,r) 乘以 b ∈ GR(p^s,kr)
                    // = 对 b 的 k 个分量各乘 c → 结果仍在 GR(p^s,kr)
                    const auto& coeffs = *all_coeffs[s];

                    // coeffs 存储为 base ring 下的 ZZ_pE, 但当前 context 是 ext ring
                    // 需要把每个 coeff 的 ZZ_pX rep 嵌入到 ext ring 中
                    clear(expected);
                    for (long i = 0; i < num_rows; i++) {
                        // coeff 的 ZZ_pX 表示 (来自 base ring, 度 < base_r)
                        // 直接 to_ZZ_pE 即可: 在 ext ring 中低次多项式自然嵌入
                        ZZ_pE coeff_ext = to_ZZ_pE(rep(coeffs[i]));
                        expected += coeff_ext * items[i];
                    }
                }
            } else {
                expected = items[0];
            }
            if (expected != encoded_combined[s][col]) {
                std::cout << "VERIFY FAIL: column consistency for row set "
                          << s << " at column " << col << "\n";
                return false;
            }
        }
    }

    return true;
}
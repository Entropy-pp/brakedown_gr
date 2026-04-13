#include <brakedown_pcs_gr.h>
#include <gr.h>
#include <cassert>
#include <cstdlib>
#include <chrono>
#include <iostream>

using namespace NTL;

// ============================================================
// Inner product over GR
// ============================================================
// static ZZ_pE inner_product_gr(const ZZ_pE* a, const ZZ_pE* b, long len) {
//     ZZ_pE result;
//     clear(result);
//     for (long i = 0; i < len; i++) {
//         result += a[i] * b[i];
//     }
//     return result;
// }

static ZZ_pE inner_product_gr(const ZZ_pE* a, const ZZ_pE* b, long len) {
    ZZ_pX sum_poly;
    ZZ_pX tmp_poly;
    clear(sum_poly);

    for (long i = 0; i < len; i++) {
        // rep() 获取 ZZ_pE 底层的 ZZ_pX 多项式表示
        // mul() 和 add() 自动处理系数在 ZZ_p (模 2^k) 的运算
        // 这里避开了每次相乘后昂贵的多项式取模开销
        mul(tmp_poly, rep(a[i]), rep(b[i])); 
        add(sum_poly, sum_poly, tmp_poly);
    }

    // 循环结束后，仅进行一次针对定义多项式 P(x) 的模约减
    return to_ZZ_pE(sum_poly);
}

static ZZ_pE inner_product_gr(const std::vector<ZZ_pE>& a,
                               const std::vector<ZZ_pE>& b) {
    assert(a.size() == b.size());
    return inner_product_gr(a.data(), b.data(), (long)a.size());
}

// ============================================================
// 辅助: 计算 Small Ring 下��际的 proximity test 次数
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
// 辅助: 将 NTL context 切换到 base ring / ext ring
// 注意: 只切 ZZ_pE::init(), 不动 ZZ_p!
// ============================================================
static void switch_to_base_ring(long base_degree) {
    ZZ_pX base_mod = primitiveIrredPoly(base_degree);
    ZZ_pE::init(base_mod);
}

static void switch_to_ext_ring(long ext_degree) {
    ZZ_pX ext_mod = primitiveIrredPoly(ext_degree);
    ZZ_pE::init(ext_mod);
}

// ============================================================
// 辅助: 将一��� base ring 元素 pack 为 ext ring 向量
// 调用前必须在 ext ring context
// ============================================================
static std::vector<ZZ_pE> pack_row(
    const std::vector<ZZ_pX>& row_polys,
    long row_len, long packed_len, long pf, long base_r)
{
    std::vector<ZZ_pE> packed(packed_len);
    for (long j = 0; j < packed_len; j++) {
        std::vector<ZZ_pX> group(pf);
        for (long t = 0; t < pf; t++) {
            long src = j * pf + t;
            if (src < row_len) group[t] = row_polys[src];
        }
        packed[j] = pack_elements_pub(group, pf, base_r);
    }
    return packed;
}

// ============================================================
// 计时辅助宏
// ============================================================
using hrc = std::chrono::high_resolution_clock;

static double ms_between(hrc::time_point start, hrc::time_point end) {
    return std::chrono::duration<double, std::milli>(end - start).count();
}

// ============================================================
// COMMIT — 带细分计时版本
// ============================================================
BrakedownCommitment brakedown_commit(
    const BrakedownCodeGR& code,
    const ZZ_pE* poly_coeffs,
    long n,
    long num_rows,
    CommitTiming& timing)
{
    auto t_total_start = hrc::now();

    BrakedownCommitment comm;
    comm.num_rows = num_rows;
    comm.codeword_len = code.codeword_len;
    long row_len = code.row_len;

    assert(n <= num_rows * row_len);

    // 初始化 timing
    timing.fill_rows_ms   = 0;
    timing.pack_rows_ms   = 0;
    timing.encode_rows_ms = 0;
    timing.column_hash_ms = 0;
    timing.merkle_build_ms = 0;
    timing.total_ms       = 0;

    if (!code.is_small_ring) {
        // ====== Large Ring ======
        comm.encoded_rows.resize((long)num_rows * code.codeword_len);

        // --- Step 1: Fill rows ---
        auto t1 = hrc::now();
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
        }
        auto t2 = hrc::now();
        timing.fill_rows_ms = ms_between(t1, t2);

        // --- Step 2: Encode rows ---
        auto t3 = hrc::now();
        // for (long i = 0; i < num_rows; i++) {
        //     brakedown_encode(code, &comm.encoded_rows[i * code.codeword_len]);
        // }

        // 在 commit 中累加每行的 encode 细分计时
    EncodeTiming enc_timing_total = {0, 0, 0, 0, 0};

    for (long i = 0; i < num_rows; i++) {
        EncodeTiming enc_t;
        brakedown_encode(code, &comm.encoded_rows[i * code.codeword_len], enc_t);
        enc_timing_total.forward_pass_ms  += enc_t.forward_pass_ms;
        enc_timing_total.bottom_rs_ms     += enc_t.bottom_rs_ms;
        enc_timing_total.backward_pass_ms += enc_t.backward_pass_ms;
        enc_timing_total.assemble_ms      += enc_t.assemble_ms;
        enc_timing_total.total_ms         += enc_t.total_ms;
    }

    std::cout << "  Encode breakdown (all rows summed):" << std::endl;
    std::cout << "    Forward pass:  " << enc_timing_total.forward_pass_ms  << " ms" << std::endl;
    std::cout << "    Bottom RS:     " << enc_timing_total.bottom_rs_ms     << " ms" << std::endl;
    std::cout << "    Backward pass: " << enc_timing_total.backward_pass_ms << " ms" << std::endl;
    std::cout << "    Assemble:      " << enc_timing_total.assemble_ms      << " ms" << std::endl;
    std::cout << "    Total:         " << enc_timing_total.total_ms         << " ms" << std::endl;

        auto t4 = hrc::now();
        timing.encode_rows_ms = ms_between(t3, t4);

    } else {
        // ====== Small Ring: pack → encode in ext ring ======
        comm.encoded_rows.resize((long)num_rows * code.codeword_len);

        double fill_acc = 0, pack_acc = 0, encode_acc = 0;

        for (long i = 0; i < num_rows; i++) {
            // --- Fill: 提取 base ring 系数 ---
            auto tf1 = hrc::now();
            switch_to_base_ring(code.base_degree);
            std::vector<ZZ_pX> row_polys(row_len);
            for (long j = 0; j < row_len; j++) {
                long idx = i * row_len + j;
                if (idx < n)
                    row_polys[j] = rep(poly_coeffs[idx]);
            }
            auto tf2 = hrc::now();
            fill_acc += ms_between(tf1, tf2);

            // --- Pack: base ring → ext ring ---
            auto tp1 = hrc::now();
            switch_to_ext_ring(code.ext_degree);
            long pf = code.packing_factor;
            long base_r = code.base_degree;
            long packed_len = code.packed_row_len;
            long cw_len = code.codeword_len;

            ZZ_pE* cw_ptr = &comm.encoded_rows[i * cw_len];
            auto packed = pack_row(row_polys, row_len, packed_len, pf, base_r);
            for (long j = 0; j < packed_len; j++) cw_ptr[j] = packed[j];
            for (long j = packed_len; j < cw_len; j++) clear(cw_ptr[j]);
            auto tp2 = hrc::now();
            pack_acc += ms_between(tp1, tp2);

            // --- Encode ---
            auto te1 = hrc::now();
            brakedown_encode(code, cw_ptr);
            auto te2 = hrc::now();
            encode_acc += ms_between(te1, te2);
        }
        switch_to_ext_ring(code.ext_degree);

        timing.fill_rows_ms   = fill_acc;
        timing.pack_rows_ms   = pack_acc;
        timing.encode_rows_ms = encode_acc;
    }

    // --- Step 3: Column hashing ---
    auto t5 = hrc::now();
    std::vector<std::vector<unsigned char>> leaf_hashes(code.codeword_len);
    for (long col = 0; col < code.codeword_len; col++) {
        std::vector<ZZ_pE> column(num_rows);
        for (long row = 0; row < num_rows; row++) {
            column[row] = comm.encoded_rows[row * code.codeword_len + col];
        }
        leaf_hashes[col] = hash_column(column.data(), num_rows);
    }
    auto t6 = hrc::now();
    timing.column_hash_ms = ms_between(t5, t6);

    // --- Step 4: Build Merkle tree ---
    auto t7 = hrc::now();
    comm.tree = build_merkle_tree(leaf_hashes);
    comm.root = comm.tree.root;
    auto t8 = hrc::now();
    timing.merkle_build_ms = ms_between(t7, t8);

    auto t_total_end = hrc::now();
    timing.total_ms = ms_between(t_total_start, t_total_end);

    return comm;
}

// ============================================================
// COMMIT — 原有接口 (不带计时), 委托给带计时版本
// ============================================================
BrakedownCommitment brakedown_commit(
    const BrakedownCodeGR& code,
    const ZZ_pE* poly_coeffs,
    long n,
    long num_rows)
{
    CommitTiming timing;
    return brakedown_commit(code, poly_coeffs, n, num_rows, timing);
}

// ============================================================
// PROVE — 核心修正
//
// Small Ring 关键变化:
//   combined_rows 不再存 row_len 个 base ring 元素,
//   而是存 packed_row_len 个 ext ring 元素.
//
//   做法: 先把每行 pack 到 ext ring, 然后在 ext ring 中
//   做 ∑ c_i * packed_row_i (c_i 嵌入到 ext ring).
//
//   这保证了和列一致性检查的左侧计算完全一致.
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
    long num_prox = effective_num_prox(code, num_rows);

    auto t_prove_start = hrc::now();

    if (!code.is_small_ring) {
        // ====== Large Ring: 原有逻辑 ======
        auto combine_rows_large = [&](const std::vector<ZZ_pE>& coeffs) {
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

        auto t1 = hrc::now();
        proof.prox_coeffs.resize(num_prox);
        for (long t = 0; t < num_prox; t++) {
            proof.prox_coeffs[t].resize(num_rows);
            for (long i = 0; i < num_rows; i++)
                proof.prox_coeffs[t][i] = random_ZZ_pE();
            proof.combined_rows.push_back(combine_rows_large(proof.prox_coeffs[t]));
        }
        auto t2 = hrc::now();

        auto q1_combined = combine_rows_large(q1);
        auto t3 = hrc::now();

        proof.eval_value = inner_product_gr(q1_combined, q2);
        auto t4 = hrc::now();

        proof.combined_rows.push_back(q1_combined);

        std::cout << "  Prove breakdown (Large Ring):" << std::endl;
        std::cout << "    num_rows=" << num_rows << " row_len=" << row_len
                  << " num_prox=" << num_prox << " n=" << n << std::endl;
        std::cout << "    Prox combine:    " << ms_between(t1, t2) << " ms" << std::endl;
        std::cout << "    Q1 combine:      " << ms_between(t2, t3) << " ms" << std::endl;
        std::cout << "    Inner product:   " << ms_between(t3, t4) << " ms" << std::endl;

    } else {
        // ====== Small Ring: 使用分量独立的乘法 ======
        long pf = code.packing_factor;
        long base_r = code.base_degree;
        long packed_len = code.packed_row_len;

        // 获取 base ring 的模多项式（用于分量独立乘法）
        ZZ_pX base_mod = get_base_ring_modulus(base_r);

        // Step 1: 在 base ring context 下提取所有行的 ZZ_pX 系数
        switch_to_base_ring(code.base_degree);
        std::vector<std::vector<ZZ_pX>> all_row_polys(num_rows);
        for (long i = 0; i < num_rows; i++) {
            all_row_polys[i].resize(row_len);
            for (long j = 0; j < row_len; j++) {
                long idx = i * row_len + j;
                if (idx < n) all_row_polys[i][j] = rep(poly_coeffs[idx]);
            }
        }

        // Step 2: 切到 ext ring, pack 每行
        switch_to_ext_ring(code.ext_degree);
        std::vector<std::vector<ZZ_pE>> packed_rows(num_rows);
        for (long i = 0; i < num_rows; i++) {
            packed_rows[i] = pack_row(all_row_polys[i], row_len, packed_len, pf, base_r);
        }

        // Step 3: 生成 prox coeffs (在 base ring context 下)
        switch_to_base_ring(code.base_degree);
        proof.prox_coeffs.resize(num_prox);
        std::vector<std::vector<ZZ_pX>> prox_coeffs_polys(num_prox);
        for (long t = 0; t < num_prox; t++) {
            proof.prox_coeffs[t].resize(num_rows);
            prox_coeffs_polys[t].resize(num_rows);
            for (long i = 0; i < num_rows; i++) {
                proof.prox_coeffs[t][i] = randomNonZeroInExceptionalSet();
                prox_coeffs_polys[t][i] = rep(proof.prox_coeffs[t][i]);
            }
        }

        // 同样提取 q1 的 ZZ_pX
        std::vector<ZZ_pX> q1_polys(num_rows);
        for (long i = 0; i < num_rows; i++) {
            q1_polys[i] = rep(q1[i]);
        }

        // Step 4: 切到 ext ring, 使用分量独立乘法做 combine
        switch_to_ext_ring(code.ext_degree);

        // ★ 关键修改：使用 component_wise_scalar_mul_add 代替直接乘法
        auto combine_packed_correct = [&](const std::vector<ZZ_pX>& coeff_polys) {
            std::vector<ZZ_pE> result(packed_len);
            for (long j = 0; j < packed_len; j++) clear(result[j]);
            for (long i = 0; i < num_rows; i++) {
                // 使用分量独立乘法：对每个分量独立乘以系数
                component_wise_scalar_mul_add(coeff_polys[i],
                                               packed_rows[i].data(),
                                               result.data(),
                                               packed_len,
                                               pf, base_r, base_mod);
            }
            return result;
        };

        for (long t = 0; t < num_prox; t++) {
            proof.combined_rows.push_back(combine_packed_correct(prox_coeffs_polys[t]));
        }
        // q1 的 combined packed row
        proof.combined_rows.push_back(combine_packed_correct(q1_polys));

        // ★ 计算 eval_value: 从 combined_rows.back() unpack 后和 q2 内积
        //    必须和 Verify 的计算方式完全一致!
        {
            const auto& q1_packed = proof.combined_rows.back();

            // 在 ext ring context 下 unpack
            std::vector<ZZ_pX> unpacked_polys;
            for (long j = 0; j < packed_len; j++) {
                auto parts = unpack_element_pub(q1_packed[j], pf, base_r);
                for (long t = 0; t < pf; t++) {
                    unpacked_polys.push_back(parts[t]);
                }
            }

            // 切到 base ring, 和 q2 做内积
            switch_to_base_ring(code.base_degree);
            ZZ_pE eval_result;
            clear(eval_result);
            for (long j = 0; j < row_len; j++) {
                ZZ_pE elem = to_ZZ_pE(unpacked_polys[j]);
                eval_result += elem * q2[j];
            }
            proof.eval_value = eval_result;
        }
    }

    // Column openings (在 ext ring context 下)
    if (code.is_small_ring) switch_to_ext_ring(code.ext_degree);

    auto t_col_start = hrc::now();

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

    auto t_col_end = hrc::now();
    auto t_prove_end = hrc::now();

    std::cout << "    Column openings: " << ms_between(t_col_start, t_col_end) << " ms" << std::endl;
    std::cout << "    Prove total:     " << ms_between(t_prove_start, t_prove_end) << " ms" << std::endl;

    return proof;
}

// ============================================================
// VERIFY — 核心修正
//
// Small Ring 关键变化:
//   combined_rows 现在是 packed_row_len 个 ext ring 元素,
//   不再需要 pack; 直接 encode 即可.
//
//   eval consistency: 需要 unpack combined_rows.back() 回
//   base ring, 然后做 inner product with q2.
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
    long num_prox     = effective_num_prox(code, num_rows);

    if ((long)proof.combined_rows.size() != num_prox + 1) {
        std::cout << "VERIFY FAIL: wrong number of combined_rows ("
                  << proof.combined_rows.size() << " vs expected "
                  << num_prox + 1 << ")\n";
        return false;
    }

    if (!code.is_small_ring) {
        // ====== Large Ring: 原有逻辑 ======

        // 1. Eval consistency
        const auto& q1_row = proof.combined_rows.back();
        ZZ_pE computed_eval = inner_product_gr(q1_row, q2);
        if (computed_eval != proof.eval_value) {
            std::cout << "VERIFY FAIL: evaluation consistency\n";
            return false;
        }

        // 2. Build coeff lists
        std::vector<const std::vector<ZZ_pE>*> all_coeffs;
        for (long t = 0; t < num_prox; t++) all_coeffs.push_back(&proof.prox_coeffs[t]);
        all_coeffs.push_back(&q1);

        // 3. Re-encode
        std::vector<std::vector<ZZ_pE>> encoded_combined(proof.combined_rows.size());
        for (size_t t = 0; t < proof.combined_rows.size(); t++) {
            encoded_combined[t].resize(codeword_len);
            for (long j = 0; j < row_len; j++)
                encoded_combined[t][j] = proof.combined_rows[t][j];
            for (long j = row_len; j < codeword_len; j++)
                clear(encoded_combined[t][j]);
            brakedown_encode(code, encoded_combined[t].data());
        }

        // 4. Column consistency
        for (long t = 0; t < code.num_col_open; t++) {
            long col = proof.column_indices[t];
            const auto& items = proof.column_items[t];
            const auto& path  = proof.merkle_paths[t];
            if ((long)items.size() != num_rows) return false;

            auto leaf_hash = hash_column(items.data(), num_rows);
            if (!verify_merkle_path(leaf_hash, col, path, comm.root)) {
                std::cout << "VERIFY FAIL: Merkle path at column " << col << "\n";
                return false;
            }
            for (size_t s = 0; s < all_coeffs.size(); s++) {
                ZZ_pE expected = (num_rows > 1)
                    ? inner_product_gr(all_coeffs[s]->data(), items.data(), num_rows)
                    : items[0];
                if (expected != encoded_combined[s][col]) {
                    std::cout << "VERIFY FAIL: column consistency for row set "
                              << s << " at column " << col << "\n";
                    return false;
                }
            }
        }
        return true;

    } else {
        // ====== Small Ring: 使用分量独立乘法 ======
        long pf = code.packing_factor;
        long base_r = code.base_degree;
        long packed_len = code.packed_row_len;

        // 获取 base ring 的模多项式
        ZZ_pX base_mod = get_base_ring_modulus(base_r);

        // 1. Eval consistency
        switch_to_ext_ring(code.ext_degree);
        const auto& q1_packed = proof.combined_rows.back();

        std::vector<ZZ_pX> unpacked_polys;
        for (long j = 0; j < packed_len; j++) {
            auto parts = unpack_element_pub(q1_packed[j], pf, base_r);
            for (long t = 0; t < pf; t++) {
                unpacked_polys.push_back(parts[t]);
            }
        }

        switch_to_base_ring(code.base_degree);
        ZZ_pE computed_eval;
        clear(computed_eval);
        for (long j = 0; j < row_len; j++) {
            ZZ_pE elem = to_ZZ_pE(unpacked_polys[j]);
            computed_eval += elem * q2[j];
        }
        if (computed_eval != proof.eval_value) {
            std::cout << "VERIFY FAIL: evaluation consistency\n";
            return false;
        }

        // 2. 提取 prox_coeffs 和 q1 的 ZZ_pX (在 base ring context 下)
        std::vector<std::vector<ZZ_pX>> all_coeff_polys;
        for (long t = 0; t < num_prox; t++) {
            std::vector<ZZ_pX> cp(num_rows);
            for (long i = 0; i < num_rows; i++) {
                cp[i] = rep(proof.prox_coeffs[t][i]);
            }
            all_coeff_polys.push_back(cp);
        }
        {
            std::vector<ZZ_pX> qp(num_rows);
            for (long i = 0; i < num_rows; i++) {
                qp[i] = rep(q1[i]);
            }
            all_coeff_polys.push_back(qp);
        }

        // 3. Re-encode
        switch_to_ext_ring(code.ext_degree);
        std::vector<std::vector<ZZ_pE>> encoded_combined(proof.combined_rows.size());
        for (size_t t = 0; t < proof.combined_rows.size(); t++) {
            encoded_combined[t].resize(codeword_len);
            for (long j = 0; j < packed_len; j++)
                encoded_combined[t][j] = proof.combined_rows[t][j];
            for (long j = packed_len; j < codeword_len; j++)
                clear(encoded_combined[t][j]);
            brakedown_encode(code, encoded_combined[t].data());
        }

        // 4. Column consistency (在 ext ring context 下)
        // ★ 关键修改：使用分量独立乘法
        for (long t = 0; t < code.num_col_open; t++) {
            long col = proof.column_indices[t];
            const auto& items = proof.column_items[t];
            const auto& path  = proof.merkle_paths[t];
            if ((long)items.size() != num_rows) return false;

            auto leaf_hash = hash_column(items.data(), num_rows);
            if (!verify_merkle_path(leaf_hash, col, path, comm.root)) {
                std::cout << "VERIFY FAIL: Merkle path at column " << col << "\n";
                return false;
            }

            for (size_t s = 0; s < all_coeff_polys.size(); s++) {
                ZZ_pE expected;
                clear(expected);
                if (num_rows > 1) {
                    // ★ 使用分量独立乘法计算 expected
                    for (long i = 0; i < num_rows; i++) {
                        ZZ_pE product = component_wise_scalar_mul(
                            all_coeff_polys[s][i], items[i],
                            pf, base_r, base_mod);
                        expected += product;
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
}
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
// COMMIT
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

    long num_prox = (num_rows > 1) ? code.num_prox_test : 0;

    // ============ Small Ring: 增加重复次数 ============
    // 论文: Verifier 从 GR(p^s,r) 选 challenge，重复 ceil(λ/r) 次
    if (code.is_small_ring && num_rows > 1) {
        long repetitions = (code.lambda + code.base_degree - 1) / code.base_degree;
        num_prox = std::max(num_prox, repetitions);
    }

    proof.prox_coeffs.resize(num_prox);
    for (long t = 0; t < num_prox; t++) {
        proof.prox_coeffs[t].resize(num_rows);
        for (long i = 0; i < num_rows; i++) {
            if (code.is_small_ring) {
                // Small ring: challenge 来自 GR(p^s,r) 的 exceptional set
                // 确保 a_i - a_j 可逆 (exceptional set 性质)
                proof.prox_coeffs[t][i] = randomInExceptionalSet();
                // 但需要非零
                while (IsZero(proof.prox_coeffs[t][i]))
                    proof.prox_coeffs[t][i] = randomInExceptionalSet();
            } else {
                // Large ring: 从 R* 中随机选
                proof.prox_coeffs[t][i] = randomInvertible();
            }
        }
        proof.combined_rows.push_back(combine_rows(proof.prox_coeffs[t]));
    }

    // Evaluation consistency: combine with q1 (不变)
    auto q1_combined = combine_rows(q1);
    proof.eval_value = inner_product_gr(q1_combined, q2);
    proof.combined_rows.push_back(q1_combined);

    // Column openings
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

    // FIX: When num_rows == 1, skip proximity testing (matches Rust code)
    long num_prox = (num_rows > 1) ? code.num_prox_test : 0;

    // Check proof structure
    if ((long)proof.combined_rows.size() != num_prox + 1) {
        std::cout << "VERIFY FAIL: wrong number of combined_rows ("
                  << proof.combined_rows.size() << " vs expected "
                  << num_prox + 1 << ")\n";
        return false;
    }

    // 1. Evaluation consistency
    const auto& q1_row = proof.combined_rows.back();
    ZZ_pE computed_eval = inner_product_gr(q1_row, q2);
    if (computed_eval != proof.eval_value) {
        std::cout << "VERIFY FAIL: evaluation consistency\n";
        return false;
    }

    // 2. Build coefficient lists for column consistency checks
    //    Proximity tests use prox_coeffs; final row uses q1.
    std::vector<const std::vector<ZZ_pE>*> all_coeffs;
    for (long t = 0; t < num_prox; t++) {
        all_coeffs.push_back(&proof.prox_coeffs[t]);
    }
    all_coeffs.push_back(&q1);

    // 3. Re-encode all combined rows
    std::vector<std::vector<ZZ_pE>> encoded_combined(proof.combined_rows.size());
    for (size_t t = 0; t < proof.combined_rows.size(); t++) {
        encoded_combined[t].resize(codeword_len);
        for (long j = 0; j < row_len; j++) {
            encoded_combined[t][j] = proof.combined_rows[t][j];
        }
        for (long j = row_len; j < codeword_len; j++) {
            clear(encoded_combined[t][j]);
        }
        brakedown_encode(code, encoded_combined[t].data());
    }

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

        // FIX: Column consistency — use inner_product when num_rows > 1,
        //      use items[0] directly when num_rows == 1 (matches Rust code)
        for (size_t s = 0; s < all_coeffs.size(); s++) {
            ZZ_pE expected;
            if (num_rows > 1) {
                const auto& coeffs = *all_coeffs[s];
                expected = inner_product_gr(coeffs.data(), items.data(), num_rows);
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
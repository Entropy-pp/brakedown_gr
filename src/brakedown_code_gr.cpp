#include <brakedown_code_gr.h>
#include <gr.h>
#include <cassert>
#include <algorithm>

using namespace NTL;

// ============================================================
// Reed-Solomon over GR: evaluate polynomial defined by
// coefficients input[0..input_len-1] at exceptional set
// points {1, 2, ..., output_len}.
// ============================================================
static void reed_solomon_gr(const ZZ_pE* input, long input_len,
                            ZZ_pE* output, long output_len)
{
    // 确保求值点不会超出 exceptional set 大小（否则会重复）
    assert(output_len < (1L << ZZ_pE::degree()));

    ZZ_pEX poly;
    for (long i = 0; i < input_len; i++) {
        SetCoeff(poly, i, input[i]);
    }
    for (long i = 0; i < output_len; i++) {
        ZZ_pE pt = indexedElementInExceptionalSet(i + 1);
        output[i] = eval(poly, pt);
    }
}

BrakedownCodeGR brakedown_code_setup(long row_len, long k, long degree) {
    BrakedownCodeGR code;
    code.spec = getDefaultSpec();
    code.row_len = row_len;

    long log2_q = k * degree;
    long n_0 = std::min(20L, row_len - 1);
    if (n_0 < 1) n_0 = 1;

    // Compute dimensions
    std::vector<LevelDim> a_dims, b_dims;
    brakedown_dimensions(code.spec, log2_q, row_len, n_0, a_dims, b_dims);

    // If no recursive levels possible, use plain RS
    if (a_dims.empty()) {
        code.codeword_len = iceil(row_len * code.spec.R);
        if (code.codeword_len < row_len + 1) code.codeword_len = row_len + 1;
        code.num_col_open = brakedown_num_column_opening(code.spec);
        code.num_prox_test = 0;
        return code;
    }

    code.codeword_len = brakedown_codeword_len(code.spec, log2_q, row_len, n_0);
    code.num_col_open = brakedown_num_column_opening(code.spec);
    code.num_prox_test = brakedown_num_proximity_testing(code.spec, log2_q, row_len, n_0);

    // Create random sparse matrices
    code.a_mats.resize(a_dims.size());
    code.b_mats.resize(b_dims.size());
    for (size_t i = 0; i < a_dims.size(); i++) {
        code.a_mats[i] = createRandomSparseMatrix(a_dims[i]);
    }
    for (size_t i = 0; i < b_dims.size(); i++) {
        code.b_mats[i] = createRandomSparseMatrix(b_dims[i]);
    }

    return code;
}

// ============================================================
// Encode: codeword[0..row_len-1] already holds the message.
//
// Layout of the codeword after encoding:
//
//   [message(row_len)] [A0_out(a0.m)] [A1_out(a1.m)] ... [RS_out] [BL_out] ... [B0_out]
//
// Forward pass:  A_i compresses from the block that starts at
//                offset_in into the block right after it.
//
// Bottom:        The smallest A output gets RS-encoded.
//
// Backward pass: B_i expands from the RS+previous-B output
//                region into the final tail blocks.
//
// For simplicity & correctness, we use a scratch buffer
// approach: forward pass produces intermediate vectors,
// RS encodes the last one, backward pass expands,
// then we assemble the final codeword.
// ============================================================
void brakedown_encode(const BrakedownCodeGR& code, ZZ_pE* codeword) {
    long row_len = code.row_len;
    long cw_len  = code.codeword_len;

    // ---- Trivial case: no recursive levels, pure RS ----
    if (code.a_mats.empty()) {
        long rs_len = cw_len - row_len;
        if (rs_len > 0) {
            reed_solomon_gr(codeword, row_len, codeword + row_len, rs_len);
        }
        return;
    }

    long L = (long)code.a_mats.size();

    // ---- Forward pass: produce intermediate compressions ----
    // intermediates[i] = A_i * intermediates[i-1]
    // intermediates[0] = A_0 * message
    std::vector<std::vector<ZZ_pE>> intermediates(L);
    {
        const ZZ_pE* src = codeword; // message
        for (long i = 0; i < L; i++) {
            const auto& A = code.a_mats[i];
            intermediates[i].resize(A.n_cols);
            for (long j = 0; j < A.n_cols; j++) clear(intermediates[i][j]);
            sparseMatVecMul(A, src, intermediates[i].data());
            src = intermediates[i].data();
        }
    }

    // ---- Bottom: RS-encode the last intermediate ----
    const auto& last_inter = intermediates[L - 1];
    long rs_input_len  = (long)last_inter.size();
    long rs_output_len = iceil(rs_input_len * code.spec.R) - rs_input_len;
    if (rs_output_len < 1) rs_output_len = 1;
    std::vector<ZZ_pE> rs_output(rs_output_len);
    reed_solomon_gr(last_inter.data(), rs_input_len,
                    rs_output.data(), rs_output_len);

    // ---- Backward pass: B_i expands ----
    // B_L-1 takes input = [rs_input | rs_output] (or just rs_output)
    // B_i   takes input = [intermediates[i] | B_{i+1} output]
    std::vector<std::vector<ZZ_pE>> b_outputs(L);
    {
        // For the last B matrix, input is the RS output
        const auto& B = code.b_mats[L - 1];
        b_outputs[L - 1].resize(B.n_cols);
        for (long j = 0; j < B.n_cols; j++) clear(b_outputs[L - 1][j]);

        // B input: combine rs_input (intermediates[L-1]) + rs_output
        long b_input_len = B.n_rows;
        std::vector<ZZ_pE> b_input(b_input_len);
        long copy1 = std::min(rs_input_len, b_input_len);
        for (long j = 0; j < copy1; j++) b_input[j] = intermediates[L - 1][j];
        long copy2 = std::min(rs_output_len, b_input_len - copy1);
        for (long j = 0; j < copy2; j++) b_input[copy1 + j] = rs_output[j];

        sparseMatVecMul(B, b_input.data(), b_outputs[L - 1].data());

        // Remaining B matrices (backward)
        for (long i = L - 2; i >= 0; i--) {
            const auto& Bi = code.b_mats[i];
            b_outputs[i].resize(Bi.n_cols);
            for (long j = 0; j < Bi.n_cols; j++) clear(b_outputs[i][j]);

            long bi_input_len = Bi.n_rows;
            std::vector<ZZ_pE> bi_input(bi_input_len);
            // Input: [intermediates[i] | b_outputs[i+1]]
            long c1 = std::min((long)intermediates[i].size(), bi_input_len);
            for (long j = 0; j < c1; j++) bi_input[j] = intermediates[i][j];
            long c2 = std::min((long)b_outputs[i + 1].size(), bi_input_len - c1);
            for (long j = 0; j < c2; j++) bi_input[c1 + j] = b_outputs[i + 1][j];

            sparseMatVecMul(Bi, bi_input.data(), b_outputs[i].data());
        }
    }

    // ---- Assemble the final codeword ----
    // Layout: [message] [intermediates[0]] [intermediates[1]] ... [rs_output] [b_outputs[L-1]] ... [b_outputs[0]]
    long offset = row_len;

    // Intermediates (except last, which is embedded in RS input)
    for (long i = 0; i < L - 1; i++) {
        long sz = (long)intermediates[i].size();
        long copy_len = std::min(sz, cw_len - offset);
        for (long j = 0; j < copy_len; j++) {
            codeword[offset + j] = intermediates[i][j];
        }
        offset += copy_len;
    }

    // Last intermediate (= RS input)
    {
        long sz = rs_input_len;
        long copy_len = std::min(sz, cw_len - offset);
        for (long j = 0; j < copy_len; j++) {
            codeword[offset + j] = intermediates[L - 1][j];
        }
        offset += copy_len;
    }

    // RS output
    {
        long copy_len = std::min(rs_output_len, cw_len - offset);
        for (long j = 0; j < copy_len; j++) {
            codeword[offset + j] = rs_output[j];
        }
        offset += copy_len;
    }

    // B outputs (from last to first)
    for (long i = L - 1; i >= 0; i--) {
        long sz = (long)b_outputs[i].size();
        long copy_len = std::min(sz, cw_len - offset);
        for (long j = 0; j < copy_len; j++) {
            codeword[offset + j] = b_outputs[i][j];
        }
        offset += copy_len;
    }

    // Zero any remaining space (shouldn't happen if dimensions are correct)
    for (long j = offset; j < cw_len; j++) {
        clear(codeword[j]);
    }
}
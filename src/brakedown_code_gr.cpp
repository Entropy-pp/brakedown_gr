#include <brakedown_code_gr.h>
#include <gr.h>
#include <cassert>
#include <algorithm>
#include <chrono>

using namespace NTL;

// ============================================================
// 计时辅助
// ============================================================
using hrc = std::chrono::high_resolution_clock;

static double ms_between(hrc::time_point start, hrc::time_point end) {
    return std::chrono::duration<double, std::milli>(end - start).count();
}

// ============================================================
// Reed-Solomon over GR
// ============================================================
static void reed_solomon_gr(const ZZ_pE* input, long input_len,
                            ZZ_pE* output, long output_len)
{
    long deg = ZZ_pE::degree();
    if (deg < 63) {
        assert(output_len < (1L << deg));
    }
    ZZ_pEX poly;
    for (long i = 0; i < input_len; i++) {
        SetCoeff(poly, i, input[i]);
    }
    for (long i = 0; i < output_len; i++) {
        ZZ_pE pt = indexedElementInExceptionalSet(i + 1);
        output[i] = eval(poly, pt);
    }
}

// ============================================================
// brakedown_code_setup (Large Ring)
// ============================================================
BrakedownCodeGR brakedown_code_setup(long row_len, long k, long degree) {
    BrakedownCodeGR code;
    code.spec = getDefaultSpec();
    code.row_len = row_len;
    code.is_small_ring = false;
    code.packing_factor = 1;
    code.ext_degree = degree;
    code.base_degree = degree;
    code.lambda = 128;
    code.packed_row_len = row_len;

    long log2_q = k * degree;
    long n_0 = std::min(20L, row_len - 1);
    if (n_0 < 1) n_0 = 1;

    std::vector<LevelDim> a_dims, b_dims;
    brakedown_dimensions(code.spec, log2_q, row_len, n_0, a_dims, b_dims);

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
// brakedown_code_setup_auto
// ============================================================
BrakedownCodeGR brakedown_code_setup_auto(long row_len, long k, long degree, long lambda) {
    bool is_small = !isLargeRing(degree, lambda);

    if (!is_small) {
        BrakedownCodeGR code = brakedown_code_setup(row_len, k, degree);
        code.lambda = lambda;
        return code;
    }

    long pf = smallRingPackingFactor(degree, lambda);
    long ext_deg = pf * degree;
    long packed_row = (row_len + pf - 1) / pf;

    BrakedownCodeGR code;
    code.spec = getDefaultSpec();
    code.row_len = row_len;
    code.is_small_ring = true;
    code.packing_factor = pf;
    code.ext_degree = ext_deg;
    code.base_degree = degree;
    code.lambda = lambda;
    code.packed_row_len = packed_row;

    long log2_q_ext = k * ext_deg;
    long n_0 = std::min(20L, packed_row - 1);
    if (n_0 < 1) n_0 = 1;

    std::vector<LevelDim> a_dims, b_dims;

    if (packed_row <= n_0) {
        code.codeword_len = iceil(packed_row * code.spec.R);
        if (code.codeword_len < packed_row + 1) code.codeword_len = packed_row + 1;
        code.num_col_open = brakedown_num_column_opening(code.spec);
        code.num_prox_test = 0;
        return code;
    }

    brakedown_dimensions(code.spec, log2_q_ext, packed_row, n_0, a_dims, b_dims);

    if (a_dims.empty()) {
        code.codeword_len = iceil(packed_row * code.spec.R);
        if (code.codeword_len < packed_row + 1) code.codeword_len = packed_row + 1;
        code.num_col_open = brakedown_num_column_opening(code.spec);
        code.num_prox_test = 0;
        return code;
    }

    code.codeword_len = brakedown_codeword_len(code.spec, log2_q_ext, packed_row, n_0);
    code.num_col_open = brakedown_num_column_opening(code.spec);
    code.num_prox_test = brakedown_num_proximity_testing(code.spec, log2_q_ext, packed_row, n_0);

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
// brakedown_encode — 带细分计时版本
// ============================================================
void brakedown_encode(const BrakedownCodeGR& code, ZZ_pE* codeword, EncodeTiming& timing) {
    auto t_total_start = hrc::now();

    long row_len = code.is_small_ring ? code.packed_row_len : code.row_len;
    long cw_len  = code.codeword_len;

    // 初始化
    timing.forward_pass_ms  = 0;
    timing.bottom_rs_ms     = 0;
    timing.backward_pass_ms = 0;
    timing.assemble_ms      = 0;
    timing.total_ms         = 0;

    if (code.a_mats.empty()) {
        auto t1 = hrc::now();
        long rs_len = cw_len - row_len;
        if (rs_len > 0) {
            reed_solomon_gr(codeword, row_len, codeword + row_len, rs_len);
        }
        auto t2 = hrc::now();
        timing.bottom_rs_ms = ms_between(t1, t2);
        timing.total_ms = ms_between(t_total_start, t2);
        return;
    }

    long L = (long)code.a_mats.size();

    // ---- Forward pass ----
    auto t_fwd_start = hrc::now();
    std::vector<std::vector<ZZ_pE>> intermediates(L);
    {
        const ZZ_pE* src = codeword;
        for (long i = 0; i < L; i++) {
            const auto& A = code.a_mats[i];
            intermediates[i].resize(A.n_cols);
            for (long j = 0; j < A.n_cols; j++) clear(intermediates[i][j]);
            sparseMatVecMul(A, src, intermediates[i].data());
            src = intermediates[i].data();
        }
    }
    auto t_fwd_end = hrc::now();
    timing.forward_pass_ms = ms_between(t_fwd_start, t_fwd_end);

    // ---- Bottom RS ----
    auto t_rs_start = hrc::now();
    const auto& last_inter = intermediates[L - 1];
    long rs_input_len  = (long)last_inter.size();
    long rs_output_len = iceil(rs_input_len * code.spec.R) - rs_input_len;
    if (rs_output_len < 1) rs_output_len = 1;
    std::vector<ZZ_pE> rs_output(rs_output_len);
    reed_solomon_gr(last_inter.data(), rs_input_len,
                    rs_output.data(), rs_output_len);
    auto t_rs_end = hrc::now();
    timing.bottom_rs_ms = ms_between(t_rs_start, t_rs_end);

    // ---- Backward pass ----
    auto t_bwd_start = hrc::now();
    std::vector<std::vector<ZZ_pE>> b_outputs(L);
    {
        const auto& B = code.b_mats[L - 1];
        b_outputs[L - 1].resize(B.n_cols);
        for (long j = 0; j < B.n_cols; j++) clear(b_outputs[L - 1][j]);

        long b_input_len = B.n_rows;
        std::vector<ZZ_pE> b_input(b_input_len);
        long copy1 = std::min(rs_input_len, b_input_len);
        for (long j = 0; j < copy1; j++) b_input[j] = intermediates[L - 1][j];
        long copy2 = std::min(rs_output_len, b_input_len - copy1);
        for (long j = 0; j < copy2; j++) b_input[copy1 + j] = rs_output[j];

        sparseMatVecMul(B, b_input.data(), b_outputs[L - 1].data());

        for (long i = L - 2; i >= 0; i--) {
            const auto& Bi = code.b_mats[i];
            b_outputs[i].resize(Bi.n_cols);
            for (long j = 0; j < Bi.n_cols; j++) clear(b_outputs[i][j]);

            long bi_input_len = Bi.n_rows;
            std::vector<ZZ_pE> bi_input(bi_input_len);
            long c1 = std::min((long)intermediates[i].size(), bi_input_len);
            for (long j = 0; j < c1; j++) bi_input[j] = intermediates[i][j];
            long c2 = std::min((long)b_outputs[i + 1].size(), bi_input_len - c1);
            for (long j = 0; j < c2; j++) bi_input[c1 + j] = b_outputs[i + 1][j];

            sparseMatVecMul(Bi, bi_input.data(), b_outputs[i].data());
        }
    }
    auto t_bwd_end = hrc::now();
    timing.backward_pass_ms = ms_between(t_bwd_start, t_bwd_end);

    // ---- Assemble ----
    auto t_asm_start = hrc::now();
    long offset = row_len;
    for (long i = 0; i < L - 1; i++) {
        long sz = (long)intermediates[i].size();
        long copy_len = std::min(sz, cw_len - offset);
        for (long j = 0; j < copy_len; j++) codeword[offset + j] = intermediates[i][j];
        offset += copy_len;
    }
    {
        long sz = rs_input_len;
        long copy_len = std::min(sz, cw_len - offset);
        for (long j = 0; j < copy_len; j++) codeword[offset + j] = intermediates[L - 1][j];
        offset += copy_len;
    }
    {
        long copy_len = std::min(rs_output_len, cw_len - offset);
        for (long j = 0; j < copy_len; j++) codeword[offset + j] = rs_output[j];
        offset += copy_len;
    }
    for (long i = L - 1; i >= 0; i--) {
        long sz = (long)b_outputs[i].size();
        long copy_len = std::min(sz, cw_len - offset);
        for (long j = 0; j < copy_len; j++) codeword[offset + j] = b_outputs[i][j];
        offset += copy_len;
    }
    for (long j = offset; j < cw_len; j++) clear(codeword[j]);
    auto t_asm_end = hrc::now();
    timing.assemble_ms = ms_between(t_asm_start, t_asm_end);

    timing.total_ms = ms_between(t_total_start, t_asm_end);
}

// ============================================================
// brakedown_encode — 原有接口 (不带计时), 委托给带计时版本
// ============================================================
void brakedown_encode(const BrakedownCodeGR& code, ZZ_pE* codeword) {
    EncodeTiming timing;
    brakedown_encode(code, codeword, timing);
}

// ============================================================
// Pack/Unpack 辅助函数
// ============================================================
static ZZ_pE pack_elements(const std::vector<ZZ_pX>& base_coeffs,
                            long pf, long base_r) {
    ZZ_pX packed;
    for (long j = 0; j < pf; j++) {
        for (long i = 0; i < base_r; i++) {
            ZZ_p c;
            if (j < (long)base_coeffs.size() && deg(base_coeffs[j]) >= i) {
                c = coeff(base_coeffs[j], i);
            } else {
                clear(c);
            }
            if (!IsZero(c)) {
                SetCoeff(packed, j * base_r + i, c);
            }
        }
    }
    return to_ZZ_pE(packed);
}

static std::vector<ZZ_pX> unpack_element(const ZZ_pE& ext_elem,
                                          long pf, long base_r) {
    std::vector<ZZ_pX> result(pf);
    const ZZ_pX& poly = rep(ext_elem);
    for (long j = 0; j < pf; j++) {
        for (long i = 0; i < base_r; i++) {
            ZZ_p c;
            long idx = j * base_r + i;
            if (deg(poly) >= idx) {
                c = coeff(poly, idx);
            } else {
                clear(c);
            }
            if (!IsZero(c)) {
                SetCoeff(result[j], i, c);
            }
        }
    }
    return result;
}

ZZ_pE pack_elements_pub(const std::vector<ZZ_pX>& base_coeffs,
                         long pf, long base_r) {
    return pack_elements(base_coeffs, pf, base_r);
}

std::vector<ZZ_pX> unpack_element_pub(const ZZ_pE& ext_elem,
                                       long pf, long base_r) {
    return unpack_element(ext_elem, pf, base_r);
}

// ============================================================
// Small Ring 编码: Algorithm 2 (Enc')
// ============================================================
void brakedown_encode_small_ring(const BrakedownCodeGR& code,
                                  const ZZ_pE* input_small,
                                  long input_len,
                                  ZZ_pE* output_ext)
{
    assert(code.is_small_ring);

    long pf     = code.packing_factor;
    long base_r = code.base_degree;
    long ext_r  = code.ext_degree;
    long packed_len = code.packed_row_len;
    long cw_len = code.codeword_len;

    std::vector<ZZ_pX> input_polys(input_len);
    for (long i = 0; i < input_len; i++) {
        input_polys[i] = rep(input_small[i]);
    }

    ZZ ps = ZZ_p::modulus();
    ZZ_pX ext_mod = primitiveIrredPoly(ext_r);
    ZZ_pE::init(ext_mod);

    for (long j = 0; j < packed_len; j++) {
        std::vector<ZZ_pX> group(pf);
        for (long t = 0; t < pf; t++) {
            long src_idx = j * pf + t;
            if (src_idx < input_len) {
                group[t] = input_polys[src_idx];
            }
        }
        output_ext[j] = pack_elements(group, pf, base_r);
    }

    for (long j = packed_len; j < cw_len; j++) {
        clear(output_ext[j]);
    }

    brakedown_encode(code, output_ext);
}
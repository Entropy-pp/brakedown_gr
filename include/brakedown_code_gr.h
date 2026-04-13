#ifndef BRAKEDOWN_CODE_GR_H
#define BRAKEDOWN_CODE_GR_H

#include <NTL/ZZ_pE.h>
#include <NTL/ZZ_pEX.h>
#include <vector>
#include <brakedown_params.h>
#include <sparse_matrix_gr.h>

using namespace NTL;

struct BrakedownCodeGR {
    BrakedownSpec spec;
    long row_len;
    long codeword_len;
    long num_col_open;
    long num_prox_test;

    // Large ring: 使用 ext ring 稀疏矩阵
    std::vector<SparseMatrixGR> a_mats;
    std::vector<SparseMatrixGR> b_mats;

    // Small ring: 使用 base ring 系数的稀疏矩阵
    std::vector<SparseMatrixSmallRing> a_mats_small;
    std::vector<SparseMatrixSmallRing> b_mats_small;

    bool is_small_ring;
    long packing_factor;
    long ext_degree;
    long base_degree;
    long lambda;
    long packed_row_len;
};

// Encode 阶段细分计时结果 (单位: 毫秒)
struct EncodeTiming {
    double forward_pass_ms;    // Forward pass (所有 A 矩阵稀疏乘法)
    double bottom_rs_ms;       // Bottom RS 编码
    double backward_pass_ms;   // Backward pass (所有 B 矩阵稀疏乘法)
    double assemble_ms;        // Assemble (拼装码字)
    double total_ms;           // 总计
};

BrakedownCodeGR brakedown_code_setup(long row_len, long k, long degree);
BrakedownCodeGR brakedown_code_setup_auto(long row_len, long k, long degree, long lambda = 128);

// 原有接口 (不带计时)
void brakedown_encode(const BrakedownCodeGR& code, ZZ_pE* codeword);

// 带细分计时的接口
void brakedown_encode(const BrakedownCodeGR& code, ZZ_pE* codeword, EncodeTiming& timing);

void brakedown_encode_small_ring(const BrakedownCodeGR& code,
                                  const ZZ_pE* input_small,
                                  long input_len,
                                  ZZ_pE* output_ext);

ZZ_pE pack_elements_pub(const std::vector<ZZ_pX>& base_coeffs,
                         long pf, long base_r);

std::vector<ZZ_pX> unpack_element_pub(const ZZ_pE& ext_elem,
                                       long pf, long base_r);

// ============================================================
// 小环正确乘法：将 base ring 元素乘以 packed 元素
// 关键：对每个分量独立相乘，避免分量混合
// ============================================================

// 将 base ring 元素 (ZZ_pX, degree < base_r) 乘以 packed 元素
// 返回新的 packed 元素，每个分量都乘以了 scalar
// 调用前必须在 ext ring context
ZZ_pE component_wise_scalar_mul(const ZZ_pX& scalar_poly,
                                 const ZZ_pE& packed_elem,
                                 long pf, long base_r,
                                 const ZZ_pX& base_mod);

// 向量版本：对 packed 向量的每个元素执行 component_wise_scalar_mul
void component_wise_scalar_mul_vec(const ZZ_pX& scalar_poly,
                                    const ZZ_pE* packed_vec,
                                    ZZ_pE* result,
                                    long vec_len,
                                    long pf, long base_r,
                                    const ZZ_pX& base_mod);

// 累加版本：result[j] += scalar * packed_vec[j]（分量独立）
void component_wise_scalar_mul_add(const ZZ_pX& scalar_poly,
                                    const ZZ_pE* packed_vec,
                                    ZZ_pE* result,
                                    long vec_len,
                                    long pf, long base_r,
                                    const ZZ_pX& base_mod);

// 获取 base ring 的模多项式
ZZ_pX get_base_ring_modulus(long base_r);

#endif
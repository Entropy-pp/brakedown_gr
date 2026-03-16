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

    std::vector<SparseMatrixGR> a_mats;
    std::vector<SparseMatrixGR> b_mats;

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

#endif
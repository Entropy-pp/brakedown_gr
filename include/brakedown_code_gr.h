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
    long row_len;        // 原始消息行长度（在 small_ring 上的元素数）
    long codeword_len;   // 编码后总长度 (对 small ring: GR(p^s,kr) 上的元素数)
    long num_col_open;
    long num_prox_test;

    std::vector<SparseMatrixGR> a_mats;
    std::vector<SparseMatrixGR> b_mats;

    // ============ Small Ring 支持 ============
    bool is_small_ring;     // 是否为 small ring
    long packing_factor;    // k: 将 k 个 GR(p^s,r) 元素打包为 GR(p^s,kr) 的 1 个元素
    long ext_degree;        // 扩展后的度数 kr (仅 small ring 使用)
    long base_degree;       // 原始度数 r
    long lambda;            // 安全参数
    long packed_row_len;    // = ceil(row_len / packing_factor)
};

// Setup: create the code for a given message length (row_len)
BrakedownCodeGR brakedown_code_setup(long row_len, long k, long degree);

// 带安全参数的 setup，自动判断 large/small ring
BrakedownCodeGR brakedown_code_setup_auto(long row_len, long k, long degree, long lambda = 128);

// Encode in-place (large ring): codeword[0..row_len-1] = message
void brakedown_encode(const BrakedownCodeGR& code, ZZ_pE* codeword);

// Small ring 编码 (Algorithm 2: Enc')
// 输入: row_len 个 GR(p^s,r) 元素 (base ring context)
// 输出: codeword_len 个 GR(p^s,kr) 元素 (ext ring context)
// 内部负责 pack → encode → 输出
// 调用者需保证:
//   - 调用前 NTL context 为 base ring GR(p^s,r)
//   - input_small 的 ZZ_pX 系数已在 base ring 下提取好
//   - output_ext 长度至少 code.codeword_len
void brakedown_encode_small_ring(const BrakedownCodeGR& code,
                                  const ZZ_pE* input_small,
                                  long input_len,
                                  ZZ_pE* output_ext);

// Pack 辅助: 将 pf 个 base ring 系数 ZZ_pX 打包为 ext ring 的 1 个 ZZ_pE
// (需在 ext ring context 下调用)
ZZ_pE pack_elements_pub(const std::vector<ZZ_pX>& base_coeffs,
                         long pf, long base_r);

// Unpack 辅助: 从 ext ring 的 1 个 ZZ_pE 拆为 pf 个 base ring 系数 ZZ_pX
std::vector<ZZ_pX> unpack_element_pub(const ZZ_pE& ext_elem,
                                       long pf, long base_r);

#endif
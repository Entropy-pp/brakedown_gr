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
    long codeword_len;   // 编码后总长度
    long num_col_open;
    long num_prox_test;

    std::vector<SparseMatrixGR> a_mats;
    std::vector<SparseMatrixGR> b_mats;

    // ============ 新增: Small Ring 支持 ============
    bool is_small_ring;     // 是否为 small ring
    long packing_factor;    // k: 将 k 个 GR(p^s,r) 元素打包为 GR(p^s,kr) 的 1 个元素
    long ext_degree;        // 扩展后的度数 kr (仅 small ring 使用)
    long base_degree;       // 原始度数 r
    long lambda;            // 安全参数
    // 对于 small ring, row_len 仍表示 GR(p^s,r) 上的元素数
    // 编码时实际在 GR(p^s,kr) 上操作的行长 = ceil(row_len / packing_factor)
    long packed_row_len;    // = ceil(row_len / packing_factor)
};

// Setup: create the code for a given message length (row_len)
BrakedownCodeGR brakedown_code_setup(long row_len, long k, long degree);

// 新增: 带安全参数的 setup，自动判断 large/small ring
BrakedownCodeGR brakedown_code_setup_auto(long row_len, long k, long degree, long lambda = 128);

// Encode in-place: codeword[0..row_len-1] = message, fills codeword[row_len..codeword_len-1]
// codeword must have length = code.codeword_len
void brakedown_encode(const BrakedownCodeGR& code, ZZ_pE* codeword);

/ 新增: Small ring 编码 (Algorithm 2: Enc')
// 输入: kn 个 GR(p^s,r) 元素
// 输出: 在 GR(p^s,kr) 上编码后的 tn 个元素
// 注意: 调用前需要先切换 NTL context 到 GR(p^s,kr)
void brakedown_encode_small_ring(const BrakedownCodeGR& code,
                                  const ZZ_pE* input_small,  // GR(p^s,r) 上的 kn 元素
                                  long input_len,
                                  ZZ_pE* output_ext);        // GR(p^s,kr) 上的结果

#endif
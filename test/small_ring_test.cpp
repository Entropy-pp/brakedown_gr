#include <iostream>
#include <NTL/ZZ_pE.h>
#include <NTL/ZZ_pX.h>
#include <gr.h>
#include <brakedown_code_gr.h>
#include <brakedown_pcs_gr.h>

using namespace NTL;
using namespace std;

// ============================================================
// 辅助: 切换 NTL context
// ============================================================
static void init_ring(long ps_exp, long degree) {
    ZZ_p::init(ZZ(1) << ps_exp);
    ZZ_pX mod = primitiveIrredPoly(degree);
    ZZ_pE::init(mod);
}

// ============================================================
// Small Ring PCS 端到端测试
//
// 选择参数使得测试有意义:
//   GR(2^2, 8): p^s = 4, base degree r = 8
//   lambda = 32 (降低安全参数以便测试)
//   packing_factor = ceil(32/8) = 4
//   ext_degree = 4 * 8 = 32
//
//   多项式变量数 = 12, 系数数 = 2^12 = 4096
//   num_rows = 2^6 = 64, row_len = 2^6 = 64
//   packed_row_len = ceil(64/4) = 16
//   这保证 packed_row_len > n_0=15, 会触发递归编码
// ============================================================
void test_small_ring_pcs() {
    cout << "====== Small Ring PCS Test ======" << endl;

    long s = 2;       // p^s = 4
    long base_r = 8;  // small ring degree
    long lambda = 32; // 降低安全参数以便测试

    // 初始化 base ring GR(4, 8)
    init_ring(s, base_r);

    cout << "Base ring: GR(2^" << s << ", " << base_r << ")" << endl;
    cout << "isLargeRing = " << isLargeRing(base_r, lambda) << endl;

    // 多项式: 12 个变量, 2^12 = 4096 个系数
    long num_vars = 12;
    long n = 1L << num_vars;        // 4096
    long half = num_vars / 2;       // 6
    long num_rows = 1L << half;     // 64
    long row_len  = 1L << (num_vars - half); // 64

    cout << "n = " << n << ", num_rows = " << num_rows
         << ", row_len = " << row_len << endl;

    // 随机多项式系数 (在 base ring 上)
    vector<ZZ_pE> poly(n);
    for (long i = 0; i < n; i++) {
        poly[i] = random_ZZ_pE();
    }

    // 用 auto setup (注意: 这会生成稀疏矩阵，需要在扩展环 context 下)
    // brakedown_code_setup_auto 内部会调用 createRandomSparseMatrix,
    // 需要 NTL context 设置正确。对 small ring, 矩阵在 ext ring 上操作。
    // 所以先切换到 ext ring 来做 setup:
    long pf = smallRingPackingFactor(base_r, lambda);
    long ext_r = pf * base_r;
    cout << "packing_factor = " << pf << ", ext_degree = " << ext_r << endl;

    // 切换到 ext ring 做 code setup (稀疏矩阵需要在 ext ring context 下生成)
    init_ring(s, ext_r);

    BrakedownCodeGR code = brakedown_code_setup_auto(row_len, s, base_r, lambda);
    cout << "is_small_ring = " << code.is_small_ring << endl;
    cout << "packed_row_len = " << code.packed_row_len << endl;
    cout << "codeword_len = " << code.codeword_len << endl;
    cout << "num_col_open = " << code.num_col_open << endl;
    cout << "num_prox_test = " << code.num_prox_test << endl;

    // 切回 base ring 来做 commit (内部会自行切换)
    init_ring(s, base_r);

    // Commit
    auto comm = brakedown_commit(code, poly.data(), n, num_rows);
    cout << "Commit done. Root size: " << comm.root.size() << " bytes" << endl;

    // 构造评估点 (在 base ring 上)
    init_ring(s, base_r);

    vector<ZZ_pE> eval_point(num_vars);
    for (long i = 0; i < num_vars; i++) {
        eval_point[i] = randomNonZeroInExceptionalSet();
    }

    // 构造 q1 (tensor product over first half variables)
    vector<ZZ_pE> q1(num_rows);
    {
        ZZ_pE one;
        set(one); // one = 1
        q1[0] = one;
        for (long i = 0; i < half; i++) {
            long cur_len = 1L << i;
            ZZ_pE one_minus_r = one - eval_point[i];
            for (long j = cur_len - 1; j >= 0; j--) {
                q1[2 * j + 1] = q1[j] * eval_point[i];
                q1[2 * j]     = q1[j] * one_minus_r;
            }
        }
    }

    // 构造 q2 (tensor product over second half variables)
    vector<ZZ_pE> q2(row_len);
    {
        ZZ_pE one;
        set(one);
        q2[0] = one;
        for (long i = 0; i < (num_vars - half); i++) {
            long cur_len = 1L << i;
            ZZ_pE one_minus_r = one - eval_point[half + i];
            for (long j = cur_len - 1; j >= 0; j--) {
                q2[2 * j + 1] = q2[j] * eval_point[half + i];
                q2[2 * j]     = q2[j] * one_minus_r;
            }
        }
    }

    // 直接计算期望值: f(r) = q1^T * S * q2
    ZZ_pE expected_eval;
    clear(expected_eval);
    for (long i = 0; i < num_rows; i++) {
        ZZ_pE row_dot;
        clear(row_dot);
        for (long j = 0; j < row_len; j++) {
            row_dot += poly[i * row_len + j] * q2[j];
        }
        expected_eval += q1[i] * row_dot;
    }
    cout << "Expected eval computed." << endl;

    // Prove
    auto proof = brakedown_prove(code, comm, poly.data(), n, q1, q2);
    cout << "Prove done." << endl;

    // 切回 base ring 比较 eval value
    init_ring(s, base_r);
    if (proof.eval_value == expected_eval) {
        cout << "[PASS] Evaluation value matches!" << endl;
    } else {
        cout << "[FAIL] Evaluation value mismatch!" << endl;
    }

    // Verify
    bool ok = brakedown_verify(code, comm, proof, q1, q2);
    cout << "Verify result: " << (ok ? "PASS" : "FAIL") << endl;
    cout << "=================================" << endl;
}

int main() {
    test_small_ring_pcs();
    return 0;
}
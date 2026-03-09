#include <iostream>
#include <NTL/ZZ_pE.h>
#include <NTL/ZZ_pX.h>
#include <gr.h>
#include <brakedown_code_gr.h>
#include <brakedown_pcs_gr.h>

using namespace NTL;
using namespace std;

// ============================================================
// 辅助: 只切换 ZZ_pE context (不动 ZZ_p！)
// ZZ_p::init() 在整个程序中只调一次！
// ============================================================
static void switch_ring(long degree) {
    ZZ_pX mod = primitiveIrredPoly(degree);
    ZZ_pE::init(mod);
}

void test_small_ring_pcs() {
    cout << "====== Small Ring PCS Test ======" << endl;

    long s = 2;       // p^s = 4
    long base_r = 8;  // small ring degree
    long lambda = 32; // 降低安全参数以便测试

    // ★ ZZ_p::init 只调一次！
    ZZ_p::init(ZZ(1) << s);

    // 初始化 base ring GR(4, 8)
    switch_ring(base_r);

    cout << "Base ring: GR(2^" << s << ", " << base_r << ")" << endl;
    cout << "isLargeRing = " << isLargeRing(base_r, lambda) << endl;

    long num_vars = 12;
    long n = 1L << num_vars;
    long half = num_vars / 2;
    long num_rows = 1L << half;
    long row_len  = 1L << (num_vars - half);

    cout << "n = " << n << ", num_rows = " << num_rows
         << ", row_len = " << row_len << endl;

    // 随机多项式系数 (在 base ring 上)
    vector<ZZ_pE> poly(n);
    for (long i = 0; i < n; i++) {
        poly[i] = random_ZZ_pE();
    }

    long pf = smallRingPackingFactor(base_r, lambda);
    long ext_r = pf * base_r;
    cout << "packing_factor = " << pf << ", ext_degree = " << ext_r << endl;

    // 切换到 ext ring 做 code setup (只切 ZZ_pE, 不动 ZZ_p)
    switch_ring(ext_r);

    BrakedownCodeGR code = brakedown_code_setup_auto(row_len, s, base_r, lambda);
    cout << "is_small_ring = " << code.is_small_ring << endl;
    cout << "packed_row_len = " << code.packed_row_len << endl;
    cout << "codeword_len = " << code.codeword_len << endl;
    cout << "num_col_open = " << code.num_col_open << endl;
    cout << "num_prox_test = " << code.num_prox_test << endl;

    // 切回 base ring (只切 ZZ_pE)
    switch_ring(base_r);

    // Commit
    auto comm = brakedown_commit(code, poly.data(), n, num_rows);
    cout << "Commit done. Root size: " << comm.root.size() << " bytes" << endl;

    // 构造评估点 (确保在 base ring)
    switch_ring(base_r);

    vector<ZZ_pE> eval_point(num_vars);
    for (long i = 0; i < num_vars; i++) {
        eval_point[i] = randomNonZeroInExceptionalSet();
    }

    // 构造 q1
    vector<ZZ_pE> q1(num_rows);
    {
        ZZ_pE one;
        set(one);
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

    // 构造 q2
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

    // 直接计算期望值
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

    // ★ 不要调 init_ring! 只切 ZZ_pE
    switch_ring(base_r);
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
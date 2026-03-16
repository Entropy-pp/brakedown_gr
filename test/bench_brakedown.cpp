#include <iostream>
#include <iomanip>
#include <chrono>
#include <vector>
#include <string>
#include <NTL/ZZ_pE.h>
#include <NTL/ZZ_pX.h>
#include <gr.h>
#include <brakedown_code_gr.h>
#include <brakedown_pcs_gr.h>

using namespace NTL;
using namespace std;

// ============================================================
// 计时工具
// ============================================================
struct Timer {
    using Clock = chrono::high_resolution_clock;
    Clock::time_point t0;
    void start() { t0 = Clock::now(); }
    double elapsed_ms() const {
        return chrono::duration<double, milli>(Clock::now() - t0).count();
    }
};

// 只切 ZZ_pE, 不动 ZZ_p
static void switch_ring(long degree) {
    ZZ_pX mod = primitiveIrredPoly(degree);
    ZZ_pE::init(mod);
}

// ============================================================
// 构造 multilinear evaluation tensor product
// ============================================================
static vector<ZZ_pE> build_tensor(const vector<ZZ_pE>& points, long start, long count) {
    long len = 1L << count;
    vector<ZZ_pE> q(len);
    ZZ_pE one; set(one);
    q[0] = one;
    for (long i = 0; i < count; i++) {
        long cur = 1L << i;
        ZZ_pE one_minus_r = one - points[start + i];
        for (long j = cur - 1; j >= 0; j--) {
            q[2 * j + 1] = q[j] * points[start + i];
            q[2 * j]     = q[j] * one_minus_r;
        }
    }
    return q;
}

// ============================================================
// 计算 expected eval (和 Prove/Verify 一致的方式)
// ============================================================
static ZZ_pE compute_expected_eval(
    const BrakedownCodeGR& code,
    const vector<ZZ_pE>& poly,
    const vector<ZZ_pE>& q1,
    const vector<ZZ_pE>& q2,
    long num_rows, long row_len,
    long base_r, long ext_r)
{
    if (!code.is_small_ring) {
        // Large ring: 直接 base ring 内积
        ZZ_pE result; clear(result);
        for (long i = 0; i < num_rows; i++) {
            ZZ_pE row_dot; clear(row_dot);
            for (long j = 0; j < row_len; j++)
                row_dot += poly[i * row_len + j] * q2[j];
            result += q1[i] * row_dot;
        }
        return result;
    }

    // Small ring: pack → ext ring combine → unpack → base ring 内积
    long pf = code.packing_factor;
    long packed_len = code.packed_row_len;

    switch_ring(base_r);
    vector<vector<ZZ_pX>> all_row_polys(num_rows);
    for (long i = 0; i < num_rows; i++) {
        all_row_polys[i].resize(row_len);
        for (long j = 0; j < row_len; j++)
            all_row_polys[i][j] = rep(poly[i * row_len + j]);
    }
    vector<ZZ_pX> q1_polys(num_rows);
    for (long i = 0; i < num_rows; i++)
        q1_polys[i] = rep(q1[i]);

    switch_ring(ext_r);
    vector<ZZ_pE> combined(packed_len);
    for (long j = 0; j < packed_len; j++) clear(combined[j]);
    for (long i = 0; i < num_rows; i++) {
        ZZ_pE c_ext = to_ZZ_pE(q1_polys[i]);
        for (long j = 0; j < packed_len; j++) {
            vector<ZZ_pX> group(pf);
            for (long t = 0; t < pf; t++) {
                long src = j * pf + t;
                group[t] = (src < row_len) ? all_row_polys[i][src] : ZZ_pX();
            }
            ZZ_pE packed_elem = pack_elements_pub(group, pf, base_r);
            combined[j] += c_ext * packed_elem;
        }
    }

    vector<ZZ_pX> unpacked;
    for (long j = 0; j < packed_len; j++) {
        auto parts = unpack_element_pub(combined[j], pf, base_r);
        for (long t = 0; t < pf; t++)
            unpacked.push_back(parts[t]);
    }

    switch_ring(base_r);
    ZZ_pE result; clear(result);
    for (long j = 0; j < row_len; j++)
        result += to_ZZ_pE(unpacked[j]) * q2[j];
    return result;
}

// ============================================================
// 单次 Benchmark
// ============================================================
struct BenchResult {
    string label;
    long num_vars;
    long n;
    long num_rows;
    long row_len;
    bool is_small_ring;
    long base_degree;
    long ext_degree;
    long packed_row_len;
    long codeword_len;
    double setup_ms;
    double commit_ms;
    double prove_ms;
    double verify_ms;
    bool correct;
    // sizes
    long proof_combined_rows;
    long proof_column_opens;
};

static BenchResult run_benchmark(
    const string& label,
    long s,         // p^s exponent (ZZ_p modulus = 2^s)
    long base_r,    // base ring degree
    long lambda,    // security parameter
    long num_vars)  // number of variables → n = 2^num_vars
{
    BenchResult res;
    res.label = label;
    res.num_vars = num_vars;
    res.n = 1L << num_vars;
    long half = num_vars / 2;
    res.num_rows = 1L << half;
    res.row_len = 1L << (num_vars - half);

    Timer timer;
    long n = res.n;
    long num_rows = res.num_rows;
    long row_len = res.row_len;

    bool is_small = !isLargeRing(base_r, lambda);
    long pf = smallRingPackingFactor(base_r, lambda);
    long ext_r = pf * base_r;

    res.is_small_ring = is_small;
    res.base_degree = base_r;
    res.ext_degree = ext_r;

    // === Setup ===
    timer.start();
    if (is_small) {
        switch_ring(ext_r);
    } else {
        switch_ring(base_r);
    }
    BrakedownCodeGR code = brakedown_code_setup_auto(row_len, s, base_r, lambda);
    res.setup_ms = timer.elapsed_ms();

    res.packed_row_len = code.packed_row_len;
    res.codeword_len = code.codeword_len;

    // 生成随机多项式 (在 base ring)
    switch_ring(base_r);
    vector<ZZ_pE> poly(n);
    for (long i = 0; i < n; i++)
        poly[i] = random_ZZ_pE();

    // 生成评估点和 tensor 向量
    vector<ZZ_pE> eval_point(num_vars);
    for (long i = 0; i < num_vars; i++)
        eval_point[i] = randomNonZeroInExceptionalSet();

    vector<ZZ_pE> q1 = build_tensor(eval_point, 0, half);
    vector<ZZ_pE> q2 = build_tensor(eval_point, half, num_vars - half);

    // 期望值
    ZZ_pE expected_eval = compute_expected_eval(
        code, poly, q1, q2, num_rows, row_len, base_r, ext_r);

    // === Commit ===
    switch_ring(base_r);
    timer.start();
    CommitTiming timing;
    auto comm = brakedown_commit(code, poly.data(), n, num_rows, timing);
    res.commit_ms = timer.elapsed_ms();

    std::cout << "  Commit breakdown:" << std::endl;
    std::cout << "    Fill rows:     " << timing.fill_rows_ms    << " ms" << std::endl;
    if (code.is_small_ring) {
        std::cout << "    Pack rows:     " << timing.pack_rows_ms << " ms" << std::endl;
    }
    std::cout << "    Encode rows:   " << timing.encode_rows_ms  << " ms" << std::endl;
    std::cout << "    Column hash:   " << timing.column_hash_ms  << " ms" << std::endl;
    std::cout << "    Merkle build:  " << timing.merkle_build_ms << " ms" << std::endl;
    std::cout << "    Total:         " << timing.total_ms        << " ms" << std::endl;

    // === Prove ===
    switch_ring(base_r);
    timer.start();
    auto proof = brakedown_prove(code, comm, poly.data(), n, q1, q2);
    res.prove_ms = timer.elapsed_ms();

    res.proof_combined_rows = (long)proof.combined_rows.size();
    res.proof_column_opens = (long)proof.column_indices.size();

    // === Verify ===
    switch_ring(base_r);
    timer.start();
    bool ok = brakedown_verify(code, comm, proof, q1, q2);
    res.verify_ms = timer.elapsed_ms();

    // 检查 eval value
    switch_ring(base_r);
    bool eval_ok = (proof.eval_value == expected_eval);
    res.correct = ok && eval_ok;

    return res;
}

// ============================================================
// 打印结果
// ============================================================
static void print_header() {
    cout << string(120, '=') << endl;
    cout << left
         << setw(28) << "Configuration"
         << setw(10) << "n"
         << setw(8)  << "rows"
         << setw(8)  << "rowlen"
         << setw(6)  << "cwlen"
         << setw(12) << "Setup(ms)"
         << setw(12) << "Commit(ms)"
         << setw(12) << "Prove(ms)"
         << setw(12) << "Verify(ms)"
         << setw(8)  << "Status"
         << endl;
    cout << string(120, '-') << endl;
}

static void print_result(const BenchResult& r) {
    cout << left
         << setw(28) << r.label
         << setw(10) << r.n
         << setw(8)  << r.num_rows
         << setw(8)  << r.row_len
         << setw(6)  << r.codeword_len
         << setw(12) << fixed << setprecision(1) << r.setup_ms
         << setw(12) << fixed << setprecision(1) << r.commit_ms
         << setw(12) << fixed << setprecision(1) << r.prove_ms
         << setw(12) << fixed << setprecision(1) << r.verify_ms
         << setw(8)  << (r.correct ? "PASS" : "FAIL")
         << endl;
}

static void print_detail(const BenchResult& r) {
    cout << "  Ring: GR(2^" << 2 << ", " << r.base_degree << ")"
         << (r.is_small_ring ? " [Small Ring]" : " [Large Ring]")
         << "  ext_degree=" << r.ext_degree
         << "  packed_row=" << r.packed_row_len
         << "  proof_rows=" << r.proof_combined_rows
         << "  col_opens=" << r.proof_column_opens
         << endl;
}

// ============================================================
// Main
// ============================================================
int main() {
    long s = 1; // p^s = 4, 即 GR(2^2, d)

    // ★ ZZ_p::init 全程只调一次
    ZZ_p::init(ZZ(1) << s);

    cout << endl;
    cout << "╔══════════════════════════════════════════════════════════════╗" << endl;
    cout << "║    Brakedown PCS over Galois Ring — Full Benchmark           ║" << endl;
    cout << "║    Base modulus: p^s = 2^" << s << " = " << (1L << s) << setw(34) << "║" << endl;
    cout << "╚══════════════════════════════════════════════════════════════╝" << endl;
    cout << endl;

    vector<BenchResult> results;

    // ============================================================
    // Part 1: Large Ring 测试
    //   degree >= lambda → isLargeRing = true, packing_factor = 1
    //   使用 degree = 128, lambda = 128
    // ============================================================
    cout << ">>> Part 1: Large Ring  GR(2^2, 128), lambda=128" << endl;
    print_header();
    {
        long base_r = 162;
        long lambda = 128;
        for (long nv : {13, 14, 15, 16}) {
            string label = "Large d=128 nv=" + to_string(nv);
            auto r = run_benchmark(label, s, base_r, lambda, nv);
            print_result(r);
            print_detail(r);
            results.push_back(r);
        }
    }
    cout << endl;

    // // ============================================================
    // // Part 2: Small Ring 测试 (不同 base degree)
    // //   degree < lambda → isLargeRing = false
    // //   lambda = 32 (测试用低安全参数, 加速)
    // // ============================================================
    // cout << ">>> Part 2: Small Ring  lambda=32, 不同 base degree" << endl;
    // print_header();
    // {
    //     long lambda = 32;
    //     for (long base_r : {8, 16}) {
    //         for (long nv : {10, 12, 14}) {
    //             long pf = smallRingPackingFactor(base_r, lambda);
    //             string label = "Small d=" + to_string(base_r)
    //                 + " k=" + to_string(pf) + " nv=" + to_string(nv);
    //             auto r = run_benchmark(label, s, base_r, lambda, nv);
    //             print_result(r);
    //             print_detail(r);
    //             results.push_back(r);
    //         }
    //     }
    // }
    // cout << endl;

    // ============================================================
    // Part 3: Small Ring 高安全参数
    //   lambda = 128, base_r = 8 → pf = 16, ext = 128
    // ============================================================
    cout << ">>> Part 2: Small Ring  GR(2^2, 8), lambda=128 (高安全参数)" << endl;
    print_header();
    {
        long base_r = 8;
        long lambda = 128;
        for (long nv : {10, 12, 14}) {
            long pf = smallRingPackingFactor(base_r, lambda);
            string label = "Small d=8 k=" + to_string(pf) + " L128 nv=" + to_string(nv);
            auto r = run_benchmark(label, s, base_r, lambda, nv);
            print_result(r);
            print_detail(r);
            results.push_back(r);
        }
    }
    cout << endl;

    // // ============================================================
    // // Part 4: 固定 num_vars, 对比 Large vs Small
    // // ============================================================
    // cout << ">>> Part 4: Large vs Small 对比 (num_vars=14)" << endl;
    // print_header();
    // {
    //     long nv = 14;
    //     // Large: d=128, lambda=128
    //     {
    //         auto r = run_benchmark("Large d=128 L=128", s, 128, 128, nv);
    //         print_result(r);
    //         print_detail(r);
    //         results.push_back(r);
    //     }
    //     // Small: d=8, lambda=32
    //     {
    //         auto r = run_benchmark("Small d=8 k=4 L=32", s, 8, 32, nv);
    //         print_result(r);
    //         print_detail(r);
    //         results.push_back(r);
    //     }
    //     // Small: d=16, lambda=32
    //     {
    //         auto r = run_benchmark("Small d=16 k=2 L=32", s, 16, 32, nv);
    //         print_result(r);
    //         print_detail(r);
    //         results.push_back(r);
    //     }
    //     // Small: d=8, lambda=128
    //     {
    //         auto r = run_benchmark("Small d=8 k=16 L=128", s, 8, 128, nv);
    //         print_result(r);
    //         print_detail(r);
    //         results.push_back(r);
    //     }
    // }
    // cout << endl;

    // ============================================================
    // 汇总
    // ============================================================
    cout << "╔══════════════════════════════════════════════════════════════╗" << endl;
    cout << "║                         Summary                              ║" << endl;
    cout << "╚══════════════════════════════════════════════════════════════╝" << endl;
    int pass = 0, fail = 0;
    for (auto& r : results) {
        if (r.correct) pass++; else fail++;
    }
    cout << "Total: " << results.size()
         << "  PASS: " << pass
         << "  FAIL: " << fail << endl;

    if (fail > 0) {
        cout << "FAILED cases:" << endl;
        for (auto& r : results) {
            if (!r.correct) cout << "  - " << r.label << endl;
        }
    }
    cout << endl;

    return 0;
}
/**
 * bench_pcs.cpp — Brakedown PCS over Galois Ring 综合性能 Benchmark
 *
 * 编译:
 *   g++ -O2 -std=c++11 -I./include test/bench_pcs.cpp \
 *       src/gr.cpp src/sparse_matrix_gr.cpp src/brakedown_code_gr.cpp \
 *       src/merkle.cpp src/brakedown_pcs_gr.cpp \
 *       -o out/bench_pcs -pthread -lntl -lgmp -lm
 *
 * 运行:
 *   ./out/bench_pcs
 */

#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pE.h>
#include <iostream>
#include <iomanip>
#include <vector>
#include <ctime>
#include <cmath>
#include <cassert>
#include <chrono>
#include <string>
#include <sstream>
#include <numeric>

#include <gr.h>
#include <brakedown_params.h>
#include <brakedown_code_gr.h>
#include <brakedown_pcs_gr.h>

using namespace std;
using namespace NTL;

// ============================================================
// 高精度计时器
// ============================================================
using hclock = chrono::high_resolution_clock;
using time_point = chrono::time_point<hclock>;

static inline time_point now() { return hclock::now(); }
static inline double elapsed_ms(time_point start, time_point end) {
    return chrono::duration<double, milli>(end - start).count();
}

// ============================================================
// 初始化 GR(2^k, d)
// ============================================================
static void initGR(long k, long degree) {
    ZZ modulus = ZZ(1) << k;
    ZZ_p::init(modulus);
    ZZ_pX P = primitiveIrredPoly(degree);
    ZZ_pE::init(P);
}

// ============================================================
// 选择合理的 row_len: 取 sqrt(n) 附近的 2 的幂
// ============================================================
static long choose_row_len(long n) {
    long target = (long)std::ceil(std::sqrt((double)n));
    long row_len = 1;
    while (row_len * 2 <= target) row_len *= 2;
    if (row_len < 2) row_len = 2;
    while (n / row_len < 1 && row_len > 2) row_len /= 2;
    return row_len;
}

// ============================================================
// 人类可读的字节大小
// ============================================================
static string human_bytes(double bytes) {
    ostringstream oss;
    if (bytes < 1024.0) {
        oss << fixed << setprecision(0) << bytes << " B";
    } else if (bytes < 1024.0 * 1024.0) {
        oss << fixed << setprecision(1) << bytes / 1024.0 << " KB";
    } else {
        oss << fixed << setprecision(2) << bytes / (1024.0 * 1024.0) << " MB";
    }
    return oss.str();
}

// ============================================================
// 估算 proof 大小 (字节)
// ============================================================
static double estimate_proof_size(const BrakedownCodeGR& code,
                                   long num_rows, long k, long degree) {
    double elem_bytes = (double)k * degree / 8.0;
    long num_prox = (num_rows > 1) ? code.num_prox_test : 0;
    long row_len = code.row_len;
    long num_open = code.num_col_open;
    long depth = (long)std::ceil(std::log2((double)code.codeword_len));

    double size = 0;
    size += (double)(num_prox + 1) * row_len * elem_bytes;   // combined_rows
    size += (double)num_prox * num_rows * elem_bytes;         // prox_coeffs
    size += (double)num_open * num_rows * elem_bytes;         // column_items
    size += (double)num_open * 8;                             // column_indices
    size += (double)num_open * depth * 32;                    // merkle_paths
    size += elem_bytes;                                       // eval_value
    return size;
}

// ============================================================
// 单次 benchmark 结果
// ============================================================
struct BenchResult {
    long n;
    long k;
    long degree;
    long row_len;
    long num_rows;
    long codeword_len;
    long num_col_open;
    long num_prox_test;
    bool is_small_ring;
    long packing_factor;
    double setup_ms;
    double commit_ms;
    double prove_ms;
    double verify_ms;
    double total_ms;
    double proof_bytes;
    bool verify_ok;
};

// ============================================================
// 执行单次 benchmark
//
// ★ 关键修复: Small Ring 路径需要在扩展环上下文中运行
//   brakedown_code_setup_auto 调用 createRandomSparseMatrix,
//   后者依赖 random_ZZ_pE() 生成扩展环上的元素。
//   如果 NTL context 仍在 base ring, 生成的元素 degree 不匹配
//   → 后续 encode/commit 中操作不一致 → 段错误。
//
//   因此: 当 use_auto=true 且 degree < lambda 时,
//   我们必须先切换到 ext_deg = packing_factor * degree 的扩展环。
// ============================================================
BenchResult run_benchmark(long n, long k, long degree, long lambda, bool use_auto) {
    BenchResult res = {};
    res.n = n;
    res.k = k;
    res.degree = degree;

    // --- 计算实际使用的环度数 ---
    long actual_degree = degree;
    if (use_auto && !isLargeRing(degree, lambda)) {
        // Small ring: 需要在扩展环上操作
        long pf = smallRingPackingFactor(degree, lambda);
        actual_degree = pf * degree;
    }

    // --- 初始化到实际使用的环 ---
    initGR(k, actual_degree);

    // --- 生成随机多项式 (在实际使用的环上) ---
    vector<ZZ_pE> poly(n);
    for (long i = 0; i < n; i++) poly[i] = random_ZZ_pE();

    long row_len = choose_row_len(n);
    long num_rows = n / row_len;
    if (num_rows < 1) { num_rows = 1; row_len = n; }
    long n_aligned = num_rows * row_len;

    res.row_len = row_len;
    res.num_rows = num_rows;

    // --- Setup ---
    time_point t0 = now();
    BrakedownCodeGR code;
    if (use_auto) {
        code = brakedown_code_setup_auto(row_len, k, degree, lambda);
    } else {
        code = brakedown_code_setup(row_len, k, degree);
    }
    time_point t1 = now();
    res.setup_ms = elapsed_ms(t0, t1);
    res.codeword_len = code.codeword_len;
    res.num_col_open = code.num_col_open;
    res.num_prox_test = code.num_prox_test;
    res.is_small_ring = code.is_small_ring;
    res.packing_factor = code.packing_factor;

    // --- Commit ---
    time_point t2 = now();
    BrakedownCommitment comm = brakedown_commit(code, poly.data(), n_aligned, num_rows);
    time_point t3 = now();
    res.commit_ms = elapsed_ms(t2, t3);

    // --- 生成 evaluation query ---
    vector<ZZ_pE> q1(num_rows), q2(row_len);
    if (num_rows == 1) {
        set(q1[0]);
    } else {
        for (long i = 0; i < num_rows; i++) q1[i] = random_ZZ_pE();
    }
    for (long j = 0; j < row_len; j++) q2[j] = random_ZZ_pE();

    // --- Prove ---
    time_point t4 = now();
    BrakedownEvalProof proof = brakedown_prove(code, comm, poly.data(), n_aligned, q1, q2);
    time_point t5 = now();
    res.prove_ms = elapsed_ms(t4, t5);

    // --- 正确性检查 ---
    ZZ_pE expected_eval;
    clear(expected_eval);
    for (long i = 0; i < num_rows; i++) {
        for (long j = 0; j < row_len; j++) {
            long idx = i * row_len + j;
            if (idx < n_aligned)
                expected_eval += q1[i] * poly[idx] * q2[j];
        }
    }
    if (expected_eval != proof.eval_value) {
        cerr << "  [ERROR] Evaluation mismatch for n=" << n
             << " k=" << k << " d=" << degree << endl;
    }

    // --- Verify ---
    time_point t6 = now();
    res.verify_ok = brakedown_verify(code, comm, proof, q1, q2);
    time_point t7 = now();
    res.verify_ms = elapsed_ms(t6, t7);

    res.total_ms = res.setup_ms + res.commit_ms + res.prove_ms + res.verify_ms;

    // 对于 small ring, proof size 仍按 actual_degree 计算 (因为元素是扩展环上的)
    long degree_for_size = code.is_small_ring ? actual_degree : degree;
    res.proof_bytes = estimate_proof_size(code, num_rows, k, degree_for_size);

    return res;
}

// ============================================================
// 打印分隔线
// ============================================================
static void print_sep(int width = 148) {
    cout << string(width, '-') << "\n";
}

// ============================================================
// Benchmark 1: 固定环参数, 变多项式规模
// ============================================================
void bench_scaling(long k, long degree) {
    cout << "\n";
    cout << "================================================================"
            "===========================================================\n";
    cout << " Benchmark 1: Scaling with polynomial size n\n";
    cout << " Ring: GR(2^" << k << ", " << degree << ")  |  "
         << "log2|R| = " << k * degree << "  |  "
         << "Element size = " << k * degree / 8 << " bytes\n";
    cout << "================================================================"
            "===========================================================\n\n";

    cout << left
         << setw(10) << "n"
         << setw(10) << "row_len"
         << setw(10) << "rows"
         << setw(10) << "cw_len"
         << setw(12) << "Setup(ms)"
         << setw(12) << "Commit(ms)"
         << setw(12) << "Prove(ms)"
         << setw(12) << "Verify(ms)"
         << setw(12) << "Total(ms)"
         << setw(12) << "Proof"
         << setw(8)  << "OK?"
         << "\n";
    print_sep();

    long sizes[] = {256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536, 131072};

    for (int i = 0; i < 10; i++) {
        long n = sizes[i];
        BenchResult r = run_benchmark(n, k, degree, 128, false);

        cout << left << fixed
             << setw(10) << r.n
             << setw(10) << r.row_len
             << setw(10) << r.num_rows
             << setw(10) << r.codeword_len
             << setw(12) << setprecision(2) << r.setup_ms
             << setw(12) << setprecision(2) << r.commit_ms
             << setw(12) << setprecision(2) << r.prove_ms
             << setw(12) << setprecision(2) << r.verify_ms
             << setw(12) << setprecision(2) << r.total_ms
             << setw(12) << human_bytes(r.proof_bytes)
             << setw(8)  << (r.verify_ok ? "PASS" : "FAIL")
             << "\n";
    }
    print_sep();
}

// ============================================================
// Benchmark 2: 固定多项式规模, 变环度数 d
// ============================================================
void bench_varying_degree(long k, long n) {
    cout << "\n";
    cout << "================================================================"
            "=======================================================================\n";
    cout << " Benchmark 2: Varying extension degree d (fixed n=" << n << ", k=" << k << ")\n";
    cout << "================================================================"
            "=======================================================================\n\n";

    cout << left
         << setw(6)  << "d"
         << setw(10) << "log2|R|"
         << setw(10) << "row_len"
         << setw(10) << "rows"
         << setw(10) << "cw_len"
         << setw(10) << "#colOpen"
         << setw(10) << "#proxTst"
         << setw(12) << "Setup(ms)"
         << setw(12) << "Commit(ms)"
         << setw(12) << "Prove(ms)"
         << setw(12) << "Verify(ms)"
         << setw(12) << "Total(ms)"
         << setw(12) << "Proof"
         << setw(6)  << "OK?"
         << "\n";
    print_sep(154);

    long degrees[] = {4, 8, 16, 32, 40, 60, 80};

    for (int i = 0; i < 7; i++) {
        long degree = degrees[i];
        BenchResult r = run_benchmark(n, k, degree, 128, false);

        cout << left << fixed
             << setw(6)  << r.degree
             << setw(10) << (k * r.degree)
             << setw(10) << r.row_len
             << setw(10) << r.num_rows
             << setw(10) << r.codeword_len
             << setw(10) << r.num_col_open
             << setw(10) << r.num_prox_test
             << setw(12) << setprecision(2) << r.setup_ms
             << setw(12) << setprecision(2) << r.commit_ms
             << setw(12) << setprecision(2) << r.prove_ms
             << setw(12) << setprecision(2) << r.verify_ms
             << setw(12) << setprecision(2) << r.total_ms
             << setw(12) << human_bytes(r.proof_bytes)
             << setw(6)  << (r.verify_ok ? "PASS" : "FAIL")
             << "\n";
    }
    print_sep(154);
}

// ============================================================
// Benchmark 3: Large Ring vs Small Ring
// ★ 这是修复 segfault 的关键 benchmark
// ============================================================
void bench_large_vs_small(long k) {
    cout << "\n";
    cout << "================================================================"
            "==========================================================================\n";
    cout << " Benchmark 3: Large Ring vs Small Ring (auto mode, k=" << k << ")\n";
    cout << " Large Ring: d >= lambda (直接编码)\n";
    cout << " Small Ring: d < lambda  (需要 packing, 在 GR(2^k, pf*d) 上编码)\n";
    cout << "================================================================"
            "==========================================================================\n\n";

    // --- λ = 128 ---
    cout << " ──── Security Parameter lambda = 128 ────\n\n";

    cout << left
         << setw(6)  << "d"
         << setw(7)  << "type"
         << setw(6)  << "pack"
         << setw(9)  << "ext_d"
         << setw(8)  << "n"
         << setw(10) << "row_len"
         << setw(10) << "rows"
         << setw(12) << "Setup(ms)"
         << setw(12) << "Commit(ms)"
         << setw(12) << "Prove(ms)"
         << setw(12) << "Verify(ms)"
         << setw(12) << "Total(ms)"
         << setw(12) << "Proof"
         << setw(6)  << "OK?"
         << "\n";
    print_sep(144);

    struct Config { long degree; long n; };

    // λ=128: 只有 d >= 128 是 large ring
    // Small ring 的 packing 会使 ext_deg 很大, 运行时间可能很长
    // 因此减小 n 以保持可行性
    Config configs_128[] = {
        {32,  1024},   // small: pf=4, ext_d=128
        {64,  1024},   // small: pf=2, ext_d=128
        {128, 1024},   // large: pf=1, ext_d=128
    };

    for (int i = 0; i < 3; i++) {
        long degree = configs_128[i].degree;
        long n = configs_128[i].n;

        cout << "  Running d=" << degree << " n=" << n << " ... " << flush;
        BenchResult r = run_benchmark(n, k, degree, 128, true);

        string ring_type = r.is_small_ring ? "SMALL" : "LARGE";
        long ext_d = r.is_small_ring ? (r.packing_factor * degree) : degree;

        cout << "\r" << left << fixed
             << setw(6)  << degree
             << setw(7)  << ring_type
             << setw(6)  << r.packing_factor
             << setw(9)  << ext_d
             << setw(8)  << r.n
             << setw(10) << r.row_len
             << setw(10) << r.num_rows
             << setw(12) << setprecision(2) << r.setup_ms
             << setw(12) << setprecision(2) << r.commit_ms
             << setw(12) << setprecision(2) << r.prove_ms
             << setw(12) << setprecision(2) << r.verify_ms
             << setw(12) << setprecision(2) << r.total_ms
             << setw(12) << human_bytes(r.proof_bytes)
             << setw(6)  << (r.verify_ok ? "PASS" : "FAIL")
             << "\n";
    }
    print_sep(144);
}

// ============================================================
// Benchmark 4: 多次运行取平均
// ============================================================
void bench_repeated(long n, long k, long degree, int num_runs) {
    cout << "\n";
    cout << "================================================================\n";
    cout << " Benchmark 4: Statistical stability (" << num_runs << " runs)\n";
    cout << " n=" << n << "  GR(2^" << k << ", " << degree << ")\n";
    cout << "================================================================\n\n";

    vector<double> t_setup, t_commit, t_prove, t_verify;

    for (int run = 0; run < num_runs; run++) {
        BenchResult r = run_benchmark(n, k, degree, 128, false);
        t_setup.push_back(r.setup_ms);
        t_commit.push_back(r.commit_ms);
        t_prove.push_back(r.prove_ms);
        t_verify.push_back(r.verify_ms);
        if (!r.verify_ok) {
            cerr << "  [WARN] Run " << (run + 1) << " verification FAILED!\n";
        }
    }

    auto avg = [](const vector<double>& v) {
        return accumulate(v.begin(), v.end(), 0.0) / (double)v.size();
    };
    auto min_val = [](const vector<double>& v) {
        return *min_element(v.begin(), v.end());
    };
    auto max_val = [](const vector<double>& v) {
        return *max_element(v.begin(), v.end());
    };
    auto stddev_fn = [&avg](const vector<double>& v) {
        double m = avg(v);
        double sq_sum = 0;
        for (double x : v) sq_sum += (x - m) * (x - m);
        return std::sqrt(sq_sum / (double)v.size());
    };

    cout << left
         << setw(12) << "Phase"
         << setw(14) << "Avg (ms)"
         << setw(14) << "Min (ms)"
         << setw(14) << "Max (ms)"
         << setw(14) << "Stddev (ms)"
         << "\n";
    print_sep(68);

    auto print_row = [&](const string& name, const vector<double>& v) {
        cout << left << fixed << setprecision(2)
             << setw(12) << name
             << setw(14) << avg(v)
             << setw(14) << min_val(v)
             << setw(14) << max_val(v)
             << setw(14) << stddev_fn(v)
             << "\n";
    };

    print_row("Setup",  t_setup);
    print_row("Commit", t_commit);
    print_row("Prove",  t_prove);
    print_row("Verify", t_verify);

    vector<double> t_total((size_t)num_runs);
    for (int i = 0; i < num_runs; i++)
        t_total[i] = t_setup[i] + t_commit[i] + t_prove[i] + t_verify[i];
    print_row("Total",  t_total);
    print_sep(68);
}

// ============================================================
// Benchmark 5: Commit 吞吐率
// ============================================================
void bench_throughput(long k, long degree) {
    cout << "\n";
    cout << "================================================================\n";
    cout << " Benchmark 5: Commit throughput (coefficients/second)\n";
    cout << " Ring: GR(2^" << k << ", " << degree << ")\n";
    cout << "================================================================\n\n";

    cout << left
         << setw(10) << "n"
         << setw(14) << "Commit(ms)"
         << setw(18) << "Throughput(coeff/s)"
         << setw(18) << "Throughput(KB/s)"
         << "\n";
    print_sep(60);

    double elem_bytes = (double)k * degree / 8.0;
    long sizes[] = {256, 1024, 4096, 16384, 65536};

    for (int i = 0; i < 5; i++) {
        long n = sizes[i];
        BenchResult r = run_benchmark(n, k, degree, 128, false);

        double coeff_per_s = (r.commit_ms > 0.001) ? (double)n / (r.commit_ms / 1000.0) : 0;
        double kb_per_s = coeff_per_s * elem_bytes / 1024.0;

        cout << left << fixed
             << setw(10) << n
             << setw(14) << setprecision(2) << r.commit_ms
             << setw(18) << setprecision(0) << coeff_per_s
             << setw(18) << setprecision(1) << kb_per_s
             << "\n";
    }
    print_sep(60);
}

// ============================================================
// Benchmark 6: 参数摘要
// ============================================================
void bench_params_summary(long k) {
    cout << "\n";
    cout << "================================================================\n";
    cout << " Parameter Summary: BrakedownSpec from GLSTW21 Figure 2\n";
    cout << "================================================================\n\n";

    BrakedownSpec spec = getDefaultSpec();
    cout << "  LAMBDA = " << spec.LAMBDA << "\n"
         << "  ALPHA  = " << spec.ALPHA  << "\n"
         << "  BETA   = " << spec.BETA   << "\n"
         << "  R      = " << spec.R      << "\n"
         << "  delta  = " << spec.delta() << "\n"
         << "  mu     = " << spec.mu()    << "\n"
         << "  nu     = " << spec.nu()    << "\n\n";

    cout << left
         << setw(6)  << "d"
         << setw(10) << "log2|R|"
         << setw(12) << "elem_bytes"
         << setw(10) << "large128"
         << setw(10) << "large256"
         << setw(10) << "pack128"
         << setw(10) << "pack256"
         << "\n";
    print_sep(68);

    long degrees[] = {2, 4, 8, 16, 32, 64, 80, 128, 160, 256};
    for (int i = 0; i < 10; i++) {
        long d = degrees[i];
        long log2R = k * d;
        long elem_b = log2R / 8;
        bool large_128 = isLargeRing(d, 128);
        bool large_256 = isLargeRing(d, 256);
        long pack_128 = large_128 ? 1 : smallRingPackingFactor(d, 128);
        long pack_256 = large_256 ? 1 : smallRingPackingFactor(d, 256);

        cout << left
             << setw(6)  << d
             << setw(10) << log2R
             << setw(12) << elem_b
             << setw(10) << (large_128 ? "YES" : "no")
             << setw(10) << (large_256 ? "YES" : "no")
             << setw(10) << pack_128
             << setw(10) << pack_256
             << "\n";
    }
    print_sep(68);
}

// ============================================================
// MAIN
// ============================================================
int main() {
    cout << "##############################################################\n";
    cout << "#                                                            #\n";
    cout << "#  Brakedown PCS over Galois Ring — Performance Benchmark    #\n";
    cout << "#                                                            #\n";
    cout << "##############################################################\n";

    long k = 64;

    // --- 0. 参数总览 ---
    bench_params_summary(k);

    // --- 1. 扩展性: 固定 d=8, 变 n ---
    bench_scaling(k, 8);

    // --- 2. 度数影响: 固定 n=4096, 变 d ---
    bench_varying_degree(k, 4096);

    // --- 3. Large Ring vs Small Ring ---
    bench_large_vs_small(k);

    // --- 4. 统计稳定性: n=4096, d=8, 5 次 ---
    bench_repeated(4096, k, 8, 5);

    // --- 5. 吞吐率 ---
    bench_throughput(k, 8);

    cout << "\n##############################################################\n";
    cout << "#  Benchmark complete.                                       #\n";
    cout << "##############################################################\n";

    return 0;
}
/**
 * bench_distributed_v2.cpp - Benchmark Protocol 2: Distributed Merkle Tree
 *
 * Compile:
 *   g++ -O3 -std=c++17 -pthread -I./include \
 *       test/bench_distributed_v2.cpp \
 *       src/gr.cpp src/sparse_matrix_gr.cpp src/brakedown_code_gr.cpp \
 *       src/merkle.cpp src/brakedown_pcs_gr.cpp src/brakedown_distributed.cpp \
 *       -o out/bench_distributed_v2 -lntl -lgmp -lm -lssl -lcrypto
 *
 * Run:
 *   ./out/bench_distributed_v2 [num_workers]
 */

#include <iostream>
#include <iomanip>
#include <vector>
#include <chrono>
#include <thread>
#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pE.h>

#include <gr.h>
#include <brakedown_params.h>
#include <brakedown_code_gr.h>
#include <brakedown_pcs_gr.h>
#include <brakedown_distributed.h>

using namespace std;
using namespace NTL;
using hrc = chrono::high_resolution_clock;

static double ms_between(hrc::time_point start, hrc::time_point end) {
    return chrono::duration<double, milli>(end - start).count();
}

static void initGR(long k, long degree) {
    ZZ modulus = ZZ(1) << k;
    ZZ_p::init(modulus);
    ZZ_pX P = primitiveIrredPoly(degree);
    ZZ_pE::init(P);
}

// ============================================================
// Single-threaded benchmark (baseline)
// ============================================================
struct SingleResult {
    double commit_ms;
    double prove_ms;
    double verify_ms;
    bool verified;
};

SingleResult run_single(
    const BrakedownCodeGR& code,
    const ZZ_pE* poly,
    long n,
    long num_rows,
    const vector<ZZ_pE>& q1,
    const vector<ZZ_pE>& q2)
{
    SingleResult res;

    auto t1 = hrc::now();
    auto comm = brakedown_commit(code, poly, n, num_rows);
    auto t2 = hrc::now();
    res.commit_ms = ms_between(t1, t2);

    auto t3 = hrc::now();
    auto proof = brakedown_prove(code, comm, poly, n, q1, q2);
    auto t4 = hrc::now();
    res.prove_ms = ms_between(t3, t4);

    auto t5 = hrc::now();
    res.verified = brakedown_verify(code, comm, proof, q1, q2);
    auto t6 = hrc::now();
    res.verify_ms = ms_between(t5, t6);

    return res;
}

// ============================================================
// Protocol 2 (V2) benchmark
// ============================================================
struct V2Result {
    int num_workers;
    DistributedCommitTimingV2 commit_timing;
    DistributedProveTiming prove_timing;
    double verify_ms;
    bool verified;
    vector<WorkerStatsV2> commit_stats;
    vector<WorkerStats> prove_stats;
};

V2Result run_v2(
    const BrakedownCodeGR& code,
    const ZZ_pE* poly,
    long n,
    long num_rows,
    const vector<ZZ_pE>& q1,
    const vector<ZZ_pE>& q2,
    int num_workers)
{
    V2Result res;
    res.num_workers = num_workers;

    auto comm = distributed_commit_v2(code, poly, n, num_rows, num_workers,
                                      res.commit_timing, res.commit_stats);

    auto proof = distributed_prove_v2(code, comm, poly, n, q1, q2, num_workers,
                                      res.prove_timing, res.prove_stats);

    auto t1 = hrc::now();
    res.verified = distributed_verify_v2(code, comm, proof, q1, q2);
    auto t2 = hrc::now();
    res.verify_ms = ms_between(t1, t2);

    return res;
}

// ============================================================
// Print results
// ============================================================
void print_comparison(const SingleResult& single, const V2Result& v2) {
    cout << "\n";
    cout << "============================================================\n";
    cout << "  Single-threaded vs Protocol 2 (" << v2.num_workers << " workers)\n";
    cout << "============================================================\n";
    cout << "\n";

    cout << left << fixed << setprecision(2);
    cout << setw(20) << "Phase"
         << setw(15) << "Single(ms)"
         << setw(15) << "V2(ms)"
         << setw(12) << "Speedup"
         << "\n";
    cout << string(62, '-') << "\n";

    double commit_speedup = single.commit_ms / v2.commit_timing.total_ms;
    double prove_speedup = single.prove_ms / v2.prove_timing.total_ms;

    cout << setw(20) << "Commit"
         << setw(15) << single.commit_ms
         << setw(15) << v2.commit_timing.total_ms
         << setw(12) << (to_string(commit_speedup).substr(0,4) + "x")
         << "\n";

    cout << setw(20) << "Prove"
         << setw(15) << single.prove_ms
         << setw(15) << v2.prove_timing.total_ms
         << setw(12) << (to_string(prove_speedup).substr(0,4) + "x")
         << "\n";

    cout << setw(20) << "Verify"
         << setw(15) << single.verify_ms
         << setw(15) << v2.verify_ms
         << setw(12) << "-"
         << "\n";

    cout << string(62, '-') << "\n";

    double total_single = single.commit_ms + single.prove_ms + single.verify_ms;
    double total_v2 = v2.commit_timing.total_ms + v2.prove_timing.total_ms + v2.verify_ms;
    double total_speedup = total_single / total_v2;

    cout << setw(20) << "TOTAL"
         << setw(15) << total_single
         << setw(15) << total_v2
         << setw(12) << (to_string(total_speedup).substr(0,4) + "x")
         << "\n";

    cout << "\n";
    cout << "Status: Single=" << (single.verified ? "PASS" : "FAIL")
         << "  V2=" << (v2.verified ? "PASS" : "FAIL") << "\n";
}

void print_v2_breakdown(const V2Result& res) {
    cout << "\n";
    cout << "--- Protocol 2 Commit Breakdown ---\n";
    cout << "  Local encode:     " << fixed << setprecision(2)
         << res.commit_timing.local_encode_ms << " ms\n";
    cout << "  Local hash:       " << res.commit_timing.local_hash_ms << " ms\n";
    cout << "  Hash exchange:    " << res.commit_timing.hash_exchange_ms << " ms (simulated)\n";
    cout << "  Global H_M:       " << res.commit_timing.global_hash_ms << " ms\n";
    cout << "  Local Merkle:     " << res.commit_timing.local_merkle_ms << " ms\n";
    cout << "  Global commit:    " << res.commit_timing.global_commit_ms << " ms\n";
    cout << "  Total:            " << res.commit_timing.total_ms << " ms\n";

    cout << "\n";
    cout << "--- Protocol 2 Prove Breakdown ---\n";
    cout << "  Distribute data:  " << res.prove_timing.distribute_ms << " ms\n";
    cout << "  Parallel combine: " << res.prove_timing.combine_ms << " ms\n";
    cout << "  Collect/aggregate:" << res.prove_timing.collect_ms << " ms\n";
    cout << "  Column openings:  " << res.prove_timing.column_open_ms << " ms\n";
    cout << "  Total:            " << res.prove_timing.total_ms << " ms\n";

    cout << "\n";
    cout << "--- Worker Stats (Commit) ---\n";
    cout << "  " << left << setw(8) << "Worker"
         << setw(8) << "Rows"
         << setw(8) << "Cols"
         << setw(12) << "Encode(ms)"
         << setw(12) << "Hash(ms)"
         << setw(12) << "H_M(ms)"
         << setw(12) << "Merkle(ms)"
         << "\n";
    for (const auto& ws : res.commit_stats) {
        cout << "  " << setw(8) << ws.worker_id
             << setw(8) << ws.rows_processed
             << setw(8) << ws.cols_processed
             << setw(12) << ws.encode_ms
             << setw(12) << ws.local_hash_ms
             << setw(12) << ws.global_hash_ms
             << setw(12) << ws.local_merkle_ms
             << "\n";
    }
}

// ============================================================
// Main
// ============================================================
int main(int argc, char* argv[]) {
    int num_workers = thread::hardware_concurrency();
    if (argc > 1) {
        num_workers = atoi(argv[1]);
    }
    if (num_workers < 1) num_workers = 1;

    cout << "\n";
    cout << "##############################################################\n";
    cout << "  Protocol 2: Distributed Merkle Tree Benchmark\n";
    cout << "  Hardware threads: " << thread::hardware_concurrency() << "\n";
    cout << "  Using workers:    " << num_workers << "\n";
    cout << "##############################################################\n";

    // Parameters
    long s = 2;        // GR(2^2, d)
    long base_r = 128; // Large ring
    long lambda = 128;

    ZZ_p::init(ZZ(1) << s);

    // Test configuration
    long nv = 14;

    cout << "\n";
    cout << "##############################################################\n";
    cout << "  Test: num_vars=" << nv << " (n=" << (1L << nv) << ")\n";
    cout << "##############################################################\n";

    long n = 1L << nv;
    long half = nv / 2;
    long num_rows = 1L << half;
    long row_len = 1L << (nv - half);

    cout << "  num_rows=" << num_rows << ", row_len=" << row_len << "\n";

    // Initialize ring
    initGR(s, base_r);

    // Setup code
    BrakedownCodeGR code = brakedown_code_setup(row_len, s, base_r);
    cout << "  codeword_len=" << code.codeword_len << "\n";

    // Generate polynomial
    vector<ZZ_pE> poly(n);
    for (long i = 0; i < n; i++) {
        poly[i] = random_ZZ_pE();
    }

    // Generate evaluation point
    vector<ZZ_pE> q1(num_rows), q2(row_len);
    for (long i = 0; i < num_rows; i++) q1[i] = random_ZZ_pE();
    for (long j = 0; j < row_len; j++) q2[j] = random_ZZ_pE();

    // Run single-threaded baseline
    cout << "\n  Running single-threaded baseline...\n";
    auto single_res = run_single(code, poly.data(), n, num_rows, q1, q2);

    // Run Protocol 2
    cout << "  Running Protocol 2 (" << num_workers << " workers)...\n";
    auto v2_res = run_v2(code, poly.data(), n, num_rows, q1, q2, num_workers);

    // Print comparison
    print_comparison(single_res, v2_res);
    print_v2_breakdown(v2_res);

    // Scalability test
    cout << "\n";
    cout << "##############################################################\n";
    cout << "  Scalability Test: n=2^" << nv << ", varying workers\n";
    cout << "##############################################################\n";

    cout << "\n";
    cout << left << fixed << setprecision(2);
    cout << setw(10) << "Workers"
         << setw(15) << "Commit(ms)"
         << setw(15) << "Prove(ms)"
         << setw(15) << "Verify(ms)"
         << setw(15) << "Total(ms)"
         << setw(12) << "Speedup"
         << "\n";
    cout << string(82, '-') << "\n";

    double single_total = single_res.commit_ms + single_res.prove_ms;
    cout << setw(10) << "1 (base)"
         << setw(15) << single_res.commit_ms
         << setw(15) << single_res.prove_ms
         << setw(15) << single_res.verify_ms
         << setw(15) << single_total
         << setw(12) << "1.00x"
         << "\n";

    int max_workers = min((int)thread::hardware_concurrency(), 16);
    for (int w = 2; w <= max_workers; w *= 2) {
        auto v2 = run_v2(code, poly.data(), n, num_rows, q1, q2, w);
        double total = v2.commit_timing.total_ms + v2.prove_timing.total_ms;
        double speedup = single_total / total;

        cout << setw(10) << w
             << setw(15) << v2.commit_timing.total_ms
             << setw(15) << v2.prove_timing.total_ms
             << setw(15) << v2.verify_ms
             << setw(15) << total
             << setw(12) << (to_string(speedup).substr(0,4) + "x")
             << "\n";
    }

    cout << "\n";
    cout << "##############################################################\n";
    cout << "  Benchmark complete.\n";
    cout << "##############################################################\n";

    return 0;
}

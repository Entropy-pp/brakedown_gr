/**
 * bench_final.cpp - Final Protocol 2 Scalability Benchmark
 *
 * Compile:
 *   make bench_final
 *
 * Run:
 *   ./out/bench_final [max_workers]
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

int main(int argc, char* argv[]) {
    int max_workers = thread::hardware_concurrency();
    if (argc > 1) {
        max_workers = atoi(argv[1]);
    }
    if (max_workers < 1) max_workers = 1;

    // Parameters
    long s = 2;
    long base_r = 128;
    long nv = 16;

    ZZ_p::init(ZZ(1) << s);

    long n = 1L << nv;
    long half = nv / 2;
    long num_rows = 1L << half;
    long row_len = 1L << (nv - half);

    initGR(s, base_r);
    BrakedownCodeGR code = brakedown_code_setup(row_len, s, base_r);

    // Generate polynomial and evaluation point
    vector<ZZ_pE> poly(n);
    for (long i = 0; i < n; i++) poly[i] = random_ZZ_pE();

    vector<ZZ_pE> q1(num_rows), q2(row_len);
    for (long i = 0; i < num_rows; i++) q1[i] = random_ZZ_pE();
    for (long j = 0; j < row_len; j++) q2[j] = random_ZZ_pE();

    // Print header
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

    // V2 with 1 worker as baseline (fair comparison)
    double baseline_total = 0;

    // Test with varying workers (starting from 1)
    for (int w = 1; w <= max_workers; w *= 2) {
        DistributedCommitTimingV2 commit_timing;
        DistributedProveTiming prove_timing;
        vector<WorkerStatsV2> commit_stats;
        vector<WorkerStats> prove_stats;

        auto comm = distributed_commit_v2(code, poly.data(), n, num_rows, w,
                                          commit_timing, commit_stats);

        auto proof = distributed_prove_v2(code, comm, poly.data(), n, q1, q2, w,
                                          prove_timing, prove_stats);

        auto tv1 = hrc::now();
        bool ok = distributed_verify_v2(code, comm, proof, q1, q2);
        auto tv2 = hrc::now();
        double verify_ms = ms_between(tv1, tv2);

        double total = commit_timing.total_ms + prove_timing.total_ms;

        if (w == 1) {
            baseline_total = total;
        }
        double speedup = baseline_total / total;

        string label = (w == 1) ? "1 (base)" : to_string(w);
        cout << setw(10) << label
             << setw(15) << commit_timing.total_ms
             << setw(15) << prove_timing.total_ms
             << setw(15) << verify_ms
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

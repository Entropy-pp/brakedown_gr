/**
 * bench_distributed.cpp - Benchmark distributed vs single-threaded Brakedown PCS
 *
 * Compile:
 *   g++ -O3 -std=c++17 -pthread -I./include \
 *       test/bench_distributed.cpp \
 *       src/gr.cpp src/sparse_matrix_gr.cpp src/brakedown_code_gr.cpp \
 *       src/merkle.cpp src/brakedown_pcs_gr.cpp src/brakedown_distributed.cpp \
 *       -o out/bench_distributed -lntl -lgmp -lm -lssl -lcrypto
 *
 * Run:
 *   ./out/bench_distributed [num_workers]
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

// ============================================================
// Communication and Proof Size Measurement
// ============================================================

// 计算单个 ZZ_pE 元素的字节大小
static long get_ZZ_pE_byte_size() {
    long k = NumBits(ZZ_p::modulus());  // 模数的位数
    long degree = ZZ_pE::degree();       // 扩展次数
    return degree * ((k + 7) / 8);       // 每个系数的字节数 * 系数个数
}

// 哈希值大小 (SHA-256 = 32 bytes)
static const long HASH_SIZE = 32;

// 通信大小统计结构
struct CommunicationStats {
    int num_workers;               // 节点数量

    // === 每个节点的通信开销 ===
    // Commit阶段 (每节点)
    long per_node_hash_send_bytes;      // 每节点发送的哈希量
    long per_node_hash_recv_bytes;      // 每节点接收的哈希量
    long per_node_commit_bytes;         // 每节点Commit阶段总通信量

    // Prove阶段 (每节点)
    long per_node_combine_bytes;        // 每节点线性组合通信量
    long per_node_column_bytes;         // 每节点列数据通信量
    long per_node_prove_bytes;          // 每节点Prove阶段总通信量

    // 每节点总通信量
    long per_node_total_bytes;

    // === 系统总通信量 ===
    long total_commit_bytes;       // Commit阶段系统总通信量
    long total_prove_bytes;        // Prove阶段系统总通信量
    long total_bytes;              // 系统总通信量
};

// 证明大小统计结构
struct ProofSizeStats {
    long commitment_bytes;         // 承诺大小 (Merkle根)
    long combined_rows_bytes;      // combined_rows 大小
    long prox_coeffs_bytes;        // prox_coeffs 大小
    long column_items_bytes;       // column_items 大小
    long merkle_paths_bytes;       // merkle_paths 大小
    long eval_value_bytes;         // eval_value 大小
    long total_proof_bytes;        // 证明总大小
};

// 计算分布式协议的通信量
CommunicationStats compute_communication_stats(
    const BrakedownCodeGR& code,
    long num_rows,
    int num_workers)
{
    CommunicationStats stats;
    stats.num_workers = num_workers;
    long elem_size = get_ZZ_pE_byte_size();
    long cw_len = code.codeword_len;
    long row_len = code.row_len;
    long num_prox = code.num_prox_test;
    long num_col_open = code.num_col_open;
    long rows_per_worker = (num_rows + num_workers - 1) / num_workers;
    long cols_per_worker = (cw_len + num_workers - 1) / num_workers;

    // ============ Commit阶段通信 ============
    // 哈希交换：每个节点向其他N-1个节点各发送 cw_len/N 个哈希
    // 每节点发送量 = (N-1) * (cw_len/N) * HASH_SIZE
    // 每节点接收量 = (N-1) * (cw_len/N) * HASH_SIZE (接收其他节点发来的哈希)
    stats.per_node_hash_send_bytes = (num_workers - 1) * cols_per_worker * HASH_SIZE;
    stats.per_node_hash_recv_bytes = (num_workers - 1) * cols_per_worker * HASH_SIZE;

    // 局部Merkle根发送给P_1：每个非主导节点发送1个哈希
    // 平均每节点: (N-1)/N * HASH_SIZE ≈ HASH_SIZE (对于大N)
    long root_send_per_node = (num_workers > 1) ? HASH_SIZE : 0;  // 非P_1节点发送

    stats.per_node_commit_bytes = stats.per_node_hash_send_bytes +
                                   stats.per_node_hash_recv_bytes + root_send_per_node;

    // ============ Prove阶段通信 ============
    // 线性组合树状规约：每个节点最多发送1次，接收log2(N)次
    // 每次发送/接收 (num_prox + 1) * row_len 个元素
    long combine_elements = (num_prox + 1) * row_len;
    long combine_msg_size = combine_elements * elem_size;
    // 平均每节点发送: (N-1)/N 的节点发送一次 ≈ combine_msg_size
    stats.per_node_combine_bytes = combine_msg_size;  // 简化：每节点约发送一次

    // 列数据发送给P_1：每个非P_1节点发送 num_col_open * rows_per_worker 个元素
    stats.per_node_column_bytes = num_col_open * rows_per_worker * elem_size;

    stats.per_node_prove_bytes = stats.per_node_combine_bytes + stats.per_node_column_bytes;

    // 每节点总通信量
    stats.per_node_total_bytes = stats.per_node_commit_bytes + stats.per_node_prove_bytes;

    // ============ 系统总通信量 ============
    // Commit: 哈希交换 + 局部根汇总
    // 哈希交换总量 = N * (N-1) * (cw_len/N) * HASH_SIZE = (N-1) * cw_len * HASH_SIZE
    stats.total_commit_bytes = (num_workers - 1) * cw_len * HASH_SIZE +
                               (num_workers - 1) * HASH_SIZE;

    // Prove: 树状规约 + 列数据汇总
    // 树状规约总量 = (N-1) * combine_msg_size
    // 列数据总量 = (N-1) * num_col_open * rows_per_worker * elem_size
    stats.total_prove_bytes = (num_workers - 1) * combine_msg_size +
                              (num_workers - 1) * num_col_open * rows_per_worker * elem_size;

    stats.total_bytes = stats.total_commit_bytes + stats.total_prove_bytes;

    return stats;
}

// 计算证明大小
ProofSizeStats compute_proof_size(
    const BrakedownCodeGR& code,
    const BrakedownEvalProof& proof,
    long num_rows)
{
    ProofSizeStats stats;
    long elem_size = get_ZZ_pE_byte_size();

    // 承诺大小 (Merkle根)
    stats.commitment_bytes = HASH_SIZE;

    // combined_rows: (num_prox + 1) 个向量，每个长度为 row_len
    stats.combined_rows_bytes = 0;
    for (const auto& row : proof.combined_rows) {
        stats.combined_rows_bytes += row.size() * elem_size;
    }

    // prox_coeffs: num_prox 个向量，每个长度为 num_rows
    stats.prox_coeffs_bytes = 0;
    for (const auto& coeffs : proof.prox_coeffs) {
        stats.prox_coeffs_bytes += coeffs.size() * elem_size;
    }

    // column_items: num_col_open 个向量，每个长度为 num_rows
    stats.column_items_bytes = 0;
    for (const auto& col : proof.column_items) {
        stats.column_items_bytes += col.size() * elem_size;
    }

    // merkle_paths: num_col_open 个路径，每个路径包含 depth 个哈希
    stats.merkle_paths_bytes = 0;
    for (const auto& path : proof.merkle_paths) {
        stats.merkle_paths_bytes += path.size() * HASH_SIZE;
    }

    // eval_value: 1个 ZZ_pE 元素
    stats.eval_value_bytes = elem_size;

    // 总证明大小
    stats.total_proof_bytes = stats.commitment_bytes + stats.combined_rows_bytes +
                              stats.prox_coeffs_bytes + stats.column_items_bytes +
                              stats.merkle_paths_bytes + stats.eval_value_bytes;

    return stats;
}

// 格式化字节大小输出
string format_bytes(long bytes) {
    if (bytes < 1024) {
        return to_string(bytes) + " B";
    } else if (bytes < 1024 * 1024) {
        return to_string(bytes / 1024) + "." + to_string((bytes % 1024) * 10 / 1024) + " KB";
    } else {
        return to_string(bytes / (1024 * 1024)) + "." +
               to_string((bytes % (1024 * 1024)) * 10 / (1024 * 1024)) + " MB";
    }
}

// 打印通信统计
void print_communication_stats(const CommunicationStats& stats, int num_workers) {
    cout << "\n";
    cout << "--- Communication Statistics (" << num_workers << " workers) ---\n";

    cout << "\n  [Per-Node Communication]\n";
    cout << "    Commit phase:\n";
    cout << "      Hash send:        " << format_bytes(stats.per_node_hash_send_bytes) << "\n";
    cout << "      Hash receive:     " << format_bytes(stats.per_node_hash_recv_bytes) << "\n";
    cout << "      Commit subtotal:  " << format_bytes(stats.per_node_commit_bytes) << "\n";
    cout << "    Prove phase:\n";
    cout << "      Linear combine:   " << format_bytes(stats.per_node_combine_bytes) << "\n";
    cout << "      Column data:      " << format_bytes(stats.per_node_column_bytes) << "\n";
    cout << "      Prove subtotal:   " << format_bytes(stats.per_node_prove_bytes) << "\n";
    cout << "    Per-node total:     " << format_bytes(stats.per_node_total_bytes) << "\n";

    cout << "\n  [System Total Communication]\n";
    cout << "    Commit phase:       " << format_bytes(stats.total_commit_bytes) << "\n";
    cout << "    Prove phase:        " << format_bytes(stats.total_prove_bytes) << "\n";
    cout << "    System total:       " << format_bytes(stats.total_bytes) << "\n";
}

// 打印证明大小统计
void print_proof_size_stats(const ProofSizeStats& stats) {
    cout << "\n";
    cout << "--- Proof Size to Verifier (independent of worker count) ---\n";
    cout << "  Commitment (root):    " << format_bytes(stats.commitment_bytes) << "\n";
    cout << "  Combined rows:        " << format_bytes(stats.combined_rows_bytes) << "\n";
    cout << "  Prox coefficients:    " << format_bytes(stats.prox_coeffs_bytes) << "\n";
    cout << "  Column items:         " << format_bytes(stats.column_items_bytes) << "\n";
    cout << "  Merkle paths:         " << format_bytes(stats.merkle_paths_bytes) << "\n";
    cout << "  Eval value:           " << format_bytes(stats.eval_value_bytes) << "\n";
    cout << "  -----------------------------------------\n";
    cout << "  [Total Proof Size]    " << format_bytes(stats.total_proof_bytes) << "\n";
    cout << "\n  Note: Proof size is the same regardless of the number of\n";
    cout << "        prover nodes (verifier transparency).\n";
}

static void initGR(long k, long degree) {
    ZZ modulus = ZZ(1) << k;
    ZZ_p::init(modulus);
    ZZ_pX P = primitiveIrredPoly(degree);
    ZZ_pE::init(P);
}

// ============================================================
// Single-threaded benchmark
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
// Distributed benchmark
// ============================================================
struct DistributedResult {
    int num_workers;
    DistributedCommitTiming commit_timing;
    DistributedProveTiming prove_timing;
    double verify_ms;
    bool verified;
    vector<WorkerStats> commit_stats;
    vector<WorkerStats> prove_stats;
    CommunicationStats comm_stats;
    ProofSizeStats proof_stats;
};

DistributedResult run_distributed(
    const BrakedownCodeGR& code,
    const ZZ_pE* poly,
    long n,
    long num_rows,
    const vector<ZZ_pE>& q1,
    const vector<ZZ_pE>& q2,
    int num_workers)
{
    DistributedResult res;
    res.num_workers = num_workers;

    auto comm = distributed_commit(code, poly, n, num_rows, num_workers,
                                   res.commit_timing, res.commit_stats);

    auto proof = distributed_prove(code, comm, poly, n, q1, q2, num_workers,
                                   res.prove_timing, res.prove_stats);

    auto t1 = hrc::now();
    res.verified = brakedown_verify(code, comm, proof, q1, q2);
    auto t2 = hrc::now();
    res.verify_ms = ms_between(t1, t2);

    // 计算通信统计和证明大小
    res.comm_stats = compute_communication_stats(code, num_rows, num_workers);
    res.proof_stats = compute_proof_size(code, proof, num_rows);

    return res;
}

// ============================================================
// Print results
// ============================================================
void print_comparison(const SingleResult& single, const DistributedResult& dist) {
    cout << "\n";
    cout << "============================================================\n";
    cout << "  Single-threaded vs Distributed (" << dist.num_workers << " workers)\n";
    cout << "============================================================\n";
    cout << "\n";

    cout << left << fixed << setprecision(2);
    cout << setw(20) << "Phase"
         << setw(15) << "Single(ms)"
         << setw(15) << "Dist(ms)"
         << setw(12) << "Speedup"
         << "\n";
    cout << string(62, '-') << "\n";

    double commit_speedup = single.commit_ms / dist.commit_timing.total_ms;
    double prove_speedup = single.prove_ms / dist.prove_timing.total_ms;

    cout << setw(20) << "Commit"
         << setw(15) << single.commit_ms
         << setw(15) << dist.commit_timing.total_ms
         << setw(12) << (to_string(commit_speedup).substr(0,4) + "x")
         << "\n";

    cout << setw(20) << "Prove"
         << setw(15) << single.prove_ms
         << setw(15) << dist.prove_timing.total_ms
         << setw(12) << (to_string(prove_speedup).substr(0,4) + "x")
         << "\n";

    cout << setw(20) << "Verify"
         << setw(15) << single.verify_ms
         << setw(15) << dist.verify_ms
         << setw(12) << "-"
         << "\n";

    cout << string(62, '-') << "\n";

    double total_single = single.commit_ms + single.prove_ms + single.verify_ms;
    double total_dist = dist.commit_timing.total_ms + dist.prove_timing.total_ms + dist.verify_ms;
    double total_speedup = total_single / total_dist;

    cout << setw(20) << "TOTAL"
         << setw(15) << total_single
         << setw(15) << total_dist
         << setw(12) << (to_string(total_speedup).substr(0,4) + "x")
         << "\n";

    cout << "\n";
    cout << "Status: Single=" << (single.verified ? "PASS" : "FAIL")
         << "  Distributed=" << (dist.verified ? "PASS" : "FAIL") << "\n";
}

void print_distributed_breakdown(const DistributedResult& res) {
    cout << "\n";
    cout << "--- Distributed Commit Breakdown ---\n";
    cout << "  Distribute rows:  " << fixed << setprecision(2)
         << res.commit_timing.distribute_ms << " ms\n";
    cout << "  Parallel encode:  " << res.commit_timing.encode_ms << " ms\n";
    cout << "  Collect results:  " << res.commit_timing.collect_ms << " ms\n";
    cout << "  Column hashing:   " << res.commit_timing.column_hash_ms << " ms\n";
    cout << "  Merkle tree:      " << res.commit_timing.merkle_build_ms << " ms\n";
    cout << "  Total:            " << res.commit_timing.total_ms << " ms\n";

    cout << "\n";
    cout << "--- Distributed Prove Breakdown ---\n";
    cout << "  Distribute data:  " << res.prove_timing.distribute_ms << " ms\n";
    cout << "  Parallel combine: " << res.prove_timing.combine_ms << " ms\n";
    cout << "  Collect/aggregate:" << res.prove_timing.collect_ms << " ms\n";
    cout << "  Column openings:  " << res.prove_timing.column_open_ms << " ms\n";
    cout << "  Total:            " << res.prove_timing.total_ms << " ms\n";

    cout << "\n";
    cout << "--- Worker Stats (Commit) ---\n";
    cout << "  " << left << setw(8) << "Worker"
         << setw(8) << "Rows"
         << setw(12) << "Start(ms)"
         << setw(12) << "End(ms)"
         << setw(12) << "Elapsed(ms)"
         << "\n";
    for (const auto& ws : res.commit_stats) {
        cout << "  " << setw(8) << ws.worker_id
             << setw(8) << ws.rows_processed
             << setw(12) << ws.encode_start_ms
             << setw(12) << ws.encode_end_ms
             << setw(12) << ws.encode_ms
             << "\n";
    }

    cout << "\n";
    cout << "--- Worker Stats (Prove) ---\n";
    cout << "  " << left << setw(8) << "Worker"
         << setw(8) << "Rows"
         << setw(12) << "Start(ms)"
         << setw(12) << "End(ms)"
         << setw(12) << "Elapsed(ms)"
         << "\n";
    for (const auto& ws : res.prove_stats) {
        cout << "  " << setw(8) << ws.worker_id
             << setw(8) << ws.rows_processed
             << setw(12) << ws.combine_start_ms
             << setw(12) << ws.combine_end_ms
             << setw(12) << ws.combine_ms
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
    cout << "  Distributed Brakedown PCS Benchmark\n";
    cout << "  Hardware threads: " << thread::hardware_concurrency() << "\n";
    cout << "  Using workers:    " << num_workers << "\n";
    cout << "##############################################################\n";

    // Parameters
    long s = 2;        // GR(2^2, d)
    long base_r = 128; // Large ring
    long lambda = 128;

    ZZ_p::init(ZZ(1) << s);

    // Test configurations
    // vector<long> test_nv = {14,15,16};

    // for (long nv : test_nv) {
    //     cout << "\n";
    //     cout << "##############################################################\n";
    //     cout << "  Test: num_vars=" << nv << " (n=" << (1L << nv) << ")\n";
    //     cout << "##############################################################\n";

    //     long n = 1L << nv;
    //     long half = nv / 2;
    //     long num_rows = 1L << half;
    //     long row_len = 1L << (nv - half);

    //     cout << "  num_rows=" << num_rows << ", row_len=" << row_len << "\n";

    //     // Initialize ring
    //     initGR(s, base_r);

    //     // Setup code
    //     BrakedownCodeGR code = brakedown_code_setup(row_len, s, base_r);
    //     cout << "  codeword_len=" << code.codeword_len << "\n";

    //     // Generate polynomial
    //     vector<ZZ_pE> poly(n);
    //     for (long i = 0; i < n; i++) {
    //         poly[i] = random_ZZ_pE();
    //     }

    //     // Generate evaluation point
    //     vector<ZZ_pE> q1(num_rows), q2(row_len);
    //     for (long i = 0; i < num_rows; i++) q1[i] = random_ZZ_pE();
    //     for (long j = 0; j < row_len; j++) q2[j] = random_ZZ_pE();

    //     // Run single-threaded
    //     cout << "\n  Running single-threaded...\n";
    //     auto single_res = run_single(code, poly.data(), n, num_rows, q1, q2);

    //     // Run distributed
    //     cout << "  Running distributed (" << num_workers << " workers)...\n";
    //     auto dist_res = run_distributed(code, poly.data(), n, num_rows, q1, q2, num_workers);

    //     // Print comparison
    //     print_comparison(single_res, dist_res);
    //     print_distributed_breakdown(dist_res);
    // }

    // Scalability test: vary number of workers

    {
        long nv = 14;

        cout << "\n";
        cout << "##############################################################\n";
        cout << "  Scalability Test: n=2^" << nv << ", varying workers\n";
        cout << "##############################################################\n";

        long n = 1L << nv;
        long half = nv / 2;
        long num_rows = 1L << half;
        long row_len = 1L << (nv - half);

        initGR(s, base_r);
        BrakedownCodeGR code = brakedown_code_setup(row_len, s, base_r);

        vector<ZZ_pE> poly(n);
        for (long i = 0; i < n; i++) poly[i] = random_ZZ_pE();

        vector<ZZ_pE> q1(num_rows), q2(row_len);
        for (long i = 0; i < num_rows; i++) q1[i] = random_ZZ_pE();
        for (long j = 0; j < row_len; j++) q2[j] = random_ZZ_pE();

        // Single-threaded baseline
        auto single_res = run_single(code, poly.data(), n, num_rows, q1, q2);

        cout << "\n";
        cout << left << fixed << setprecision(2);
        cout << setw(10) << "Workers"
             << setw(15) << "Commit(ms)"
             << setw(15) << "Prove(ms)"
             << setw(15) << "Verify(ms)"
             << setw(15) << "Total(ms)"
             << setw(12) << "Speedup"
             << "\n";
        cout << string(67, '-') << "\n";

        // double single_total = single_res.commit_ms + single_res.prove_ms + single_res.verify_ms;
        double single_total = single_res.commit_ms + single_res.prove_ms;
        cout << setw(10) << "1 (base)"
             << setw(15) << single_res.commit_ms
             << setw(15) << single_res.prove_ms
             << setw(15) << single_res.verify_ms
             << setw(15) << single_total
             << setw(12) << "1.00x"
             << "\n";

        int max_workers = min((int)thread::hardware_concurrency(), 16);
        DistributedResult last_dist;
        for (int w = 2; w <= max_workers; w *= 2) {
            auto dist = run_distributed(code, poly.data(), n, num_rows, q1, q2, w);
            //  double total = dist.commit_timing.total_ms + dist.prove_timing.total_ms + dist.verify_ms;
            double total = dist.commit_timing.total_ms + dist.prove_timing.total_ms;
            double speedup = single_total / total;

            cout << setw(10) << w
                 << setw(15) << dist.commit_timing.total_ms
                 << setw(15) << dist.prove_timing.total_ms
                 << setw(15) << dist.verify_ms
                 << setw(15) << total
                 << setw(12) << (to_string(speedup).substr(0,4) + "x")
                 << "\n";

            last_dist = dist;
        }

        // 打印通信统计和证明大小 (使用最后一次分布式运行的结果)
        cout << "\n";
        cout << "##############################################################\n";
        cout << "  Communication and Proof Size Analysis\n";
        cout << "##############################################################\n";

        // 打印不同节点数下的通信量 (每节点)
        cout << "\n--- Per-Node Communication by Worker Count ---\n";
        cout << left << fixed << setprecision(2);
        cout << setw(10) << "Workers"
             << setw(18) << "Commit/Node"
             << setw(18) << "Prove/Node"
             << setw(18) << "Total/Node"
             << "\n";
        cout << string(64, '-') << "\n";

        for (int w = 2; w <= max_workers; w *= 2) {
            auto comm_stats = compute_communication_stats(code, num_rows, w);
            cout << setw(10) << w
                 << setw(18) << format_bytes(comm_stats.per_node_commit_bytes)
                 << setw(18) << format_bytes(comm_stats.per_node_prove_bytes)
                 << setw(18) << format_bytes(comm_stats.per_node_total_bytes)
                 << "\n";
        }

        // 打印不同节点数下的系统总通信量
        cout << "\n--- System Total Communication by Worker Count ---\n";
        cout << left << fixed << setprecision(2);
        cout << setw(10) << "Workers"
             << setw(18) << "Commit Total"
             << setw(18) << "Prove Total"
             << setw(18) << "System Total"
             << "\n";
        cout << string(64, '-') << "\n";

        for (int w = 2; w <= max_workers; w *= 2) {
            auto comm_stats = compute_communication_stats(code, num_rows, w);
            cout << setw(10) << w
                 << setw(18) << format_bytes(comm_stats.total_commit_bytes)
                 << setw(18) << format_bytes(comm_stats.total_prove_bytes)
                 << setw(18) << format_bytes(comm_stats.total_bytes)
                 << "\n";
        }

        // 打印详细的通信统计 (以最大worker数为例)
        print_communication_stats(last_dist.comm_stats, max_workers);

        // 打印证明大小统计
        print_proof_size_stats(last_dist.proof_stats);
    }

    cout << "\n";
    cout << "##############################################################\n";
    cout << "  Benchmark complete.\n";
    cout << "##############################################################\n";

    return 0;
}

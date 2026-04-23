#include <brakedown_distributed.h>
#include <brakedown_code_gr.h>
#include <brakedown_pcs_gr.h>
#include <merkle.h>
#include <gr.h>
#include <iostream>
#include <chrono>
#include <cstring>
#include <atomic>

using namespace NTL;
using hrc = std::chrono::high_resolution_clock;

static double ms_between(hrc::time_point start, hrc::time_point end) {
    return std::chrono::duration<double, std::milli>(end - start).count();
}

// ============================================================
// Thread-local NTL context management
// NTL's ZZ_p and ZZ_pE use thread-local storage, so each thread
// needs to initialize its own context.
// primitiveIrredPoly may not be thread-safe, so we protect it.
// ============================================================

static std::mutex ntl_init_mutex;

static void init_thread_context(long k, long degree) {
    ZZ modulus = ZZ(1) << k;
    ZZ_p::init(modulus);

    ZZ_pX P;
    {
        // Protect primitiveIrredPoly call with mutex
        std::lock_guard<std::mutex> lock(ntl_init_mutex);
        P = primitiveIrredPoly(degree);
    }
    ZZ_pE::init(P);
}

// ============================================================
// Worker task: encode assigned rows
// ============================================================

struct EncodeTask {
    const BrakedownCodeGR* code;
    const ZZ_pE* input_rows;     // Pointer to start of this worker's input
    ZZ_pE* output_rows;          // Pointer to start of this worker's output
    long start_row;              // Global row index start
    long end_row;                // Global row index end (exclusive)
    long row_len;
    long codeword_len;
    long k;                      // For NTL context
    long degree;                 // For NTL context

    // Timing: use atomic to safely record timestamps from worker thread
    std::atomic<double> start_offset_ms{0};  // Time from global_start to worker start
    std::atomic<double> end_offset_ms{0};    // Time from global_start to worker end
    hrc::time_point* global_start;           // Shared reference time point
};

static void worker_encode(EncodeTask* task) {
    // Record worker start time relative to global start
    auto worker_start = hrc::now();
    task->start_offset_ms.store(ms_between(*task->global_start, worker_start));

    // Initialize NTL context for this thread
    init_thread_context(task->k, task->degree);

    long row_len = task->row_len;
    long cw_len = task->codeword_len;

    for (long i = task->start_row; i < task->end_row; i++) {
        long local_i = i - task->start_row;

        // Copy input to output buffer (message portion)
        for (long j = 0; j < row_len; j++) {
            task->output_rows[local_i * cw_len + j] = task->input_rows[local_i * row_len + j];
        }
        // Zero out parity portion
        for (long j = row_len; j < cw_len; j++) {
            clear(task->output_rows[local_i * cw_len + j]);
        }
        // Encode in place
        brakedown_encode(*task->code, &task->output_rows[local_i * cw_len]);
    }

    // Record worker end time relative to global start
    auto worker_end = hrc::now();
    task->end_offset_ms.store(ms_between(*task->global_start, worker_end));
}

// ============================================================
// Distributed Commit Implementation
// ============================================================

BrakedownCommitment distributed_commit(
    const BrakedownCodeGR& code,
    const ZZ_pE* poly_coeffs,
    long n,
    long num_rows,
    int num_workers,
    DistributedCommitTiming& timing,
    std::vector<WorkerStats>& worker_stats)
{
    auto t_total_start = hrc::now();

    BrakedownCommitment comm;
    comm.num_rows = num_rows;
    comm.codeword_len = code.codeword_len;
    long row_len = code.row_len;
    long cw_len = code.codeword_len;

    // Get current NTL context parameters
    long k = NumBits(ZZ_p::modulus()) - 1;  // ZZ_p::modulus() = 2^k
    long degree = ZZ_pE::degree();

    // Allocate output buffer
    comm.encoded_rows.resize(num_rows * cw_len);

    // ============================================================
    // Phase 1: Distribute rows to workers (simulated by partitioning)
    // ============================================================
    auto t1 = hrc::now();

    // Prepare input matrix (rows of polynomial coefficients)
    std::vector<ZZ_pE> input_matrix(num_rows * row_len);
    for (long i = 0; i < num_rows; i++) {
        for (long j = 0; j < row_len; j++) {
            long idx = i * row_len + j;
            if (idx < n) {
                input_matrix[i * row_len + j] = poly_coeffs[idx];
            } else {
                clear(input_matrix[i * row_len + j]);
            }
        }
    }

    // Partition rows among workers
    std::vector<EncodeTask> tasks(num_workers);
    std::vector<std::vector<ZZ_pE>> worker_outputs(num_workers);

    long rows_per_worker = (num_rows + num_workers - 1) / num_workers;

    for (int w = 0; w < num_workers; w++) {
        long start = w * rows_per_worker;
        long end = std::min(start + rows_per_worker, num_rows);
        if (start >= num_rows) {
            start = end = num_rows;  // Empty task
        }

        long worker_rows = end - start;
        worker_outputs[w].resize(worker_rows * cw_len);

        tasks[w].code = &code;
        tasks[w].input_rows = &input_matrix[start * row_len];
        tasks[w].output_rows = worker_outputs[w].data();
        tasks[w].start_row = start;
        tasks[w].end_row = end;
        tasks[w].row_len = row_len;
        tasks[w].codeword_len = cw_len;
        tasks[w].k = k;
        tasks[w].degree = degree;
        // global_start will be set just before launching threads
    }

    auto t2 = hrc::now();
    timing.distribute_ms = ms_between(t1, t2);

    // ============================================================
    // Phase 2: Workers encode in parallel
    // ============================================================
    auto t3 = hrc::now();

    // Set global start time for all workers
    for (int w = 0; w < num_workers; w++) {
        tasks[w].global_start = &t3;
    }

    std::vector<std::thread> threads;
    for (int w = 0; w < num_workers; w++) {
        if (tasks[w].start_row < tasks[w].end_row) {
            threads.emplace_back(worker_encode, &tasks[w]);
        }
    }

    // Wait for all workers to finish
    for (auto& t : threads) {
        t.join();
    }

    auto t4 = hrc::now();
    timing.encode_ms = ms_between(t3, t4);

    // ============================================================
    // Phase 3: Collect encoded rows from workers
    // ============================================================
    auto t5 = hrc::now();

    for (int w = 0; w < num_workers; w++) {
        long start = tasks[w].start_row;
        long end = tasks[w].end_row;
        for (long i = start; i < end; i++) {
            long local_i = i - start;
            for (long j = 0; j < cw_len; j++) {
                comm.encoded_rows[i * cw_len + j] = worker_outputs[w][local_i * cw_len + j];
            }
        }
    }

    auto t6 = hrc::now();
    timing.collect_ms = ms_between(t5, t6);

    // ============================================================
    // Phase 4: Master computes column hashes
    // ============================================================
    auto t7 = hrc::now();

    std::vector<std::vector<unsigned char>> leaf_hashes(cw_len);
    for (long col = 0; col < cw_len; col++) {
        std::vector<ZZ_pE> column(num_rows);
        for (long row = 0; row < num_rows; row++) {
            column[row] = comm.encoded_rows[row * cw_len + col];
        }
        leaf_hashes[col] = hash_column(column.data(), num_rows);
    }

    auto t8 = hrc::now();
    timing.column_hash_ms = ms_between(t7, t8);

    // ============================================================
    // Phase 5: Master builds Merkle tree
    // ============================================================
    auto t9 = hrc::now();

    comm.tree = build_merkle_tree(leaf_hashes);
    comm.root = comm.tree.root;

    auto t10 = hrc::now();
    timing.merkle_build_ms = ms_between(t9, t10);

    auto t_total_end = hrc::now();
    timing.total_ms = ms_between(t_total_start, t_total_end);

    // Collect worker stats (compute elapsed from start/end offsets)
    worker_stats.resize(num_workers);
    for (int w = 0; w < num_workers; w++) {
        worker_stats[w].worker_id = w;
        worker_stats[w].rows_processed = tasks[w].end_row - tasks[w].start_row;
        worker_stats[w].encode_start_ms = tasks[w].start_offset_ms.load();
        worker_stats[w].encode_end_ms = tasks[w].end_offset_ms.load();
        // Elapsed = end_offset - start_offset (both relative to same global_start)
        worker_stats[w].encode_ms = worker_stats[w].encode_end_ms - worker_stats[w].encode_start_ms;
        worker_stats[w].combine_ms = 0;
        worker_stats[w].combine_start_ms = 0;
        worker_stats[w].combine_end_ms = 0;
    }

    return comm;
}

// Convenience overload
BrakedownCommitment distributed_commit(
    const BrakedownCodeGR& code,
    const ZZ_pE* poly_coeffs,
    long n,
    long num_rows,
    int num_workers)
{
    DistributedCommitTiming timing;
    std::vector<WorkerStats> stats;
    return distributed_commit(code, poly_coeffs, n, num_rows, num_workers, timing, stats);
}

// ============================================================
// Worker task: compute ALL partial combines (batched)
// ============================================================

struct BatchCombineTask {
    const ZZ_pE* poly_rows;      // Pointer to this worker's polynomial rows
    const std::vector<const ZZ_pE*>* all_coeffs;  // All coefficient vectors
    std::vector<std::vector<ZZ_pE>>* all_partial_results;  // Output: all partial sums
    long start_row;
    long end_row;
    long row_len;
    long n;
    long k;
    long degree;
    int worker_id;
    int num_workers;             // Total number of workers (for indexing)
    long num_combines;           // Number of combined rows to compute

    std::atomic<double> start_offset_ms{0};
    std::atomic<double> end_offset_ms{0};
    hrc::time_point* global_start;
};

static void worker_batch_combine(BatchCombineTask* task) {
    auto worker_start = hrc::now();
    task->start_offset_ms.store(ms_between(*task->global_start, worker_start));

    // Initialize NTL context ONCE for this thread
    init_thread_context(task->k, task->degree);

    long row_len = task->row_len;
    int w = task->worker_id;
    int num_workers = task->num_workers;

    // Process all combined rows in one go
    for (long c = 0; c < task->num_combines; c++) {
        // Index: combine_idx * num_workers + worker_id
        // Partial result is already pre-allocated and zeroed
        std::vector<ZZ_pE>& partial = (*task->all_partial_results)[c * num_workers + w];

        const ZZ_pE* coeffs = (*task->all_coeffs)[c];

        // Compute partial sum (accumulate into pre-zeroed vector)
        for (long i = task->start_row; i < task->end_row; i++) {
            long local_i = i - task->start_row;
            for (long j = 0; j < row_len; j++) {
                long idx = i * row_len + j;
                ZZ_pE val;
                if (idx < task->n) {
                    val = task->poly_rows[local_i * row_len + j];
                } else {
                    clear(val);
                }
                partial[j] += coeffs[i] * val;
            }
        }
    }

    auto worker_end = hrc::now();
    task->end_offset_ms.store(ms_between(*task->global_start, worker_end));
}

// ============================================================
// Distributed Prove Implementation
// ============================================================

// Helper: compute effective num_prox (copied from brakedown_pcs_gr.cpp)
static long effective_num_prox_dist(const BrakedownCodeGR& code, long num_rows) {
    if (num_rows <= 1) return 0;
    long num_prox = code.num_prox_test;
    if (code.is_small_ring) {
        long repetitions = (code.lambda + code.base_degree - 1) / code.base_degree;
        num_prox = std::max(num_prox, repetitions);
    }
    return num_prox;
}

// Helper: inner product
static ZZ_pE inner_product_dist(const ZZ_pE* a, const ZZ_pE* b, long len) {
    ZZ_pX sum_poly;
    ZZ_pX tmp_poly;
    clear(sum_poly);
    for (long i = 0; i < len; i++) {
        mul(tmp_poly, rep(a[i]), rep(b[i]));
        add(sum_poly, sum_poly, tmp_poly);
    }
    return to_ZZ_pE(sum_poly);
}

BrakedownEvalProof distributed_prove(
    const BrakedownCodeGR& code,
    const BrakedownCommitment& comm,
    const ZZ_pE* poly_coeffs,
    long n,
    const std::vector<ZZ_pE>& q1,
    const std::vector<ZZ_pE>& q2,
    int num_workers,
    DistributedProveTiming& timing,
    std::vector<WorkerStats>& worker_stats)
{
    auto t_total_start = hrc::now();

    BrakedownEvalProof proof;
    long row_len = code.row_len;
    long num_rows = comm.num_rows;
    long num_prox = effective_num_prox_dist(code, num_rows);

    // Get current NTL context parameters
    long k = NumBits(ZZ_p::modulus()) - 1;
    long degree = ZZ_pE::degree();

    // Only handle large ring for now
    if (code.is_small_ring) {
        std::cerr << "distributed_prove: small ring not yet supported" << std::endl;
        return proof;
    }

    // ============================================================
    // Phase 1: Distribute data to workers
    // ============================================================
    auto t1 = hrc::now();

    // Prepare input matrix
    std::vector<ZZ_pE> input_matrix(num_rows * row_len);
    for (long i = 0; i < num_rows; i++) {
        for (long j = 0; j < row_len; j++) {
            long idx = i * row_len + j;
            if (idx < n) {
                input_matrix[i * row_len + j] = poly_coeffs[idx];
            } else {
                clear(input_matrix[i * row_len + j]);
            }
        }
    }

    // Partition rows among workers
    long rows_per_worker = (num_rows + num_workers - 1) / num_workers;

    // Generate prox coefficients
    proof.prox_coeffs.resize(num_prox);
    for (long t = 0; t < num_prox; t++) {
        proof.prox_coeffs[t].resize(num_rows);
        for (long i = 0; i < num_rows; i++) {
            proof.prox_coeffs[t][i] = random_ZZ_pE();
        }
    }

    auto t2 = hrc::now();
    timing.distribute_ms = ms_between(t1, t2);

    // ============================================================
    // Phase 2: Workers compute ALL partial combines in ONE batch
    // This avoids creating threads multiple times
    // ============================================================
    auto t3 = hrc::now();

    long num_combines = num_prox + 1;
    std::vector<std::vector<ZZ_pE>> all_combined(num_combines);

    // Collect all coefficient vectors
    std::vector<const ZZ_pE*> all_coeffs(num_combines);
    for (long t = 0; t < num_prox; t++) {
        all_coeffs[t] = proof.prox_coeffs[t].data();
    }
    all_coeffs[num_prox] = q1.data();

    // Prepare batch tasks - each worker processes ALL combines for its rows
    std::vector<BatchCombineTask> tasks(num_workers);

    // Pre-allocate storage for all partial results to avoid resize in threads
    // Layout: [combine_idx * num_workers + worker_id]
    std::vector<std::vector<ZZ_pE>> all_partial_results(num_combines * num_workers);
    for (long c = 0; c < num_combines; c++) {
        for (int w = 0; w < num_workers; w++) {
            all_partial_results[c * num_workers + w].resize(row_len);
            // Initialize to zero
            for (long j = 0; j < row_len; j++) {
                clear(all_partial_results[c * num_workers + w][j]);
            }
        }
    }

    for (int w = 0; w < num_workers; w++) {
        long start = w * rows_per_worker;
        long end = std::min(start + rows_per_worker, num_rows);
        if (start >= num_rows) {
            start = end = num_rows;
        }

        tasks[w].poly_rows = &input_matrix[start * row_len];
        tasks[w].all_coeffs = &all_coeffs;
        tasks[w].all_partial_results = &all_partial_results;
        tasks[w].start_row = start;
        tasks[w].end_row = end;
        tasks[w].row_len = row_len;
        tasks[w].n = n;
        tasks[w].k = k;
        tasks[w].degree = degree;
        tasks[w].worker_id = w;
        tasks[w].num_workers = num_workers;
        tasks[w].num_combines = num_combines;
        tasks[w].global_start = &t3;
    }

    // Launch workers ONCE for all combines
    std::vector<std::thread> threads;
    for (int w = 0; w < num_workers; w++) {
        if (tasks[w].start_row < tasks[w].end_row) {
            threads.emplace_back(worker_batch_combine, &tasks[w]);
        }
    }

    for (auto& t : threads) {
        t.join();
    }

    // Aggregate all partial results
    for (long c = 0; c < num_combines; c++) {
        all_combined[c].resize(row_len);
        for (long j = 0; j < row_len; j++) {
            clear(all_combined[c][j]);
        }
        for (int w = 0; w < num_workers; w++) {
            const auto& partial = all_partial_results[c * num_workers + w];
            for (long j = 0; j < row_len; j++) {
                all_combined[c][j] += partial[j];
            }
        }
    }

    // Collect worker stats
    worker_stats.resize(num_workers);
    for (int w = 0; w < num_workers; w++) {
        worker_stats[w].worker_id = w;
        worker_stats[w].rows_processed = tasks[w].end_row - tasks[w].start_row;
        worker_stats[w].encode_ms = 0;
        worker_stats[w].encode_start_ms = 0;
        worker_stats[w].encode_end_ms = 0;
        worker_stats[w].combine_start_ms = tasks[w].start_offset_ms.load();
        worker_stats[w].combine_end_ms = tasks[w].end_offset_ms.load();
        worker_stats[w].combine_ms = worker_stats[w].combine_end_ms - worker_stats[w].combine_start_ms;
    }

    auto t4 = hrc::now();
    timing.combine_ms = ms_between(t3, t4);

    // ============================================================
    // Phase 3: Master collects and finalizes
    // ============================================================
    auto t5 = hrc::now();

    // Store combined rows in proof
    for (long c = 0; c < num_prox; c++) {
        proof.combined_rows.push_back(all_combined[c]);
    }

    // Compute eval value
    proof.eval_value = inner_product_dist(all_combined[num_prox].data(), q2.data(), row_len);
    proof.combined_rows.push_back(all_combined[num_prox]);

    auto t6 = hrc::now();
    timing.collect_ms = ms_between(t5, t6);
    timing.tree_reduce_ms = 0;  // V1 doesn't use tree reduction

    // ============================================================
    // Phase 4: Column openings
    // ============================================================
    auto t7 = hrc::now();

    long num_open = code.num_col_open;
    proof.column_indices.resize(num_open);
    proof.column_items.resize(num_open);
    proof.merkle_paths.resize(num_open);

    for (long t = 0; t < num_open; t++) {
        long col_idx = RandomBnd(code.codeword_len);
        proof.column_indices[t] = col_idx;

        std::vector<ZZ_pE> column(num_rows);
        for (long row = 0; row < num_rows; row++) {
            column[row] = comm.encoded_rows[row * comm.codeword_len + col_idx];
        }
        proof.column_items[t] = column;
        proof.merkle_paths[t] = get_merkle_path(comm.tree, col_idx);
    }

    auto t8 = hrc::now();
    timing.column_open_ms = ms_between(t7, t8);

    auto t_total_end = hrc::now();
    timing.total_ms = ms_between(t_total_start, t_total_end);

    return proof;
}

// Convenience overload
BrakedownEvalProof distributed_prove(
    const BrakedownCodeGR& code,
    const BrakedownCommitment& comm,
    const ZZ_pE* poly_coeffs,
    long n,
    const std::vector<ZZ_pE>& q1,
    const std::vector<ZZ_pE>& q2,
    int num_workers)
{
    DistributedProveTiming timing;
    std::vector<WorkerStats> stats;
    return distributed_prove(code, comm, poly_coeffs, n, q1, q2, num_workers, timing, stats);
}

// ============================================================
// Protocol 2: Distributed Commit with Distributed Merkle Tree
// ============================================================

// Worker task for V2: encode rows AND compute local column hashes
struct EncodeAndHashTask {
    const BrakedownCodeGR* code;
    const ZZ_pE* input_rows;
    ZZ_pE* output_rows;                          // Encoded rows output
    std::vector<std::vector<unsigned char>>* local_hashes;  // h^(i)[k] for all columns k
    long start_row;
    long end_row;
    long row_len;
    long codeword_len;
    long k;
    long degree;
    int worker_id;

    std::atomic<double> encode_end_ms{0};
    std::atomic<double> hash_end_ms{0};
    hrc::time_point* global_start;
};

static void worker_encode_and_hash(EncodeAndHashTask* task) {
    auto worker_start = hrc::now();

    // Initialize NTL context
    init_thread_context(task->k, task->degree);

    long row_len = task->row_len;
    long cw_len = task->codeword_len;
    long num_local_rows = task->end_row - task->start_row;

    // Phase 1: Encode all rows
    for (long i = task->start_row; i < task->end_row; i++) {
        long local_i = i - task->start_row;

        // Copy input to output buffer
        for (long j = 0; j < row_len; j++) {
            task->output_rows[local_i * cw_len + j] = task->input_rows[local_i * row_len + j];
        }
        for (long j = row_len; j < cw_len; j++) {
            clear(task->output_rows[local_i * cw_len + j]);
        }
        brakedown_encode(*task->code, &task->output_rows[local_i * cw_len]);
    }

    auto encode_end = hrc::now();
    task->encode_end_ms.store(ms_between(*task->global_start, encode_end));

    // Phase 2: Compute local column hashes h^(i)[k] for all columns k
    // h^(i)[k] = H(Û^(i)_{*,k}) - hash of worker i's portion of column k
    task->local_hashes->resize(cw_len);
    for (long col = 0; col < cw_len; col++) {
        // Extract this worker's portion of column col
        std::vector<ZZ_pE> col_portion(num_local_rows);
        for (long local_i = 0; local_i < num_local_rows; local_i++) {
            col_portion[local_i] = task->output_rows[local_i * cw_len + col];
        }
        (*task->local_hashes)[col] = hash_column(col_portion.data(), num_local_rows);
    }

    auto hash_end = hrc::now();
    task->hash_end_ms.store(ms_between(*task->global_start, hash_end));
}

// Worker task for V2 Phase 2: Global hash aggregation + local Merkle tree
struct AggregateAndMerkleTask {
    const std::vector<std::vector<std::vector<unsigned char>>>* all_worker_hashes;  // [worker][col] -> hash
    std::vector<std::vector<unsigned char>>* global_hashes;  // Output: h[k] for assigned columns
    MerkleTree* local_tree;                                   // Output: local Merkle tree
    long col_start;
    long col_end;
    int num_workers;
    int worker_id;
    long k;
    long degree;

    std::atomic<double> aggregate_end_ms{0};
    std::atomic<double> merkle_end_ms{0};
    hrc::time_point* global_start;
};

static void worker_aggregate_and_merkle(AggregateAndMerkleTask* task) {
    // Initialize NTL context (may be needed for some operations)
    init_thread_context(task->k, task->degree);

    long num_cols = task->col_end - task->col_start;

    // Phase 1: Aggregate hashes with H_M
    // h[k] = H_M(h^(0)[k], h^(1)[k], ..., h^(M-1)[k])
    task->global_hashes->resize(num_cols);
    for (long col = task->col_start; col < task->col_end; col++) {
        long local_col = col - task->col_start;

        // Collect all workers' hashes for this column
        std::vector<std::vector<unsigned char>> col_hashes(task->num_workers);
        for (int w = 0; w < task->num_workers; w++) {
            col_hashes[w] = (*task->all_worker_hashes)[w][col];
        }

        // H_M aggregation
        (*task->global_hashes)[local_col] = hash_M(col_hashes);
    }

    auto aggregate_end = hrc::now();
    task->aggregate_end_ms.store(ms_between(*task->global_start, aggregate_end));

    // Phase 2: Build local Merkle tree over assigned global hashes
    *task->local_tree = build_merkle_tree(*task->global_hashes);

    auto merkle_end = hrc::now();
    task->merkle_end_ms.store(ms_between(*task->global_start, merkle_end));
}

DistributedCommitmentV2 distributed_commit_v2(
    const BrakedownCodeGR& code,
    const ZZ_pE* poly_coeffs,
    long n,
    long num_rows,
    int num_workers,
    DistributedCommitTimingV2& timing,
    std::vector<WorkerStatsV2>& worker_stats)
{
    auto t_total_start = hrc::now();

    DistributedCommitmentV2 comm;
    comm.num_rows = num_rows;
    comm.codeword_len = code.codeword_len;
    comm.num_workers = num_workers;
    long row_len = code.row_len;
    long cw_len = code.codeword_len;

    long k = NumBits(ZZ_p::modulus()) - 1;
    long degree = ZZ_pE::degree();

    // Prepare input matrix
    std::vector<ZZ_pE> input_matrix(num_rows * row_len);
    for (long i = 0; i < num_rows; i++) {
        for (long j = 0; j < row_len; j++) {
            long idx = i * row_len + j;
            if (idx < n) {
                input_matrix[i * row_len + j] = poly_coeffs[idx];
            } else {
                clear(input_matrix[i * row_len + j]);
            }
        }
    }

    // ============================================================
    // Phase 1: Workers encode rows AND compute local column hashes
    // ============================================================
    auto t1 = hrc::now();

    long rows_per_worker = (num_rows + num_workers - 1) / num_workers;

    std::vector<EncodeAndHashTask> encode_tasks(num_workers);
    std::vector<std::vector<ZZ_pE>> worker_outputs(num_workers);
    std::vector<std::vector<std::vector<unsigned char>>> all_worker_hashes(num_workers);

    for (int w = 0; w < num_workers; w++) {
        long start = w * rows_per_worker;
        long end = std::min(start + rows_per_worker, num_rows);
        if (start >= num_rows) {
            start = end = num_rows;
        }

        long worker_rows = end - start;
        worker_outputs[w].resize(worker_rows * cw_len);

        encode_tasks[w].code = &code;
        encode_tasks[w].input_rows = &input_matrix[start * row_len];
        encode_tasks[w].output_rows = worker_outputs[w].data();
        encode_tasks[w].local_hashes = &all_worker_hashes[w];
        encode_tasks[w].start_row = start;
        encode_tasks[w].end_row = end;
        encode_tasks[w].row_len = row_len;
        encode_tasks[w].codeword_len = cw_len;
        encode_tasks[w].k = k;
        encode_tasks[w].degree = degree;
        encode_tasks[w].worker_id = w;
        encode_tasks[w].global_start = &t1;
    }

    std::vector<std::thread> threads;
    for (int w = 0; w < num_workers; w++) {
        if (encode_tasks[w].start_row < encode_tasks[w].end_row) {
            threads.emplace_back(worker_encode_and_hash, &encode_tasks[w]);
        }
    }
    for (auto& t : threads) {
        t.join();
    }

    auto t2 = hrc::now();
    timing.local_encode_ms = 0;
    timing.local_hash_ms = 0;
    for (int w = 0; w < num_workers; w++) {
        timing.local_encode_ms = std::max(timing.local_encode_ms, encode_tasks[w].encode_end_ms.load());
        timing.local_hash_ms = std::max(timing.local_hash_ms,
            encode_tasks[w].hash_end_ms.load() - encode_tasks[w].encode_end_ms.load());
    }

    // Collect encoded rows into comm (still needed for prove phase)
    comm.encoded_rows.resize(num_rows * cw_len);
    for (int w = 0; w < num_workers; w++) {
        long start = encode_tasks[w].start_row;
        long end = encode_tasks[w].end_row;
        for (long i = start; i < end; i++) {
            long local_i = i - start;
            for (long j = 0; j < cw_len; j++) {
                comm.encoded_rows[i * cw_len + j] = worker_outputs[w][local_i * cw_len + j];
            }
        }
    }

    // ============================================================
    // Phase 2: Hash exchange (simulated - in real distributed system,
    // workers would exchange hashes for their assigned columns)
    // ============================================================
    auto t3 = hrc::now();
    timing.hash_exchange_ms = ms_between(t2, t3);  // Simulated as zero-time in shared memory

    // ============================================================
    // Phase 3: Global hash aggregation + local Merkle tree construction
    // Each worker is responsible for cw_len/num_workers columns
    // ============================================================
    auto t4 = hrc::now();

    long cols_per_worker = (cw_len + num_workers - 1) / num_workers;

    std::vector<AggregateAndMerkleTask> agg_tasks(num_workers);
    std::vector<std::vector<std::vector<unsigned char>>> worker_global_hashes(num_workers);
    comm.local_trees.resize(num_workers);

    for (int w = 0; w < num_workers; w++) {
        long col_start = w * cols_per_worker;
        long col_end = std::min(col_start + cols_per_worker, cw_len);
        if (col_start >= cw_len) {
            col_start = col_end = cw_len;
        }

        agg_tasks[w].all_worker_hashes = &all_worker_hashes;
        agg_tasks[w].global_hashes = &worker_global_hashes[w];
        agg_tasks[w].local_tree = &comm.local_trees[w];
        agg_tasks[w].col_start = col_start;
        agg_tasks[w].col_end = col_end;
        agg_tasks[w].num_workers = num_workers;
        agg_tasks[w].worker_id = w;
        agg_tasks[w].k = k;
        agg_tasks[w].degree = degree;
        agg_tasks[w].global_start = &t4;
    }

    threads.clear();
    for (int w = 0; w < num_workers; w++) {
        if (agg_tasks[w].col_start < agg_tasks[w].col_end) {
            threads.emplace_back(worker_aggregate_and_merkle, &agg_tasks[w]);
        }
    }
    for (auto& t : threads) {
        t.join();
    }

    auto t5 = hrc::now();
    timing.global_hash_ms = 0;
    timing.local_merkle_ms = 0;
    for (int w = 0; w < num_workers; w++) {
        timing.global_hash_ms = std::max(timing.global_hash_ms, agg_tasks[w].aggregate_end_ms.load());
        timing.local_merkle_ms = std::max(timing.local_merkle_ms,
            agg_tasks[w].merkle_end_ms.load() - agg_tasks[w].aggregate_end_ms.load());
    }

    // Collect global column hashes (needed for proof verification)
    comm.global_column_hashes.resize(cw_len);
    for (int w = 0; w < num_workers; w++) {
        long col_start = w * cols_per_worker;
        long col_end = std::min(col_start + cols_per_worker, cw_len);
        for (long col = col_start; col < col_end; col++) {
            comm.global_column_hashes[col] = worker_global_hashes[w][col - col_start];
        }
    }

    // ============================================================
    // Phase 4: Global commitment
    // MC(Û) = H(r_0, r_1, ..., r_{M-1})
    // ============================================================
    auto t6 = hrc::now();

    comm.local_roots.resize(num_workers);
    for (int w = 0; w < num_workers; w++) {
        comm.local_roots[w] = comm.local_trees[w].root;
    }

    // Final root = H(r_0, r_1, ..., r_{M-1})
    comm.root = hash_M(comm.local_roots);

    auto t7 = hrc::now();
    timing.global_commit_ms = ms_between(t6, t7);

    auto t_total_end = hrc::now();
    timing.total_ms = ms_between(t_total_start, t_total_end);

    // Collect worker stats
    worker_stats.resize(num_workers);
    for (int w = 0; w < num_workers; w++) {
        worker_stats[w].worker_id = w;
        worker_stats[w].rows_processed = encode_tasks[w].end_row - encode_tasks[w].start_row;
        worker_stats[w].cols_processed = agg_tasks[w].col_end - agg_tasks[w].col_start;
        worker_stats[w].encode_ms = encode_tasks[w].encode_end_ms.load();
        worker_stats[w].local_hash_ms = encode_tasks[w].hash_end_ms.load() - encode_tasks[w].encode_end_ms.load();
        worker_stats[w].global_hash_ms = agg_tasks[w].aggregate_end_ms.load();
        worker_stats[w].local_merkle_ms = agg_tasks[w].merkle_end_ms.load() - agg_tasks[w].aggregate_end_ms.load();
    }

    return comm;
}

// Convenience overload
DistributedCommitmentV2 distributed_commit_v2(
    const BrakedownCodeGR& code,
    const ZZ_pE* poly_coeffs,
    long n,
    long num_rows,
    int num_workers)
{
    DistributedCommitTimingV2 timing;
    std::vector<WorkerStatsV2> stats;
    return distributed_commit_v2(code, poly_coeffs, n, num_rows, num_workers, timing, stats);
}

// ============================================================
// Tree Reduction for aggregating partial results
// ============================================================

// Task for one round of tree reduction
struct TreeReduceTask {
    std::vector<std::vector<ZZ_pE>>* partial_results;  // [combine_idx * num_workers + worker_id]
    int src_worker;      // Worker to receive from
    int dst_worker;      // Worker to accumulate into
    int num_workers;
    long num_combines;
    long row_len;
    long k;
    long degree;
    hrc::time_point* global_start;
};

static void worker_tree_reduce(TreeReduceTask* task) {
    // Initialize NTL context
    init_thread_context(task->k, task->degree);

    int src = task->src_worker;
    int dst = task->dst_worker;
    long row_len = task->row_len;

    // For each combine, add src's result to dst's result
    for (long c = 0; c < task->num_combines; c++) {
        auto& dst_partial = (*task->partial_results)[c * task->num_workers + dst];
        const auto& src_partial = (*task->partial_results)[c * task->num_workers + src];

        for (long j = 0; j < row_len; j++) {
            dst_partial[j] += src_partial[j];
        }
    }
}

// Perform tree reduction on partial results
// Returns timing for the reduction phase
static double tree_reduce_partial_results(
    std::vector<std::vector<ZZ_pE>>& all_partial_results,
    int num_workers,
    long num_combines,
    long row_len,
    long k,
    long degree)
{
    auto t_start = hrc::now();

    // Find e such that 2^e >= num_workers
    int e = 0;
    int padded_workers = 1;
    while (padded_workers < num_workers) {
        padded_workers <<= 1;
        e++;
    }

    // Tree reduction: e rounds
    // Round k (1-indexed): for i < 2^{e-k}, P_i += P_{i + 2^{e-k}}
    for (int round = 1; round <= e; round++) {
        int stride = 1 << (e - round);  // 2^{e-k}
        int active_workers = 1 << (e - round);  // Number of receiving workers

        std::vector<TreeReduceTask> tasks;
        std::vector<std::thread> threads;

        for (int i = 0; i < active_workers; i++) {
            int dst = i;
            int src = i + stride;

            // Skip if src is beyond actual workers
            if (src >= num_workers) continue;
            // Skip if dst is beyond actual workers (shouldn't happen)
            if (dst >= num_workers) continue;

            TreeReduceTask task;
            task.partial_results = &all_partial_results;
            task.src_worker = src;
            task.dst_worker = dst;
            task.num_workers = num_workers;
            task.num_combines = num_combines;
            task.row_len = row_len;
            task.k = k;
            task.degree = degree;
            task.global_start = &t_start;
            tasks.push_back(task);
        }

        // Launch threads for this round
        for (auto& task : tasks) {
            threads.emplace_back(worker_tree_reduce, &task);
        }
        for (auto& t : threads) {
            t.join();
        }
    }

    auto t_end = hrc::now();
    return ms_between(t_start, t_end);
}

// ============================================================
// Protocol 2: Distributed Prove Implementation
// ============================================================

DistributedEvalProofV2 distributed_prove_v2(
    const BrakedownCodeGR& code,
    const DistributedCommitmentV2& comm,
    const ZZ_pE* poly_coeffs,
    long n,
    const std::vector<ZZ_pE>& q1,
    const std::vector<ZZ_pE>& q2,
    int num_workers,
    DistributedProveTiming& timing,
    std::vector<WorkerStats>& worker_stats)
{
    auto t_total_start = hrc::now();

    DistributedEvalProofV2 proof;
    long row_len = code.row_len;
    long num_rows = comm.num_rows;
    long cw_len = comm.codeword_len;
    long num_prox = effective_num_prox_dist(code, num_rows);

    long k = NumBits(ZZ_p::modulus()) - 1;
    long degree = ZZ_pE::degree();

    if (code.is_small_ring) {
        std::cerr << "distributed_prove_v2: small ring not yet supported" << std::endl;
        return proof;
    }

    // ============================================================
    // Phase 1: Distribute data to workers (same as V1)
    // ============================================================
    auto t1 = hrc::now();

    std::vector<ZZ_pE> input_matrix(num_rows * row_len);
    for (long i = 0; i < num_rows; i++) {
        for (long j = 0; j < row_len; j++) {
            long idx = i * row_len + j;
            if (idx < n) {
                input_matrix[i * row_len + j] = poly_coeffs[idx];
            } else {
                clear(input_matrix[i * row_len + j]);
            }
        }
    }

    long rows_per_worker = (num_rows + num_workers - 1) / num_workers;

    proof.prox_coeffs.resize(num_prox);
    for (long t = 0; t < num_prox; t++) {
        proof.prox_coeffs[t].resize(num_rows);
        for (long i = 0; i < num_rows; i++) {
            proof.prox_coeffs[t][i] = random_ZZ_pE();
        }
    }

    auto t2 = hrc::now();
    timing.distribute_ms = ms_between(t1, t2);

    // ============================================================
    // Phase 2: Workers compute partial combines (same as V1)
    // ============================================================
    auto t3 = hrc::now();

    long num_combines = num_prox + 1;
    std::vector<std::vector<ZZ_pE>> all_combined(num_combines);

    std::vector<const ZZ_pE*> all_coeffs(num_combines);
    for (long t = 0; t < num_prox; t++) {
        all_coeffs[t] = proof.prox_coeffs[t].data();
    }
    all_coeffs[num_prox] = q1.data();

    std::vector<BatchCombineTask> tasks(num_workers);
    std::vector<std::vector<ZZ_pE>> all_partial_results(num_combines * num_workers);

    for (long c = 0; c < num_combines; c++) {
        for (int w = 0; w < num_workers; w++) {
            all_partial_results[c * num_workers + w].resize(row_len);
            for (long j = 0; j < row_len; j++) {
                clear(all_partial_results[c * num_workers + w][j]);
            }
        }
    }

    for (int w = 0; w < num_workers; w++) {
        long start = w * rows_per_worker;
        long end = std::min(start + rows_per_worker, num_rows);
        if (start >= num_rows) {
            start = end = num_rows;
        }

        tasks[w].poly_rows = &input_matrix[start * row_len];
        tasks[w].all_coeffs = &all_coeffs;
        tasks[w].all_partial_results = &all_partial_results;
        tasks[w].start_row = start;
        tasks[w].end_row = end;
        tasks[w].row_len = row_len;
        tasks[w].n = n;
        tasks[w].k = k;
        tasks[w].degree = degree;
        tasks[w].worker_id = w;
        tasks[w].num_workers = num_workers;
        tasks[w].num_combines = num_combines;
        tasks[w].global_start = &t3;
    }

    std::vector<std::thread> threads;
    for (int w = 0; w < num_workers; w++) {
        if (tasks[w].start_row < tasks[w].end_row) {
            threads.emplace_back(worker_batch_combine, &tasks[w]);
        }
    }
    for (auto& t : threads) {
        t.join();
    }

    auto t4 = hrc::now();
    timing.combine_ms = ms_between(t3, t4);

    // ============================================================
    // Phase 2.5: Tree Reduction for aggregating partial results
    // Instead of master collecting all, use tree-based reduction
    // ============================================================
    timing.tree_reduce_ms = tree_reduce_partial_results(
        all_partial_results, num_workers, num_combines, row_len, k, degree);

    // After tree reduction, results are in worker 0's slots
    for (long c = 0; c < num_combines; c++) {
        all_combined[c] = std::move(all_partial_results[c * num_workers + 0]);
    }

    worker_stats.resize(num_workers);
    for (int w = 0; w < num_workers; w++) {
        worker_stats[w].worker_id = w;
        worker_stats[w].rows_processed = tasks[w].end_row - tasks[w].start_row;
        worker_stats[w].encode_ms = 0;
        worker_stats[w].encode_start_ms = 0;
        worker_stats[w].encode_end_ms = 0;
        worker_stats[w].combine_start_ms = tasks[w].start_offset_ms.load();
        worker_stats[w].combine_end_ms = tasks[w].end_offset_ms.load();
        worker_stats[w].combine_ms = worker_stats[w].combine_end_ms - worker_stats[w].combine_start_ms;
    }

    // ============================================================
    // Phase 3: Collect combined rows and compute eval
    // ============================================================
    auto t5 = hrc::now();

    for (long c = 0; c < num_prox; c++) {
        proof.combined_rows.push_back(all_combined[c]);
    }
    proof.eval_value = inner_product_dist(all_combined[num_prox].data(), q2.data(), row_len);
    proof.combined_rows.push_back(all_combined[num_prox]);

    auto t6 = hrc::now();
    timing.collect_ms = ms_between(t5, t6);

    // ============================================================
    // Phase 4: Column openings with distributed Merkle proof
    // For Protocol 2, we need:
    //   - Column values from all workers
    //   - Merkle path within the responsible worker's local tree
    //   - All local roots r_j for final verification
    // ============================================================
    auto t7 = hrc::now();

    long num_open = code.num_col_open;
    long cols_per_worker = (cw_len + comm.num_workers - 1) / comm.num_workers;

    proof.column_indices.resize(num_open);
    proof.column_items.resize(num_open);
    proof.merkle_paths.resize(num_open);

    // Store all local roots (needed for verification)
    proof.all_local_roots = comm.local_roots;

    for (long t = 0; t < num_open; t++) {
        long col_idx = RandomBnd(cw_len);
        proof.column_indices[t] = col_idx;

        // Collect full column from encoded_rows
        std::vector<ZZ_pE> column(num_rows);
        for (long row = 0; row < num_rows; row++) {
            column[row] = comm.encoded_rows[row * cw_len + col_idx];
        }
        proof.column_items[t] = column;

        // Find which worker owns this column
        int owner_worker = col_idx / cols_per_worker;
        if (owner_worker >= comm.num_workers) {
            owner_worker = comm.num_workers - 1;
        }

        // Get local index within that worker's tree
        long local_col_idx = col_idx - owner_worker * cols_per_worker;

        // Get Merkle path within that worker's local tree
        proof.merkle_paths[t] = get_merkle_path(comm.local_trees[owner_worker], local_col_idx);
    }

    auto t8 = hrc::now();
    timing.column_open_ms = ms_between(t7, t8);

    auto t_total_end = hrc::now();
    timing.total_ms = ms_between(t_total_start, t_total_end);

    return proof;
}

// Convenience overload
DistributedEvalProofV2 distributed_prove_v2(
    const BrakedownCodeGR& code,
    const DistributedCommitmentV2& comm,
    const ZZ_pE* poly_coeffs,
    long n,
    const std::vector<ZZ_pE>& q1,
    const std::vector<ZZ_pE>& q2,
    int num_workers)
{
    DistributedProveTiming timing;
    std::vector<WorkerStats> stats;
    return distributed_prove_v2(code, comm, poly_coeffs, n, q1, q2, num_workers, timing, stats);
}

// ============================================================
// Protocol 2: Distributed Verify Implementation
// ============================================================

bool distributed_verify_v2(
    const BrakedownCodeGR& code,
    const DistributedCommitmentV2& comm,
    const DistributedEvalProofV2& proof,
    const std::vector<ZZ_pE>& q1,
    const std::vector<ZZ_pE>& q2)
{
    long row_len = code.row_len;
    long cw_len = code.codeword_len;
    long num_rows = comm.num_rows;
    long num_prox = effective_num_prox_dist(code, num_rows);
    int num_workers = comm.num_workers;
    long rows_per_worker = (num_rows + num_workers - 1) / num_workers;
    long cols_per_worker = (cw_len + num_workers - 1) / num_workers;

    if (code.is_small_ring) {
        std::cerr << "distributed_verify_v2: small ring not yet supported" << std::endl;
        return false;
    }

    // ============================================================
    // Step 1: Verify final root = H(r_0, ..., r_{M-1})
    // ============================================================
    std::vector<unsigned char> computed_root = hash_M(proof.all_local_roots);
    if (computed_root != comm.root) {
        std::cerr << "distributed_verify_v2: root mismatch" << std::endl;
        return false;
    }

    // ============================================================
    // Step 2: Verify evaluation
    // ============================================================
    const auto& combined_q1 = proof.combined_rows.back();
    ZZ_pE computed_eval = inner_product_dist(combined_q1.data(), q2.data(), row_len);
    if (computed_eval != proof.eval_value) {
        std::cerr << "distributed_verify_v2: eval mismatch" << std::endl;
        return false;
    }

    // ============================================================
    // Step 3: Re-encode combined rows
    // ============================================================
    std::vector<std::vector<ZZ_pE>> encoded_combined(proof.combined_rows.size());
    for (size_t t = 0; t < proof.combined_rows.size(); t++) {
        encoded_combined[t].resize(cw_len);
        for (long j = 0; j < row_len; j++) {
            encoded_combined[t][j] = proof.combined_rows[t][j];
        }
        for (long j = row_len; j < cw_len; j++) {
            clear(encoded_combined[t][j]);
        }
        brakedown_encode(code, encoded_combined[t].data());
    }

    // ============================================================
    // Step 4: Verify column openings with distributed Merkle proof
    // ============================================================

    // Build coefficient vectors for proximity tests
    std::vector<std::vector<ZZ_pE>> coeffs(proof.combined_rows.size());
    for (long t = 0; t < num_prox; t++) {
        coeffs[t] = proof.prox_coeffs[t];
    }
    coeffs[num_prox] = q1;

    for (size_t idx = 0; idx < proof.column_indices.size(); idx++) {
        long col = proof.column_indices[idx];
        const auto& column = proof.column_items[idx];
        const auto& merkle_path = proof.merkle_paths[idx];

        // Step 4a: Recompute local hashes h^(i)[col] for each worker
        std::vector<std::vector<unsigned char>> local_hashes(num_workers);
        for (int w = 0; w < num_workers; w++) {
            long start = w * rows_per_worker;
            long end = std::min(start + rows_per_worker, num_rows);
            long worker_rows = end - start;

            // Extract this worker's portion of the column
            std::vector<ZZ_pE> col_portion(worker_rows);
            for (long i = 0; i < worker_rows; i++) {
                col_portion[i] = column[start + i];
            }
            local_hashes[w] = hash_column(col_portion.data(), worker_rows);
        }

        // Step 4b: Aggregate h[col] = H_M(h^(0)[col], ..., h^(M-1)[col])
        std::vector<unsigned char> global_hash = hash_M(local_hashes);

        // Step 4c: Find owner worker and verify Merkle path
        int owner_worker = col / cols_per_worker;
        if (owner_worker >= num_workers) {
            owner_worker = num_workers - 1;
        }
        long local_col_idx = col - owner_worker * cols_per_worker;

        // Verify path from h[col] to r_{owner_worker}
        if (!verify_merkle_path(global_hash, local_col_idx, merkle_path,
                                proof.all_local_roots[owner_worker])) {
            std::cerr << "distributed_verify_v2: Merkle path failed for column " << col << std::endl;
            return false;
        }

        // Step 4d: Column consistency check
        // Check: encoded_combined[s][col] == sum of coeffs[s][i] * column[i]
        for (size_t s = 0; s < proof.combined_rows.size(); s++) {
            ZZ_pE expected;
            clear(expected);
            for (long i = 0; i < num_rows; i++) {
                expected += column[i] * coeffs[s][i];
            }
            if (encoded_combined[s][col] != expected) {
                std::cerr << "distributed_verify_v2: column consistency failed" << std::endl;
                return false;
            }
        }
    }

    return true;
}

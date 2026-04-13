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

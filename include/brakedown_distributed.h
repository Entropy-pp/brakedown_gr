#ifndef BRAKEDOWN_DISTRIBUTED_H
#define BRAKEDOWN_DISTRIBUTED_H

#include <NTL/ZZ_pE.h>
#include <vector>
#include <thread>
#include <mutex>
#include <atomic>
#include <brakedown_code_gr.h>
#include <brakedown_pcs_gr.h>

using namespace NTL;

// ============================================================
// Distributed Brakedown PCS - Timing structures
// ============================================================

struct DistributedCommitTiming {
    double distribute_ms;      // Master: distribute rows to workers
    double encode_ms;          // Workers: parallel encoding (max across workers)
    double collect_ms;         // Master: collect encoded rows
    double column_hash_ms;     // Master: compute column hashes
    double merkle_build_ms;    // Master: build Merkle tree
    double total_ms;
};

struct DistributedProveTiming {
    double distribute_ms;      // Master: distribute data to workers
    double combine_ms;         // Workers: parallel combine (max across workers)
    double collect_ms;         // Master: collect and aggregate
    double column_open_ms;     // Master: column openings
    double total_ms;
};

struct WorkerStats {
    int worker_id;
    long rows_processed;
    double encode_ms;          // Actual encoding time (end - start)
    double combine_ms;         // Actual combine time (end - start)
    double encode_start_ms;    // When worker started (relative to global start)
    double encode_end_ms;      // When worker finished (relative to global start)
    double combine_start_ms;
    double combine_end_ms;
};

// ============================================================
// Distributed Commit
// ============================================================

// Distributed commit: workers encode rows in parallel, master builds Merkle tree
BrakedownCommitment distributed_commit(
    const BrakedownCodeGR& code,
    const ZZ_pE* poly_coeffs,
    long n,
    long num_rows,
    int num_workers,
    DistributedCommitTiming& timing,
    std::vector<WorkerStats>& worker_stats);

// Convenience overload without timing
BrakedownCommitment distributed_commit(
    const BrakedownCodeGR& code,
    const ZZ_pE* poly_coeffs,
    long n,
    long num_rows,
    int num_workers);

// ============================================================
// Distributed Prove
// ============================================================

// Distributed prove: workers compute partial combines in parallel
BrakedownEvalProof distributed_prove(
    const BrakedownCodeGR& code,
    const BrakedownCommitment& comm,
    const ZZ_pE* poly_coeffs,
    long n,
    const std::vector<ZZ_pE>& q1,
    const std::vector<ZZ_pE>& q2,
    int num_workers,
    DistributedProveTiming& timing,
    std::vector<WorkerStats>& worker_stats);

// Convenience overload without timing
BrakedownEvalProof distributed_prove(
    const BrakedownCodeGR& code,
    const BrakedownCommitment& comm,
    const ZZ_pE* poly_coeffs,
    long n,
    const std::vector<ZZ_pE>& q1,
    const std::vector<ZZ_pE>& q2,
    int num_workers);

#endif

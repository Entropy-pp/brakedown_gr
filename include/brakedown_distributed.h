#ifndef BRAKEDOWN_DISTRIBUTED_H
#define BRAKEDOWN_DISTRIBUTED_H

#include <NTL/ZZ_pE.h>
#include <vector>
#include <thread>
#include <mutex>
#include <atomic>
#include <brakedown_code_gr.h>
#include <brakedown_pcs_gr.h>
#include <merkle.h>

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
    double tree_reduce_ms;     // Tree reduction aggregation
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

// ============================================================
// Protocol 2: Distributed Commitment with Distributed Merkle Tree
// ============================================================

// Timing for new distributed commit protocol
struct DistributedCommitTimingV2 {
    double local_encode_ms;        // Workers: parallel encoding
    double local_hash_ms;          // Workers: local column hashing
    double hash_exchange_ms;       // All-to-all hash exchange
    double global_hash_ms;         // Workers: H_M aggregation
    double local_merkle_ms;        // Workers: local Merkle tree construction
    double global_commit_ms;       // Master: final root aggregation
    double total_ms;
};

// Commitment structure for Protocol 2
// Stores local roots and metadata for distributed verification
struct DistributedCommitmentV2 {
    long num_rows;
    long codeword_len;
    int num_workers;

    // Each worker's local Merkle root
    std::vector<std::vector<unsigned char>> local_roots;  // r_0, r_1, ..., r_{M-1}

    // Final commitment: MC(Û) = H(r_0, r_1, ..., r_{M-1})
    std::vector<unsigned char> root;

    // Local Merkle trees (one per worker, for generating proofs)
    std::vector<MerkleTree> local_trees;

    // Global column hashes h[k] (for proof generation)
    std::vector<std::vector<unsigned char>> global_column_hashes;

    // Encoded rows (still needed for column openings in prove phase)
    std::vector<ZZ_pE> encoded_rows;
};

// Proof structure for Protocol 2
struct DistributedEvalProofV2 {
    // Evaluation value
    ZZ_pE eval_value;

    // Combined rows for proximity test
    std::vector<std::vector<ZZ_pE>> combined_rows;
    std::vector<std::vector<ZZ_pE>> prox_coeffs;

    // Column openings
    std::vector<long> column_indices;
    std::vector<std::vector<ZZ_pE>> column_items;          // Full column values
    std::vector<std::vector<std::vector<unsigned char>>> merkle_paths;  // Path within local tree
    std::vector<std::vector<unsigned char>> all_local_roots;  // All r_j for verification
};

// Worker statistics for V2
struct WorkerStatsV2 {
    int worker_id;
    long rows_processed;
    long cols_processed;           // Columns this worker is responsible for
    double encode_ms;
    double local_hash_ms;
    double global_hash_ms;
    double local_merkle_ms;
};

// ============================================================
// Protocol 2: Distributed Commit
// ============================================================

DistributedCommitmentV2 distributed_commit_v2(
    const BrakedownCodeGR& code,
    const ZZ_pE* poly_coeffs,
    long n,
    long num_rows,
    int num_workers,
    DistributedCommitTimingV2& timing,
    std::vector<WorkerStatsV2>& worker_stats);

// Convenience overload without timing
DistributedCommitmentV2 distributed_commit_v2(
    const BrakedownCodeGR& code,
    const ZZ_pE* poly_coeffs,
    long n,
    long num_rows,
    int num_workers);

// ============================================================
// Protocol 2: Distributed Prove
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
    std::vector<WorkerStats>& worker_stats);

// Convenience overload without timing
DistributedEvalProofV2 distributed_prove_v2(
    const BrakedownCodeGR& code,
    const DistributedCommitmentV2& comm,
    const ZZ_pE* poly_coeffs,
    long n,
    const std::vector<ZZ_pE>& q1,
    const std::vector<ZZ_pE>& q2,
    int num_workers);

// ============================================================
// Protocol 2: Verify
// ============================================================

bool distributed_verify_v2(
    const BrakedownCodeGR& code,
    const DistributedCommitmentV2& comm,
    const DistributedEvalProofV2& proof,
    const std::vector<ZZ_pE>& q1,
    const std::vector<ZZ_pE>& q2);

#endif

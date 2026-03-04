#ifndef BRAKEDOWN_PCS_GR_H
#define BRAKEDOWN_PCS_GR_H

#include <NTL/ZZ_pE.h>
#include <NTL/ZZ_pEX.h>
#include <vector>
#include <brakedown_code_gr.h>
#include <merkle.h>

using namespace NTL;

// ============================================================
// Brakedown Polynomial Commitment Scheme over Galois Ring
// Protocol: [GLSTW21] adapted to GR(2^k, d)
// ============================================================

struct BrakedownCommitment {
    std::vector<unsigned char> root;
    std::vector<ZZ_pE> encoded_rows; // num_rows * codeword_len flat
    MerkleTree tree;
    long num_rows;
    long codeword_len;
};

struct BrakedownEvalProof {
    // Combined rows: first num_prox_test from proximity testing,
    // then 1 from evaluation consistency (q1-combined).
    // Each vector has length = row_len.
    std::vector<std::vector<ZZ_pE>> combined_rows;

    // Proximity test coefficients so verifier can reproduce the checks.
    // prox_coeffs[t] has length = num_rows.
    std::vector<std::vector<ZZ_pE>> prox_coeffs;

    // Column openings
    std::vector<long> column_indices;
    std::vector<std::vector<ZZ_pE>> column_items;
    std::vector<std::vector<std::vector<unsigned char>>> merkle_paths;

    // Evaluation value
    ZZ_pE eval_value;
};

BrakedownCommitment brakedown_commit(
    const BrakedownCodeGR& code,
    const ZZ_pE* poly_coeffs,
    long n,
    long num_rows);

BrakedownEvalProof brakedown_prove(
    const BrakedownCodeGR& code,
    const BrakedownCommitment& comm,
    const ZZ_pE* poly_coeffs,
    long n,
    const std::vector<ZZ_pE>& q1,
    const std::vector<ZZ_pE>& q2);

bool brakedown_verify(
    const BrakedownCodeGR& code,
    const BrakedownCommitment& comm,
    const BrakedownEvalProof& proof,
    const std::vector<ZZ_pE>& q1,
    const std::vector<ZZ_pE>& q2);

#endif
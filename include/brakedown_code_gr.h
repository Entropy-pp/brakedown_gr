#ifndef BRAKEDOWN_CODE_GR_H
#define BRAKEDOWN_CODE_GR_H

#include <NTL/ZZ_pE.h>
#include <NTL/ZZ_pEX.h>
#include <vector>
#include <brakedown_params.h>
#include <sparse_matrix_gr.h>

using namespace NTL;

struct BrakedownCodeGR {
    long row_len;       // message length (= number of columns in the matrix view)
    long codeword_len;  // total codeword length after encoding
    long num_col_open;  // number of column openings (Θ(λ))
    long num_prox_test; // number of proximity tests
    BrakedownSpec spec;
    std::vector<SparseMatrixGR> a_mats; // A_1, ..., A_L
    std::vector<SparseMatrixGR> b_mats; // B_1, ..., B_L
};

// Setup: create the code for a given message length (row_len)
BrakedownCodeGR brakedown_code_setup(long row_len, long k, long degree);

// Encode in-place: codeword[0..row_len-1] = message, fills codeword[row_len..codeword_len-1]
// codeword must have length = code.codeword_len
void brakedown_encode(const BrakedownCodeGR& code, ZZ_pE* codeword);

#endif
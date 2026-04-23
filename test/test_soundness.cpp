/**
 * test_soundness.cpp - Test soundness of Brakedown PCS
 *
 * Verifies that:
 * 1. Correct evaluation passes verification
 * 2. Incorrect evaluation fails verification
 *
 * Compile:
 *   make test_soundness
 *
 * Run:
 *   ./out/test_soundness
 */

#include <iostream>
#include <iomanip>
#include <vector>
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

static void initGR(long k, long degree) {
    ZZ modulus = ZZ(1) << k;
    ZZ_p::init(modulus);
    ZZ_pX P = primitiveIrredPoly(degree);
    ZZ_pE::init(P);
}

// Compute correct evaluation: sum of poly[i] * tensor(q1, q2)[i]
static ZZ_pE compute_correct_eval(
    const vector<ZZ_pE>& poly,
    const vector<ZZ_pE>& q1,
    const vector<ZZ_pE>& q2,
    long num_rows,
    long row_len)
{
    ZZ_pE result;
    clear(result);

    for (long i = 0; i < num_rows; i++) {
        for (long j = 0; j < row_len; j++) {
            long idx = i * row_len + j;
            if (idx < (long)poly.size()) {
                result += poly[idx] * q1[i] * q2[j];
            }
        }
    }
    return result;
}

int main() {
    cout << "\n";
    cout << "##############################################################\n";
    cout << "  Soundness Test: Brakedown PCS\n";
    cout << "##############################################################\n";

    // Parameters
    long s = 2;
    long base_r = 128;
    long nv = 10;  // Small for quick test

    ZZ_p::init(ZZ(1) << s);

    long n = 1L << nv;
    long half = nv / 2;
    long num_rows = 1L << half;
    long row_len = 1L << (nv - half);

    cout << "\n  Parameters:\n";
    cout << "    n = 2^" << nv << " = " << n << "\n";
    cout << "    num_rows = " << num_rows << ", row_len = " << row_len << "\n";

    initGR(s, base_r);
    BrakedownCodeGR code = brakedown_code_setup(row_len, s, base_r);

    // Generate random polynomial
    vector<ZZ_pE> poly(n);
    for (long i = 0; i < n; i++) {
        poly[i] = random_ZZ_pE();
    }

    // Generate random evaluation point
    vector<ZZ_pE> q1(num_rows), q2(row_len);
    for (long i = 0; i < num_rows; i++) q1[i] = random_ZZ_pE();
    for (long j = 0; j < row_len; j++) q2[j] = random_ZZ_pE();

    // Compute correct evaluation
    ZZ_pE correct_eval = compute_correct_eval(poly, q1, q2, num_rows, row_len);

    cout << "\n";
    cout << "============================================================\n";
    cout << "  Test 1: Correct evaluation should PASS\n";
    cout << "============================================================\n";

    // Commit
    auto comm = brakedown_commit(code, poly.data(), n, num_rows);
    cout << "  Commitment generated.\n";

    // Prove with correct evaluation
    auto proof = brakedown_prove(code, comm, poly.data(), n, q1, q2);
    cout << "  Proof generated.\n";
    cout << "  Claimed eval = " << proof.eval_value << "\n";
    cout << "  Correct eval = " << correct_eval << "\n";

    // Verify
    bool result1 = brakedown_verify(code, comm, proof, q1, q2);
    cout << "  Verification: " << (result1 ? "PASS" : "FAIL") << "\n";

    cout << "\n";
    cout << "============================================================\n";
    cout << "  Test 2: Wrong evaluation should FAIL\n";
    cout << "============================================================\n";

    // Create a proof with wrong evaluation
    BrakedownEvalProof bad_proof = proof;
    bad_proof.eval_value = proof.eval_value + to_ZZ_pE(1);  // Add 1 to make it wrong

    cout << "  Original eval = " << proof.eval_value << "\n";
    cout << "  Tampered eval = " << bad_proof.eval_value << "\n";

    // Verify should fail
    bool result2 = brakedown_verify(code, comm, bad_proof, q1, q2);
    cout << "  Verification: " << (result2 ? "PASS (BUG!)" : "FAIL (expected)") << "\n";

    cout << "\n";
    cout << "============================================================\n";
    cout << "  Test 3: Wrong combined_row should FAIL\n";
    cout << "============================================================\n";

    // Create a proof with tampered combined_row
    BrakedownEvalProof bad_proof2 = proof;
    if (!bad_proof2.combined_rows.empty() && !bad_proof2.combined_rows.back().empty()) {
        bad_proof2.combined_rows.back()[0] += to_ZZ_pE(1);  // Tamper first element
    }

    cout << "  Tampered combined_row[0]\n";

    bool result3 = brakedown_verify(code, comm, bad_proof2, q1, q2);
    cout << "  Verification: " << (result3 ? "PASS (BUG!)" : "FAIL (expected)") << "\n";

    cout << "\n";
    cout << "============================================================\n";
    cout << "  Test 4: Wrong column opening should FAIL\n";
    cout << "============================================================\n";

    // Create a proof with tampered column value
    BrakedownEvalProof bad_proof3 = proof;
    if (!bad_proof3.column_items.empty() && !bad_proof3.column_items[0].empty()) {
        bad_proof3.column_items[0][0] += to_ZZ_pE(1);  // Tamper first column item
    }

    cout << "  Tampered column_items[0][0]\n";

    bool result4 = brakedown_verify(code, comm, bad_proof3, q1, q2);
    cout << "  Verification: " << (result4 ? "PASS (BUG!)" : "FAIL (expected)") << "\n";

    cout << "\n";
    cout << "##############################################################\n";
    cout << "  Summary:\n";
    cout << "    Test 1 (correct eval):    " << (result1 ? "PASS" : "FAIL") << "\n";
    cout << "    Test 2 (wrong eval):      " << (!result2 ? "PASS" : "FAIL") << "\n";
    cout << "    Test 3 (wrong combined):  " << (!result3 ? "PASS" : "FAIL") << "\n";
    cout << "    Test 4 (wrong column):    " << (!result4 ? "PASS" : "FAIL") << "\n";
    cout << "##############################################################\n";

    bool all_pass = result1 && !result2 && !result3 && !result4;
    cout << "\n  Overall: " << (all_pass ? "ALL TESTS PASSED" : "SOME TESTS FAILED") << "\n\n";

    return all_pass ? 0 : 1;
}

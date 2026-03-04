/**
 * brakedown_gr_test.cpp — Fixed version
 */

#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pE.h>
#include <iostream>
#include <iomanip>
#include <ctime>
#include <cassert>
#include <vector>

#include <gr.h>
#include <brakedown_params.h>
#include <brakedown_code_gr.h>
#include <brakedown_pcs_gr.h>

using namespace std;
using namespace NTL;

void initGR(long k, long degree) {
    ZZ modulus = ZZ(1) << k;
    ZZ_p::init(modulus);
    ZZ_pX P = primitiveIrredPoly(degree);
    ZZ_pE::init(P);
}

// ============================================================
// Test 1: Encode linearity test
// ============================================================
bool test_linearity(long row_len, long k, long degree) {
    cout << "  [Linearity Test] row_len=" << row_len << " k=" << k << " d=" << degree << "...\n";
    initGR(k, degree);

    BrakedownCodeGR code = brakedown_code_setup(row_len, k, degree);
    long cw_len = code.codeword_len;

    vector<ZZ_pE> u(cw_len), v(cw_len), uv(cw_len);
    for (long i = 0; i < row_len; i++) {
        u[i] = random_ZZ_pE();
        v[i] = random_ZZ_pE();
    }
    for (long i = row_len; i < cw_len; i++) { clear(u[i]); clear(v[i]); }

    ZZ_pE alpha = randomInvertible();
    for (long i = 0; i < row_len; i++) uv[i] = alpha * u[i] + v[i];
    for (long i = row_len; i < cw_len; i++) clear(uv[i]);

    brakedown_encode(code, u.data());
    brakedown_encode(code, v.data());
    brakedown_encode(code, uv.data());

    bool ok = true;
    for (long i = 0; i < cw_len; i++) {
        if (uv[i] != alpha * u[i] + v[i]) {
            cout << "    FAIL at position " << i << "\n";
            ok = false;
            break;
        }
    }
    if (ok) cout << "    PASS\n";
    return ok;
}

// ============================================================
// Test 2: Full PCS test (Commit + Prove + Verify)
// FIX: When num_rows == 1, q1 must be [1] (matching Rust behavior)
// ============================================================
bool test_full_pcs(long n, long k, long degree) {
    cout << "  [Full PCS Test] n=" << n << " k=" << k << " d=" << degree << "...\n";
    initGR(k, degree);

    vector<ZZ_pE> poly(n);
    for (long i = 0; i < n; i++) poly[i] = random_ZZ_pE();

    long row_len = 1;
    while (row_len * 2 <= n && row_len * 2 <= 256) row_len *= 2;
    long num_rows = n / row_len;

    BrakedownCodeGR code = brakedown_code_setup(row_len, k, degree);

    clock_t t = clock();
    BrakedownCommitment comm = brakedown_commit(code, poly.data(), n, num_rows);
    double commit_time = double(clock() - t) / CLOCKS_PER_SEC;

    // FIX: Build q1, q2 correctly.
    //   When num_rows == 1: q1 = [1], q2 = random (length row_len)
    //   When num_rows >  1: both random
    vector<ZZ_pE> q1(num_rows), q2(row_len);
    if (num_rows == 1) {
        set(q1[0]); // q1[0] = 1
    } else {
        for (long i = 0; i < num_rows; i++) q1[i] = random_ZZ_pE();
    }
    for (long j = 0; j < row_len; j++) q2[j] = random_ZZ_pE();

    t = clock();
    BrakedownEvalProof proof = brakedown_prove(code, comm, poly.data(), n, q1, q2);
    double prove_time = double(clock() - t) / CLOCKS_PER_SEC;

    // Compute expected evaluation manually
    ZZ_pE expected_eval;
    clear(expected_eval);
    for (long i = 0; i < num_rows; i++) {
        for (long j = 0; j < row_len; j++) {
            long idx = i * row_len + j;
            if (idx < n) expected_eval += q1[i] * poly[idx] * q2[j];
        }
    }

    if (expected_eval != proof.eval_value) {
        cout << "    FAIL: evaluation mismatch\n";
        return false;
    }

    t = clock();
    bool ok = brakedown_verify(code, comm, proof, q1, q2);
    double verify_time = double(clock() - t) / CLOCKS_PER_SEC;

    cout << "    Commit: " << fixed << setprecision(4) << commit_time << "s"
         << "  Prove: " << prove_time << "s"
         << "  Verify: " << verify_time << "s"
         << "  " << (ok ? "PASS" : "FAIL") << "\n";
    return ok;
}

// ============================================================
// Benchmark
// ============================================================
void benchmark(long k, long degree) {
    cout << "\n============================================================\n";
    cout << " Brakedown over GR(2^" << k << ", " << degree << ") Benchmark\n";
    cout << "============================================================\n";

    initGR(k, degree);

    cout << "\n" << left
         << setw(8) << "n"
         << setw(10) << "rows"
         << setw(10) << "row_len"
         << setw(10) << "cw_len"
         << setw(12) << "Commit(s)"
         << setw(12) << "Prove(s)"
         << setw(12) << "Verify(s)"
         << setw(12) << "Total(s)"
         << "\n";
    cout << string(96, '-') << "\n";

    long sizes[] = {64, 128, 256, 512, 1024, 2048, 4096};
    for (long n : sizes) {
        long row_len = 1;
        while (row_len * 2 <= n && row_len * 2 <= 512) row_len *= 2;
        long num_rows = n / row_len;
        if (num_rows < 1) num_rows = 1;

        initGR(k, degree);

        vector<ZZ_pE> poly(n);
        for (long i = 0; i < n; i++) poly[i] = random_ZZ_pE();

        BrakedownCodeGR code = brakedown_code_setup(row_len, k, degree);

        clock_t t = clock();
        BrakedownCommitment comm = brakedown_commit(code, poly.data(), n, num_rows);
        double commit_time = double(clock() - t) / CLOCKS_PER_SEC;

        // FIX: q1 = [1] when num_rows == 1
        vector<ZZ_pE> q1(num_rows), q2(row_len);
        if (num_rows == 1) {
            set(q1[0]);
        } else {
            for (long i = 0; i < num_rows; i++) q1[i] = random_ZZ_pE();
        }
        for (long j = 0; j < row_len; j++) q2[j] = random_ZZ_pE();

        t = clock();
        BrakedownEvalProof proof = brakedown_prove(code, comm, poly.data(), n, q1, q2);
        double prove_time = double(clock() - t) / CLOCKS_PER_SEC;

        t = clock();
        bool ok = brakedown_verify(code, comm, proof, q1, q2);
        double verify_time = double(clock() - t) / CLOCKS_PER_SEC;

        double total = commit_time + prove_time + verify_time;

        cout << left << fixed << setprecision(4)
             << setw(8) << n
             << setw(10) << num_rows
             << setw(10) << row_len
             << setw(10) << code.codeword_len
             << setw(12) << commit_time
             << setw(12) << prove_time
             << setw(12) << verify_time
             << setw(12) << total
             << (ok ? "" : " FAIL!")
             << "\n";
    }
}

// ============================================================
// MAIN
// ============================================================
int main() {
    cout << "##############################################################\n";
    cout << "  Brakedown PCS over Galois Ring — Test Suite\n";
    cout << "##############################################################\n";

    long k = 64;
    long degree = 8;

    cout << "\n>>> Test 1: Encode Linearity <<<\n";
    test_linearity(32, k, degree);
    test_linearity(64, k, degree);
    test_linearity(128, k, degree);

    cout << "\n>>> Test 2: Full PCS (Commit + Prove + Verify) <<<\n";
    test_full_pcs(64, k, degree);
    test_full_pcs(256, k, degree);
    test_full_pcs(1024, k, degree);

    benchmark(k, degree);

    cout << "\n============================================================\n";
    cout << " Comparison with Rinocchio (from previous benchmark)\n";
    cout << "============================================================\n";
    cout << "\n";
    cout << " Rinocchio 1024 gates: Setup=9.09s  Prove=7.11s  Verify=0.01s  Total=34.76s\n";
    cout << " (Brakedown results shown above — compare Commit+Prove vs Setup+CRS+Prove)\n";
    cout << "\n Key advantages of Brakedown over GR:\n";
    cout << "   - NO trusted setup (transparent)\n";
    cout << "   - Hash-based (no modular exponentiation)\n";
    cout << "   - Prover should be 10-100x faster\n";
    cout << "   - Verify may be slightly slower (linear-time re-encoding)\n";

    cout << "\n##############################################################\n";
    cout << "  All tests complete.\n";
    cout << "##############################################################\n";

    return 0;
}
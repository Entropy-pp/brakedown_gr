/**
 * fair_benchmark.cpp
 *
 * 公平对比 Brakedown PCS vs Rinocchio：
 *   - Rinocchio: 使用 n×n×n 矩阵乘法电路（与原始 rinocchio_test.cpp 一致）
 *   - Brakedown: 对相同规模的多项式做 Commit/Prove/Verify
 *   - 两者在同一进程中运行，使用相同的 GR 参数，计时公平
 *   - 每个规模运行 NUM_TRIALS 次取平均值，数据更客观
 *
 * 编译:
 *   g++ -O2 -std=c++11 -I./include test/fair_benchmark.cpp src/*.cpp \
 *       -o out/fair_benchmark -pthread -lntl -lgmp -lm
 *
 * 运行:
 *   ./out/fair_benchmark
 */

#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pX.h>
#include <NTL/ZZ_pE.h>
#include <NTL/ZZ_pEX.h>
#include <NTL/vec_ZZ_pE.h>
#include <vector>
#include <iostream>
#include <iomanip>
#include <cassert>
#include <ctime>
#include <cmath>
#include <algorithm>

// Rinocchio 相关
#include <gr.h>
#include <circuit.h>
#include <qrp.h>
#include <setup.h>
#include <rinocchio.h>

// Brakedown 相关
#include <brakedown_params.h>
#include <brakedown_code_gr.h>
#include <brakedown_pcs_gr.h>

using namespace std;
using namespace NTL;

// ============================================================
//  每个规模的运行次数（取平均值）
// ============================================================
static const int NUM_TRIALS = 10;

// ============================================================
//  初始化 GR(2^k, d)
// ============================================================
void initGR(long k, long degree) {
    ZZ modulus = ZZ(1) << k;
    ZZ_p::init(modulus);
    ZZ_pX P = primitiveIrredPoly(degree);
    ZZ_pE::init(P);
}

// ============================================================
//  自动扩展次数（与 rinocchio_test.cpp 完全一致）
// ============================================================
long autoExtensionDegree(long numMultGates) {
    long d = 2;
    if (numMultGates > 4)     d = 3;
    if (numMultGates > 8)     d = 4;
    if (numMultGates > 16)    d = 5;
    if (numMultGates > 32)    d = 6;
    if (numMultGates > 64)    d = 7;
    if (numMultGates > 128)   d = 8;
    if (numMultGates > 256)   d = 9;
    if (numMultGates > 512)   d = 10;
    if (numMultGates > 1024)  d = 11;
    if (numMultGates > 2048)  d = 12;
    if (numMultGates > 4096)  d = 13;
    if (numMultGates > 8192)  d = 14;
    if (numMultGates > 16384) d = 15;
    if (numMultGates > 32768) d = 16;
    return d;
}

// ============================================================
//  构造 n×n×n 矩阵乘法电路
//  C[i][j] = Σ_k A[i][k] * B[k][j]
//  乘法门数 = n^3
//  这是 Rinocchio 原始 benchmark 使用的电路
// ============================================================
Circuit buildMatrixMultCircuit(int n) {
    Circuit c;
    long numInput = 2 * n * n;
    long numOutput = n * n;
    long numMultGates = n * n * n;
    long numMid = numMultGates - numOutput;
    long numWires = numInput + numMid + numOutput;

    c.numberOfWires = numWires;
    c.numberOfInputWires = numInput;
    c.numberOfMidWires = numMid;
    c.numberOfOutputWires = numOutput;
    c.numberOfMultiplicationGates = numMultGates;

    long midStart = numInput;

    c.gates.resize(numMultGates);
    int gateIdx = 0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            for (int k = 0; k < n; k++) {
                Gate g;
                long leftA = i * n + k;
                long rightB = n * n + k * n + j;

                if (k == 0) {
                    g.leftInputs.push_back(leftA);
                    g.rightInputs.push_back(rightB);
                } else {
                    long prevWire = numInput + (gateIdx - 1);
                    g.leftInputs.push_back(prevWire);
                    g.leftInputs.push_back(leftA);
                    sort(g.leftInputs.begin(), g.leftInputs.end());
                    g.rightInputs.push_back(rightB);
                }
                c.gates[gateIdx] = g;
                gateIdx++;
            }
        }
    }
    return c;
}

// ============================================================
//  Rinocchio benchmark 结果
// ============================================================
struct RinocchioResult {
    long numMultGates;
    long degree;
    double qrpTime;
    double setupTime;
    double crsTime;
    double computeHTime;
    double proveTime;
    double verifyTime;
};

RinocchioResult benchRinocchio(const Circuit& circuit, long k) {
    RinocchioResult r;
    r.numMultGates = circuit.numberOfMultiplicationGates;
    long degree = autoExtensionDegree(r.numMultGates);
    r.degree = degree;

    initGR(k, degree);

    // QRP
    clock_t t = clock();
    QRP qrp = getQRP(circuit, k, degree);
    r.qrpTime = double(clock() - t) / CLOCKS_PER_SEC;

    // Setup
    t = clock();
    SecretState secret = setup(qrp, 512, 64);
    r.setupTime = double(clock() - t) / CLOCKS_PER_SEC;

    // CRS
    t = clock();
    CRS crs = getCRS(qrp, secret);
    r.crsTime = double(clock() - t) / CLOCKS_PER_SEC;

    // Eval
    Vec<ZZ_p> input;
    input.SetLength(circuit.numberOfInputWires);
    for (long i = 0; i < circuit.numberOfInputWires; i++) {
        input[i] = to_ZZ_p(RandomBits_ZZ(k));
    }
    Vec<ZZ_p> allWireValues = eval(circuit, input);

    Vec<ZZ_p> output;
    output.SetLength(circuit.numberOfOutputWires);
    for (long i = 0; i < circuit.numberOfOutputWires; i++) {
        output[i] = allWireValues[i + qrp.outOffset];
    }

    // Compute H
    t = clock();
    ZZ_pEX H = proverComputeH(qrp, allWireValues);
    r.computeHTime = double(clock() - t) / CLOCKS_PER_SEC;

    // Prove
    t = clock();
    Proof pi = prove(H, crs, allWireValues, qrp.midOffset, qrp.outOffset);
    r.proveTime = double(clock() - t) / CLOCKS_PER_SEC;

    // Verify
    t = clock();
    bool ok = verify(secret, crs, pi, input, output);
    r.verifyTime = double(clock() - t) / CLOCKS_PER_SEC;

    assert(ok && "Rinocchio verification failed!");
    return r;
}

// ============================================================
//  Brakedown benchmark 结果
// ============================================================
struct BrakedownResult {
    long n;           // polynomial size (= numMultGates)
    long degree;
    double setupTime; // code_setup time (included in commit for fairness)
    double commitTime;
    double proveTime;
    double verifyTime;
};

BrakedownResult benchBrakedown(long n, long k, long degree) {
    BrakedownResult r;
    r.n = n;
    r.degree = degree;

    initGR(k, degree);

    // Random polynomial of size n
    vector<ZZ_pE> poly(n);
    for (long i = 0; i < n; i++) poly[i] = random_ZZ_pE();

    // Choose row_len (square-root split, capped at 512)
    long row_len = 1;
    while (row_len * 2 <= n && row_len * 2 <= 512) row_len *= 2;
    long num_rows = (n + row_len - 1) / row_len;

    // Code setup (this is Brakedown's "preprocessing", should be counted)
    clock_t t = clock();
    BrakedownCodeGR code = brakedown_code_setup(row_len, k, degree);
    r.setupTime = double(clock() - t) / CLOCKS_PER_SEC;

    // Commit
    t = clock();
    BrakedownCommitment comm = brakedown_commit(code, poly.data(), n, num_rows);
    r.commitTime = double(clock() - t) / CLOCKS_PER_SEC;

    // Prepare evaluation query
    vector<ZZ_pE> q1(num_rows), q2(row_len);
    if (num_rows == 1) {
        set(q1[0]);
    } else {
        for (long i = 0; i < num_rows; i++) q1[i] = random_ZZ_pE();
    }
    for (long j = 0; j < row_len; j++) q2[j] = random_ZZ_pE();

    // Prove
    t = clock();
    BrakedownEvalProof proof = brakedown_prove(code, comm, poly.data(), n, q1, q2);
    r.proveTime = double(clock() - t) / CLOCKS_PER_SEC;

    // Verify
    t = clock();
    bool ok = brakedown_verify(code, comm, proof, q1, q2);
    r.verifyTime = double(clock() - t) / CLOCKS_PER_SEC;

    assert(ok && "Brakedown verification failed!");
    return r;
}

// ============================================================
//  对多次运行的 RinocchioResult 取平均值
// ============================================================
RinocchioResult averageRinocchio(const vector<RinocchioResult>& trials) {
    RinocchioResult avg;
    avg.numMultGates = trials[0].numMultGates;
    avg.degree = trials[0].degree;
    avg.qrpTime = avg.setupTime = avg.crsTime = 0;
    avg.computeHTime = avg.proveTime = avg.verifyTime = 0;
    for (const auto& t : trials) {
        avg.qrpTime      += t.qrpTime;
        avg.setupTime    += t.setupTime;
        avg.crsTime      += t.crsTime;
        avg.computeHTime += t.computeHTime;
        avg.proveTime    += t.proveTime;
        avg.verifyTime   += t.verifyTime;
    }
    double n = (double)trials.size();
    avg.qrpTime      /= n;
    avg.setupTime    /= n;
    avg.crsTime      /= n;
    avg.computeHTime /= n;
    avg.proveTime    /= n;
    avg.verifyTime   /= n;
    return avg;
}

// ============================================================
//  对多次运行的 BrakedownResult 取平均值
// ============================================================
BrakedownResult averageBrakedown(const vector<BrakedownResult>& trials) {
    BrakedownResult avg;
    avg.n = trials[0].n;
    avg.degree = trials[0].degree;
    avg.setupTime = avg.commitTime = avg.proveTime = avg.verifyTime = 0;
    for (const auto& t : trials) {
        avg.setupTime  += t.setupTime;
        avg.commitTime += t.commitTime;
        avg.proveTime  += t.proveTime;
        avg.verifyTime += t.verifyTime;
    }
    double n = (double)trials.size();
    avg.setupTime  /= n;
    avg.commitTime /= n;
    avg.proveTime  /= n;
    avg.verifyTime /= n;
    return avg;
}

// ============================================================
//  打印分隔线
// ============================================================
void printSep(int width = 140) {
    cout << string(width, '-') << "\n";
}

// ============================================================
//  MAIN
// ============================================================
int main() {
    long k = 64;

    cout << "##############################################################\n";
    cout << "  Fair Benchmark: Brakedown PCS vs Rinocchio\n";
    cout << "  Circuit: n×n×n matrix multiplication (same as original Rinocchio)\n";
    cout << "  Both schemes run in the SAME process with SAME GR parameters\n";
    cout << "  Each size runs " << NUM_TRIALS << " times, results are averaged\n";
    cout << "##############################################################\n\n";

    // 矩阵大小序列: n=2..8 对应乘法门数 8..512
    int matSizes[] = {2, 3, 4, 5, 6, 7, 8};
    int numTests = sizeof(matSizes) / sizeof(matSizes[0]);

    // ============================================================
    //  Phase 1: 逐个运行 Rinocchio (矩阵乘法电路), 每个规模 NUM_TRIALS 次
    // ============================================================
    cout << ">>> Phase 1: Running Rinocchio on matrix multiplication circuits <<<\n";
    cout << "    (" << NUM_TRIALS << " trials per size, showing averages)\n\n";

    cout << left
         << setw(6)  << "MatN"
         << setw(8)  << "Gates"
         << setw(5)  << "d"
         << setw(10) << "QRP(s)"
         << setw(10) << "Setup(s)"
         << setw(10) << "CRS(s)"
         << setw(12) << "ComputeH"
         << setw(10) << "Prove(s)"
         << setw(10) << "Verify(s)"
         << setw(12) << "Prover*"
         << "\n";
    printSep(103);

    vector<RinocchioResult> rResults;
    for (int ti = 0; ti < numTests; ti++) {
        int matN = matSizes[ti];
        long numGates = (long)matN * matN * matN;

        cout << "  Building " << matN << "x" << matN << "x" << matN
             << " circuit (gates=" << numGates << ")..." << flush;

        Circuit c = buildMatrixMultCircuit(matN);

        cout << " running Rinocchio x" << NUM_TRIALS << "..." << flush;

        vector<RinocchioResult> trials(NUM_TRIALS);
        for (int trial = 0; trial < NUM_TRIALS; trial++) {
            trials[trial] = benchRinocchio(c, k);
        }
        RinocchioResult rr = averageRinocchio(trials);
        rResults.push_back(rr);

        double proverTotal = rr.computeHTime + rr.proveTime;
        cout << " done.\n";

        cout << left << fixed << setprecision(4)
             << setw(6)  << matN
             << setw(8)  << rr.numMultGates
             << setw(5)  << rr.degree
             << setw(10) << rr.qrpTime
             << setw(10) << rr.setupTime
             << setw(10) << rr.crsTime
             << setw(12) << rr.computeHTime
             << setw(10) << rr.proveTime
             << setw(10) << rr.verifyTime
             << setw(12) << proverTotal
             << "\n";
    }
    printSep(103);
    cout << "  * Prover = ComputeH + Prove (fair prover total)\n";
    cout << "  * All values are averages over " << NUM_TRIALS << " runs\n\n";

    // ============================================================
    //  Phase 2: 逐个运行 Brakedown (相同多项式规模), 每个规模 NUM_TRIALS 次
    // ============================================================
    cout << ">>> Phase 2: Running Brakedown on same polynomial sizes <<<\n";
    cout << "    (" << NUM_TRIALS << " trials per size, showing averages)\n\n";

    cout << left
         << setw(6)  << "MatN"
         << setw(8)  << "n"
         << setw(5)  << "d"
         << setw(12) << "CodeSetup"
         << setw(10) << "Commit(s)"
         << setw(10) << "Prove(s)"
         << setw(10) << "Verify(s)"
         << setw(12) << "Prover*"
         << "\n";
    printSep(73);

    vector<BrakedownResult> bResults;
    for (int ti = 0; ti < numTests; ti++) {
        int matN = matSizes[ti];
        long numGates = (long)matN * matN * matN;
        long degree = autoExtensionDegree(numGates);

        cout << "  Brakedown n=" << numGates << " d=" << degree
             << " x" << NUM_TRIALS << "..." << flush;

        vector<BrakedownResult> trials(NUM_TRIALS);
        for (int trial = 0; trial < NUM_TRIALS; trial++) {
            trials[trial] = benchBrakedown(numGates, k, degree);
        }
        BrakedownResult br = averageBrakedown(trials);
        bResults.push_back(br);

        double proverTotal = br.setupTime + br.commitTime + br.proveTime;
        cout << " done.\n";

        cout << left << fixed << setprecision(4)
             << setw(6)  << matN
             << setw(8)  << br.n
             << setw(5)  << br.degree
             << setw(12) << br.setupTime
             << setw(10) << br.commitTime
             << setw(10) << br.proveTime
             << setw(10) << br.verifyTime
             << setw(12) << proverTotal
             << "\n";
    }
    printSep(73);
    cout << "  * Prover = CodeSetup + Commit + Prove\n";
    cout << "  * All values are averages over " << NUM_TRIALS << " runs\n\n";

    // ============================================================
    //  Phase 3: 对比表格
    // ============================================================
    cout << ">>> Phase 3: Head-to-head comparison (averaged over "
         << NUM_TRIALS << " runs) <<<\n\n";

    cout << left
         << setw(6)  << "MatN"
         << setw(8)  << "Gates"
         << setw(5)  << "d"
         << " | "
         << setw(12) << "R.Setup†"
         << setw(12) << "R.CRS"
         << setw(12) << "R.Prover"
         << setw(12) << "R.Verify"
         << " | "
         << setw(12) << "B.Setup"
         << setw(12) << "B.Commit"
         << setw(12) << "B.Prover"
         << setw(12) << "B.Verify"
         << " | "
         << setw(10) << "Speedup"
         << "\n";
    printSep(145);

    for (int ti = 0; ti < numTests; ti++) {
        const RinocchioResult& rr = rResults[ti];
        const BrakedownResult& br = bResults[ti];

        double r_prover = rr.crsTime + rr.computeHTime + rr.proveTime;
        double b_prover = br.setupTime + br.commitTime + br.proveTime;
        double speedup = (b_prover > 1e-6) ? r_prover / b_prover : 0.0;

        cout << left << fixed << setprecision(4)
             << setw(6)  << matSizes[ti]
             << setw(8)  << rr.numMultGates
             << setw(5)  << rr.degree
             << " | "
             << setw(12) << rr.setupTime
             << setw(12) << rr.crsTime
             << setw(12) << (rr.computeHTime + rr.proveTime)
             << setw(12) << rr.verifyTime
             << " | "
             << setw(12) << br.setupTime
             << setw(12) << br.commitTime
             << setw(12) << br.proveTime
             << setw(12) << br.verifyTime
             << " | "
             << setprecision(1) << setw(6) << speedup << "x"
             << "\n";
    }
    printSep(145);

    cout << "\n";
    cout << "  注释:\n";
    cout << "    所有数据均为 " << NUM_TRIALS << " 次运行的平均值\n";
    cout << "    R.Setup†  = Rinocchio 可信设置 (一次性, Brakedown 无需)\n";
    cout << "    R.Prover  = ComputeH + Prove (多项式除法 + 同态加密)\n";
    cout << "    B.Prover  = Prove (稀疏矩阵编码 + Merkle 哈希)\n";
    cout << "    Speedup   = R.(CRS+Prover) / B.(Setup+Commit+Prove)\n";
    cout << "    电路类型: n×n×n 矩阵乘法 (与原始 Rinocchio 论文一致)\n";
    cout << "\n";
    cout << "  Rinocchio 的 QRP 构造时间未纳入 Speedup 计算,\n";
    cout << "  因为它是电路编译步骤, 不属于多项式承诺操作.\n";

    cout << "\n##############################################################\n";
    cout << "  Benchmark complete.\n";
    cout << "##############################################################\n";

    return 0;
}
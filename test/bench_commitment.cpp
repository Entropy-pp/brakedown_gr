/**
 * bench_commitment.cpp
 * 
 * 独立测试 Rinocchio 中多项式承诺相关操作的性能
 * 不依赖任何外部电路文件，在代码内直接构造任意规模的电路
 * 
 * 测试的三个阶段 (对应多项式承诺):
 *   1. Setup  = setup() + getCRS()   → 相当于 Commit 预处理
 *   2. Prove  = proverComputeH() + prove()  → 相当于 Eval/Open
 *   3. Verify = verify()              → 相当于 Verify
 * 
 * 编译:
 *   g++ -O2 -std=c++11 -I./include test/bench_commitment.cpp src/*.cpp \
 *       -o out/bench_commitment -pthread -lntl -lgmp -lm
 * 
 * 运行:
 *   ./out/bench_commitment
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
#include <gr.h>
#include <circuit.h>
#include <qrp.h>
#include <setup.h>
#include <rinocchio.h>

using namespace std;
using namespace NTL;

// ============================================================
//  辅助函数：初始化 Galois Ring GR(2^k, degree)
// ============================================================
void initGR(long k, long degree) {
    ZZ modulus = ZZ(1) << k;
    ZZ_p::init(modulus);
    ZZ_pX P = primitiveIrredPoly(degree);
    ZZ_pE::init(P);
}

// ============================================================
//  辅助函数：根据乘法门数量自动选扩展次数
//  与 rinocchio_test.cpp 中 testFile() 完全一致
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
//  辅助函数：构造 n×n×n 矩阵乘法电路 (不依赖外部文件)
//
//  与 rinocchio_test.cpp 里 circuitFromFile 读入的电路等价
//  电路结构:
//    input:  2*n*n 个输入线 (两个 n×n 矩阵)
//    output: n*n 个输出线 (结果矩阵)
//    每个输出 C[i][j] = sum_{k} A[i][k] * B[k][j]
//    需要 n*n*n 个乘法门, n*n*(n-1) 个加法门(隐式)
//
//  乘法门数 = n^3
//  wires = 2n^2 (input) + n^3 (mid, 各乘法门的输出) + n^2 (output)
//       但这里的简化模型: 每个 C[i][j] 用一条链, 中间节点 = n-1
//       总 midWires = n*n*(n-1), 但 Rinocchio 把非 input/output 都算 mid
// ============================================================
Circuit buildMatrixMultCircuit(int n) {
    Circuit c;
    long numInput = 2 * n * n;
    long numOutput = n * n;
    long numMultGates = n * n * n;
    // mid wires = 每条乘积链的中间累加结果: n*n*(n-1) 个
    // 但第一个乘法直接产生一个 mid wire, 后续 n-1 个也是 mid
    // 总 mid = numMultGates - numOutput = n^3 - n^2
    long numMid = numMultGates - numOutput;
    long numWires = numInput + numMid + numOutput;

    c.numberOfWires = numWires;
    c.numberOfInputWires = numInput;
    c.numberOfMidWires = numMid;
    c.numberOfOutputWires = numOutput;
    c.numberOfMultiplicationGates = numMultGates;

    // 线编号:
    //   [0, 2n^2 - 1]                 = input (A 占前 n^2, B 占后 n^2)
    //   [2n^2, 2n^2 + numMid - 1]     = mid wires
    //   [numWires - n^2, numWires - 1] = output wires

    // A[i][k] 对应 wire: i*n + k
    // B[k][j] 对应 wire: n*n + k*n + j

    long midStart = numInput;
    long outStart = numWires - numOutput;

    c.gates.resize(numMultGates);

    int gateIdx = 0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            // C[i][j] = A[i][0]*B[0][j] + A[i][1]*B[1][j] + ... + A[i][n-1]*B[n-1][j]
            // 第一个乘法门: A[i][0] * B[0][j]
            // 后续: (前一个结果 + A[i][k]) * B[k][j]
            // 但 Rinocchio 电路模型: 每个 gate 是一个乘法门
            // leftInputs 和 rightInputs 可以有多条线 (隐式加法)
            //
            // 简化方案 (与标准做法一致):
            //   gate k=0: left={A[i][0]}, right={B[0][j]}  → wire = 某 mid wire
            //   gate k=1: left={上一个mid, A[i][1]}, right={B[1][j]} → wire = 某 mid wire
            //   ...
            //   gate k=n-1: left={上一个mid, A[i][n-1]}, right={B[n-1][j]} → output wire

            for (int k = 0; k < n; k++) {
                Gate g;
                long leftA = i * n + k;      // A[i][k]
                long rightB = n * n + k * n + j; // B[k][j]

                if (k == 0) {
                    // 第一个乘法: A[i][0] * B[0][j]
                    g.leftInputs.push_back(leftA);
                    g.rightInputs.push_back(rightB);
                } else {
                    // 累加: (prev_result + A[i][k]) * B[k][j]
                    // leftInputs = {prev_wire, A[i][k]}
                    long prevWire;
                    if (k == 1) {
                        // 前一个门的输出 wire
                        prevWire = midStart + (gateIdx - 1) - (i * n + j) * 0;
                        // 更简单: 直接用 gateIdx 算
                        prevWire = numInput + (gateIdx - 1);
                    } else {
                        prevWire = numInput + (gateIdx - 1);
                    }
                    g.leftInputs.push_back(prevWire);
                    g.leftInputs.push_back(leftA);
                    // 排序 (Rinocchio 代码中 binary_search 要求有序)
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
//  辅助函数：构造简单的 numGates 个独立乘法门的电路
//  更纯粹的 benchmark：每个门是独立的 a_i * b_i = c_i
//  input = 2*numGates, output = numGates, mid = 0
// ============================================================
Circuit buildSimpleMultCircuit(long numGates) {
    Circuit c;
    c.numberOfInputWires = 2 * numGates;
    c.numberOfOutputWires = numGates;
    c.numberOfMidWires = 0;
    c.numberOfWires = c.numberOfInputWires + c.numberOfMidWires + c.numberOfOutputWires;
    c.numberOfMultiplicationGates = numGates;

    c.gates.resize(numGates);
    for (long i = 0; i < numGates; i++) {
        Gate g;
        g.leftInputs.push_back(2 * i);      // a_i
        g.rightInputs.push_back(2 * i + 1);  // b_i
        c.gates[i] = g;
    }

    return c;
}

// ============================================================
//  核心性能测试函数
// ============================================================
struct BenchResult {
    double qrpTime;       // 电路 → QRP 的时间
    double setupTime;     // setup() 时间
    double crsTime;       // getCRS() 时间
    double evalTime;      // 电路求值时间
    double computeHTime;  // proverComputeH() 时间
    double proveTime;     // prove() 时间
    double verifyTime;    // verify() 时间
    double totalProverTime;  // computeH + prove
    double totalTime;     // 全部时间
    long numMultGates;
    long numWires;
    long extensionDegree;
};

BenchResult benchCircuit(const Circuit& circuit, long k = 64) {
    BenchResult result;
    result.numMultGates = circuit.numberOfMultiplicationGates;
    result.numWires = circuit.numberOfWires;

    // 选择扩展次数
    long degree = autoExtensionDegree(circuit.numberOfMultiplicationGates);
    result.extensionDegree = degree;

    // 初始化 Galois Ring
    initGR(k, degree);

    // ============ Phase 1: Circuit → QRP ============
    clock_t t = clock();
    QRP qrp = getQRP(circuit, k, degree);
    result.qrpTime = double(clock() - t) / CLOCKS_PER_SEC;

    // ============ Phase 2: Setup (Trusted Setup) ============
    // 这一步相当于多项式承诺的 KeyGen + Commit 预处理
    t = clock();
    SecretState secret = setup(qrp, 512, 64);
    result.setupTime = double(clock() - t) / CLOCKS_PER_SEC;

    t = clock();
    CRS crs = getCRS(qrp, secret);
    result.crsTime = double(clock() - t) / CLOCKS_PER_SEC;

    // ============ Phase 3: 准备输入、求值电路 ============
    Vec<ZZ_p> input;
    input.SetLength(circuit.numberOfInputWires);
    for (long i = 0; i < circuit.numberOfInputWires; i++) {
        input[i] = to_ZZ_p(RandomBits_ZZ(k));
    }

    t = clock();
    Vec<ZZ_p> allWireValues = eval(circuit, input);
    result.evalTime = double(clock() - t) / CLOCKS_PER_SEC;

    Vec<ZZ_p> output;
    output.SetLength(circuit.numberOfOutputWires);
    for (long i = 0; i < circuit.numberOfOutputWires; i++) {
        output[i] = allWireValues[i + qrp.outOffset];
    }

    // ============ Phase 4: Prove ============
    // 4a: 计算 H(x) = (V·W - Y) / t —— 多项式除法，最耗时的部分之一
    t = clock();
    ZZ_pEX H = proverComputeH(qrp, allWireValues);
    result.computeHTime = double(clock() - t) / CLOCKS_PER_SEC;

    // 4b: 用 CRS 中的加密值构造证明 —— 同态加密操作
    t = clock();
    Proof pi = prove(H, crs, allWireValues, qrp.midOffset, qrp.outOffset);
    result.proveTime = double(clock() - t) / CLOCKS_PER_SEC;

    result.totalProverTime = result.computeHTime + result.proveTime;

    // ============ Phase 5: Verify ============
    t = clock();
    bool ok = verify(secret, crs, pi, input, output);
    result.verifyTime = double(clock() - t) / CLOCKS_PER_SEC;

    assert(ok && "Verification failed!");

    result.totalTime = result.qrpTime + result.setupTime + result.crsTime
                     + result.evalTime + result.totalProverTime + result.verifyTime;

    return result;
}

// ============================================================
//  打印结果
// ============================================================
void printHeader() {
    cout << "\n";
    cout << string(120, '=') << "\n";
    cout << left
         << setw(10) << "MultGates"
         << setw(10) << "Wires"
         << setw(6)  << "Deg"
         << setw(12) << "QRP(s)"
         << setw(12) << "Setup(s)"
         << setw(12) << "CRS(s)"
         << setw(12) << "ComputeH"
         << setw(12) << "Prove(s)"
         << setw(12) << "Verify(s)"
         << setw(14) << "Prover Tot"
         << setw(12) << "Total(s)"
         << "\n";
    cout << string(120, '-') << "\n";
}

void printResult(const BenchResult& r) {
    cout << left << fixed << setprecision(4)
         << setw(10) << r.numMultGates
         << setw(10) << r.numWires
         << setw(6)  << r.extensionDegree
         << setw(12) << r.qrpTime
         << setw(12) << r.setupTime
         << setw(12) << r.crsTime
         << setw(12) << r.computeHTime
         << setw(12) << r.proveTime
         << setw(12) << r.verifyTime
         << setw(14) << r.totalProverTime
         << setw(12) << r.totalTime
         << "\n";
}

// ============================================================
//  main
// ============================================================
int main() {

    cout << "##############################################################\n";
    cout << "  Rinocchio Polynomial Commitment Performance Benchmark\n";
    cout << "##############################################################\n";

    // ==========================================
    //  测试 1: 简单独立乘法门 (纯测承诺开销)
    //  乘法门数: 4, 8, 16, 32, 64, 128, 256, 512, 1024
    // ==========================================
    cout << "\n>>> Test 1: Simple independent multiplication gates <<<\n";
    printHeader();

    long simpleSizes[] = {4, 8, 16, 32, 64, 128, 256, 512, 1024};
    for (long numGates : simpleSizes) {
        Circuit c = buildSimpleMultCircuit(numGates);
        BenchResult r = benchCircuit(c);
        printResult(r);
    }

    // ==========================================
    //  测试 2: n×n×n 矩阵乘法电路 (与原测试一致)
    //  n = 2,3,4,5,6,7,8,9,10
    //  乘法门数 = n^3
    // ==========================================
    cout << "\n>>> Test 2: n×n×n Matrix multiplication circuits <<<\n";
    printHeader();

    for (int n = 2; n <= 10; n++) {
        cout << "  Building " << n << "x" << n << "x" << n
             << " matrix mult circuit (gates=" << n*n*n << ")..." << flush;
        Circuit c = buildMatrixMultCircuit(n);
        cout << " done.\n";
        BenchResult r = benchCircuit(c);
        printResult(r);
    }

    // ==========================================
    //  测试 3: 分阶段占比分析 (取一个中等规模)
    // ==========================================
    cout << "\n>>> Test 3: Detailed breakdown for 6x6x6 matrix mult <<<\n";
    {
        Circuit c = buildMatrixMultCircuit(6);
        BenchResult r = benchCircuit(c);

        cout << "\n";
        cout << "  Circuit: " << r.numMultGates << " mult gates, "
             << r.numWires << " wires, degree=" << r.extensionDegree << "\n\n";

        double total = r.totalTime;
        cout << fixed << setprecision(4);
        cout << "  QRP construction:     " << setw(10) << r.qrpTime
             << " s  (" << setw(5) << setprecision(1) << 100.0*r.qrpTime/total << "%)\n";
        cout << "  Setup (keygen):       " << setw(10) << setprecision(4) << r.setupTime
             << " s  (" << setw(5) << setprecision(1) << 100.0*r.setupTime/total << "%)\n";
        cout << "  CRS generation:       " << setw(10) << setprecision(4) << r.crsTime
             << " s  (" << setw(5) << setprecision(1) << 100.0*r.crsTime/total << "%)\n";
        cout << "  Circuit evaluation:   " << setw(10) << setprecision(4) << r.evalTime
             << " s  (" << setw(5) << setprecision(1) << 100.0*r.evalTime/total << "%)\n";
        cout << "  Compute H (poly div): " << setw(10) << setprecision(4) << r.computeHTime
             << " s  (" << setw(5) << setprecision(1) << 100.0*r.computeHTime/total << "%)\n";
        cout << "  Prove (hom. enc.):    " << setw(10) << setprecision(4) << r.proveTime
             << " s  (" << setw(5) << setprecision(1) << 100.0*r.proveTime/total << "%)\n";
        cout << "  Verify (decrypt):     " << setw(10) << setprecision(4) << r.verifyTime
             << " s  (" << setw(5) << setprecision(1) << 100.0*r.verifyTime/total << "%)\n";
        cout << "  ────────────────────────────────────────\n";
        cout << "  Total:                " << setw(10) << setprecision(4) << total << " s\n";
        cout << "  Prover total (H+Prove): " << setw(10) << r.totalProverTime << " s\n";

        cout << "\n  [承诺相关操作]:\n";
        double commitRelated = r.setupTime + r.crsTime;
        double openRelated = r.totalProverTime;
        double verifyRelated = r.verifyTime;
        cout << "    Commit (Setup+CRS): " << setw(10) << commitRelated << " s\n";
        cout << "    Open   (Prove):     " << setw(10) << openRelated << " s\n";
        cout << "    Verify:             " << setw(10) << verifyRelated << " s\n";
    }

    cout << "\n##############################################################\n";
    cout << "  Benchmark complete.\n";
    cout << "##############################################################\n";

    return 0;
}
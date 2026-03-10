#ifndef BRAKEDOWN_PARAMS_H
#define BRAKEDOWN_PARAMS_H

#include <cmath>
#include <algorithm>
#include <vector>
#include <cassert>

// ============================================================
// Brakedown code parameters from [GLSTW21] Figure 2
// Ported from: SuccinctPaul/di-spartan-brakedown brakedown.rs
// ============================================================

struct BrakedownSpec {
    double LAMBDA; // security parameter
    double ALPHA;  // compression ratio per level
    double BETA;   // sparsity of A matrices
    double R;      // code rate (codeword_len / message_len)

    double delta() const { return BETA / R; }
    double mu()    const { return R - 1.0 - R * ALPHA; }
    double nu()    const { return BETA + ALPHA * BETA + 0.03; }
};

// Binary entropy function H(p)
inline double binary_entropy(double p) {
    assert(p > 0.0 && p < 1.0);
    return -p * std::log2(p) - (1.0 - p) * std::log2(1.0 - p);
}

inline long iceil(double v) { return (long)std::ceil(v); }

// Sparsity c_n for A matrices (number of non-zeros per row)
inline long brakedown_c_n(const BrakedownSpec& spec, long n) {
    double alpha = spec.ALPHA;
    double beta  = spec.BETA;
    double nf    = (double)n;
    long opt1 = std::max(iceil(1.28 * beta * nf), iceil(beta * nf) + 4);
    long opt2 = iceil(
        ((110.0 / nf) + binary_entropy(beta) + alpha * binary_entropy(1.28 * beta / alpha))
        / (beta * std::log2(alpha / (1.28 * beta)))
    );
    return std::min(opt1, opt2);
}

// Sparsity d_n for B matrices
inline long brakedown_d_n(const BrakedownSpec& spec, long log2_q, long n) {
    double alpha = spec.ALPHA;
    double beta  = spec.BETA;
    double r     = spec.R;
    double mu_v  = spec.mu();
    double nu_v  = spec.nu();
    double lq    = (double)log2_q;
    double nf    = (double)n;
    long opt1 = iceil((2.0 * beta + ((r - 1.0) + 110.0 / nf) / lq) * nf);
    long opt2 = iceil(
        (r * alpha * binary_entropy(beta / r) + mu_v * binary_entropy(nu_v / mu_v) + 110.0 / nf)
        / (alpha * beta * std::log2(mu_v / nu_v))
    );
    return std::min(opt1, opt2);
}

inline long brakedown_num_column_opening(const BrakedownSpec& spec) {
    // return iceil(-spec.LAMBDA / std::log2(1.0 - spec.delta() / 3.0));
    return spec.LAMBDA;
}

// Dimension descriptor for one level of the recursive code
struct LevelDim {
    long n; // input dimension
    long m; // output dimension
    long d; // sparsity (non-zeros per row)
};

// Compute dimensions for all levels of the recursive code
inline void brakedown_dimensions(
    const BrakedownSpec& spec, long log2_q, long n, long n_0,
    std::vector<LevelDim>& a_dims, std::vector<LevelDim>& b_dims)
{
    assert(n > n_0);
    a_dims.clear();
    b_dims.clear();

    // Build A-chain: n → ceil(n*α) → ceil(ceil(n*α)*α) → ... until ≤ n_0
    std::vector<long> sizes;
    sizes.push_back(n);
    while (true) {
        long prev = sizes.back();
        long next = iceil(prev * spec.ALPHA);
        if (next >= prev) break;
        sizes.push_back(next);
        if (next <= n_0) break;
    }

    for (size_t i = 0; i + 1 < sizes.size(); i++) {
        if (sizes[i] <= n_0) break;
        long nn = sizes[i];
        long mm = sizes[i + 1];
        long dd = std::min(brakedown_c_n(spec, nn), mm);
        a_dims.push_back({nn, mm, dd});
    }

    // Build B-chain
    for (size_t i = 0; i < a_dims.size(); i++) {
        long n_prime = iceil(a_dims[i].m * spec.R);
        long m_prime = iceil(a_dims[i].n * spec.R) - a_dims[i].n - n_prime;
        if (m_prime <= 0) m_prime = 1;
        long d_prime = std::min(brakedown_d_n(spec, log2_q, a_dims[i].n), m_prime);
        b_dims.push_back({n_prime, m_prime, d_prime});
    }
}

// Total codeword length
inline long brakedown_codeword_len(
    const BrakedownSpec& spec, long log2_q, long n, long n_0)
{
    std::vector<LevelDim> a, b;
    brakedown_dimensions(spec, log2_q, n, n_0, a, b);
    if (a.empty()) return iceil(n * spec.R);

    // 与 brakedown_encode 的 assemble 布局完全对齐：
    // [message] [inter[0]] ... [inter[L-2]] [inter[L-1]] [rs_output] [b_out[L-1]] ... [b_out[0]]

    long total = a[0].n;  // message

    // 所有 A 中间结果 (inter[0] .. inter[L-1])
    for (size_t i = 0; i < a.size(); i++) total += a[i].m;

    // RS 额外输出 (rs_output)，输入是最后一个 intermediate
    long rs_input  = a.back().m;
    long rs_output = iceil(rs_input * spec.R) - rs_input;
    if (rs_output < 1) rs_output = 1;
    total += rs_output;

    // 所有 B 输出
    for (size_t i = 0; i < b.size(); i++) total += b[i].m;

    return total;
}

inline long brakedown_num_proximity_testing(
    const BrakedownSpec& spec, long log2_q, long n, long n_0)
{
    long cw_len = brakedown_codeword_len(spec, log2_q, n, n_0);
    double denom = (double)log2_q - std::log2((double)cw_len);
    assert(denom > 0.0 && "log2_q must be larger than log2(codeword_len)");
    return iceil(spec.LAMBDA / denom);
}

// Pre-configured spec from GLSTW21 Figure 2 (Spec 6 — best for large fields)
inline BrakedownSpec getDefaultSpec() {
    return {128.0, 0.2380, 0.1205, 1.720};
}

#endif
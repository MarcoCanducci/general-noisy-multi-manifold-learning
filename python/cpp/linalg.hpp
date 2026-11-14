// Small dense linear algebra helpers (symmetric eigen-decomposition, covariance).
// Part of manydream (AGPL-3.0-or-later).
#pragma once
#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>

namespace manydream {

// Cyclic Jacobi eigen-decomposition of a symmetric d x d matrix A (row-major, destroyed).
// evals are returned in descending order; evecs holds the matching unit eigenvectors as rows.
inline void sym_eig(double* A, int d, double* evals, double* evecs) {
    std::vector<double> V(static_cast<size_t>(d) * d, 0.0);
    for (int i = 0; i < d; ++i) V[i * d + i] = 1.0;
    for (int sweep = 0; sweep < 64; ++sweep) {
        double off = 0.0, diag = 0.0;
        for (int i = 0; i < d; ++i) {
            diag += A[i * d + i] * A[i * d + i];
            for (int j = i + 1; j < d; ++j) off += A[i * d + j] * A[i * d + j];
        }
        if (off <= 1e-30 * (diag + off) || off == 0.0) break;
        for (int p = 0; p < d - 1; ++p) {
            for (int q = p + 1; q < d; ++q) {
                const double apq = A[p * d + q];
                if (apq == 0.0) continue;
                const double theta = (A[q * d + q] - A[p * d + p]) / (2.0 * apq);
                const double t = (theta >= 0 ? 1.0 : -1.0) / (std::fabs(theta) + std::sqrt(theta * theta + 1.0));
                const double c = 1.0 / std::sqrt(t * t + 1.0), s = t * c;
                for (int k = 0; k < d; ++k) {
                    const double akp = A[k * d + p], akq = A[k * d + q];
                    A[k * d + p] = c * akp - s * akq;
                    A[k * d + q] = s * akp + c * akq;
                }
                for (int k = 0; k < d; ++k) {
                    const double apk = A[p * d + k], aqk = A[q * d + k];
                    A[p * d + k] = c * apk - s * aqk;
                    A[q * d + k] = s * apk + c * aqk;
                }
                for (int k = 0; k < d; ++k) {
                    const double vkp = V[k * d + p], vkq = V[k * d + q];
                    V[k * d + p] = c * vkp - s * vkq;
                    V[k * d + q] = s * vkp + c * vkq;
                }
            }
        }
    }
    std::vector<int> order(d);
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](int a, int b) { return A[a * d + a] > A[b * d + b]; });
    for (int r = 0; r < d; ++r) {
        const int c = order[r];
        evals[r] = A[c * d + c];
        if (evecs)
            for (int k = 0; k < d; ++k) evecs[r * d + k] = V[k * d + c];
    }
}

// Mean and (n-1)-normalised covariance of the rows listed in pts.
template <class GetRow>
inline void covariance(GetRow&& row, size_t count, int d, double* mean, double* cov) {
    std::fill(mean, mean + d, 0.0);
    std::fill(cov, cov + d * d, 0.0);
    for (size_t i = 0; i < count; ++i) {
        const double* x = row(i);
        for (int k = 0; k < d; ++k) mean[k] += x[k];
    }
    for (int k = 0; k < d; ++k) mean[k] /= static_cast<double>(count);
    for (size_t i = 0; i < count; ++i) {
        const double* x = row(i);
        for (int a = 0; a < d; ++a) {
            const double da = x[a] - mean[a];
            for (int b = a; b < d; ++b) cov[a * d + b] += da * (x[b] - mean[b]);
        }
    }
    const double norm = count > 1 ? 1.0 / static_cast<double>(count - 1) : 1.0;
    for (int a = 0; a < d; ++a)
        for (int b = a; b < d; ++b) cov[a * d + b] = cov[b * d + a] = cov[a * d + b] * norm;
}

} // namespace manydream

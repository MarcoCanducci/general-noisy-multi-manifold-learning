// Minimal static KD-tree for double precision points in any dimension.
// Part of manydream (AGPL-3.0-or-later).
#pragma once
#include <algorithm>
#include <cstdint>
#include <limits>
#include <numeric>
#include <vector>

namespace manydream {

class KDTree {
public:
    struct Node {
        int64_t lo, hi;      // range [lo, hi) in tree order
        int32_t left, right; // children, -1 for a leaf
        int32_t dim;
        double split;
    };

    int64_t n = 0;
    int d = 0;
    std::vector<double> data;   // points in tree order, row-major
    std::vector<int64_t> index; // tree position -> original row
    std::vector<Node> nodes;

    KDTree() = default;
    KDTree(const double* X, int64_t n_, int d_, int leaf_size = 16) : n(n_), d(d_) {
        index.resize(n);
        std::iota(index.begin(), index.end(), int64_t(0));
        if (n > 0) {
            nodes.reserve(static_cast<size_t>(2 * n / leaf_size + 16));
            build(X, 0, n, leaf_size);
        }
        data.resize(static_cast<size_t>(n) * d);
        for (int64_t p = 0; p < n; ++p)
            std::copy(X + index[p] * d, X + index[p] * d + d, data.begin() + p * d);
    }

    // Calls f(original_index, squared_distance) for every point within sqrt(r2) of q.
    template <class F>
    void radius(const double* q, double r2, F&& f) const {
        if (n > 0) radius_rec(0, q, r2, f);
    }

    // Nearest point accepted by pred(original_index) with squared distance <= max_d2.
    // Returns -1 if there is none.
    template <class Pred>
    int64_t nearest(const double* q, double max_d2, double& best_d2, Pred&& pred) const {
        int64_t best = -1;
        best_d2 = max_d2;
        if (n > 0) nearest_rec(0, q, best, best_d2, pred);
        return best;
    }

private:
    int32_t build(const double* X, int64_t lo, int64_t hi, int leaf_size) {
        const int32_t id = static_cast<int32_t>(nodes.size());
        nodes.push_back(Node{lo, hi, -1, -1, 0, 0.0});
        if (hi - lo <= leaf_size) return id;
        int best_dim = 0;
        double best_spread = -1.0;
        for (int k = 0; k < d; ++k) {
            double mn = std::numeric_limits<double>::infinity(), mx = -mn;
            for (int64_t p = lo; p < hi; ++p) {
                const double v = X[index[p] * d + k];
                mn = std::min(mn, v);
                mx = std::max(mx, v);
            }
            if (mx - mn > best_spread) { best_spread = mx - mn; best_dim = k; }
        }
        if (!(best_spread > 0.0)) return id; // all points identical: keep as a leaf
        const int64_t mid = lo + (hi - lo) / 2;
        std::nth_element(index.begin() + lo, index.begin() + mid, index.begin() + hi,
                         [&](int64_t a, int64_t b) { return X[a * d + best_dim] < X[b * d + best_dim]; });
        const double split = X[index[mid] * d + best_dim];
        const int32_t l = build(X, lo, mid, leaf_size);
        const int32_t r = build(X, mid, hi, leaf_size);
        nodes[id].left = l;
        nodes[id].right = r;
        nodes[id].dim = best_dim;
        nodes[id].split = split;
        return id;
    }

    inline double dist2(const double* q, int64_t pos) const {
        const double* p = data.data() + pos * d;
        double s = 0.0;
        for (int k = 0; k < d; ++k) { const double t = q[k] - p[k]; s += t * t; }
        return s;
    }

    template <class F>
    void radius_rec(int32_t id, const double* q, double r2, F& f) const {
        const Node& nd = nodes[id];
        if (nd.left < 0) {
            for (int64_t p = nd.lo; p < nd.hi; ++p) {
                const double d2 = dist2(q, p);
                if (d2 <= r2) f(index[p], d2);
            }
            return;
        }
        const double diff = q[nd.dim] - nd.split;
        radius_rec(diff < 0 ? nd.left : nd.right, q, r2, f);
        if (diff * diff <= r2) radius_rec(diff < 0 ? nd.right : nd.left, q, r2, f);
    }

    template <class Pred>
    void nearest_rec(int32_t id, const double* q, int64_t& best, double& best_d2, Pred& pred) const {
        const Node& nd = nodes[id];
        if (nd.left < 0) {
            for (int64_t p = nd.lo; p < nd.hi; ++p) {
                const double d2 = dist2(q, p);
                if ((d2 < best_d2 || (best < 0 && d2 <= best_d2)) && pred(index[p])) {
                    best_d2 = d2;
                    best = index[p];
                }
            }
            return;
        }
        const double diff = q[nd.dim] - nd.split;
        nearest_rec(diff < 0 ? nd.left : nd.right, q, best, best_d2, pred);
        if (diff * diff <= best_d2) nearest_rec(diff < 0 ? nd.right : nd.left, q, best, best_d2, pred);
    }
};

} // namespace manydream

// examm._core: compiled kernels for the multi-manifold pipeline.
// Part of examm (AGPL-3.0-or-later).
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <deque>
#include <memory>
#include <stdexcept>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "kdtree.hpp"
#include "linalg.hpp"

namespace py = pybind11;
using examm::KDTree;

using DArr = py::array_t<double, py::array::c_style | py::array::forcecast>;
using IArr = py::array_t<int64_t, py::array::c_style | py::array::forcecast>;
using BArr = py::array_t<uint8_t, py::array::c_style | py::array::forcecast>;

namespace {

int g_threads = 0; // 0 = OpenMP default

int n_threads() {
#ifdef _OPENMP
    return g_threads > 0 ? g_threads : omp_get_max_threads();
#else
    return 1;
#endif
}

void check_2d(const DArr& a, int d, const char* name) {
    if (a.ndim() != 2) throw std::invalid_argument(std::string(name) + " must be a 2-D array");
    if (d > 0 && a.shape(1) != d)
        throw std::invalid_argument(std::string(name) + " has the wrong number of columns");
}

inline double dot(const double* a, const double* b, int d) {
    double s = 0.0;
    for (int k = 0; k < d; ++k) s += a[k] * b[k];
    return s;
}

struct Tree {
    KDTree t;
    std::vector<double> X; // original order
    int64_t n;
    int d;

    explicit Tree(DArr x, int leaf_size) {
        check_2d(x, 0, "X");
        n = x.shape(0);
        d = static_cast<int>(x.shape(1));
        if (d < 1) throw std::invalid_argument("X must have at least one column");
        X.assign(x.data(), x.data() + n * d);
        for (double v : X)
            if (!std::isfinite(v)) throw std::invalid_argument("X contains NaN or infinite values");
        t = KDTree(X.data(), n, d, leaf_size);
    }
    const double* row(int64_t i) const { return X.data() + i * d; }

    // ---------------------------------------------------------------- queries
    IArr radius_count(DArr Q, double r) const {
        check_2d(Q, d, "Q");
        const int64_t m = Q.shape(0);
        IArr out(m);
        int64_t* o = out.mutable_data();
        const double* q = Q.data();
        const double r2 = r * r;
        {
            py::gil_scoped_release nogil;
#pragma omp parallel for schedule(dynamic, 256) num_threads(n_threads())
            for (int64_t i = 0; i < m; ++i) {
                int64_t c = 0;
                t.radius(q + i * d, r2, [&](int64_t, double) { ++c; });
                o[i] = c;
            }
        }
        return out;
    }

    py::tuple radius_neighbors(DArr Q, double r, bool sort) const {
        check_2d(Q, d, "Q");
        const int64_t m = Q.shape(0);
        const double* q = Q.data();
        const double r2 = r * r;
        std::vector<std::vector<std::pair<double, int64_t>>> res(m);
        {
            py::gil_scoped_release nogil;
#pragma omp parallel for schedule(dynamic, 256) num_threads(n_threads())
            for (int64_t i = 0; i < m; ++i) {
                auto& v = res[i];
                t.radius(q + i * d, r2, [&](int64_t j, double d2) { v.emplace_back(d2, j); });
                if (sort) std::sort(v.begin(), v.end());
                else std::sort(v.begin(), v.end(), [](const auto& a, const auto& b) { return a.second < b.second; });
            }
        }
        IArr indptr(m + 1);
        int64_t* ip = indptr.mutable_data();
        ip[0] = 0;
        for (int64_t i = 0; i < m; ++i) ip[i + 1] = ip[i] + static_cast<int64_t>(res[i].size());
        IArr indices(ip[m]);
        py::array_t<double> dists(ip[m]);
        int64_t* id = indices.mutable_data();
        double* dd = dists.mutable_data();
        for (int64_t i = 0; i < m; ++i) {
            int64_t p = ip[i];
            for (const auto& e : res[i]) { id[p] = e.second; dd[p] = std::sqrt(e.first); ++p; }
        }
        return py::make_tuple(indptr, indices, dists);
    }

    // For each query: is any point with flags != 0 within r?
    BArr radius_any(DArr Q, double r, BArr flags) const {
        check_2d(Q, d, "Q");
        if (flags.size() != n) throw std::invalid_argument("flags must have one entry per tree point");
        const int64_t m = Q.shape(0);
        BArr out(m);
        uint8_t* o = out.mutable_data();
        const uint8_t* f = flags.data();
        const double* q = Q.data();
        const double r2 = r * r;
        {
            py::gil_scoped_release nogil;
#pragma omp parallel for schedule(dynamic, 256) num_threads(n_threads())
            for (int64_t i = 0; i < m; ++i) {
                double best;
                o[i] = t.nearest(q + i * d, r2, best, [&](int64_t j) { return f[j] != 0; }) >= 0;
            }
        }
        return out;
    }

    // Mask of the tree points lying within r of at least one query.
    BArr radius_union(DArr Q, double r) const {
        check_2d(Q, d, "Q");
        const int64_t m = Q.shape(0);
        BArr out(n);
        uint8_t* o = out.mutable_data();
        std::fill(o, o + n, uint8_t(0));
        const double* q = Q.data();
        const double r2 = r * r;
        {
            py::gil_scoped_release nogil;
#pragma omp parallel for schedule(dynamic, 64) num_threads(n_threads())
            for (int64_t i = 0; i < m; ++i)
                t.radius(q + i * d, r2, [&](int64_t j, double) {
                    if (!o[j]) o[j] = 1; // benign race: every writer stores 1
                });
        }
        return out;
    }

    py::tuple nearest(DArr Q, double max_dist, py::object alive_obj) const {
        check_2d(Q, d, "Q");
        const int64_t m = Q.shape(0);
        BArr alive;
        const uint8_t* al = nullptr;
        if (!alive_obj.is_none()) {
            alive = alive_obj.cast<BArr>();
            if (alive.size() != n) throw std::invalid_argument("alive must have one entry per tree point");
            al = alive.data();
        }
        IArr idx(m);
        py::array_t<double> dist(m);
        int64_t* oi = idx.mutable_data();
        double* od = dist.mutable_data();
        const double* q = Q.data();
        const double max2 = std::isfinite(max_dist) ? max_dist * max_dist : std::numeric_limits<double>::max();
        {
            py::gil_scoped_release nogil;
#pragma omp parallel for schedule(dynamic, 256) num_threads(n_threads())
            for (int64_t i = 0; i < m; ++i) {
                double b;
                oi[i] = t.nearest(q + i * d, max2, b, [&](int64_t j) { return !al || al[j]; });
                od[i] = oi[i] >= 0 ? std::sqrt(b) : std::numeric_limits<double>::infinity();
            }
        }
        return py::make_tuple(idx, dist);
    }

    // ------------------------------------------------------------- local PCA
    bool local_frame(const double* q, double r2, int n_vec, size_t min_pts, std::vector<int64_t>& nb,
                     double* mean, double* cov, double* evals, double* evecs_full, size_t* count = nullptr) const {
        nb.clear();
        t.radius(q, r2, [&](int64_t j, double) { nb.push_back(j); });
        if (count) *count = nb.size();
        if (nb.size() < std::max<size_t>(min_pts, 2)) return false;
        examm::covariance([&](size_t i) { return row(nb[i]); }, nb.size(), d, mean, cov);
        examm::sym_eig(cov, d, evals, n_vec > 0 ? evecs_full : nullptr);
        return true;
    }

    py::tuple local_pca(DArr Q, double r, int n_vec, int min_pts) const {
        check_2d(Q, d, "Q");
        if (n_vec < 0 || n_vec > d) throw std::invalid_argument("n_vec must be between 0 and d");
        const int64_t m = Q.shape(0);
        IArr counts(m);
        py::array_t<double> evals({m, static_cast<int64_t>(d)});
        py::array_t<double> evecs({m, static_cast<int64_t>(n_vec), static_cast<int64_t>(d)});
        int64_t* oc = counts.mutable_data();
        double* oe = evals.mutable_data();
        double* ov = evecs.mutable_data();
        const double* q = Q.data();
        const double r2 = r * r;
        {
            py::gil_scoped_release nogil;
#pragma omp parallel num_threads(n_threads())
            {
                std::vector<int64_t> nb;
                std::vector<double> mean(d), cov(d * d), ev(d), vec(d * d);
#pragma omp for schedule(dynamic, 128)
                for (int64_t i = 0; i < m; ++i) {
                    size_t c = 0;
                    const bool ok = local_frame(q + i * d, r2, n_vec, static_cast<size_t>(min_pts), nb,
                                                mean.data(), cov.data(), ev.data(), vec.data(), &c);
                    oc[i] = static_cast<int64_t>(c);
                    for (int k = 0; k < d; ++k) oe[i * d + k] = ok ? std::max(ev[k], 0.0) : 0.0;
                    for (int k = 0; k < n_vec * d; ++k) ov[i * n_vec * d + k] = ok ? vec[k] : 0.0;
                }
            }
        }
        return py::make_tuple(counts, evals, evecs);
    }

    // Gaussian-weighted average of the rows of P over the valid neighbours of each tree point.
    py::array_t<double> kernel_smooth(DArr P, double r, double sigma, BArr valid) const {
        if (P.ndim() != 2 || P.shape(0) != n) throw std::invalid_argument("P must have one row per tree point");
        if (valid.size() != n) throw std::invalid_argument("valid must have one entry per tree point");
        const int64_t c = P.shape(1);
        py::array_t<double> out({n, c});
        double* o = out.mutable_data();
        const double* p = P.data();
        const uint8_t* va = valid.data();
        const double r2 = r * r, inv = -1.0 / (2.0 * sigma * sigma);
        {
            py::gil_scoped_release nogil;
#pragma omp parallel num_threads(n_threads())
            {
                std::vector<double> acc(c);
#pragma omp for schedule(dynamic, 128)
                for (int64_t i = 0; i < n; ++i) {
                    std::fill(acc.begin(), acc.end(), 0.0);
                    double ws = 0.0;
                    if (va[i])
                        t.radius(row(i), r2, [&](int64_t j, double d2) {
                            if (!va[j]) return;
                            const double w = std::exp(d2 * inv);
                            ws += w;
                            for (int64_t k = 0; k < c; ++k) acc[k] += w * p[j * c + k];
                        });
                    for (int64_t k = 0; k < c; ++k) o[i * c + k] = ws > 0 ? acc[k] / ws : 0.0;
                }
            }
        }
        return out;
    }

    // -------------------------------------------------------------------- SAF
    // One Jacobi iteration of the structure-aware filter; the tree holds the original points.
    py::array_t<double> saf_step(DArr S, double h, double mu, double support) const {
        check_2d(S, d, "samples");
        const int64_t m = S.shape(0);
        py::array_t<double> out({m, static_cast<int64_t>(d)});
        double* o = out.mutable_data();
        const double* s = S.data();
        const double r2 = support * support, inv = -1.0 / (2.0 * h * h);
        {
            py::gil_scoped_release nogil;
            KDTree st(s, m, d);
#pragma omp parallel num_threads(n_threads())
            {
                std::vector<double> avg(d), rep(d);
#pragma omp for schedule(dynamic, 128)
                for (int64_t i = 0; i < m; ++i) {
                    const double* v = s + i * d;
                    std::fill(avg.begin(), avg.end(), 0.0);
                    std::fill(rep.begin(), rep.end(), 0.0);
                    double wa = 0.0, wr = 0.0;
                    t.radius(v, r2, [&](int64_t j, double d2) {
                        const double w = std::exp(d2 * inv);
                        wa += w;
                        const double* x = row(j);
                        for (int k = 0; k < d; ++k) avg[k] += w * x[k];
                    });
                    if (mu != 0.0)
                        st.radius(v, r2, [&](int64_t j, double d2) {
                            if (j == i) return;
                            const double w = std::exp(d2 * inv) / std::max(std::sqrt(d2), 1e-8);
                            wr += w;
                            const double* x = s + j * d;
                            for (int k = 0; k < d; ++k) rep[k] += w * (v[k] - x[k]);
                        });
                    for (int k = 0; k < d; ++k) {
                        double val = wa > 1e-20 ? avg[k] / wa : v[k];
                        if (wr > 1e-20) val += mu * rep[k] / wr;
                        o[i * d + k] = val;
                    }
                }
            }
        }
        return out;
    }

    // -------------------------------------------------- kernel covariances
    // C_i = sum_n k_in (t_n - c_i)(t_n - c_i)^T / sum_n k_in, k = exp(-|t - c|^2 / (div * scale_i^2)),
    // truncated at support * scale_i. Returns (C, kernel mass).
    py::tuple weighted_cov(DArr C, DArr scales, double div, double support, py::object weights_obj) const {
        DArr pw;
        const double* wt = nullptr;
        if (!weights_obj.is_none()) {
            pw = weights_obj.cast<DArr>();
            if (pw.size() != n) throw std::invalid_argument("weights must have one entry per tree point");
            wt = pw.data();
        }
        check_2d(C, d, "centres");
        const int64_t m = C.shape(0);
        if (scales.size() != m) throw std::invalid_argument("scales must have one entry per centre");
        py::array_t<double> out({m, static_cast<int64_t>(d), static_cast<int64_t>(d)});
        py::array_t<double> mass(m);
        double* o = out.mutable_data();
        double* om = mass.mutable_data();
        const double* c = C.data();
        const double* sc = scales.data();
        {
            py::gil_scoped_release nogil;
#pragma omp parallel for schedule(dynamic, 32) num_threads(n_threads())
            for (int64_t i = 0; i < m; ++i) {
                double* M = o + i * d * d;
                std::fill(M, M + d * d, 0.0);
                const double l2 = sc[i] * sc[i];
                const double inv = -1.0 / (div * l2);
                const double* ci = c + i * d;
                double ws = 0.0;
                t.radius(ci, support * support * l2, [&](int64_t j, double d2) {
                    const double w = std::exp(d2 * inv) * (wt ? wt[j] : 1.0);
                    ws += w;
                    const double* x = row(j);
                    for (int a = 0; a < d; ++a) {
                        const double da = w * (x[a] - ci[a]);
                        for (int b = a; b < d; ++b) M[a * d + b] += da * (x[b] - ci[b]);
                    }
                });
                om[i] = ws;
                for (int a = 0; a < d; ++a)
                    for (int b = a; b < d; ++b) M[a * d + b] = M[b * d + a] = ws > 0 ? M[a * d + b] / ws : 0.0;
            }
        }
        return py::make_tuple(out, mass);
    }

    // ------------------------------------------------------------ ball cover
    // Greedy cover by balls of radius c, visiting points in the given order.
    // Returns (label per point, index of the seed point of each ball).
    py::tuple ball_cover(double c, IArr order) const {
        if (order.size() != n) throw std::invalid_argument("order must be a permutation of the tree points");
        IArr labels(n);
        int64_t* lab = labels.mutable_data();
        std::fill(lab, lab + n, int64_t(-1));
        std::vector<int64_t> seeds;
        const int64_t* ord = order.data();
        const double c2 = c * c;
        {
            py::gil_scoped_release nogil;
            for (int64_t k = 0; k < n; ++k) {
                const int64_t i = ord[k];
                if (i < 0 || i >= n) throw std::invalid_argument("order contains an invalid index");
                if (lab[i] >= 0) continue;
                const int64_t id = static_cast<int64_t>(seeds.size());
                seeds.push_back(i);
                t.radius(row(i), c2, [&](int64_t j, double) { if (lab[j] < 0) lab[j] = id; });
            }
        }
        return py::make_tuple(labels, IArr(seeds.size(), seeds.data()));
    }

    // -------------------------------------------------------------- crawling
    // Frame transport: project the parent frame Up (j x d, rows) onto span(G) and re-orthonormalise.
    void transport(const double* G, const double* Up, int j, double* U) const {
        std::vector<double> B(j * j), S(j * j), ev(j), V(j * j), R(j * j);
        for (int b = 0; b < j; ++b)
            for (int a = 0; a < j; ++a) B[b * j + a] = dot(G + a * d, Up + b * d, d); // U' = B G
        for (int a = 0; a < j; ++a)
            for (int b = 0; b < j; ++b) {
                double s = 0.0;
                for (int k = 0; k < j; ++k) s += B[a * j + k] * B[b * j + k];
                S[a * j + b] = s;
            }
        examm::sym_eig(S.data(), j, ev.data(), V.data());
        if (ev[j - 1] < 1e-8) { // degenerate projection: fall back to the PCA frame
            std::copy(G, G + j * d, U);
            return;
        }
        for (int a = 0; a < j; ++a)
            for (int b = 0; b < j; ++b) {
                double s = 0.0;
                for (int k = 0; k < j; ++k) s += V[k * j + a] * V[k * j + b] / std::sqrt(ev[k]);
                R[a * j + b] = s; // (B B^T)^(-1/2)
            }
        std::vector<double> RB(j * j);
        for (int a = 0; a < j; ++a)
            for (int b = 0; b < j; ++b) {
                double s = 0.0;
                for (int k = 0; k < j; ++k) s += R[a * j + k] * B[k * j + b];
                RB[a * j + b] = s;
            }
        for (int a = 0; a < j; ++a)
            for (int k = 0; k < d; ++k) {
                double s = 0.0;
                for (int b = 0; b < j; ++b) s += RB[a * j + b] * G[b * d + k];
                U[a * d + k] = s;
            }
    }

    // Crawl one manifold of intrinsic dimension j starting at point `seed`.
    // alive: candidate points (nodes are cleared in place); node_of: workspace of -1, restored on exit.
    py::tuple crawl(int64_t seed, py::array_t<uint8_t, py::array::c_style> alive_arr,
                    py::array_t<int64_t, py::array::c_style> node_of_arr, double r, int j, double eta,
                    double beta, int min_pts, int64_t max_nodes) const {
        if (j < 1 || j > d) throw std::invalid_argument("intrinsic dimension must satisfy 1 <= j <= d");
        if (seed < 0 || seed >= n) throw std::invalid_argument("seed out of range");
        if (alive_arr.size() != n || node_of_arr.size() != n)
            throw std::invalid_argument("alive and node_of must have one entry per tree point");
        uint8_t* alive = alive_arr.mutable_data();
        int64_t* node_of = node_of_arr.mutable_data();
        const int ns = 2 * j;
        const double r2 = r * r, br2 = beta * beta * r2;
        const size_t minp = static_cast<size_t>(std::max(min_pts, j + 1));

        std::vector<int64_t> node_pt, slots, edges;
        std::vector<double> frames;
        {
            py::gil_scoped_release nogil;
            std::vector<int64_t> nb;
            std::vector<double> mean(d), cov(d * d), ev(d), G(d * d), U(j * d), z(d), w(d);
            std::deque<int64_t> queue;

            auto add_node = [&](int64_t pt, const double* frame) {
                const int64_t id = static_cast<int64_t>(node_pt.size());
                node_pt.push_back(pt);
                node_of[pt] = id;
                alive[pt] = 0;
                frames.insert(frames.end(), frame, frame + j * d);
                slots.insert(slots.end(), ns, int64_t(-1));
                queue.push_back(id);
                return id;
            };
            // Free slot of node a best aligned with direction dir; -1 if none points forward.
            auto best_slot = [&](int64_t a, const double* dir) {
                int best = -1;
                double bv = 0.0;
                for (int s = 0; s < ns; ++s) {
                    if (slots[a * ns + s] >= 0) continue;
                    const double v = (s < j ? 1.0 : -1.0) * dot(frames.data() + (a * j + s % j) * d, dir, d);
                    if (v > bv) { bv = v; best = s; }
                }
                return best;
            };

            alive[seed] = 0;
            if (local_frame(row(seed), r2, j, minp, nb, mean.data(), cov.data(), ev.data(), G.data()) &&
                ev[j - 1] > 0.0) {
                add_node(seed, G.data());
                while (!queue.empty() && static_cast<int64_t>(node_pt.size()) < max_nodes) {
                    const int64_t i = queue.front();
                    queue.pop_front();
                    const double* v = row(node_pt[i]);
                    for (int s = 0; s < ns; ++s) {
                        if (slots[i * ns + s] >= 0) continue;
                        const double sign = s < j ? 1.0 : -1.0;
                        const double* u = frames.data() + (i * j + s % j) * d;
                        for (int k = 0; k < d; ++k) z[k] = v[k] + sign * eta * r * u[k];
                        double bd2;
                        const int64_t tp = t.nearest(z.data(), r2, bd2, [&](int64_t q) { return alive[q] != 0; });
                        if (tp < 0) continue;
                        const double* tc = row(tp);
                        for (int k = 0; k < d; ++k) w[k] = tc[k] - v[k];
                        if (sign * dot(w.data(), u, d) <= 0.0) continue; // snapped behind the node
                        // contraction: nearest existing node within beta*r of the candidate
                        int64_t e = -1;
                        double ed2 = br2;
                        t.radius(tc, br2, [&](int64_t q, double d2) {
                            if (node_of[q] >= 0 && d2 <= ed2) { ed2 = d2; e = node_of[q]; }
                        });
                        if (e >= 0) {
                            if (e == i) continue;
                            const double* ve = row(node_pt[e]);
                            double de2 = 0.0;
                            for (int k = 0; k < d; ++k) { w[k] = v[k] - ve[k]; de2 += w[k] * w[k]; }
                            if (de2 > r2) continue;
                            bool dup = false;
                            for (int q = 0; q < ns; ++q) dup = dup || slots[i * ns + q] == e;
                            if (dup) continue;
                            const int se = best_slot(e, w.data());
                            if (se < 0) continue; // e is saturated in that direction
                            slots[i * ns + s] = e;
                            slots[e * ns + se] = i;
                            edges.push_back(i);
                            edges.push_back(e);
                        } else {
                            if (dot(w.data(), w.data(), d) > r2) continue;
                            if (!local_frame(tc, r2, j, minp, nb, mean.data(), cov.data(), ev.data(), G.data()) ||
                                !(ev[j - 1] > 0.0)) {
                                alive[tp] = 0; // border point: never a node
                                continue;
                            }
                            transport(G.data(), frames.data() + i * j * d, j, U.data());
                            const int64_t x = add_node(tp, U.data());
                            for (int k = 0; k < d; ++k) w[k] = -w[k];
                            int sx = best_slot(x, w.data());
                            if (sx < 0) sx = (s + j) % ns;
                            slots[i * ns + s] = x;
                            slots[x * ns + sx] = i;
                            edges.push_back(i);
                            edges.push_back(x);
                        }
                    }
                }
            }
            for (int64_t pt : node_pt) node_of[pt] = -1;
        }
        const int64_t K = static_cast<int64_t>(node_pt.size());
        IArr nodes(K, node_pt.data());
        IArr E({static_cast<int64_t>(edges.size() / 2), int64_t(2)});
        std::copy(edges.begin(), edges.end(), E.mutable_data());
        py::array_t<double> F({K, static_cast<int64_t>(j), static_cast<int64_t>(d)});
        std::copy(frames.begin(), frames.end(), F.mutable_data());
        return py::make_tuple(nodes, E, F);
    }
};

// ------------------------------------------------------------------ AGTM E-step
// Mixture with component log-constants logc_k = log(pi_k) - 0.5*logdet_k - 0.5*d*log(2*pi) and
// whitening matrices Wh_k (lower-triangular inverse Cholesky factors): D2 = |Wh_k (t - mu_k)|^2.
// Returns per-point log density and, if stats is true, the sufficient statistics
//   G_k = sum_n R_nk,  S_k = sum_n R_nk t_n,  C_k = sum_n R_nk t_n t_n^T,  Vs_k = sum_n R_nk values_n.
// Optional per-point weights multiply every point's contribution to the statistics (not logp).
// If radii (K,) is given, component k is only evaluated for points within radii[k] of its mean
// (points with no component in reach fall back to the full sum).
py::tuple agtm_estep(DArr T, DArr means, DArr Wh, DArr logc, py::object values_obj, bool stats,
                     py::object radii_obj, py::object weights_obj) {
    check_2d(T, 0, "T");
    const int64_t N = T.shape(0);
    const int d = static_cast<int>(T.shape(1));
    check_2d(means, d, "means");
    const int64_t K = means.shape(0);
    if (Wh.ndim() != 3 || Wh.shape(0) != K || Wh.shape(1) != d || Wh.shape(2) != d)
        throw std::invalid_argument("whitening matrices must have shape (K, d, d)");
    if (logc.size() != K) throw std::invalid_argument("logc must have one entry per component");
    DArr values, radii;
    int64_t q = 0;
    const double* val = nullptr;
    if (!values_obj.is_none()) {
        values = values_obj.cast<DArr>();
        if (values.ndim() != 2 || values.shape(0) != N)
            throw std::invalid_argument("values must have one row per data point");
        q = values.shape(1);
        val = values.data();
    }
    DArr pweights;
    const double* wt = nullptr;
    if (!weights_obj.is_none()) {
        pweights = weights_obj.cast<DArr>();
        if (pweights.size() != T.shape(0)) throw std::invalid_argument("weights must have one entry per data point");
        wt = pweights.data();
    }
    const double* rad = nullptr;
    if (!radii_obj.is_none()) {
        radii = radii_obj.cast<DArr>();
        if (radii.size() != K) throw std::invalid_argument("radii must have one entry per component");
        rad = radii.data();
    }
    py::array_t<double> logp(N);
    py::array_t<double> G(stats ? K : 0);
    py::array_t<double> S({stats ? K : int64_t(0), static_cast<int64_t>(d)});
    py::array_t<double> C({stats ? K : int64_t(0), static_cast<int64_t>(d), static_cast<int64_t>(d)});
    py::array_t<double> Vs({stats ? K : int64_t(0), q});
    double* lp = logp.mutable_data();
    const double *t = T.data(), *mu = means.data(), *wh = Wh.data(), *lc = logc.data();
    const int64_t dd = static_cast<int64_t>(d) * d;
    std::vector<double> gG(stats ? K : 0, 0.0), gS(stats ? K * d : 0, 0.0), gC(stats ? K * dd : 0, 0.0),
        gV(stats ? K * q : 0, 0.0);
    {
        py::gil_scoped_release nogil;
        // components are binned by reach (geometric bins); each bin gets its own tree over the
        // means and is queried with the bin's largest reach. Tiny bins are always evaluated.
        struct Bin { KDTree tree; std::vector<int64_t> ids; std::vector<double> mu; double r2; };
        std::vector<Bin> bins;
        std::vector<int64_t> big;
        if (rad && K > 64) {
            std::vector<int64_t> order(K);
            std::iota(order.begin(), order.end(), int64_t(0));
            std::sort(order.begin(), order.end(), [&](int64_t a, int64_t b) { return rad[a] < rad[b]; });
            int64_t i = 0;
            while (i < K) {
                int64_t j = i;
                const double top = std::max(rad[order[i]], 1e-300) * 1.5;
                while (j < K && rad[order[j]] <= top) ++j;
                if (j - i < 16 || !std::isfinite(rad[order[j - 1]])) {
                    big.insert(big.end(), order.begin() + i, order.begin() + j);
                } else {
                    Bin b;
                    b.ids.assign(order.begin() + i, order.begin() + j);
                    for (int64_t k : b.ids) b.mu.insert(b.mu.end(), mu + k * d, mu + k * d + d);
                    b.r2 = rad[order[j - 1]] * rad[order[j - 1]];
                    bins.push_back(std::move(b));
                }
                i = j;
            }
            for (auto& b : bins) b.tree = KDTree(b.mu.data(), static_cast<int64_t>(b.ids.size()), d);
        } else {
            rad = nullptr;
        }
#pragma omp parallel num_threads(n_threads())
        {
            std::vector<double> la, diff(d);
            std::vector<int64_t> cand;
            std::vector<double> lG(gG.size(), 0.0), lS(gS.size(), 0.0), lC(gC.size(), 0.0), lV(gV.size(), 0.0);
#pragma omp for schedule(dynamic, 256)
            for (int64_t n = 0; n < N; ++n) {
                const double* x = t + n * d;
                cand.clear();
                if (rad) {
                    cand = big;
                    for (const auto& b : bins)
                        b.tree.radius(x, b.r2, [&](int64_t i, double d2) {
                            const int64_t k = b.ids[i];
                            if (d2 <= rad[k] * rad[k]) cand.push_back(k);
                        });
                }
                const bool all = cand.empty();
                const int64_t nc = all ? K : static_cast<int64_t>(cand.size());
                la.resize(nc);
                double mx = -std::numeric_limits<double>::infinity();
                for (int64_t c = 0; c < nc; ++c) {
                    const int64_t k = all ? c : cand[c];
                    const double* m = mu + k * d;
                    const double* w = wh + k * dd;
                    for (int a = 0; a < d; ++a) diff[a] = x[a] - m[a];
                    double d2 = 0.0;
                    for (int a = 0; a < d; ++a) {
                        double s = 0.0;
                        for (int b = 0; b <= a; ++b) s += w[a * d + b] * diff[b];
                        d2 += s * s;
                    }
                    la[c] = lc[k] - 0.5 * d2;
                    mx = std::max(mx, la[c]);
                }
                double sum = 0.0;
                for (int64_t c = 0; c < nc; ++c) {
                    la[c] = la[c] - mx > -40.0 ? std::exp(la[c] - mx) : 0.0;
                    sum += la[c];
                }
                lp[n] = mx + std::log(sum);
                if (!stats) continue;
                for (int64_t c = 0; c < nc; ++c) {
                    const double rnk = la[c] / sum * (wt ? wt[n] : 1.0);
                    if (rnk < 1e-14) continue;
                    const int64_t k = all ? c : cand[c];
                    lG[k] += rnk;
                    for (int a = 0; a < d; ++a) {
                        const double ra = rnk * x[a];
                        lS[k * d + a] += ra;
                        for (int b = a; b < d; ++b) lC[k * dd + a * d + b] += ra * x[b];
                    }
                    for (int64_t v = 0; v < q; ++v) lV[k * q + v] += rnk * val[n * q + v];
                }
            }
#pragma omp critical
            {
                for (size_t i = 0; i < gG.size(); ++i) gG[i] += lG[i];
                for (size_t i = 0; i < gS.size(); ++i) gS[i] += lS[i];
                for (size_t i = 0; i < gC.size(); ++i) gC[i] += lC[i];
                for (size_t i = 0; i < gV.size(); ++i) gV[i] += lV[i];
            }
        }
        if (stats)
            for (int64_t k = 0; k < K; ++k)
                for (int a = 0; a < d; ++a)
                    for (int b = 0; b < a; ++b) gC[k * dd + a * d + b] = gC[k * dd + b * d + a];
    }
    std::copy(gG.begin(), gG.end(), G.mutable_data());
    std::copy(gS.begin(), gS.end(), S.mutable_data());
    std::copy(gC.begin(), gC.end(), C.mutable_data());
    std::copy(gV.begin(), gV.end(), Vs.mutable_data());
    return py::make_tuple(logp, G, S, C, Vs);
}

// Dense responsibilities (N x K); intended for moderate N*K.
py::array_t<double> agtm_resp(DArr T, DArr means, DArr Wh, DArr logc) {
    check_2d(T, 0, "T");
    const int64_t N = T.shape(0);
    const int d = static_cast<int>(T.shape(1));
    check_2d(means, d, "means");
    const int64_t K = means.shape(0);
    if (Wh.ndim() != 3 || Wh.shape(0) != K || Wh.shape(1) != d || Wh.shape(2) != d)
        throw std::invalid_argument("whitening matrices must have shape (K, d, d)");
    if (logc.size() != K) throw std::invalid_argument("logc must have one entry per component");
    py::array_t<double> R({N, K});
    double* o = R.mutable_data();
    const double *t = T.data(), *mu = means.data(), *wh = Wh.data(), *lc = logc.data();
    {
        py::gil_scoped_release nogil;
#pragma omp parallel for schedule(static) num_threads(n_threads())
        for (int64_t n = 0; n < N; ++n) {
            double* la = o + n * K;
            double mx = -std::numeric_limits<double>::infinity();
            for (int64_t k = 0; k < K; ++k) {
                double d2 = 0.0;
                for (int a = 0; a < d; ++a) {
                    double s = 0.0;
                    for (int b = 0; b <= a; ++b) s += wh[(k * d + a) * d + b] * (t[n * d + b] - mu[k * d + b]);
                    d2 += s * s;
                }
                la[k] = lc[k] - 0.5 * d2;
                mx = std::max(mx, la[k]);
            }
            double sum = 0.0;
            for (int64_t k = 0; k < K; ++k) { la[k] = std::exp(la[k] - mx); sum += la[k]; }
            for (int64_t k = 0; k < K; ++k) la[k] /= sum;
        }
    }
    return R;
}

} // namespace

PYBIND11_MODULE(_core, m) {
    m.doc() = "Compiled kernels of examm";
    m.def("set_num_threads", [](int n) { g_threads = n > 0 ? n : 0; },
          "Number of OpenMP threads used by the kernels (0 = all available).");
    m.def("get_num_threads", [] { return n_threads(); });
    m.def("has_openmp", [] {
#ifdef _OPENMP
        return true;
#else
        return false;
#endif
    });

    py::class_<Tree>(m, "KDTree")
        .def(py::init<DArr, int>(), py::arg("X"), py::arg("leaf_size") = 16)
        .def_readonly("n", &Tree::n)
        .def_readonly("d", &Tree::d)
        .def("radius_count", &Tree::radius_count, py::arg("Q"), py::arg("r"))
        .def("radius_neighbors", &Tree::radius_neighbors, py::arg("Q"), py::arg("r"), py::arg("sort") = false)
        .def("radius_any", &Tree::radius_any, py::arg("Q"), py::arg("r"), py::arg("flags"))
        .def("radius_union", &Tree::radius_union, py::arg("Q"), py::arg("r"))
        .def("nearest", &Tree::nearest, py::arg("Q"), py::arg("max_dist") = INFINITY,
             py::arg("alive") = py::none())
        .def("local_pca", &Tree::local_pca, py::arg("Q"), py::arg("r"), py::arg("n_vec") = 0,
             py::arg("min_pts") = 2)
        .def("kernel_smooth", &Tree::kernel_smooth, py::arg("P"), py::arg("r"), py::arg("sigma"),
             py::arg("valid"))
        .def("saf_step", &Tree::saf_step, py::arg("samples"), py::arg("h"), py::arg("mu"), py::arg("support"))
        .def("weighted_cov", &Tree::weighted_cov, py::arg("centres"), py::arg("scales"),
             py::arg("div") = 3.0, py::arg("support") = 5.0, py::arg("weights") = py::none())
        .def("ball_cover", &Tree::ball_cover, py::arg("c"), py::arg("order"))
        .def("crawl", &Tree::crawl, py::arg("seed"), py::arg("alive"), py::arg("node_of"), py::arg("r"),
             py::arg("j"), py::arg("eta") = 0.75, py::arg("beta") = 0.4, py::arg("min_pts") = 5,
             py::arg("max_nodes") = INT64_MAX);

    m.def("agtm_estep", &agtm_estep, py::arg("T"), py::arg("means"), py::arg("whiten"), py::arg("logc"),
          py::arg("values") = py::none(), py::arg("stats") = true, py::arg("radii") = py::none(),
          py::arg("weights") = py::none());
    m.def("agtm_resp", &agtm_resp, py::arg("T"), py::arg("means"), py::arg("whiten"), py::arg("logc"));
}

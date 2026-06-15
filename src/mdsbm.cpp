#include <Rcpp.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using namespace Rcpp;

namespace {

const double NEG_INF = -std::numeric_limits<double>::infinity();

int argmax_first(const std::vector<double>& values) {
  int best = 0;
  double best_value = values.empty() ? NEG_INF : values[0];
  for (int i = 1; i < static_cast<int>(values.size()); ++i) {
    if (values[i] > best_value) {
      best = i;
      best_value = values[i];
    }
  }
  return best;
}

double clip_prob(double value, double eps) {
  if (value > 1.0 - eps) return 1.0 - eps;
  if (value < eps) return eps;
  return value;
}

std::vector<int> cumsum_starts(const std::vector<int>& sizes) {
  std::vector<int> starts(sizes.size() + 1, 0);
  for (int i = 0; i < static_cast<int>(sizes.size()); ++i) {
    starts[i + 1] = starts[i] + sizes[i];
  }
  return starts;
}

std::vector<int> to_std_int(IntegerVector x, const char* name) {
  std::vector<int> out(x.size());
  for (int i = 0; i < x.size(); ++i) {
    if (IntegerVector::is_na(x[i])) {
      stop(std::string(name) + " contains NA.");
    }
    out[i] = x[i];
  }
  return out;
}

struct DenseIntMatrix {
  int nrow = 0;
  int ncol = 0;
  std::vector<int> value;

  DenseIntMatrix() = default;

  DenseIntMatrix(int rows, int cols) : nrow(rows), ncol(cols), value(rows * cols, 0) {}

  int& operator()(int row, int col) {
    return value[static_cast<std::size_t>(row) * ncol + col];
  }

  int operator()(int row, int col) const {
    return value[static_cast<std::size_t>(row) * ncol + col];
  }
};

struct Link {
  int from = -1;
  int to = -1;
};

struct Blockmodel {
  DenseIntMatrix x;
  std::vector<int> sets;
  std::vector<int> extsets;
  std::vector<int> startsets;
  std::vector<int> nem;
  std::vector<int> nvan;
  std::vector<std::vector<Link> > links;
  std::vector<std::vector<int> > left_link;
  std::vector<std::vector<int> > right_link;
  bool symmetric = true;
  int number_of_sets = 0;

  Blockmodel(const IntegerMatrix& input, const std::vector<int>& sets_, bool symmetric_)
      : sets(sets_), symmetric(symmetric_), number_of_sets(static_cast<int>(sets_.size())) {
    if (input.nrow() != input.ncol()) {
      stop("x must be square.");
    }
    DenseIntMatrix original(input.nrow(), input.ncol());
    for (int col = 0; col < input.ncol(); ++col) {
      for (int row = 0; row < input.nrow(); ++row) {
        original(row, col) = input(row, col);
      }
    }
    build(original);
  }

  void build(const DenseIntMatrix& original) {
    std::vector<int> original_starts = cumsum_starts(sets);
    std::vector<std::vector<Link> > original_links = linked_units(original, sets, original_starts);

    nem.assign(number_of_sets, 0);
    nvan.assign(number_of_sets, 0);
    std::vector<std::vector<int> > emerging(std::max(0, number_of_sets - 1));
    std::vector<std::vector<int> > vanishing(std::max(0, number_of_sets - 1));

    for (int s = 0; s < number_of_sets - 1; ++s) {
      std::vector<int> seen_col(sets[s + 1], 0);
      std::vector<int> seen_row(sets[s], 0);
      for (const Link& link : original_links[s]) {
        if (link.to >= 0 && link.to < sets[s + 1]) seen_col[link.to] = 1;
        if (link.from >= 0 && link.from < sets[s]) seen_row[link.from] = 1;
      }
      for (int i = 0; i < sets[s + 1]; ++i) {
        if (!seen_col[i]) emerging[s].push_back(i);
      }
      for (int i = 0; i < sets[s]; ++i) {
        if (!seen_row[i]) vanishing[s].push_back(i);
      }
      nem[s] = static_cast<int>(emerging[s].size());
      nvan[s + 1] = static_cast<int>(vanishing[s].size());
    }

    extsets.resize(number_of_sets);
    for (int s = 0; s < number_of_sets; ++s) {
      extsets[s] = sets[s] + nem[s] + nvan[s];
    }
    startsets = cumsum_starts(extsets);
    x = DenseIntMatrix(startsets.back(), startsets.back());

    for (int sr = 0; sr < number_of_sets; ++sr) {
      for (int i = 0; i < sets[sr]; ++i) {
        const int old_row = original_starts[sr] + i;
        const int new_row = startsets[sr] + i;
        for (int sc = 0; sc < number_of_sets; ++sc) {
          for (int j = 0; j < sets[sc]; ++j) {
            const int old_col = original_starts[sc] + j;
            const int new_col = startsets[sc] + j;
            x(new_row, new_col) = original(old_row, old_col);
          }
        }
      }
    }

    for (int s = 0; s < number_of_sets - 1; ++s) {
      const int row_start = startsets[s] + sets[s];
      const int col_start = startsets[s + 1];
      for (int i = 0; i < nem[s]; ++i) {
        x(row_start + i, col_start + emerging[s][i]) = 1;
      }
    }
    for (int s = 1; s < number_of_sets; ++s) {
      const int row_start = startsets[s - 1];
      const int col_start = startsets[s] + sets[s] + nem[s];
      for (int i = 0; i < nvan[s]; ++i) {
        x(row_start + vanishing[s - 1][i], col_start + i) = 1;
      }
    }

    links = linked_units(x, extsets, startsets);
    precompute_link_lookup();
  }

  static std::vector<std::vector<Link> > linked_units(const DenseIntMatrix& mat,
                                                       const std::vector<int>& sizes,
                                                       const std::vector<int>& starts) {
    const int S = static_cast<int>(sizes.size());
    std::vector<std::vector<Link> > out(std::max(0, S - 1));
    for (int s = 0; s < S - 1; ++s) {
      for (int i = 0; i < sizes[s]; ++i) {
        const int row = starts[s] + i;
        for (int j = 0; j < sizes[s + 1]; ++j) {
          const int col = starts[s + 1] + j;
          if (mat(row, col) == 1) {
            out[s].push_back(Link{i, j});
          }
        }
      }
    }
    return out;
  }

  void precompute_link_lookup() {
    left_link.assign(number_of_sets, std::vector<int>());
    right_link.assign(number_of_sets, std::vector<int>());
    for (int s = 0; s < number_of_sets; ++s) {
      left_link[s].assign(extsets[s], -1);
      right_link[s].assign(extsets[s], -1);
    }
    for (int s = 0; s < number_of_sets - 1; ++s) {
      for (const Link& link : links[s]) {
        if (right_link[s][link.from] == -1) right_link[s][link.from] = link.to;
        if (left_link[s + 1][link.to] == -1) left_link[s + 1][link.to] = link.from;
      }
    }
  }
};

struct MatrixD {
  int nrow = 0;
  int ncol = 0;
  std::vector<double> value;

  MatrixD() = default;

  MatrixD(int rows, int cols) : nrow(rows), ncol(cols), value(rows * cols, 0.0) {}

  double& operator()(int row, int col) {
    return value[static_cast<std::size_t>(row) * ncol + col];
  }

  double operator()(int row, int col) const {
    return value[static_cast<std::size_t>(row) * ncol + col];
  }
};

struct Partition {
  const Blockmodel* bm = nullptr;
  std::vector<int> cl;
  std::vector<std::vector<int> > z;
  std::vector<std::vector<int> > cls;
  std::vector<MatrixD> gamma;
  std::vector<MatrixD> log_gamma;
  std::vector<MatrixD> log_1_gamma;
  std::vector<MatrixD> tau;
  std::vector<MatrixD> tau_reversed;
  std::vector<MatrixD> log_tau;
  std::vector<MatrixD> log_tau_reversed;
  std::vector<std::vector<double> > pi;
  std::vector<std::vector<double> > log_pi;
  std::vector<MatrixD> block_counts_diagonal;
  std::vector<MatrixD> block_sizes_diagonal;
  std::vector<MatrixD> block_counts_off;
  std::vector<std::vector<std::vector<double> > > loglike;
  double ll = 0.0;
  double margin_ll = 0.0;
  double inter_ll = 0.0;
  double intra_ll = 0.0;
  double total_ll = NEG_INF;
  double ICL = NEG_INF;

  Partition() = default;

  Partition(const Blockmodel& bm_, const std::vector<int>& cl_, double epsilon, double epsilonTrans)
      : bm(&bm_), cl(cl_) {
    random_partition();
    recompute_statistics(epsilon, epsilonTrans);
  }

  Partition(const Blockmodel& bm_, const std::vector<int>& cl_, const List& clu,
            double epsilon, double epsilonTrans)
      : bm(&bm_), cl(cl_) {
    set_partition(clu);
    recompute_statistics(epsilon, epsilonTrans);
  }

  void set_partition(const List& clu) {
    const int S = bm->number_of_sets;
    z.assign(S, std::vector<int>());
    if (clu.size() != S) stop("clu must contain one vector per set.");
    for (int s = 0; s < S; ++s) {
      IntegerVector zi = as<IntegerVector>(clu[s]);
      if (zi.size() != bm->sets[s]) {
        stop("Each supplied partition vector must match its original set size.");
      }
      z[s].reserve(bm->extsets[s]);
      for (int i = 0; i < zi.size(); ++i) {
        if (IntegerVector::is_na(zi[i]) || zi[i] < 0 || zi[i] >= cl[s]) {
          stop("Partition labels must be zero-based and smaller than k.");
        }
        z[s].push_back(zi[i]);
      }
      for (int i = 0; i < bm->nem[s]; ++i) z[s].push_back(cl[s]);
      for (int i = 0; i < bm->nvan[s]; ++i) z[s].push_back(cl[s] + 1);
    }
  }

  void random_partition() {
    const int S = bm->number_of_sets;
    z.assign(S, std::vector<int>());
    const int initial_mincluster = 3;
    for (int s = 0; s < S; ++s) {
      if (bm->sets[s] <= initial_mincluster * cl[s]) {
        stop("Cannot initialize random partition: set size is too small for k and the Python-compatible minimum cluster rule.");
      }
      std::vector<int> current(bm->sets[s], 0);
      bool ok = false;
      int attempts = 0;
      while (!ok) {
        std::vector<int> counts(cl[s], 0);
        for (int i = 0; i < bm->sets[s]; ++i) {
          int label = static_cast<int>(std::floor(R::runif(0.0, static_cast<double>(cl[s]))));
          if (label >= cl[s]) label = cl[s] - 1;
          current[i] = label;
          counts[label]++;
        }
        int min_count = *std::min_element(counts.begin(), counts.end());
        ok = min_count > initial_mincluster;
        if (++attempts > 100000) {
          stop("Could not initialize a valid random partition after many attempts.");
        }
      }
      z[s].reserve(bm->extsets[s]);
      z[s].insert(z[s].end(), current.begin(), current.end());
      for (int i = 0; i < bm->nem[s]; ++i) z[s].push_back(cl[s]);
      for (int i = 0; i < bm->nvan[s]; ++i) z[s].push_back(cl[s] + 1);
    }
  }

  void recompute_statistics(double epsilon, double epsilonTrans) {
    compute_cluster_sizes();
    compute_densities(epsilon, epsilonTrans);
    compute_block_frequencies();
  }

  void compute_cluster_sizes() {
    const int S = bm->number_of_sets;
    cls.assign(S, std::vector<int>());
    for (int s = 0; s < S; ++s) {
      cls[s].assign(cl[s] + 2, 0);
      for (int label : z[s]) {
        if (label >= 0 && label < cl[s] + 2) cls[s][label]++;
      }
    }
  }

  void compute_densities(double epsilon, double epsilonTrans) {
    const int S = bm->number_of_sets;
    gamma.assign(S, MatrixD());
    log_gamma.assign(S, MatrixD());
    log_1_gamma.assign(S, MatrixD());
    tau.assign(std::max(0, S - 1), MatrixD());
    tau_reversed.assign(std::max(0, S - 1), MatrixD());
    log_tau.assign(std::max(0, S - 1), MatrixD());
    log_tau_reversed.assign(std::max(0, S - 1), MatrixD());
    pi.assign(S, std::vector<double>());
    log_pi.assign(S, std::vector<double>());

    for (int s = 0; s < S; ++s) {
      const int K = cl[s] + 2;
      MatrixD counts(K, K);
      for (int i = 0; i < bm->extsets[s]; ++i) {
        const int ci = z[s][i];
        const int row = bm->startsets[s] + i;
        for (int j = 0; j < bm->extsets[s]; ++j) {
          const int cj = z[s][j];
          const int col = bm->startsets[s] + j;
          counts(ci, cj) += bm->x(row, col);
        }
      }
      gamma[s] = MatrixD(K, K);
      log_gamma[s] = MatrixD(K, K);
      log_1_gamma[s] = MatrixD(K, K);
      for (int c = 0; c < K; ++c) {
        for (int d = 0; d < K; ++d) {
          double density = epsilon;
          const double product = static_cast<double>(cls[s][c]) * cls[s][d];
          if (product > 1.0) {
            double denom = product;
            if (c == d) denom -= cls[s][c];
            if (denom > 0.0) {
              density = counts(c, d) / denom;
            } else {
              density = 0.5;
            }
          }
          density = clip_prob(density, epsilon);
          gamma[s](c, d) = density;
          log_gamma[s](c, d) = std::log(density);
          log_1_gamma[s](c, d) = std::log(1.0 - density);
        }
      }

      pi[s].assign(cl[s], 0.0);
      log_pi[s].assign(cl[s], NEG_INF);
      double original_total = 0.0;
      for (int c = 0; c < cl[s]; ++c) original_total += cls[s][c];
      for (int c = 0; c < cl[s]; ++c) {
        pi[s][c] = cls[s][c] / original_total;
        log_pi[s][c] = std::log(pi[s][c]);
      }

      if (s < S - 1) {
        const int KR = cl[s] + 2;
        const int KC = cl[s + 1] + 2;
        MatrixD off_counts(KR, KC);
        for (int i = 0; i < bm->extsets[s]; ++i) {
          const int ci = z[s][i];
          const int row = bm->startsets[s] + i;
          for (int j = 0; j < bm->extsets[s + 1]; ++j) {
            const int cj = z[s + 1][j];
            const int col = bm->startsets[s + 1] + j;
            off_counts(ci, cj) += bm->x(row, col);
          }
        }

        tau[s] = MatrixD(KR, KC);
        tau_reversed[s] = MatrixD(KR, KC);
        log_tau[s] = MatrixD(KR, KC);
        log_tau_reversed[s] = MatrixD(KR, KC);
        for (int c = 0; c < KR; ++c) {
          for (int d = 0; d < KC; ++d) {
            double density = 0.0;
            const double product = static_cast<double>(cls[s][c]) * cls[s + 1][d];
            if (product > 1.0) {
              density = off_counts(c, d) / product;
            }
            const double next_size = std::max(cls[s + 1][d], 1);
            const double prev_size = std::max(cls[s][c], 1);
            const double t = clip_prob(density * next_size, epsilonTrans);
            const double tr = clip_prob(prev_size * density, epsilonTrans);
            tau[s](c, d) = t;
            tau_reversed[s](c, d) = tr;
            log_tau[s](c, d) = std::log(t);
            log_tau_reversed[s](c, d) = std::log(tr);
          }
        }
      }
    }
  }

  void compute_block_frequencies() {
    const int S = bm->number_of_sets;
    block_counts_diagonal.assign(S, MatrixD());
    block_sizes_diagonal.assign(S, MatrixD());
    block_counts_off.assign(std::max(0, S - 1), MatrixD());

    for (int s = 0; s < S; ++s) {
      const int Q = cl[s];
      block_counts_diagonal[s] = MatrixD(Q, Q);
      block_sizes_diagonal[s] = MatrixD(Q, Q);
      for (int i = 0; i < bm->extsets[s]; ++i) {
        const int ci = z[s][i];
        if (ci >= Q) continue;
        const int row = bm->startsets[s] + i;
        for (int j = 0; j < bm->extsets[s]; ++j) {
          const int cj = z[s][j];
          if (cj >= Q) continue;
          const int col = bm->startsets[s] + j;
          block_counts_diagonal[s](ci, cj) += bm->x(row, col);
        }
      }
      for (int c = 0; c < Q; ++c) {
        for (int d = 0; d < Q; ++d) {
          double size = static_cast<double>(cls[s][c]) * cls[s][d];
          if (c == d) size -= cls[s][c];
          block_sizes_diagonal[s](c, d) = size;
        }
      }

      if (s < S - 1) {
        const int QN = cl[s + 1];
        block_counts_off[s] = MatrixD(Q, QN);
        for (int i = 0; i < bm->extsets[s]; ++i) {
          const int ci = z[s][i];
          if (ci >= Q) continue;
          const int row = bm->startsets[s] + i;
          for (int j = 0; j < bm->extsets[s + 1]; ++j) {
            const int cj = z[s + 1][j];
            if (cj >= QN) continue;
            const int col = bm->startsets[s + 1] + j;
            block_counts_off[s](ci, cj) += bm->x(row, col);
          }
        }
      }
    }
  }

  std::vector<double> cll(int s, int ind) const {
    const int Q = cl[s];
    std::vector<int> c1(Q, 0), c0(Q, 0), r1(Q, 0), r0(Q, 0);
    const int unit_global = bm->startsets[s] + ind;
    for (int j = 0; j < bm->extsets[s]; ++j) {
      const int label = z[s][j];
      if (label >= Q) continue;
      const int other_global = bm->startsets[s] + j;
      if (bm->x(other_global, unit_global) == 1) {
        c1[label]++;
      } else {
        c0[label]++;
      }
      if (bm->x(unit_global, other_global) == 1) {
        r1[label]++;
      } else {
        r0[label]++;
      }
    }

    std::vector<double> out(Q, 0.0);
    for (int cand = 0; cand < Q; ++cand) {
      double value = 0.0;
      if (s == 0) value += log_pi[s][cand];
      for (int other = 0; other < Q; ++other) {
        value += log_gamma[s](cand, other) * r1[other];
        value += log_1_gamma[s](cand, other) * r0[other];
        if (!bm->symmetric) {
          value += c1[other] * log_gamma[s](other, cand);
          value += c0[other] * log_1_gamma[s](other, cand);
        }
      }
      const int left = bm->left_link[s][ind];
      if (left != -1) {
        const int left_cluster = z[s - 1][left];
        value += log_tau_reversed[s - 1](left_cluster, cand);
      }
      const int right = bm->right_link[s][ind];
      if (right != -1) {
        const int right_cluster = z[s + 1][right];
        value += log_tau[s](cand, right_cluster);
      }
      out[cand] = value;
    }
    return out;
  }

  void ll_elements() {
    ll = 0.0;
    const int S = bm->number_of_sets;
    loglike.assign(S, std::vector<std::vector<double> >());
    for (int s = 0; s < S; ++s) {
      loglike[s].assign(bm->sets[s], std::vector<double>());
      for (int i = 0; i < bm->sets[s]; ++i) {
        std::vector<double> current = cll(s, i);
        ll += current[z[s][i]];
        loglike[s][i] = current;
      }
    }
  }

  void ll_total() {
    margin_ll = 0.0;
    for (int c = 0; c < cl[0]; ++c) {
      margin_ll += cls[0][c] * log_pi[0][c];
    }

    inter_ll = 0.0;
    for (int s = 0; s < bm->number_of_sets - 1; ++s) {
      for (int c = 0; c < cl[s]; ++c) {
        for (int d = 0; d < cl[s + 1]; ++d) {
          inter_ll += block_counts_off[s](c, d) * log_tau[s](c, d);
        }
      }
    }

    intra_ll = 0.0;
    for (int s = 0; s < bm->number_of_sets; ++s) {
      for (int c = 0; c < cl[s]; ++c) {
        for (int d = 0; d < cl[s]; ++d) {
          const double count = block_counts_diagonal[s](c, d);
          const double size = block_sizes_diagonal[s](c, d);
          intra_ll += count * log_gamma[s](c, d);
          intra_ll += (size - count) * log_1_gamma[s](c, d);
        }
      }
    }

    total_ll = bm->symmetric ? margin_ll + inter_ll + intra_ll / 2.0
                             : margin_ll + inter_ll + intra_ll;
    set_icl();
  }

  int recluster(double ratio, int mincluster, double epsilon, double epsilonTrans) {
    int changes = 0;
    for (int s = 0; s < bm->number_of_sets; ++s) {
      for (int i = 0; i < bm->sets[s]; ++i) {
        if (R::runif(0.0, 1.0) < ratio) {
          const int old_cluster = z[s][i];
          const int new_cluster = argmax_first(loglike[s][i]);
          z[s][i] = new_cluster;
          if (new_cluster != old_cluster) ++changes;
        }
      }
    }

    bool ok = false;
    while (!ok) {
      compute_cluster_sizes();
      ok = true;
      for (int s = 0; s < bm->number_of_sets; ++s) {
        for (int c = 0; c < cl[s]; ++c) {
          if (cls[s][c] < mincluster) {
            ok = false;
            int largest = 0;
            for (int candidate = 1; candidate < cl[s]; ++candidate) {
              if (cls[s][candidate] > cls[s][largest]) largest = candidate;
            }
            const int need = mincluster - cls[s][c];
            std::vector<int> indices;
            for (int i = 0; i < static_cast<int>(z[s].size()); ++i) {
              if (z[s][i] == largest) indices.push_back(i);
            }
            if (static_cast<int>(indices.size()) < need) {
              stop("Not enough units in the largest cluster to enforce mincluster.");
            }
            for (int draw = 0; draw < need; ++draw) {
              const int chosen_pos = draw + static_cast<int>(
                  std::floor(R::runif(0.0, static_cast<double>(indices.size() - draw))));
              std::swap(indices[draw], indices[chosen_pos]);
              z[s][indices[draw]] = c;
            }
          }
        }
      }
    }

    recompute_statistics(epsilon, epsilonTrans);
    return changes;
  }

  void set_icl() {
    double pen = 0.5 * (cl[0] - 1.0) * std::log(static_cast<double>(bm->sets[0]));
    int q_next = 0;
    for (int s = 0; s < bm->number_of_sets; ++s) {
      const int Q = cl[s];
      if (s < bm->number_of_sets - 1) q_next = cl[s + 1];
      const double n = bm->sets[s];
      if (bm->symmetric) {
        pen += 0.5 * Q * (Q + 1.0) / 2.0 * std::log(n * (n - 1.0) / 2.0);
      } else {
        pen += 0.5 * Q * Q * std::log(n * (n - 1.0));
      }
      // Keep the original Python penalty behavior, including the final-set term.
      pen += 0.5 * (Q + 1.0) * q_next * std::log(n + bm->nem[s]);
    }
    ICL = total_ll - pen;
  }

  List result() const {
    const int S = bm->number_of_sets;
    List clu(S);
    List clu_size(S);
    for (int s = 0; s < S; ++s) {
      IntegerVector zi(bm->sets[s]);
      for (int i = 0; i < bm->sets[s]; ++i) zi[i] = z[s][i];
      clu[s] = zi;

      IntegerVector sizes(cl[s]);
      for (int c = 0; c < cl[s]; ++c) sizes[c] = cls[s][c];
      clu_size[s] = sizes;
    }

    IntegerVector sets_out(bm->sets.begin(), bm->sets.end());
    IntegerVector k_out(cl.begin(), cl.end());
    return List::create(
      _["clu"] = clu,
      _["sets"] = sets_out,
      _["k"] = k_out,
      _["cluSize"] = clu_size,
      _["ICL"] = ICL,
      _["logLik"] = total_ll,
      _["marginLogLik"] = margin_ll,
      _["interLogLik"] = inter_ll,
      _["intraLogLik"] = intra_ll
    );
  }
};

Partition optimize_partition(Partition part, int mincluster, int maxiter, double chng_ratio,
                             bool verbose, double epsilon, double epsilonTrans,
                             int maxRuns, int maxNoImp) {
  part.ll_elements();
  part.ll_total();
  int changes = 10;
  int runs = 0;
  int iter = 0;
  Partition best = part;
  double old_total = part.total_ll;
  int no_improvement = 0;

  if (verbose) Rcout << "Optimizing...\n";
  while ((changes > 0 || runs < maxRuns) && iter < maxiter && no_improvement < maxNoImp) {
    ++iter;
    if (verbose) Rcout << "ICL " << part.ICL << "\n";
    const double ratio = std::min(chng_ratio * runs + chng_ratio, 1.0);
    changes = part.recluster(ratio, mincluster, epsilon, epsilonTrans);
    if (changes < 2) {
      ++runs;
    } else {
      runs = 0;
    }
    if (verbose) Rcout << changes << " changes\n";
    part.ll_elements();
    part.ll_total();
    if (part.total_ll <= old_total) ++no_improvement;
    old_total = part.total_ll;
    if (part.total_ll > best.total_ll) {
      best = part;
      no_improvement = 0;
    }
  }
  if (verbose) Rcout << "Finished optimizing...\n";
  return best;
}

List add_diagnostics(List out, int runs, int best_run) {
  out["runs"] = runs;
  out["bestRun"] = best_run + 1;
  return out;
}

}  // namespace

// [[Rcpp::export]]
List bm_cpp_icl_partition(IntegerMatrix x, IntegerVector sets, IntegerVector k, List clu,
                          bool symmetric, double epsilon, double epsilonTrans) {
  Blockmodel bm(x, to_std_int(sets, "sets"), symmetric);
  Partition part(bm, to_std_int(k, "k"), clu, epsilon, epsilonTrans);
  part.ll_total();
  return add_diagnostics(part.result(), 1, 0);
}

// [[Rcpp::export]]
List bm_cpp_optimize_partition(IntegerMatrix x, IntegerVector sets, IntegerVector k, List clu,
                               bool symmetric, bool verbose, double epsilon,
                               double epsilonTrans, int mincluster, int maxiter,
                               double chng_ratio, int maxRuns, int maxNoImp) {
  Blockmodel bm(x, to_std_int(sets, "sets"), symmetric);
  Partition part(bm, to_std_int(k, "k"), clu, epsilon, epsilonTrans);
  Partition best = optimize_partition(part, mincluster, maxiter, chng_ratio, verbose,
                                      epsilon, epsilonTrans, maxRuns, maxNoImp);
  return add_diagnostics(best.result(), 1, 0);
}

// [[Rcpp::export]]
List bm_cpp_multiple_optimize(IntegerMatrix x, IntegerVector sets, IntegerVector k, int runs,
                              bool symmetric, bool verbose, double epsilon,
                              double epsilonTrans, int mincluster, int maxiter,
                              double chng_ratio, int maxRuns, int maxNoImp) {
  if (runs < 1) stop("runs must be at least 1.");
  Blockmodel bm(x, to_std_int(sets, "sets"), symmetric);
  const std::vector<int> cl = to_std_int(k, "k");

  if (verbose) Rcout << "optimizing partitions...\n";
  Partition best;
  int best_run = 0;
  bool have_best = false;
  for (int run = 0; run < runs; ++run) {
    Partition part(bm, cl, epsilon, epsilonTrans);
    Partition opt = optimize_partition(part, mincluster, maxiter, chng_ratio, verbose,
                                       epsilon, epsilonTrans, maxRuns, maxNoImp);
    if (!have_best || opt.total_ll > best.total_ll) {
      best = opt;
      best_run = run;
      have_best = true;
    }
  }
  return add_diagnostics(best.result(), runs, best_run);
}


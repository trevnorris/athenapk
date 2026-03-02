#include "mode_tables.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace Modes4D {
namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;

void HermitePhysicists(const int n, const double x, double *hn, double *hnm1) {
  if (n < 0) {
    throw std::invalid_argument("Hermite polynomial order must be non-negative");
  }

  if (n == 0) {
    *hn = 1.0;
    *hnm1 = 0.0;
    return;
  }

  double hm2 = 1.0;
  double hm1 = 2.0 * x;
  if (n == 1) {
    *hn = hm1;
    *hnm1 = hm2;
    return;
  }

  for (int k = 2; k <= n; ++k) {
    const double h = (2.0 * x * hm1) - (2.0 * static_cast<double>(k - 1) * hm2);
    hm2 = hm1;
    hm1 = h;
  }

  *hn = hm1;
  *hnm1 = hm2;
}

void ComputeGaussHermiteStandard(const int n, std::vector<double> *nodes,
                                 std::vector<double> *weights) {
  if (n <= 0) {
    throw std::invalid_argument("Gauss-Hermite order must be positive");
  }

  // Build Hermite (e^{-x^2}) Jacobi matrix and diagonalize it.
  // Off-diagonal entries for orthonormal Hermite functions are sqrt(k/2).
  std::vector<std::vector<double>> a(n, std::vector<double>(n, 0.0));
  std::vector<std::vector<double>> v(n, std::vector<double>(n, 0.0));
  for (int i = 0; i < n; ++i) {
    v[i][i] = 1.0;
  }
  for (int i = 0; i < n - 1; ++i) {
    const double beta = std::sqrt(static_cast<double>(i + 1) / 2.0);
    a[i][i + 1] = beta;
    a[i + 1][i] = beta;
  }

  constexpr double kTol = 1.0e-15;
  const int max_iters = 128 * n * n;
  for (int iter = 0; iter < max_iters; ++iter) {
    int p = 0;
    int q = 1;
    double max_offdiag = 0.0;
    for (int i = 0; i < n; ++i) {
      for (int j = i + 1; j < n; ++j) {
        const double val = std::abs(a[i][j]);
        if (val > max_offdiag) {
          max_offdiag = val;
          p = i;
          q = j;
        }
      }
    }

    if (max_offdiag < kTol) {
      break;
    }

    const double app = a[p][p];
    const double aqq = a[q][q];
    const double apq = a[p][q];
    const double tau = (aqq - app) / (2.0 * apq);
    const double t = ((tau >= 0.0) ? 1.0 : -1.0) /
                     (std::abs(tau) + std::sqrt(1.0 + (tau * tau)));
    const double c = 1.0 / std::sqrt(1.0 + (t * t));
    const double s = t * c;

    a[p][p] = app - (t * apq);
    a[q][q] = aqq + (t * apq);
    a[p][q] = 0.0;
    a[q][p] = 0.0;

    for (int k = 0; k < n; ++k) {
      if (k == p || k == q) {
        continue;
      }
      const double akp = a[k][p];
      const double akq = a[k][q];
      a[k][p] = (c * akp) - (s * akq);
      a[p][k] = a[k][p];
      a[k][q] = (s * akp) + (c * akq);
      a[q][k] = a[k][q];
    }

    for (int k = 0; k < n; ++k) {
      const double vkp = v[k][p];
      const double vkq = v[k][q];
      v[k][p] = (c * vkp) - (s * vkq);
      v[k][q] = (s * vkp) + (c * vkq);
    }
  }

  std::vector<std::pair<double, double>> node_weight_pairs(n);
  for (int i = 0; i < n; ++i) {
    const double node = a[i][i];
    const double weight = std::sqrt(kPi) * v[0][i] * v[0][i];
    node_weight_pairs[i] = std::make_pair(node, weight);
  }
  std::sort(node_weight_pairs.begin(), node_weight_pairs.end(),
            [](const auto &lhs, const auto &rhs) { return lhs.first < rhs.first; });

  nodes->assign(n, 0.0);
  weights->assign(n, 0.0);
  for (int i = 0; i < n; ++i) {
    (*nodes)[i] = node_weight_pairs[i].first;
    (*weights)[i] = node_weight_pairs[i].second;
  }
}

} // namespace

ModeTables::ModeTables(ModeConfig config) : config_(config) {
  if (config_.n_modes <= 0) {
    throw std::invalid_argument("n_modes must be positive");
  }
  if (config_.lambda <= 0.0) {
    throw std::invalid_argument("lambda must be positive");
  }

  if (config_.identity_mode0) {
    if (config_.n_modes != 1) {
      throw std::invalid_argument("identity_mode0 requires n_modes == 1");
    }
    // Hard controlled-limit identity basis:
    // store mode-0 directly as the physical coefficient with no quadrature
    // reconstruction/projection scaling.
    config_.n_quadrature = 1;
    nodes_ = {0.0};
    weights_ = {1.0};
    phi_ = {1.0};
    dphi_ = {0.0};
    mass_squared_ = {0.0};
    ComputeGramErrors();
    return;
  }

  if (config_.n_quadrature < config_.n_modes) {
    throw std::invalid_argument("n_quadrature must be >= n_modes");
  }

  BuildGaussHermiteTables();
  BuildBasisTables();

  mass_squared_.assign(config_.n_modes, 0.0);
  const double lambda_sq = config_.lambda * config_.lambda;
  for (int n = 0; n < config_.n_modes; ++n) {
    mass_squared_[n] = (2.0 * static_cast<double>(n)) / lambda_sq;
  }

  ComputeGramErrors();
}

void ModeTables::BuildGaussHermiteTables() {
  std::vector<double> standard_nodes;
  std::vector<double> standard_weights;
  ComputeGaussHermiteStandard(config_.n_quadrature, &standard_nodes,
                              &standard_weights);

  nodes_.assign(config_.n_quadrature, 0.0);
  weights_.assign(config_.n_quadrature, 0.0);

  for (int q = 0; q < config_.n_quadrature; ++q) {
    nodes_[q] = config_.lambda * standard_nodes[q];
    weights_[q] = config_.lambda * standard_weights[q];
  }
}

void ModeTables::BuildBasisTables() {
  const int total_size = config_.n_modes * config_.n_quadrature;
  phi_.assign(total_size, 0.0);
  dphi_.assign(total_size, 0.0);

  const double normalizer_factor = config_.lambda * std::sqrt(kPi);

  for (int n = 0; n < config_.n_modes; ++n) {
    const double two_to_n = std::ldexp(1.0, n);
    const double n_factorial = std::tgamma(static_cast<double>(n) + 1.0);
    const double normalization =
        std::sqrt(normalizer_factor * two_to_n * n_factorial);

    for (int q = 0; q < config_.n_quadrature; ++q) {
      double hn = 0.0;
      double hnm1 = 0.0;
      const double x = nodes_[q] / config_.lambda;
      HermitePhysicists(n, x, &hn, &hnm1);
      PhiAt(n, q) = hn / normalization;

      if (n == 0) {
        DPhiAt(n, q) = 0.0;
      } else {
        DPhiAt(n, q) = (std::sqrt(2.0 * static_cast<double>(n)) / config_.lambda) *
                       Phi(n - 1, q);
      }
    }
  }
}

void ModeTables::ComputeGramErrors() {
  gram_diag_max_abs_error_ = 0.0;
  gram_offdiag_max_abs_ = 0.0;

  for (int n = 0; n < config_.n_modes; ++n) {
    for (int m = 0; m < config_.n_modes; ++m) {
      double gram_nm = 0.0;
      for (int q = 0; q < config_.n_quadrature; ++q) {
        gram_nm += weights_[q] * Phi(n, q) * Phi(m, q);
      }
      if (n == m) {
        gram_diag_max_abs_error_ =
            std::max(gram_diag_max_abs_error_, std::abs(gram_nm - 1.0));
      } else {
        gram_offdiag_max_abs_ = std::max(gram_offdiag_max_abs_, std::abs(gram_nm));
      }
    }
  }
}

double ModeTables::Phi(const int n, const int q) const {
  if (n < 0 || n >= config_.n_modes) {
    throw std::out_of_range("mode index n out of range in Phi(n,q)");
  }
  if (q < 0 || q >= config_.n_quadrature) {
    throw std::out_of_range("quadrature index q out of range in Phi(n,q)");
  }
  return phi_[n * config_.n_quadrature + q];
}

double ModeTables::DPhi(const int n, const int q) const {
  if (n < 0 || n >= config_.n_modes) {
    throw std::out_of_range("mode index n out of range in DPhi(n,q)");
  }
  if (q < 0 || q >= config_.n_quadrature) {
    throw std::out_of_range("quadrature index q out of range in DPhi(n,q)");
  }
  return dphi_[n * config_.n_quadrature + q];
}

std::vector<double> ModeTables::ModesToNodes(
    const std::vector<double> &mode_values) const {
  if (static_cast<int>(mode_values.size()) != config_.n_modes) {
    throw std::invalid_argument("mode_values size must equal n_modes");
  }

  std::vector<double> node_values(config_.n_quadrature, 0.0);
  for (int q = 0; q < config_.n_quadrature; ++q) {
    double value = 0.0;
    for (int n = 0; n < config_.n_modes; ++n) {
      value += mode_values[n] * Phi(n, q);
    }
    node_values[q] = value;
  }
  return node_values;
}

std::vector<double> ModeTables::NodesToModes(
    const std::vector<double> &node_values) const {
  if (static_cast<int>(node_values.size()) != config_.n_quadrature) {
    throw std::invalid_argument("node_values size must equal n_quadrature");
  }

  std::vector<double> mode_values(config_.n_modes, 0.0);
  for (int n = 0; n < config_.n_modes; ++n) {
    double value = 0.0;
    for (int q = 0; q < config_.n_quadrature; ++q) {
      value += weights_[q] * Phi(n, q) * node_values[q];
    }
    mode_values[n] = value;
  }
  return mode_values;
}

double ModeTables::ApplyID3Raising(const std::vector<double> &mode_values,
                                   const int n) const {
  if (static_cast<int>(mode_values.size()) != config_.n_modes) {
    throw std::invalid_argument("mode_values size must equal n_modes");
  }
  if (n < 0 || n >= config_.n_modes) {
    throw std::out_of_range("mode index n out of range in ApplyID3Raising");
  }
  const int m = n + 1;
  if (m >= config_.n_modes) {
    return 0.0;
  }
  return (std::sqrt(2.0 * static_cast<double>(m)) / config_.lambda) * mode_values[m];
}

double ModeTables::ApplyID4Lowering(const std::vector<double> &mode_values,
                                    const int n) const {
  if (static_cast<int>(mode_values.size()) != config_.n_modes) {
    throw std::invalid_argument("mode_values size must equal n_modes");
  }
  if (n < 0 || n >= config_.n_modes) {
    throw std::out_of_range("mode index n out of range in ApplyID4Lowering");
  }
  const int m = n - 1;
  if (m < 0) {
    return 0.0;
  }
  return -(std::sqrt(2.0 * static_cast<double>(n)) / config_.lambda) * mode_values[m];
}

double ModeTables::ApplyID7DoubleRaising(const std::vector<double> &mode_values,
                                         const int n) const {
  if (static_cast<int>(mode_values.size()) != config_.n_modes) {
    throw std::invalid_argument("mode_values size must equal n_modes");
  }
  if (n < 0 || n >= config_.n_modes) {
    throw std::out_of_range("mode index n out of range in ApplyID7DoubleRaising");
  }
  const int m = n + 2;
  if (m >= config_.n_modes) {
    return 0.0;
  }
  const double numerator = 2.0 * std::sqrt(static_cast<double>((n + 1) * (n + 2)));
  return (numerator / (config_.lambda * config_.lambda)) * mode_values[m];
}

double &ModeTables::PhiAt(const int n, const int q) {
  return phi_[n * config_.n_quadrature + q];
}

double &ModeTables::DPhiAt(const int n, const int q) {
  return dphi_[n * config_.n_quadrature + q];
}

} // namespace Modes4D

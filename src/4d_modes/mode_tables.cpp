#include "mode_tables.hpp"

#include <cmath>
#include <stdexcept>

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

  nodes->assign(n, 0.0);
  weights->assign(n, 0.0);

  const int m = (n + 1) / 2;
  constexpr double kNewtonTolerance = 1.0e-14;
  constexpr int kNewtonMaxIters = 32;

  double z = 0.0;
  const double n_as_double = static_cast<double>(n);

  for (int i = 0; i < m; ++i) {
    if (i == 0) {
      z = std::sqrt((2.0 * n_as_double) + 1.0) -
          (1.85575 * std::pow((2.0 * n_as_double) + 1.0, -1.0 / 6.0));
    } else if (i == 1) {
      z -= 1.14 * std::pow(n_as_double, 0.426) / z;
    } else if (i == 2) {
      z = (1.86 * z) - (0.86 * (*nodes)[0]);
    } else if (i == 3) {
      z = (1.91 * z) - (0.91 * (*nodes)[1]);
    } else {
      z = (2.0 * z) - (*nodes)[i - 2];
    }

    for (int iteration = 0; iteration < kNewtonMaxIters; ++iteration) {
      double hn = 0.0;
      double hnm1 = 0.0;
      HermitePhysicists(n, z, &hn, &hnm1);
      const double derivative = 2.0 * n_as_double * hnm1;
      const double z_new = z - (hn / derivative);
      if (std::abs(z_new - z) < kNewtonTolerance) {
        z = z_new;
        break;
      }
      z = z_new;
    }

    double hn = 0.0;
    double hnm1 = 0.0;
    HermitePhysicists(n, z, &hn, &hnm1);

    const double log_numerator = ((n - 1) * std::log(2.0)) +
                                 std::lgamma(n_as_double + 1.0) +
                                 (0.5 * std::log(kPi));
    const double denominator = (n_as_double * n_as_double * hnm1 * hnm1);
    const double weight = std::exp(log_numerator) / denominator;

    (*nodes)[i] = -z;
    (*nodes)[n - 1 - i] = z;
    (*weights)[i] = weight;
    (*weights)[n - 1 - i] = weight;
  }
}

} // namespace

ModeTables::ModeTables(ModeConfig config) : config_(config) {
  if (config_.n_modes <= 0) {
    throw std::invalid_argument("n_modes must be positive");
  }
  if (config_.n_quadrature < config_.n_modes) {
    throw std::invalid_argument("n_quadrature must be >= n_modes");
  }
  if (config_.lambda <= 0.0) {
    throw std::invalid_argument("lambda must be positive");
  }

  BuildGaussHermiteTables();
  BuildBasisTables();

  mass_squared_.assign(config_.n_modes, 0.0);
  const double lambda_sq = config_.lambda * config_.lambda;
  for (int n = 0; n < config_.n_modes; ++n) {
    mass_squared_[n] = (2.0 * static_cast<double>(n)) / lambda_sq;
  }
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

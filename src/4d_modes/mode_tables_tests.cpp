#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "mode_tables.hpp"
#include "mode_transform.hpp"

namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;

bool NearlyEqual(const double a, const double b, const double atol, const double rtol) {
  return std::abs(a - b) <= (atol + rtol * std::abs(b));
}

void Require(const bool cond, const std::string &message, bool *ok) {
  if (cond) {
    return;
  }
  *ok = false;
  std::cerr << "FAILED: " << message << "\n";
}

} // namespace

int main() {
  bool ok = true;

  Modes4D::ModeConfig config;
  config.n_modes = 5;
  config.n_quadrature = 8;
  config.lambda = 1.3;

  Modes4D::ModeTables tables(config);

  const auto &w = tables.Weights();
  const auto &x = tables.Nodes();

  double z_int = 0.0;
  for (const double wi : w) {
    z_int += wi;
  }
  Require(NearlyEqual(z_int, config.lambda * std::sqrt(kPi), 1.0e-12, 1.0e-12),
          "Quadrature weights do not recover Z_int = lambda*sqrt(pi)", &ok);

  for (int n = 0; n < config.n_modes; ++n) {
    for (int m = 0; m < config.n_modes; ++m) {
      double ip = 0.0;
      for (int q = 0; q < config.n_quadrature; ++q) {
        ip += w[q] * tables.Phi(n, q) * tables.Phi(m, q);
      }
      const double expected = (n == m) ? 1.0 : 0.0;
      Require(NearlyEqual(ip, expected, 1.0e-11, 1.0e-11),
              "Orthonormality check failed for (n,m)= (" + std::to_string(n) + "," +
                  std::to_string(m) + ")",
              &ok);
    }
  }

  const std::vector<double> mode_values{0.75, -0.25, 0.15, 0.05, -0.02};
  const auto node_values = tables.ModesToNodes(mode_values);
  const auto recovered_modes = tables.NodesToModes(node_values);
  const auto wrapped_node_values = Modes4D::ForwardTransform(tables, mode_values);
  const auto wrapped_modes = Modes4D::InverseTransform(tables, wrapped_node_values);

  for (int n = 0; n < config.n_modes; ++n) {
    Require(NearlyEqual(recovered_modes[n], mode_values[n], 1.0e-11, 1.0e-11),
            "Mode->node->mode round-trip failed at n=" + std::to_string(n), &ok);
    Require(NearlyEqual(wrapped_modes[n], mode_values[n], 1.0e-11, 1.0e-11),
            "Wrapped transform round-trip failed at n=" + std::to_string(n), &ok);
  }

  const auto &m2 = tables.MassSquared();
  const double lambda_sq = config.lambda * config.lambda;
  for (int n = 0; n < config.n_modes; ++n) {
    const double expected = 2.0 * static_cast<double>(n) / lambda_sq;
    Require(NearlyEqual(m2[n], expected, 1.0e-14, 1.0e-14),
            "Mass spectrum mismatch at n=" + std::to_string(n), &ok);
  }

  std::vector<double> g(config.n_quadrature, 0.0);
  std::vector<double> gp(config.n_quadrature, 0.0);
  std::vector<double> gpp(config.n_quadrature, 0.0);
  for (int q = 0; q < config.n_quadrature; ++q) {
    for (int m = 0; m < config.n_modes; ++m) {
      g[q] += mode_values[m] * tables.Phi(m, q);
      gp[q] += mode_values[m] * tables.DPhi(m, q);
      if (m >= 2) {
        const double pref =
            (2.0 * std::sqrt(static_cast<double>(m * (m - 1)))) / lambda_sq;
        gpp[q] += mode_values[m] * pref * tables.Phi(m - 2, q);
      }
    }
  }

  for (int n = 0; n < config.n_modes; ++n) {
    double lhs_id3 = 0.0;
    double lhs_id4 = 0.0;
    double lhs_id7 = 0.0;
    for (int q = 0; q < config.n_quadrature; ++q) {
      const double phi_n = tables.Phi(n, q);
      lhs_id3 += w[q] * phi_n * gp[q];
      lhs_id4 += w[q] * phi_n * (gp[q] - (2.0 * x[q] / lambda_sq) * g[q]);
      lhs_id7 += w[q] * phi_n * gpp[q];
    }

    Require(NearlyEqual(lhs_id3, tables.ApplyID3Raising(mode_values, n), 5.0e-11, 5.0e-11),
            "ID-3 mismatch at n=" + std::to_string(n), &ok);
    Require(NearlyEqual(lhs_id4, tables.ApplyID4Lowering(mode_values, n), 5.0e-11, 5.0e-11),
            "ID-4 mismatch at n=" + std::to_string(n), &ok);
    Require(
        NearlyEqual(lhs_id7, tables.ApplyID7DoubleRaising(mode_values, n), 5.0e-10, 5.0e-10),
        "ID-7 mismatch at n=" + std::to_string(n), &ok);
  }

  if (!ok) {
    return 1;
  }

  std::cout << "modes4d_unit_tests: PASS\n";
  return 0;
}

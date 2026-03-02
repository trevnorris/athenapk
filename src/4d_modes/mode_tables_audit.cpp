#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#include "mode_tables.hpp"
#include "mode_transform.hpp"

namespace {

struct AuditCase {
  int n_modes = 1;
  int n_quadrature = 1;
  double lambda = 1.0;
  double gram_diag_max_abs = 0.0;
  double gram_offdiag_max_abs = 0.0;
  double parseval_max_abs = 0.0;
  double roundtrip_max_rel_l2 = 0.0;
  double roundtrip_rms_rel_l2 = 0.0;
  double wrapped_roundtrip_max_rel_l2 = 0.0;
  double wrapped_roundtrip_rms_rel_l2 = 0.0;
};

struct Config {
  int samples = 512;
  int roundtrips = 128;
  double gram_tol = 5.0e-10;
  double parseval_tol = 5.0e-10;
  double roundtrip_tol = 1.0e-9;
  std::string out_csv;
};

double L2Norm(const std::vector<double> &v) {
  double sum = 0.0;
  for (const double x : v) {
    sum += x * x;
  }
  return std::sqrt(sum);
}

double RelativeL2Error(const std::vector<double> &a, const std::vector<double> &b) {
  if (a.size() != b.size()) return std::numeric_limits<double>::infinity();
  std::vector<double> d(a.size(), 0.0);
  for (std::size_t i = 0; i < a.size(); ++i) {
    d[i] = a[i] - b[i];
  }
  const double den = std::max(L2Norm(a), 1.0e-30);
  return L2Norm(d) / den;
}

void WriteCsv(const std::string &path, const std::vector<AuditCase> &rows) {
  std::ofstream out(path);
  if (!out) {
    throw std::runtime_error("Failed to open csv output: " + path);
  }
  out << "n_modes,n_quadrature,lambda,gram_diag_max_abs,gram_offdiag_max_abs,"
         "parseval_max_abs,roundtrip_max_rel_l2,roundtrip_rms_rel_l2,"
         "wrapped_roundtrip_max_rel_l2,wrapped_roundtrip_rms_rel_l2\n";
  out << std::scientific << std::setprecision(9);
  for (const auto &r : rows) {
    out << r.n_modes << "," << r.n_quadrature << "," << r.lambda << ","
        << r.gram_diag_max_abs << "," << r.gram_offdiag_max_abs << ","
        << r.parseval_max_abs << "," << r.roundtrip_max_rel_l2 << ","
        << r.roundtrip_rms_rel_l2 << "," << r.wrapped_roundtrip_max_rel_l2
        << "," << r.wrapped_roundtrip_rms_rel_l2 << "\n";
  }
}

Config ParseArgs(const int argc, char **argv) {
  Config cfg;
  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    auto need_value = [&](const std::string &name) {
      if (i + 1 >= argc) {
        throw std::runtime_error("Missing value for " + name);
      }
      return std::string(argv[++i]);
    };
    if (arg == "--samples") {
      cfg.samples = std::stoi(need_value(arg));
    } else if (arg == "--roundtrips") {
      cfg.roundtrips = std::stoi(need_value(arg));
    } else if (arg == "--gram-tol") {
      cfg.gram_tol = std::stod(need_value(arg));
    } else if (arg == "--parseval-tol") {
      cfg.parseval_tol = std::stod(need_value(arg));
    } else if (arg == "--roundtrip-tol") {
      cfg.roundtrip_tol = std::stod(need_value(arg));
    } else if (arg == "--out-csv") {
      cfg.out_csv = need_value(arg);
    } else if (arg == "-h" || arg == "--help") {
      std::cout
          << "modes4d_transform_audit options:\n"
          << "  --samples N         random vectors per config (default 512)\n"
          << "  --roundtrips N      repeated mode<->node transforms per sample (default 128)\n"
          << "  --gram-tol X        max allowed Gram error (default 5e-10)\n"
          << "  --parseval-tol X    max allowed Parseval error (default 5e-10)\n"
          << "  --roundtrip-tol X   max allowed rel-L2 roundtrip error (default 1e-9)\n"
          << "  --out-csv PATH      write full per-config table\n";
      std::exit(0);
    } else {
      throw std::runtime_error("Unknown argument: " + arg);
    }
  }
  if (cfg.samples <= 0 || cfg.roundtrips <= 0) {
    throw std::runtime_error("samples and roundtrips must be positive");
  }
  return cfg;
}

} // namespace

int main(int argc, char **argv) {
  try {
    const Config cfg = ParseArgs(argc, argv);

    const std::vector<int> n_modes_grid{1, 2, 3, 4, 6, 8, 12};
    const std::vector<double> lambda_grid{0.5, 1.0, 2.0, 4.0};
    std::vector<AuditCase> rows;
    rows.reserve(n_modes_grid.size() * lambda_grid.size() * 3);

    std::mt19937_64 rng(0x4d4f444553345f44ULL);
    std::normal_distribution<double> normal(0.0, 1.0);

    for (const int n_modes : n_modes_grid) {
      const std::vector<int> nq_grid{
          n_modes, std::max(n_modes + 2, 3), std::max(2 * n_modes, n_modes + 4)};
      for (const int nq : nq_grid) {
        for (const double lambda : lambda_grid) {
          Modes4D::ModeConfig mc;
          mc.n_modes = n_modes;
          mc.n_quadrature = nq;
          mc.lambda = lambda;
          Modes4D::ModeTables tables(mc);

          AuditCase row;
          row.n_modes = n_modes;
          row.n_quadrature = nq;
          row.lambda = lambda;

          const auto &w = tables.Weights();

          for (int n = 0; n < n_modes; ++n) {
            for (int m = 0; m < n_modes; ++m) {
              double ip = 0.0;
              for (int q = 0; q < nq; ++q) {
                ip += w[q] * tables.Phi(n, q) * tables.Phi(m, q);
              }
              if (n == m) {
                row.gram_diag_max_abs =
                    std::max(row.gram_diag_max_abs, std::abs(ip - 1.0));
              } else {
                row.gram_offdiag_max_abs = std::max(row.gram_offdiag_max_abs, std::abs(ip));
              }
            }
          }

          double rt_sq_acc = 0.0;
          double wrt_sq_acc = 0.0;
          for (int s = 0; s < cfg.samples; ++s) {
            std::vector<double> mode0(static_cast<std::size_t>(n_modes), 0.0);
            for (double &v : mode0) v = normal(rng);

            auto mode = mode0;
            auto wrapped_mode = mode0;
            for (int it = 0; it < cfg.roundtrips; ++it) {
              const auto node = tables.ModesToNodes(mode);
              mode = tables.NodesToModes(node);

              const auto wrapped_node = Modes4D::ForwardTransform(tables, wrapped_mode);
              wrapped_mode = Modes4D::InverseTransform(tables, wrapped_node);
            }

            const double rt_err = RelativeL2Error(mode, mode0);
            const double wrt_err = RelativeL2Error(wrapped_mode, mode0);
            row.roundtrip_max_rel_l2 = std::max(row.roundtrip_max_rel_l2, rt_err);
            row.wrapped_roundtrip_max_rel_l2 =
                std::max(row.wrapped_roundtrip_max_rel_l2, wrt_err);
            rt_sq_acc += rt_err * rt_err;
            wrt_sq_acc += wrt_err * wrt_err;

            const auto node = tables.ModesToNodes(mode0);
            double node_energy = 0.0;
            for (int q = 0; q < nq; ++q) node_energy += w[q] * node[q] * node[q];
            double mode_energy = 0.0;
            for (const double v : mode0) mode_energy += v * v;
            row.parseval_max_abs =
                std::max(row.parseval_max_abs, std::abs(node_energy - mode_energy));
          }

          row.roundtrip_rms_rel_l2 = std::sqrt(rt_sq_acc / static_cast<double>(cfg.samples));
          row.wrapped_roundtrip_rms_rel_l2 =
              std::sqrt(wrt_sq_acc / static_cast<double>(cfg.samples));

          rows.push_back(row);
        }
      }
    }

    if (!cfg.out_csv.empty()) {
      WriteCsv(cfg.out_csv, rows);
    }

    double worst_gram = 0.0;
    double worst_parseval = 0.0;
    double worst_rt = 0.0;
    double worst_wrt = 0.0;
    for (const auto &r : rows) {
      worst_gram = std::max(worst_gram, std::max(r.gram_diag_max_abs, r.gram_offdiag_max_abs));
      worst_parseval = std::max(worst_parseval, r.parseval_max_abs);
      worst_rt = std::max(worst_rt, r.roundtrip_max_rel_l2);
      worst_wrt = std::max(worst_wrt, r.wrapped_roundtrip_max_rel_l2);
    }

    const bool pass = (worst_gram <= cfg.gram_tol) &&
                      (worst_parseval <= cfg.parseval_tol) &&
                      (worst_rt <= cfg.roundtrip_tol) &&
                      (worst_wrt <= cfg.roundtrip_tol);

    std::cout << std::scientific << std::setprecision(9);
    std::cout << "audit_cases," << rows.size() << "\n";
    std::cout << "worst_gram_error," << worst_gram << "\n";
    std::cout << "worst_parseval_error," << worst_parseval << "\n";
    std::cout << "worst_roundtrip_rel_l2," << worst_rt << "\n";
    std::cout << "worst_wrapped_roundtrip_rel_l2," << worst_wrt << "\n";
    std::cout << "gram_tol," << cfg.gram_tol << "\n";
    std::cout << "parseval_tol," << cfg.parseval_tol << "\n";
    std::cout << "roundtrip_tol," << cfg.roundtrip_tol << "\n";
    if (!cfg.out_csv.empty()) {
      std::cout << "audit_csv," << cfg.out_csv << "\n";
    }
    std::cout << "audit_status," << (pass ? "PASS" : "FAIL") << "\n";

    return pass ? 0 : 1;
  } catch (const std::exception &e) {
    std::cerr << "modes4d_transform_audit: ERROR: " << e.what() << "\n";
    return 2;
  }
}

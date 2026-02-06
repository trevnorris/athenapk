#ifndef MODES4D_MODE_TABLES_HPP_
#define MODES4D_MODE_TABLES_HPP_

#include <vector>

namespace Modes4D {

struct ModeConfig {
  int n_modes = 1;
  int n_quadrature = 1;
  double lambda = 1.0;
};

class ModeTables {
 public:
  explicit ModeTables(ModeConfig config);

  int NumModes() const { return config_.n_modes; }
  int NumQuadrature() const { return config_.n_quadrature; }
  double Lambda() const { return config_.lambda; }

  const std::vector<double> &Nodes() const { return nodes_; }
  const std::vector<double> &Weights() const { return weights_; }
  const std::vector<double> &MassSquared() const { return mass_squared_; }

  double Phi(int n, int q) const;
  double DPhi(int n, int q) const;

  std::vector<double> ModesToNodes(const std::vector<double> &mode_values) const;
  std::vector<double> NodesToModes(const std::vector<double> &node_values) const;

  // ID-3: <phi_n, d_w G>_Z = sqrt(2(n+1))/lambda * G^(n+1)
  double ApplyID3Raising(const std::vector<double> &mode_values, int n) const;
  // ID-4: int phi_n d_w(Z G) dw = -sqrt(2n)/lambda * G^(n-1)
  double ApplyID4Lowering(const std::vector<double> &mode_values, int n) const;
  // ID-7: <phi_n, d_w^2 G>_Z = 2*sqrt((n+1)(n+2))/lambda^2 * G^(n+2)
  double ApplyID7DoubleRaising(const std::vector<double> &mode_values, int n) const;

 private:
  ModeConfig config_;
  std::vector<double> nodes_;
  std::vector<double> weights_;
  std::vector<double> phi_;
  std::vector<double> dphi_;
  std::vector<double> mass_squared_;

  void BuildGaussHermiteTables();
  void BuildBasisTables();

  double &PhiAt(int n, int q);
  double &DPhiAt(int n, int q);
};

} // namespace Modes4D

#endif // MODES4D_MODE_TABLES_HPP_

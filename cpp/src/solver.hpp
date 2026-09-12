#pragma once

#include <cstddef>
#include <string>
#include <vector>

namespace wavefield {

enum class Boundary { Dirichlet, Mur };

struct Config {
  int nx = 128;
  int ny = 128;
  double dx = 1.0;
  double dy = 1.0;
  double c = 1.0;
  double cfl = 0.5;
  double dt = -1.0;  // <0 → auto from cfl
  Boundary boundary = Boundary::Mur;
  double damping = 0.0;

  double resolved_dt() const;
  double courant() const;
  bool cfl_ok() const;
};

class Simulator {
 public:
  explicit Simulator(Config cfg);

  void reset();
  void add_gaussian(double x0, double y0, double amplitude = 1.0,
                    double sigma = 3.0);
  void step();
  void run(int steps);

  const std::vector<double>& field() const { return u_; }
  int nx() const { return cfg_.nx; }
  int ny() const { return cfg_.ny; }
  int step_index() const { return step_; }
  double max_abs() const;
  double energy() const;
  const Config& config() const { return cfg_; }

 private:
  Config cfg_;
  double dt_{};
  double rx2_{};
  double ry2_{};
  double mur_x_{};
  double mur_y_{};
  int step_{0};
  std::vector<double> u_, u_prev_, u_next_;

  inline std::size_t idx(int i, int j) const {
    return static_cast<std::size_t>(j) * static_cast<std::size_t>(cfg_.nx) +
           static_cast<std::size_t>(i);
  }
};

}  // namespace wavefield

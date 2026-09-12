#include "solver.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace wavefield {

double Config::resolved_dt() const {
  if (dt > 0.0) return dt;
  const double h = std::min(dx, dy);
  return cfl * h / (c * std::sqrt(2.0));
}

double Config::courant() const { return c * resolved_dt() / dx; }

bool Config::cfl_ok() const {
  const double h = std::min(dx, dy);
  const double r = c * resolved_dt() / h;
  return r <= 1.0 / std::sqrt(2.0) + 1e-12;
}

Simulator::Simulator(Config cfg) : cfg_(cfg) {
  if (cfg_.nx < 4 || cfg_.ny < 4) {
    throw std::invalid_argument("grid must be at least 4x4");
  }
  dt_ = cfg_.resolved_dt();
  rx2_ = (cfg_.c * dt_ / cfg_.dx) * (cfg_.c * dt_ / cfg_.dx);
  ry2_ = (cfg_.c * dt_ / cfg_.dy) * (cfg_.c * dt_ / cfg_.dy);
  mur_x_ = (cfg_.c * dt_ - cfg_.dx) / (cfg_.c * dt_ + cfg_.dx);
  mur_y_ = (cfg_.c * dt_ - cfg_.dy) / (cfg_.c * dt_ + cfg_.dy);
  const std::size_t n =
      static_cast<std::size_t>(cfg_.nx) * static_cast<std::size_t>(cfg_.ny);
  u_.assign(n, 0.0);
  u_prev_.assign(n, 0.0);
  u_next_.assign(n, 0.0);
}

void Simulator::reset() {
  std::fill(u_.begin(), u_.end(), 0.0);
  std::fill(u_prev_.begin(), u_prev_.end(), 0.0);
  std::fill(u_next_.begin(), u_next_.end(), 0.0);
  step_ = 0;
}

void Simulator::add_gaussian(double x0, double y0, double amplitude,
                             double sigma) {
  const double inv = 1.0 / (2.0 * sigma * sigma);
  for (int j = 0; j < cfg_.ny; ++j) {
    for (int i = 0; i < cfg_.nx; ++i) {
      const double dx = static_cast<double>(i) - x0;
      const double dy = static_cast<double>(j) - y0;
      const double v = amplitude * std::exp(-(dx * dx + dy * dy) * inv);
      u_[idx(i, j)] += v;
      u_prev_[idx(i, j)] += v;
    }
  }
}

void Simulator::step() {
  const int nx = cfg_.nx;
  const int ny = cfg_.ny;

  for (int j = 1; j < ny - 1; ++j) {
    for (int i = 1; i < nx - 1; ++i) {
      const double uij = u_[idx(i, j)];
      const double lap =
          rx2_ * (u_[idx(i + 1, j)] + u_[idx(i - 1, j)]) +
          ry2_ * (u_[idx(i, j + 1)] + u_[idx(i, j - 1)]) -
          2.0 * (rx2_ + ry2_) * uij;
      double next = 2.0 * uij - u_prev_[idx(i, j)] + lap;
      if (cfg_.damping > 0.0) next *= (1.0 - cfg_.damping);
      u_next_[idx(i, j)] = next;
    }
  }

  if (cfg_.boundary == Boundary::Dirichlet) {
    for (int i = 0; i < nx; ++i) {
      u_next_[idx(i, 0)] = 0.0;
      u_next_[idx(i, ny - 1)] = 0.0;
    }
    for (int j = 0; j < ny; ++j) {
      u_next_[idx(0, j)] = 0.0;
      u_next_[idx(nx - 1, j)] = 0.0;
    }
  } else {
    for (int j = 1; j < ny - 1; ++j) {
      u_next_[idx(0, j)] =
          u_[idx(1, j)] + mur_x_ * (u_next_[idx(1, j)] - u_[idx(0, j)]);
      u_next_[idx(nx - 1, j)] = u_[idx(nx - 2, j)] +
                                mur_x_ * (u_next_[idx(nx - 2, j)] -
                                          u_[idx(nx - 1, j)]);
    }
    for (int i = 1; i < nx - 1; ++i) {
      u_next_[idx(i, 0)] =
          u_[idx(i, 1)] + mur_y_ * (u_next_[idx(i, 1)] - u_[idx(i, 0)]);
      u_next_[idx(i, ny - 1)] = u_[idx(i, ny - 2)] +
                                mur_y_ * (u_next_[idx(i, ny - 2)] -
                                          u_[idx(i, ny - 1)]);
    }
    u_next_[idx(0, 0)] = 0.5 * (u_next_[idx(1, 0)] + u_next_[idx(0, 1)]);
    u_next_[idx(nx - 1, 0)] =
        0.5 * (u_next_[idx(nx - 2, 0)] + u_next_[idx(nx - 1, 1)]);
    u_next_[idx(0, ny - 1)] =
        0.5 * (u_next_[idx(1, ny - 1)] + u_next_[idx(0, ny - 2)]);
    u_next_[idx(nx - 1, ny - 1)] =
        0.5 * (u_next_[idx(nx - 2, ny - 1)] + u_next_[idx(nx - 1, ny - 2)]);
  }

  u_prev_.swap(u_);
  u_.swap(u_next_);
  std::fill(u_next_.begin(), u_next_.end(), 0.0);
  ++step_;
}

void Simulator::run(int steps) {
  for (int s = 0; s < steps; ++s) step();
}

double Simulator::max_abs() const {
  double m = 0.0;
  for (double v : u_) m = std::max(m, std::abs(v));
  return m;
}

double Simulator::energy() const {
  double kinetic = 0.0;
  for (std::size_t k = 0; k < u_.size(); ++k) {
    const double d = u_[k] - u_prev_[k];
    kinetic += d * d;
  }
  double potential = 0.0;
  for (int j = 0; j < cfg_.ny; ++j) {
    for (int i = 0; i < cfg_.nx - 1; ++i) {
      const double g = u_[idx(i + 1, j)] - u_[idx(i, j)];
      potential += g * g;
    }
  }
  for (int j = 0; j < cfg_.ny - 1; ++j) {
    for (int i = 0; i < cfg_.nx; ++i) {
      const double g = u_[idx(i, j + 1)] - u_[idx(i, j)];
      potential += g * g;
    }
  }
  return kinetic + rx2_ * potential;
}

}  // namespace wavefield

#include "solver.hpp"
#include <cstdio>

#include <cstdint>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace {

void write_f32_frame(const std::string& path, const std::vector<double>& u) {
  std::ofstream out(path, std::ios::binary);
  for (double v : u) {
    float f = static_cast<float>(v);
    out.write(reinterpret_cast<const char*>(&f), sizeof(float));
  }
}

void usage(const char* argv0) {
  std::cerr
      << "Usage: " << argv0
      << " [--nx N] [--ny N] [--steps N] [--cfl R] [--boundary mur|dirichlet]\n"
      << "       [--dump-dir DIR] [--stride K]\n";
}

}  // namespace

int main(int argc, char** argv) {
  wavefield::Config cfg;
  int steps = 200;
  int stride = 4;
  std::string dump_dir;

  for (int a = 1; a < argc; ++a) {
    std::string arg = argv[a];
    auto need = [&](const char* name) -> std::string {
      if (a + 1 >= argc) {
        std::cerr << "missing value for " << name << "\n";
        std::exit(2);
      }
      return argv[++a];
    };
    if (arg == "--nx")
      cfg.nx = std::stoi(need("--nx"));
    else if (arg == "--ny")
      cfg.ny = std::stoi(need("--ny"));
    else if (arg == "--steps")
      steps = std::stoi(need("--steps"));
    else if (arg == "--cfl")
      cfg.cfl = std::stod(need("--cfl"));
    else if (arg == "--boundary") {
      auto b = need("--boundary");
      if (b == "mur")
        cfg.boundary = wavefield::Boundary::Mur;
      else if (b == "dirichlet")
        cfg.boundary = wavefield::Boundary::Dirichlet;
      else {
        std::cerr << "unknown boundary\n";
        return 2;
      }
    } else if (arg == "--dump-dir")
      dump_dir = need("--dump-dir");
    else if (arg == "--stride")
      stride = std::stoi(need("--stride"));
    else if (arg == "-h" || arg == "--help") {
      usage(argv[0]);
      return 0;
    } else {
      std::cerr << "unknown arg: " << arg << "\n";
      usage(argv[0]);
      return 2;
    }
  }

  wavefield::Simulator sim(cfg);
  sim.add_gaussian(cfg.nx * 0.5, cfg.ny * 0.5, 1.0, 4.0);

  if (!dump_dir.empty()) {
    // meta.bin: nx, ny as little-endian uint32
    std::ofstream meta(dump_dir + "/meta.bin", std::ios::binary);
    std::uint32_t nx = static_cast<std::uint32_t>(cfg.nx);
    std::uint32_t ny = static_cast<std::uint32_t>(cfg.ny);
    meta.write(reinterpret_cast<const char*>(&nx), 4);
    meta.write(reinterpret_cast<const char*>(&ny), 4);
  }

  for (int k = 0; k < steps; ++k) {
    sim.step();
    if (!dump_dir.empty() && k % stride == 0) {
      char name[64];
      std::snprintf(name, sizeof(name), "/frame_%05d.f32", k);
      write_f32_frame(dump_dir + name, sim.field());
    }
  }

  std::cout << "done steps=" << sim.step_index() << " max|u|=" << sim.max_abs()
            << " energy=" << sim.energy() << " cfl=" << cfg.courant()
            << " ok=" << (cfg.cfl_ok() ? "true" : "false") << "\n";
  return 0;
}

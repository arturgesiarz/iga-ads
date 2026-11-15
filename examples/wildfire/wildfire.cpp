// SPDX-FileCopyrightText: 2015 - 2024 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT
//
// Wildfire simulation (IGA-ADS) -- project "Pozar lasow".
// Based on examples/fire (M. Los, M. Paszynski), extended with a scenario
// reproducing the February 2024 megafire in the Valparaiso region (Chile).
//
// Usage:
//   wildfire [scenario] [steps] [dt] [out_every] [threads]
//
//   scenario   baseline | valparaiso | valparaiso_nowind   (default: valparaiso)
//   steps      number of time steps                        (default: 10000)
//   dt         time step                                   (default: 1e-3)
//   out_every  write out_/fuel_/wind_ files every N steps   (default: 100)
//   threads    number of Galois worker threads             (default: 4)

#include <cstdlib>
#include <iostream>
#include <string>

#include "wildfire.hpp"

int main(int argc, char* argv[]) {
    std::string name = argc > 1 ? argv[1] : "valparaiso";
    int nsteps = argc > 2 ? std::atoi(argv[2]) : 10000;
    double dt = argc > 3 ? std::atof(argv[3]) : 1e-3;
    int out_every = argc > 4 ? std::atoi(argv[4]) : 100;
    int threads = argc > 5 ? std::atoi(argv[5]) : 4;

    ads::problems::scenario sc;
    if (name == "baseline") {
        sc = ads::problems::scenario::baseline;
    } else if (name == "valparaiso") {
        sc = ads::problems::scenario::valparaiso;
    } else if (name == "valparaiso_nowind") {
        sc = ads::problems::scenario::valparaiso_nowind;
    } else {
        std::cerr << "Unknown scenario: " << name << "\n"
                  << "Available: baseline | valparaiso | valparaiso_nowind\n";
        return 1;
    }

    std::cout << "Scenario " << name << ", steps " << nsteps << ", dt " << dt
              << ", output every " << out_every << ", threads " << threads << std::endl;

    // 200 x 200 quadratic B-spline elements on [0,100] x [0,100]
    ads::dim_config dim{2, 200, 0.0, 100.0};
    ads::timesteps_config steps{nsteps, dt};
    int ders = 1;

    ads::config_2d c{dim, dim, steps, ders};
    ads::problems::wildfire sim{c, sc, out_every, threads};
    sim.run();
}

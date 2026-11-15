// Oil extraction by hydraulic fracturing (IGA-ADS) -- project "Wydobycie ropy".
// Based on oil2d.cpp from the course materials.
//
// Usage:
//   oil2d [scenario] [elements] [steps] [dt] [out_every] [threads]
//
//   scenario   original_fixed | geo | geo_shale | geo_noinj   (default: geo)
//   elements   number of B-spline elements along each axis     (default: 80)
//   steps      number of time steps                            (default: 40000)
//   dt         time step                                       (default: 5e-9)
//   out_every  write out_*.data every N steps                  (default: 400)
//   threads    number of Galois worker threads                 (default: 4)

#include <cstdlib>
#include <iostream>
#include <string>

#include "oil2d.hpp"

int main(int argc, char* argv[]) {
    std::string name = argc > 1 ? argv[1] : "geo";
    int elems = argc > 2 ? std::atoi(argv[2]) : 80;
    int nsteps = argc > 3 ? std::atoi(argv[3]) : 40000;
    double dt = argc > 4 ? std::atof(argv[4]) : 5e-9;
    int out_every = argc > 5 ? std::atoi(argv[5]) : 400;
    int threads = argc > 6 ? std::atoi(argv[6]) : 4;

    ads::scenario sc;
    if (name == "original_fixed") {
        sc = ads::scenario::original_fixed;
    } else if (name == "geo") {
        sc = ads::scenario::geo;
    } else if (name == "geo_shale") {
        sc = ads::scenario::geo_shale;
    } else if (name == "geo_noinj") {
        sc = ads::scenario::geo_noinj;
    } else {
        std::cerr << "Unknown scenario: " << name << "\n"
                  << "Available: original_fixed | geo | geo_shale | geo_noinj\n";
        return 1;
    }

    std::cout << "Scenario " << name << ", elements " << elems << ", steps " << nsteps << ", dt "
              << dt << ", output every " << out_every << ", threads " << threads << std::endl;

    // quadratic B-splines, `elems` elements along each axis of [0,1]^2
    ads::dim_config dim{2, elems};
    ads::timesteps_config steps{nsteps, dt};
    int ders = 1;

    ads::config_2d c{dim, dim, steps, ders};
    ads::oil2d sim{c, sc, out_every, threads};
    sim.run();
}

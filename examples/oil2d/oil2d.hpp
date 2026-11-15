// Oil extraction by hydraulic fracturing -- non-linear flow in heterogeneous media (2D).
// Based on oil2d.* from the course materials (M. Paszynski, M. Los, IGA-ADS).
//
//   du/dt - div( K(x) exp(mu u) grad u ) = h(x, t, u)    in [0,1]^2 x [0,T]
//   grad u . n = 0                                         on the boundary
//   h = sum_p pumping_strength * phi(|x - x_p|) - sum_s draining_strength * u * phi(|x - x_s|)
//
// The domain is a vertical cross-section of the formation: x - horizontal, y - height
// (y = 1 is the top of the section, y = 0 its bottom).
//
// Changes w.r.t. the original oil2d.hpp:
//   * env_permeability(x, y): own geological map (folded + faulted layered formation,
//     segmented sandstone layers and rows of fracture zones on a shale background),
//   * own placement of pumps (injection wells) and sinks (production wells),
//   * initial state derived from the permeability map (as in the lecture slides),
//   * FIX: forcing() called draining(x, y, u) although the signature is draining(u, x, y),
//     so in the original code the sinks never drained anything,
//   * scenarios selected from the command line, output of the permeability map,
//     well positions and of the integral quantities (production / injection rates).

#ifndef OIL2D_HPP
#define OIL2D_HPP

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "ads/executor/galois.hpp"
#include "ads/output_manager.hpp"
#include "ads/simulation.hpp"

namespace ads {

struct vec2d {
    double x;
    double y;
};

inline vec2d operator-(const vec2d& a, const vec2d& b) {
    return {a.x - b.x, a.y - b.y};
}

inline vec2d operator+(const vec2d& b, const vec2d& a) {
    return {a.x + b.x, a.y + b.y};
}

inline vec2d operator*(double c, const vec2d& a) {
    return {c * a.x, c * a.y};
}

inline double dot(const vec2d& a, const vec2d& b) {
    return a.x * b.x + a.y * b.y;
}

inline double len_sq(const vec2d& a) {
    return dot(a, a);
}

inline double len(const vec2d& a) {
    return std::sqrt(len_sq(a));
}

inline double falloff(double r, double R, double t) {
    if (t < r)
        return 1.0;
    if (t > R)
        return 0.0;
    double h = (t - r) / (R - r);
    return std::pow((h - 1) * (h + 1), 2);
}

// r < R in [0, 1]
inline double bump(double r, double R, double x, double y) {
    double dx = x - 0.5;
    double dy = y - 0.5;
    double t = std::sqrt(dx * dx + dy * dy);
    return falloff(r / 2, R / 2, t);
}

// ---------------------------------------------------------------------------
// Pumps (injection wells) and sinks (production wells)
// ---------------------------------------------------------------------------
struct pumps {
    std::vector<ads::vec2d> sources;  // pumps: locally increase the pressure
    std::vector<ads::vec2d> sinks;    // sinks: locally decrease the pressure (the higher, the faster)

    double radius = 0.15;
    double pumping_strength = 1;
    double draining_strength = 1e5;

    double pumping(double x, double y) const {
        ads::vec2d v{x, y};
        double p = 0;
        for (const auto& pos : sources) {
            double dist = len(v - pos);
            p += pumping_strength * ads::falloff(0, radius, dist);
        }
        return p;
    }

    double draining(double u, double x, double y) const {
        ads::vec2d v{x, y};
        double p = 0;
        for (const auto& pos : sinks) {
            double dist = len(v - pos);
            double s = draining_strength * ads::falloff(0, radius, dist);
            p += u * s;
        }
        return p;
    }

    // draining rate of a single sink (used for the per-well production statistics)
    double draining_single(std::size_t i, double u, double x, double y) const {
        ads::vec2d v{x, y};
        return u * draining_strength * ads::falloff(0, radius, len(v - sinks[i]));
    }

    auto pumping_fun() const {
        return [this](double x, double y) { return pumping(x, y); };
    }
};

// ---------------------------------------------------------------------------
// Geological map
// ---------------------------------------------------------------------------
namespace geology {

// Permeability scale: K = K_SCALE * k_rel, k_rel in [K_SHALE, 1].
// With K_SCALE = 1000 the map takes values in [1, 1000], as K_q in the lecture.
constexpr double K_SCALE = 1000.0;
constexpr double K_SHALE = 1e-3;  // background: shale, practically impermeable

// Structure of the formation: anticline (dome, oil trap) with its crest at x = 0.42
// and a normal fault at x = 0.74 (the right block is shifted down by 0.06).
constexpr double CREST_X = 0.42;
constexpr double DOME_H = 0.07;
constexpr double DOME_W = 0.28;
constexpr double FAULT_X = 0.74;
constexpr double FAULT_THROW = 0.06;
constexpr double FAULT_GAP = 0.012;  // crushed, sealing zone around the fault plane

inline double structure(double x) {
    double dome = DOME_H * std::exp(-std::pow((x - CREST_X) / DOME_W, 2));
    double fault = x > FAULT_X ? -FAULT_THROW : 0.0;
    return dome + fault;
}

// One geological horizon: a row of bars (sandstone segments) or a row of small blocks
// (fracture zones). Every bar/block gets its own permeability from [kmin, kmax].
struct layer {
    const char* name;
    double y;      // depth of the horizon (without folding)
    double thick;  // thickness
    double kmin, kmax;
    double cell;   // horizontal size of one bar / spacing of the blocks
    double piece;  // horizontal size of the permeable part inside one cell
    double shift;  // horizontal shift of the pattern
};

// clang-format off
inline const layer LAYERS[] = {
    // name                 y      thick  kmin   kmax  cell   piece  shift
    {"shallow aquifer",     0.900, 0.030, 0.10,  0.30, 0.160, 0.150, 0.03},
    {"reservoir A",         0.740, 0.040, 0.30,  1.00, 0.120, 0.112, 0.00},
    {"fractures A/B",       0.655, 0.018, 0.04,  0.10, 0.075, 0.020, 0.02},
    {"reservoir B (tight)", 0.565, 0.034, 0.02,  0.08, 0.140, 0.132, 0.05},
    {"fractures B/C",       0.475, 0.018, 0.04,  0.15, 0.075, 0.020, 0.05},
    {"reservoir C",         0.380, 0.040, 0.15,  0.50, 0.120, 0.112, 0.06},
    {"fractures C",         0.290, 0.018, 0.10,  0.30, 0.075, 0.020, 0.01},
    {"basal sandstone",     0.160, 0.034, 0.01,  0.04, 0.200, 0.190, 0.10},
};
// clang-format on
constexpr int N_LAYERS = sizeof(LAYERS) / sizeof(LAYERS[0]);

// deterministic pseudo-random number in [0, 1) for bar `i` of layer `l`
inline double hash01(int l, int i) {
    double s = std::sin(12.9898 * (i + 37) + 78.233 * (l + 11)) * 43758.5453;
    return s - std::floor(s);
}

// relative permeability of layer `l` at (x, y), or a negative value if (x, y) is outside it
inline double layer_k(int l, double x, double y) {
    const layer& L = LAYERS[l];
    double s = (x - L.shift) / L.cell;
    int i = static_cast<int>(std::floor(s));
    double xi = (s - i) * L.cell;            // position inside the cell
    double xc = L.shift + (i + 0.5) * L.cell;  // centre of the cell
    if (std::abs(xi - 0.5 * L.cell) > 0.5 * L.piece)
        return -1;
    // bars do not cross the fault plane
    if ((x > FAULT_X) != (xc > FAULT_X) || std::abs(x - FAULT_X) < FAULT_GAP)
        return -1;
    // bars are flat, the whole bar follows the structure at its centre (stepped layers)
    // plus a small individual vertical offset
    double yc = L.y + structure(xc) + 0.012 * (hash01(l + 100, i) - 0.5);
    if (std::abs(y - yc) > 0.5 * L.thick)
        return -1;
    double r = hash01(l, i);
    return L.kmin * std::pow(L.kmax / L.kmin, r);  // log-uniform in [kmin, kmax]
}

// centre of layer `l` above/below point x (used to put wells inside a layer)
inline vec2d in_layer(int l, double x) {
    const layer& L = LAYERS[l];
    double s = (x - L.shift) / L.cell;
    int i = static_cast<int>(std::floor(s));
    double xc = L.shift + (i + 0.5) * L.cell;
    double yc = L.y + structure(xc) + 0.012 * (hash01(l + 100, i) - 0.5);
    return {x, yc};
}

}  // namespace geology

// own geological map, relative permeability in [K_SHALE, 1]
inline double env_permeability_rel(double x, double y) {
    double k = geology::K_SHALE;
    for (int l = 0; l < geology::N_LAYERS; ++l) {
        k = std::max(k, geology::layer_k(l, x, y));
    }
    return k;
}

// permeability K_q(x) used by the model, values in [1, 1000]
inline double env_permeability(double x, double y) {
    return geology::K_SCALE * env_permeability_rel(x, y);
}

// ---------------------------------------------------------------------------
// Scenarios
// ---------------------------------------------------------------------------
enum class scenario {
    original_fixed,  // original setup (k = 100, 4 wells, bump), only the draining bug fixed
    geo,             // own geology + designed wells
    geo_shale,       // own geology, wells at the same x but perforated in shale
    geo_noinj,       // own geology, production wells only (no injection)
};

struct scenario_params {
    bool own_geology = true;
    pumps wells;
};

inline scenario_params make_scenario(scenario sc) {
    using geology::in_layer;
    scenario_params p;
    if (sc == scenario::original_fixed) {
        p.own_geology = false;
        p.wells.sources = {{0.25, 0.25}, {0.75, 0.75}};
        p.wells.sinks = {{0.25, 0.75}, {0.75, 0.25}};
        return p;  // radius 0.15, strengths 1 and 1e5 as in the original
    }
    // Layer indices: 1 - reservoir A, 5 - reservoir C
    constexpr int A = 1;
    constexpr int C = 5;
    p.wells.radius = 0.05;
    p.wells.pumping_strength = 2e3;
    p.wells.draining_strength = 1e5;
    // production wells (sinks) at the crest of the anticline, where the oil is trapped
    p.wells.sinks = {in_layer(A, geology::CREST_X), in_layer(C, geology::CREST_X)};
    // injection wells (pumps) on both flanks of the anticline, in the same reservoirs,
    // west of the sealing fault
    p.wells.sources = {in_layer(A, 0.10), in_layer(A, 0.66), in_layer(C, 0.10), in_layer(C, 0.66)};

    if (sc == scenario::geo_shale) {
        // the same wells, but perforated 0.06 above the reservoirs, in the shale cap rock
        for (auto& w : p.wells.sinks)
            w.y += 0.06;
        for (auto& w : p.wells.sources)
            w.y += 0.06;
    }
    if (sc == scenario::geo_noinj) {
        p.wells.sources.clear();
    }
    return p;
}

// ---------------------------------------------------------------------------
// Simulation
// ---------------------------------------------------------------------------
class oil2d : public simulation_2d {
private:
    using Base = simulation_2d;
    vector_type u, u_prev;

    galois_executor executor;

    scenario_params params;
    pumps& process;
    int out_every;
    lin::tensor<double, 4> kq;
    output_manager<2> output;
    std::ofstream stats;

    static constexpr double mi = 10;

public:
    oil2d(const config_2d& config, scenario sc, int out_every_, int threads)
    : Base{config}
    , u{shape()}
    , u_prev{shape()}
    , executor{threads}
    , params{make_scenario(sc)}
    , process{params.wells}
    , out_every{out_every_}
    , kq{{x.basis.elements, y.basis.elements, x.basis.quad_order + 1, y.basis.quad_order + 1}}
    , output{x.B, y.B, 200} { }

    // Initial state.
    // Original: small bump in the centre of the domain.
    // Own geology: the initial pressure (oil) is derived from the permeability map, as in the
    // lecture: u0 = 0.1 * K~(x) * theta_{0.2,0.3}(|x - c|), K~ = (K - 1) / (1000 - 1),
    // with c at the crest of the anticline (the trap) and the distance stretched horizontally.
    double init_state(double x, double y) const {
        if (!params.own_geology) {
            double r = 0.1;
            double R = 0.5;
            return 1e-3 * ads::bump(r, R, x, y);
        }
        // K~ is averaged over a small 5x5 stencil (size ~ one element): the projection of
        // a discontinuous function onto B-splines would otherwise produce strong undershoots
        double Kt = 0;
        constexpr int ns = 2;
        constexpr double hs = 0.006;
        for (int i = -ns; i <= ns; ++i) {
            for (int j = -ns; j <= ns; ++j) {
                double K = env_permeability(x + i * hs, y + j * hs);
                Kt += (K - 1) / (geology::K_SCALE - 1);
            }
        }
        Kt /= (2 * ns + 1) * (2 * ns + 1);
        double dx = (x - geology::CREST_X) / 1.5;
        double dy = y - 0.62;
        double d = std::sqrt(dx * dx + dy * dy);
        return 0.1 * Kt * falloff(0.2, 0.3, d);
    }

private:
    void before() override {
        fill_permeability_map();
        prepare_matrices();

        auto init = [this](double x, double y) { return init_state(x, y); };
        projection(u, init);
        solve(u);
        output.to_file(u, "out_%d.data", 0);
        write_maps();
        stats.open("stats.csv");
        stats << "iter,t,mass,energy,umin,umax,injection_rate,production_rate";
        for (std::size_t i = 0; i < process.sinks.size(); ++i)
            stats << ",sink" << i;
        stats << "\n";
        write_stats(0, 0.0);
    }

    void fill_permeability_map() {
        for (auto e : elements()) {
            for (auto q : quad_points()) {
                auto x = point(e, q);
                kq(e[0], e[1], q[0], q[1]) = params.own_geology ? env_permeability(x[0], x[1]) : 1e2;
            }
        }
    }

    // permeability map, initial state function and wells -> text files for the plots
    void write_maps() const {
        std::ofstream perm{"permeability.data"};
        int n = 400;
        for (int i = 0; i <= n; ++i) {
            for (int j = 0; j <= n; ++j) {
                double x = double(i) / n;
                double y = double(j) / n;
                double k = params.own_geology ? env_permeability(x, y) : 1e2;
                perm << x << " " << y << " " << k << " " << init_state(x, y) << "\n";
            }
        }
        std::ofstream w{"wells.data"};
        w << "# type x y radius strength\n";
        for (const auto& p : process.sources)
            w << "pump " << p.x << " " << p.y << " " << process.radius << " " << process.pumping_strength
              << "\n";
        for (const auto& p : process.sinks)
            w << "sink " << p.x << " " << p.y << " " << process.radius << " "
              << process.draining_strength << "\n";
    }

    void before_step(int /*iter*/, double /*t*/) override {
        using std::swap;
        swap(u, u_prev);
    }

    void step(int /*iter*/, double t) override {
        compute_rhs(t);
        solve(u);
    }

    // Explicit Euler:
    // (u_{t+1}, w) = (u_t + dt h, w) - dt (K e^{mu u_t} grad u_t, grad w)
    void compute_rhs(double t) {
        auto& rhs = u;

        zero(rhs);
        executor.for_each(elements(), [&](index_type e) {
            auto U = element_rhs();

            double J = jacobian(e);
            for (auto q : quad_points()) {
                double w = weight(q);
                auto x = point(e, q);

                double k = permeability(e, q);
                value_type u = eval_fun(u_prev, e, q);
                double h = forcing(x, t, u.val);

                for (auto a : dofs_on_element(e)) {
                    auto aa = dof_global_to_local(e, a);
                    value_type v = eval_basis(e, q, a);

                    double val = -k * std::exp(mi * u.val) * grad_dot(u, v) + h * v.val;
                    U(aa[0], aa[1]) += (u.val * v.val + steps.dt * val) * w * J;
                }
            }
            executor.synchronized([&] { update_global_rhs(rhs, U, e); });
        });
    }

    double energy(const vector_type& u) const {
        double E = 0;
        for (auto e : elements()) {
            double J = jacobian(e);
            for (auto q : quad_points()) {
                double w = weight(q);
                value_type a = eval_fun(u, e, q);
                E += a.val * a.val * w * J;
            }
        }
        return E;
    }

    // integral quantities: total "mass" of the pressure field, min/max, injection rate
    // int pumping dx, production rate int draining dx (total and for each sink)
    void write_stats(int iter, double t) {
        double mass = 0, E = 0, inj = 0, prod = 0;
        double umin = 1e300, umax = -1e300;
        std::vector<double> per_sink(process.sinks.size(), 0.0);
        for (auto e : elements()) {
            double J = jacobian(e);
            for (auto q : quad_points()) {
                double w = weight(q) * J;
                auto x = point(e, q);
                double a = eval_fun(u, e, q).val;
                mass += a * w;
                E += a * a * w;
                umin = std::min(umin, a);
                umax = std::max(umax, a);
                inj += process.pumping(x[0], x[1]) * w;
                for (std::size_t i = 0; i < per_sink.size(); ++i) {
                    double s = process.draining_single(i, a, x[0], x[1]) * w;
                    per_sink[i] += s;
                    prod += s;
                }
            }
        }
        char buf[512];
        std::snprintf(buf, sizeof(buf), "%d,%.6e,%.8e,%.8e,%.6e,%.6e,%.8e,%.8e", iter, t, mass, E,
                      umin, umax, inj, prod);
        stats << buf;
        for (double s : per_sink) {
            std::snprintf(buf, sizeof(buf), ",%.8e", s);
            stats << buf;
        }
        stats << "\n";
    }

    void after_step(int iter, double t) override {
        int stats_every = std::max(1, out_every / 10);
        if ((iter + 1) % stats_every == 0) {
            write_stats(iter + 1, t + steps.dt);
        }
        if ((iter + 1) % out_every == 0) {
            std::cout << "Step " << iter + 1 << ", energy: " << energy(u) << std::endl;
            output.to_file(u, "out_%d.data", iter + 1);
        }
    }

    void after() override { stats.flush(); }

    double permeability(index_type e, index_type q) const { return kq(e[0], e[1], q[0], q[1]); }

    double forcing(point_type x, double /*t*/, double u) const {
        // FIX: the original code called process.draining(x[0], x[1], u)
        return process.pumping(x[0], x[1]) - process.draining(u, x[0], x[1]);
    }
};

}  // namespace ads

#endif  // OIL2D_HPP

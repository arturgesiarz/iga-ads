// SPDX-FileCopyrightText: 2015 - 2024 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT
//
// Wildfire model of Grasso & Innocente (2020) solved with IGA-ADS.
// Modified for the project "Pozar lasow": space- and time-dependent wind
// bx(x,y,t), by(x,y,t), custom initial fire configuration (init_state) and a
// custom fuel map (fuel_init) reproducing the February 2024 Valparaiso fire.

#ifndef WILDFIRE_WILDFIRE_HPP
#define WILDFIRE_WILDFIRE_HPP

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>

#include "ads/executor/galois.hpp"
#include "ads/output_manager.hpp"
#include "ads/simulation.hpp"
#include "ads/util.hpp"

// ---------------------------------------------------------------------------
// Geometric helpers
// ---------------------------------------------------------------------------

inline double falloff(double r, double R, double t) {
    if (t < r)
        return 1.0;
    if (t > R)
        return 0.0;
    double h = (t - r) / (R - r);
    return std::pow((h - 1) * (h + 1), 2);
}

// Original bump from examples/fire, generalised to an arbitrary centre (cx, cy).
// It equals 1 for distance < r/2 and decays smoothly to 0 at distance R/2.
inline double bump(double r, double R, double x, double y, double cx = 50, double cy = 50) {
    double dx = x - cx;
    double dy = y - cy;
    double t = std::sqrt(dx * dx + dy * dy) / 100;
    return falloff(r / 200, R / 200, t);
}

// C^1 smooth step from 0 (s <= a) to 1 (s >= b)
inline double smoothstep(double a, double b, double s) {
    double h = std::clamp((s - a) / (b - a), 0.0, 1.0);
    return h * h * (3 - 2 * h);
}

// ~1 inside an ellipse with semi-axes (a, b) rotated by `angle` [deg], ~0 outside;
// `w` is the width of the smooth transition zone
inline double in_ellipse(double x, double y, double cx, double cy, double a, double b,
                         double angle, double w) {
    double c = std::cos(angle * M_PI / 180);
    double s = std::sin(angle * M_PI / 180);
    double dx = x - cx;
    double dy = y - cy;
    double u = (c * dx + s * dy) / a;
    double v = (-s * dx + c * dy) / b;
    double d = (std::sqrt(u * u + v * v) - 1) * std::min(a, b);  // ~ signed distance
    return 1 - smoothstep(-w, w, d);
}

// ~1 within distance `hw` of the segment (x0,y0)-(x1,y1), ~0 further away
inline double near_segment(double x, double y, double x0, double y0, double x1, double y1,
                           double hw, double w) {
    double ex = x1 - x0;
    double ey = y1 - y0;
    double t = std::clamp(((x - x0) * ex + (y - y0) * ey) / (ex * ex + ey * ey), 0.0, 1.0);
    double dx = x - (x0 + t * ex);
    double dy = y - (y0 + t * ey);
    double d = std::sqrt(dx * dx + dy * dy);
    return 1 - smoothstep(hw - w, hw + w, d);
}

// Same as near_segment, but for a polyline given as an array of points
template <std::size_t N>
inline double near_polyline(double x, double y, const double (&px)[N], const double (&py)[N],
                            double hw, double w) {
    double m = 0;
    for (std::size_t i = 0; i + 1 < N; ++i)
        m = std::max(m, near_segment(x, y, px[i], py[i], px[i + 1], py[i + 1], hw, w));
    return m;
}

namespace ads::problems {

enum class scenario {
    baseline,           // original example: no wind, fuel = 1, one fire in the centre
    valparaiso,         // Valparaiso / Vina del Mar, 2-3 Feb 2024
    valparaiso_nowind,  // same fuel map and ignition, but without wind (reference)
};

class wildfire : public simulation_2d {
private:
    using Base = simulation_2d;
    vector_type u, u_prev;

    vector_type fuel, fuel_prev;

    scenario sc;
    int out_every;

    double ch = 1.0;
    double Ar = 5.7e-5;
    double rho = 1.293;
    double cp = 1.0;
    double cw = 0.5;
    double kappa = 0.3;
    double xi = 2e-2;
    double sigma = 5.67e-8;
    double eps = 0.05;  // "percentage error"
    double delta_x = 3.5e-2 / eps;
    double delta_z = 1.5 * eps;
    double T0 = 300;
    double Tcomb = 1200;
    double Tig = 800;
    double Ta = 300;  // activation temperature

    double M = 2;
    double M1 = 1;

    // scenario parameters (wind and fuel calibration)
    static constexpr double WIND_SPEED = 60;  // reference |b|
    // time schedule of the wind: t = 0 at 12:00 (ignition), 1 time unit
    // ~ 3.5 h (calibrated with the arrival of the fire in El Olivar ~18:30)
    static constexpr double T_PEAK = 1.7;   // ~18:00 strongest wind (48 km/h)
    static constexpr double T_VEER0 = 0.6;  // veering from SSW to S-SSE
    static constexpr double T_VEER1 = 1.6;
    static constexpr double T_CALM0 = 2.0;  // ~19:00-21:00 wind dies down
    static constexpr double T_CALM1 = 2.6;
    static constexpr double WIND_DIR0 = 80;   // initial direction [deg] (from SSW)
    static constexpr double WIND_DIR1 = 97;   // direction after veering [deg] (from S-SSE)
    static constexpr double FUEL_SCRUB = 0.3;  // background vegetation
    static constexpr double FUEL_MOIST = 0.22;  // green, moist vegetation

    galois_executor executor;
    output_manager<2> output;

public:
    wildfire(const config_2d& config, scenario sc, int out_every, int threads)
    : Base{config}
    , u{shape()}
    , u_prev{shape()}
    , fuel{shape()}
    , fuel_prev{shape()}
    , sc{sc}
    , out_every{out_every}
    , executor{threads}
    , output{x.B, y.B, 200} { }

    // -----------------------------------------------------------------------
    // Map of the scenario.  Domain [0,100] x [0,100], 1 unit = 250 m, so the
    // domain covers 25 km x 25 km: from 6 km W to 19 km E and from 21 km S to
    // 4 km N of the centre of Vina del Mar (X = 4 (x_km + 6), Y = 4 (y_km + 21)).
    // x grows to the East, y to the North.  Geometry traced from the Landsat 9
    // image of the burn scars (5 Feb 2024) and coordinates of the places.
    // -----------------------------------------------------------------------

    // > 0 over the Pacific Ocean (NW corner of the domain)
    static double sea_dist(double x, double y) { return std::min(y - (72 + 0.7 * x), 19 - x); }

    // -----------------------------------------------------------------------
    // Wind  b = (bx, by)  [model units]
    // -----------------------------------------------------------------------
    double bx(double x, double y, double t) const {
        if (sc != scenario::valparaiso)
            return 0;
        double s, d;
        wind(x, y, t, s, d);
        return s * std::cos(d);
    }

    double by(double x, double y, double t) const {
        if (sc != scenario::valparaiso)
            return 0;
        double s, d;
        wind(x, y, t, s, d);
        return s * std::sin(d);
    }

    // Wind speed s and direction d (angle of the vector b measured from the
    // x-axis, i.e. the direction the wind blows TO; 90 deg = towards North)
    void wind(double x, double y, double t, double& s, double& d) const {
        // 2 Feb 2024 (station Rodelillo): hot and dry wind from S-SSW
        // (160-210 deg), increasing during the afternoon (33-35 km/h at
        // 16-17 h, gusts 48-50 km/h at 18 h) and calming down at night
        // (4-7 km/h).  Direction: from SSW (blowing to NNE) veering to S-SSE
        // (blowing to N-NNW) later in the afternoon.
        double dir = WIND_DIR0 + (WIND_DIR1 - WIND_DIR0) * smoothstep(T_VEER0, T_VEER1, t);
        // local meandering of the direction (+-6 deg), travelling with the flow
        dir += 6 * std::sin(2 * M_PI * (x / 40.0 - t / 0.8));

        double ramp = 0.55 + 0.45 * smoothstep(0.0, T_PEAK, t);     // afternoon build-up
        double calm = 1 - 0.8 * smoothstep(T_CALM0, T_CALM1, t);     // evening calm
        // gusts: waves travelling northwards with the flow
        double gust = 1 + 0.2 * std::sin(2 * M_PI * (t / 0.3 - y / 30.0));
        // topography: channelling along the N-S ravines (quebradas)
        double ridge = 1 + 0.12 * std::sin(2 * M_PI * x / 25) * std::cos(2 * M_PI * y / 50);
        // weaker wind close to the coast (sea breeze, buildings)
        double sea = std::exp(std::min(0.0, sea_dist(x, y)) / 8.0);

        // numerical "sponge": no wind within 1 km of the domain boundary
        // (the explicit scheme is unstable at the inflow boundary otherwise)
        double edge = smoothstep(0.0, 4.0, std::min({x, 100 - x, y, 100 - y}));

        s = WIND_SPEED * ramp * calm * gust * ridge * (1 - 0.5 * sea) * edge;
        d = dir * M_PI / 180;
    }

    // -----------------------------------------------------------------------
    // Fuel map, values in (0,1): 1 = highly flammable, 0 = nothing to burn
    // (the model burns only where fuel > 0.2)
    // -----------------------------------------------------------------------
    double fuel_init(double x, double y) const {
        if (sc == scenario::baseline)
            return 1;

        double f = FUEL_SCRUB;  // sclerophyll scrub, dry grass, farmland

        // west of Ruta 68 / Lago Penuelas: green, moister vegetation of the
        // coastal range (bright green on the Landsat image) - hardly flammable
        double west = smoothstep(-3, 3, (36 - x) - 0.45 * (y - 30));
        f = f * (1 - west) + FUEL_MOIST * west;

        // eucalyptus / pine / acacia forest between Las Tablas and Vina del
        // Mar (Reserva Lago Penuelas, fundos, Santuario El Palmar)
        f = std::max(f, in_ellipse(x, y, 48, 45, 40, 14, 92, 4));
        // ravines (quebradas El Quiteno, El Salto) with dense vegetation
        f = std::max(f, in_ellipse(x, y, 41, 66, 12, 2.5, 100, 1.5));
        f = std::max(f, in_ellipse(x, y, 52, 70, 10, 2.5, 80, 1.5));

        // wildland-urban interface: El Olivar, Villa Independencia,
        // Achupallas - self-built houses on wooded hills
        double wui = std::max(in_ellipse(x, y, 45, 86, 8, 4.5, 0, 1.5),
                              in_ellipse(x, y, 37, 92, 5, 3, 20, 1.5));
        f = f * (1 - wui) + 0.55 * wui;

        // dense urban areas: Vina del Mar centre and Renaca, the Marga Marga
        // valley (Chorrillos - Quilpue - Villa Alemana), Placilla / Curauma
        // (between Chorrillos and Quilpue the valley - El Salto, Jardin
        // Botanico Nacional - is vegetated, so the city band has a gap there)
        static const double mmx[] = {19, 27, 35};
        static const double mmy[] = {86, 83, 80};
        static const double qux[] = {58, 65, 90, 100};
        static const double quy[] = {74, 73, 75, 76};
        double city = std::max(near_polyline(x, y, mmx, mmy, 3, 1.5),
                               near_polyline(x, y, qux, quy, 3, 1.5));
        city = std::max(city, in_ellipse(x, y, 25, 92, 6, 10, 0, 1.5));
        city = std::max(city, in_ellipse(x, y, 66, 71, 7, 4, 10, 1.5));
        city = std::max(city, in_ellipse(x, y, 91, 75, 8, 4, 0, 1.5));
        city = std::max(city, in_ellipse(x, y, 22, 46, 6, 6, 0, 1.5));
        f = f * (1 - city) + 0.1 * city;

        // Estero Marga Marga (stream along the valley)
        static const double emx[] = {20, 35, 46, 58, 70, 100};
        static const double emy[] = {87, 81, 78.5, 75, 74, 75};
        f *= 1 - 0.9 * near_polyline(x, y, emx, emy, 0.25, 0.25);

        // Lago Penuelas reservoir
        f *= 1 - in_ellipse(x, y, 31, 28, 6, 2.5, -30, 1);

        // Ruta 68 (Valparaiso - Santiago) and Ruta 60-CH (Troncal Sur):
        // narrow fuel-free corridors (~70-100 m)
        static const double r68x[] = {0, 18, 25, 36, 42, 48};
        static const double r68y[] = {70, 52, 44, 30, 15, 0};
        static const double r60x[] = {20, 40, 65, 100};
        static const double r60y[] = {84, 77, 71, 73};
        f *= 1 - 0.9 * near_polyline(x, y, r68x, r68y, 0.25, 0.25);
        f *= 1 - 0.9 * near_polyline(x, y, r60x, r60y, 0.25, 0.25);

        // Pacific Ocean
        f *= 1 - smoothstep(-1, 1, sea_dist(x, y));

        return std::clamp(f, 0.0, 1.0);
    }

    // -----------------------------------------------------------------------
    // Initial temperature [K]
    // -----------------------------------------------------------------------
    double init_state(double x, double y) const {
        if (sc == scenario::baseline) {
            double r = 10;
            double R = 30;
            return T0 + Tcomb * bump(r, R, x, y);
        }
        // 2 Feb 2024, ~12:00: four fires set within minutes in Fundo Las
        // Tablas, next to Ruta 68 (km 91), south of the Penuelas reserve:
        // Las Tablas, La Isla, Animita, Las Taguas
        double b = 0;
        b = std::max(b, bump(4, 10, x, y, 40.5, 10));   // Las Tablas
        b = std::max(b, bump(3, 8, x, y, 45, 13));      // La Isla
        b = std::max(b, bump(3, 8, x, y, 37, 14));      // Animita
        b = std::max(b, bump(2, 6, x, y, 43, 7));       // Las Taguas
        return T0 + Tcomb * b;
    }

private:
    void before() override {
        prepare_matrices();

        auto init = [this](double x, double y) { return init_state(x, y); };
        projection(u, init);
        solve(u);
        output.to_file(u, "out_0.data");

        auto fuel_fun = [this](double x, double y) { return fuel_init(x, y); };
        projection(fuel, fuel_fun);
        solve(fuel);
        output.to_file(fuel, "fuel_0.data");

        write_wind(0, 0.0);
    }

    void before_step(int /*iter*/, double /*t*/) override {
        using std::swap;
        swap(u, u_prev);
        swap(fuel, fuel_prev);
    }

    void step(int /*iter*/, double t) override {
        compute_rhs(t);
        solve(u);
        solve(fuel);
    }

    void compute_rhs(double t) {
        auto& rhs = u;
        auto& rhs_fuel = fuel;

        zero(rhs);
        zero(rhs_fuel);
        executor.for_each(elements(), [&](index_type e) {
            auto U = element_rhs();
            auto F = element_rhs();

            double J = jacobian(e);
            for (auto q : quad_points()) {
                double w = weight(q);
                auto p = point(e, q);
                double px = p[0];
                double py = p[1];
                // space- and time-dependent wind
                double wx = bx(px, py, t);
                double wy = by(px, py, t);

                value_type u = eval_fun(u_prev, e, q);
                value_type fuel = eval_fun(fuel_prev, e, q);

                for (auto a : dofs_on_element(e)) {
                    auto aa = dof_global_to_local(e, a);
                    value_type v = eval_basis(e, q, a);
                    double inv = 1.0 / (rho * cp);

                    double delta = u.val > Tig && fuel.val > 0.2 ? 1.0 : 0.0;
                    double hc = -70;  // enthalpy
                    // T+ = max(T, 0): guards T^3, T^4 and exp(-Ta/T) against the
                    // (non-physical) undershoots of the explicit scheme, which
                    // otherwise make the radiative "conductivity" negative
                    double Tp = std::max(u.val, 0.0);
                    double r = delta > 0 ? Ar * Tp * std::exp(-Ta / Tp) : 0.0;
                    double Rc = -1e4 * rho * ch * hc * M / M1 * r;
                    double Qw = -rho * cw * (wx * u.dx + wy * u.dy);
                    double qc = -kappa * grad_dot(u, v);
                    double qd = 0.0;  // omitted
                    double qr = -4 * sigma * eps * delta_x * std::pow(Tp, 3) * grad_dot(u, v);
                    double Qconv = xi * (T0 - u.val);
                    double Qrz = sigma * eps / delta_z * (std::pow(T0, 4) - std::pow(Tp, 4));

                    double val = (Rc + Qw + Qconv + Qrz) * v.val + qc + qd + qr;
                    U(aa[0], aa[1]) += (u.val * v.val + steps.dt * inv * val) * w * J;

                    double fval = -delta * 3e2 * r * fuel.val * v.val;
                    F(aa[0], aa[1]) += (fuel.val * v.val + steps.dt * fval) * w * J;
                }
            }
            executor.synchronized([&] { update_global_rhs(rhs, U, e); });
            executor.synchronized([&] { update_global_rhs(rhs_fuel, F, e); });
        });
    }

    void after_step(int iter, double t) override {
        auto i = iter + 1;
        if (i % out_every == 0) {
            std::cout << "Step " << i << ", t = " << t + steps.dt << std::endl;
            output.to_file(u, "out_%d.data", i);
            output.to_file(fuel, "fuel_%d.data", i);
            write_wind(i, t + steps.dt);
        }
    }

    // Wind field sampled on a coarse 41 x 41 grid: x y bx by
    void write_wind(int i, double t) const {
        char name[64];
        std::snprintf(name, sizeof(name), "wind_%d.data", i);
        std::ofstream os{name};
        os << std::setprecision(8);
        int n = 40;
        for (int ix = 0; ix <= n; ++ix) {
            for (int iy = 0; iy <= n; ++iy) {
                double px = 100.0 * ix / n;
                double py = 100.0 * iy / n;
                os << px << ' ' << py << ' ' << bx(px, py, t) << ' ' << by(px, py, t) << '\n';
            }
            os << '\n';
        }
    }
};

}  // namespace ads::problems

#endif  // WILDFIRE_WILDFIRE_HPP

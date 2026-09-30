// Issue #461 — N4 convection scheme verification.
//
// Verifies the convection face schemes on a uniform 1-D slab grid with a
// uniform cross-section flux:
//   - boundedness: TVD (limited) face values stay inside the adjacent-cell
//     envelope on a discontinuous profile, for every registered limiter;
//   - linear exactness: all TVD limiters reproduce a linear field exactly on
//     interior faces (psi(1)=1 with an exact gradient);
//   - observed order: first-order upwind vs the TVD limiters on a smooth
//     monotone field;
//   - conservation: the discrete divergence sums to the net boundary flux.

#include "cfdx/core/numerics/convection.h"
#include "cfdx/core/mesh/mesh.h"
#include "cfdx/core/geometry/geometry_cache.h"
#include "cfdx/core/field/field.h"
#include "verification_metrics.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace cfdx::core;
using namespace cfdx::verification;

namespace {

struct Grid {
    Mesh mesh;
    GeometryCache geometry;
    std::vector<std::size_t> cross;              // cross-section face ids, planes 0..n
    std::vector<std::size_t> interior_faces;     // internal cross-section faces i in [1, n-1]
    std::vector<std::size_t> deep_faces;         // i in [2, n-2] (upwind neighbour interior)
};

Grid make_1d_slab_grid(std::size_t n)
{
    if (n < 4) throw std::invalid_argument("make_1d_slab_grid: n must be >= 4");
    const double h = 1.0 / static_cast<double>(n);

    Mesh m;
    m.points().resize(4 * (n + 1));
    const double yz[4][2] = {{0,0},{1,0},{1,1},{0,1}};
    for (std::size_t i = 0; i <= n; ++i)
        for (unsigned q = 0; q < 4; ++q)
            m.points().set(4 * i + q, i * h, yz[q][0], yz[q][1]);

    std::vector<std::vector<std::size_t>> cell_faces(n);
    const auto cross_face = [&](std::size_t i) {
        const std::size_t id = m.faces().n_faces();
        m.faces().push_face({4 * i + 0, 4 * i + 1, 4 * i + 2, 4 * i + 3});
        return id;
    };
    const auto side_face = [&](std::initializer_list<std::size_t> v) {
        const std::size_t id = m.faces().n_faces();
        m.faces().push_face(std::vector<FaceIndex>(v.begin(), v.end()));
        return id;
    };

    Grid grid;
    grid.cross.resize(n + 1);
    for (std::size_t i = 0; i <= n; ++i) grid.cross[i] = cross_face(i);

    for (std::size_t c = 0; c < n; ++c) {
        const std::size_t a = 4 * c, b = 4 * (c + 1);
        cell_faces[c] = {grid.cross[c], grid.cross[c + 1],
            side_face({a + 0, b + 0, b + 1, a + 1}),
            side_face({a + 3, a + 2, b + 2, b + 3}),
            side_face({a + 0, a + 3, b + 3, b + 0}),
            side_face({a + 1, b + 1, b + 2, a + 2})};
    }

    m.ownership().resize(m.n_faces());
    for (std::size_t i = 0; i <= n; ++i) {
        const std::size_t owner = (i == 0) ? 0 : i - 1;
        m.ownership().set_owner(grid.cross[i], owner);
        if (i == 0 || i == n) {
            m.ownership().set_neighbour(grid.cross[i], FaceOwnership::BOUNDARY);
        } else {
            m.ownership().set_neighbour(grid.cross[i], static_cast<std::int64_t>(i));
            grid.interior_faces.push_back(grid.cross[i]);
        }
        if (i >= 2 && i <= n - 2) grid.deep_faces.push_back(grid.cross[i]);
    }
    std::size_t face_id = n + 1;
    for (std::size_t c = 0; c < n; ++c)
        for (unsigned q = 0; q < 4; ++q, ++face_id) {
            m.ownership().set_owner(face_id, c);
            m.ownership().set_neighbour(face_id, FaceOwnership::BOUNDARY);
        }
    for (std::size_t c = 0; c < n; ++c) m.cells().push_cell(cell_faces[c]);

    grid.mesh = std::move(m);
    grid.geometry = make_geometry_cache(grid.mesh);
    return grid;
}

void require(bool c, const std::string& m) { if (!c) throw std::runtime_error(m); }

std::vector<double> sample_value(const Grid& grid, double (*fn)(double x)) {
    std::vector<double> v(grid.mesh.n_cells(), 0.0);
    for (std::size_t c = 0; c < grid.mesh.n_cells(); ++c)
        v[c] = fn(grid.geometry.cell_centres[c].x);
    return v;
}

Field<double, Location::CELL> cell_field(const Grid& grid, const std::vector<double>& values)
{
    Field<double, Location::CELL> f(grid.mesh.n_cells(), "phi", "1", 1);
    for (std::size_t c = 0; c < grid.mesh.n_cells(); ++c) f(c) = values[c];
    return f;
}

Field<double, Location::FACE> uniform_flux(const Grid& grid, double f)
{
    Field<double, Location::FACE> flux(grid.mesh.n_faces(), "F", "m2/s", 1);
    for (std::size_t id = 0; id < grid.mesh.n_faces(); ++id) flux(id) = 0.0;
    flux(grid.cross.front()) = f;                       // inflow at the left boundary
    flux(grid.cross.back()) = -f;                       // outflow at the right boundary
    for (std::size_t i = 1; i + 1 < grid.cross.size(); ++i)
        flux(grid.cross[i]) = f;                        // internal owner->neighbour, left->right
    return flux;
}

// L2 error of the face values against an analytic profile at the face centres,
// over the deep interior cross-section faces.
ErrorMetrics profile_face_error(const Grid& grid, const Field<double, Location::FACE>& face,
                                double (*fn)(double x))
{
    std::vector<double> got, want;
    for (const std::size_t i : grid.deep_faces) {
        const std::size_t f = grid.cross[i];
        got.push_back(face(f));
        want.push_back(fn(grid.geometry.face_centres[f].x));
    }
    return error_norms(got, want);
}

void check_boundedness()
{
    const Grid grid = make_1d_slab_grid(16);
    const auto flux = uniform_flux(grid, 1.0);
    const auto phi = cell_field(grid, sample_value(grid, [](double x) {
        return x < 0.5 ? 2.0 : 0.0;   // step
    }));
    const auto grad = compute_gradient_gauss(phi, grid.mesh, grid.geometry);

    const std::vector<LimiterType> limiters = {
        LimiterType::MINMOD, LimiterType::VANLEER, LimiterType::SUPERBEE,
        LimiterType::VAN_ALBADA, LimiterType::MC};
    for (const LimiterType limiter : limiters) {
        const auto face = interpolate_cell_to_face(
            phi, grid.mesh, grid.geometry, InterpScheme::LIMITED, &flux, limiter, &grad);
        const auto& own = grid.mesh.ownership();
        for (const std::size_t i : grid.interior_faces) {
            const std::size_t f = grid.cross[i];
            const std::size_t o = own.owner(f);
            const std::size_t nb = static_cast<std::size_t>(own.neighbour(f));
            const double lo = std::min(phi(o), phi(nb));
            const double hi = std::max(phi(o), phi(nb));
            const double v = face(f);
            require(v >= lo - 1e-12 && v <= hi + 1e-12,
                    std::string(to_string(limiter)) +
                    ": TVD face value must stay in the adjacent-cell envelope");
        }
        std::cout << "CONV_BOUNDED limiter=" << to_string(limiter) << " ok\n";
    }
}

void check_linear_exactness()
{
    const Grid grid = make_1d_slab_grid(16);
    const auto flux = uniform_flux(grid, 1.0);
    const auto phi = cell_field(grid, sample_value(grid, [](double x) { return 2.0 * x + 1.0; }));
    const auto grad = compute_gradient_gauss(phi, grid.mesh, grid.geometry);

    const std::vector<LimiterType> limiters = {
        LimiterType::MINMOD, LimiterType::VANLEER, LimiterType::SUPERBEE,
        LimiterType::VAN_ALBADA, LimiterType::MC};
    for (const LimiterType limiter : limiters) {
        const auto face = interpolate_cell_to_face(
            phi, grid.mesh, grid.geometry, InterpScheme::LIMITED, &flux, limiter, &grad);
        const auto e = profile_face_error(grid, face,
            [](double x) { return 2.0 * x + 1.0; });
        require(e.linf <= 1e-9,
                std::string(to_string(limiter)) + ": TVD must be linear-exact on deep interior faces");
        std::cout << "CONV_LINEAR limiter=" << to_string(limiter)
                  << " Linf=" << e.linf << " L2=" << e.l2 << "\n";
    }
}

void check_order()
{
    const std::vector<LimiterType> limiters = {
        LimiterType::MINMOD, LimiterType::VANLEER, LimiterType::MC};
    const std::vector<std::size_t> ns = {8u, 16u, 32u};

    // First-order upwind.
    {
        std::vector<double> errors;
        for (const std::size_t n : ns) {
            const Grid grid = make_1d_slab_grid(n);
            const auto flux = uniform_flux(grid, 1.0);
            const auto phi = cell_field(grid, sample_value(grid, [](double x) { return 1.0 + x * x; }));
            const auto face = interpolate_cell_to_face(phi, grid.mesh, grid.geometry, InterpScheme::UPWIND, &flux);
            errors.push_back(profile_face_error(grid, face,
                [](double x) { return 1.0 + x * x; }).l2);
        }
        std::cout << "CONV_ORDER scheme=upwind";
        for (std::size_t k = 0; k < errors.size(); ++k) {
            std::cout << " n=" << ns[k] << " L2=" << errors[k];
            if (k) std::cout << " order=" << observed_order(errors[k - 1], errors[k]);
        }
        std::cout << "\n";
        require_order(errors, 1.0, 0.8, "upwind on smooth monotone profile");
    }

    for (const LimiterType limiter : limiters) {
        std::vector<double> errors;
        for (const std::size_t n : ns) {
            const Grid grid = make_1d_slab_grid(n);
            const auto flux = uniform_flux(grid, 1.0);
            const auto phi = cell_field(grid, sample_value(grid, [](double x) { return 1.0 + x * x; }));
            const auto grad = compute_gradient_gauss(phi, grid.mesh, grid.geometry);
            const auto face = interpolate_cell_to_face(
                phi, grid.mesh, grid.geometry, InterpScheme::LIMITED, &flux, limiter, &grad);
            errors.push_back(profile_face_error(grid, face,
                [](double x) { return 1.0 + x * x; }).l2);
        }
        std::cout << "CONV_ORDER scheme=tvd_" << to_string(limiter);
        for (std::size_t k = 0; k < errors.size(); ++k) {
            std::cout << " n=" << ns[k] << " L2=" << errors[k];
            if (k) std::cout << " order=" << observed_order(errors[k - 1], errors[k]);
        }
        std::cout << "\n";
        require_order(errors, 1.0, 0.8, std::string("tvd_") + to_string(limiter));
    }
}

void check_conservation()
{
    const Grid grid = make_1d_slab_grid(16);
    const auto flux = uniform_flux(grid, 1.0);
    const auto phi = cell_field(grid, sample_value(grid, [](double x) { return 1.0 + x * x; }));

    const std::vector<std::pair<std::string, InterpScheme>> schemes = {
        {"upwind", InterpScheme::UPWIND},
        {"linear", InterpScheme::LINEAR},
        {"limited", InterpScheme::LIMITED}};
    for (const auto& [name, scheme] : schemes) {
        const auto conv = compute_convection(
            phi, flux, grid.mesh, scheme,
            scheme == InterpScheme::LIMITED ? LimiterType::VANLEER : LimiterType::NONE);
        double sum = 0.0;
        for (std::size_t c = 0; c < grid.mesh.n_cells(); ++c) sum += conv(c);
        // Net boundary flux: F*phi(left) - F*phi(right), boundaries take the
        // owner value under every scheme.
        const double net = phi(0) - phi(grid.mesh.n_cells() - 1);
        require(std::abs(sum - net) <= 1e-9,
                name + ": convection must be discretely conservative");
        std::cout << "CONV_CONSERVATION scheme=" << name
                  << " sum=" << sum << " net=" << net << "\n";
    }
}

} // namespace

int main()
{
    try {
        std::cout << std::setprecision(12);
        {
            const Grid g = make_1d_slab_grid(16);
            const auto& own = g.mesh.ownership();
            double mind = 1e30;
            for (std::size_t f = 0; f < g.mesh.n_faces(); ++f) {
                if (own.neighbour(f) < 0) continue;
                const std::size_t o = own.owner(f);
                const std::size_t nb = static_cast<std::size_t>(own.neighbour(f));
                const double d = (g.geometry.cell_centres[nb] - g.geometry.cell_centres[o]).mag();
                if (d < mind) mind = d;
                if (d < 1e-9) std::cout << "DEGENERATE face=" << f
                    << " owner=" << o << " nb=" << nb << " d=" << d
                    << " oC=" << g.geometry.cell_centres[o].x << " nbC=" << g.geometry.cell_centres[nb].x << "\n";
            }
            std::cout << "min internal centre distance=" << mind << "\n";
        }
        std::cout << std::setprecision(12);
        check_boundedness();
        check_linear_exactness();
        check_order();
        check_conservation();
        std::cout << "CONVECTION_SCHEME_VERIFICATION: PASS\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "CONVECTION_SCHEME_VERIFICATION: FAIL: " << e.what() << "\n";
        return 1;
    }
}
#pragma once

// Geometry-only diagnostic of the *joint* omega, wall nut and production G
// reference at cells incident to multiple wall faces.  It intentionally does
// NOT replace SST wall boundary conditions: OpenFOAM uses patch-specific
// wall-function weights/constraints and the near-wall production is not the
// cell-centred G proxy supplied here.  All three channels are audited together.

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <map>
#include <stdexcept>
#include <vector>

namespace cfdx::physics {

struct SstWallFaceTripletV13 {
    std::size_t cell = 0;
    std::size_t face = 0;
    double area = 0.0;
    double reference_omega = 0.0;
    double reference_nut = 0.0;
    double reference_G = 0.0;
    double owner_omega = 0.0;
    double owner_nut = 0.0;
    double owner_G = 0.0;
};

struct SstWallCellTripletV13 {
    std::size_t cell = 0;
    std::size_t face_count = 0;
    double total_area = 0.0;
    double reference_omega_area_weighted = 0.0;
    double reference_nut_area_weighted = 0.0;
    double reference_G_area_weighted = 0.0;
    double omega_owner_gap = 0.0;
    double nut_owner_proxy_gap = 0.0;
    double G_owner_proxy_gap = 0.0;
};

inline std::vector<SstWallCellTripletV13> audit_sst_wall_triplet_cells_v13(
    const std::vector<SstWallFaceTripletV13>& faces)
{
    if (faces.empty()) throw std::invalid_argument("V13 SST wall triplet: no wall faces");
    struct State {
        SstWallCellTripletV13 out;
        long double omega_area = 0.0L;
        long double nut_area = 0.0L;
        long double G_area = 0.0L;
        double owner_omega = 0.0;
        double owner_nut = 0.0;
        double owner_G = 0.0;
    };
    std::map<std::size_t, State> by_cell;
    std::map<std::size_t, bool> seen_face;
    for (const auto& f: faces) {
        const double values[] = {f.area, f.reference_omega, f.reference_nut,
            f.reference_G, f.owner_omega, f.owner_nut, f.owner_G};
        for (double value:values)
            if (!std::isfinite(value) || value < 0.0)
                throw std::invalid_argument("V13 SST wall triplet: negative/nonfinite measurement");
        if (!(f.area > 0.0) || !seen_face.emplace(f.face, true).second)
            throw std::invalid_argument("V13 SST wall triplet: invalid/duplicate wall face");
        auto& state = by_cell[f.cell];
        if (state.out.face_count) {
            if (state.owner_omega != f.owner_omega ||
                state.owner_nut != f.owner_nut || state.owner_G != f.owner_G)
                throw std::invalid_argument("V13 SST wall triplet: inconsistent owner state");
        } else {
            state.out.cell = f.cell;
            state.owner_omega = f.owner_omega;
            state.owner_nut = f.owner_nut;
            state.owner_G = f.owner_G;
        }
        ++state.out.face_count;
        state.out.total_area += f.area;
        state.omega_area += static_cast<long double>(f.area)*f.reference_omega;
        state.nut_area += static_cast<long double>(f.area)*f.reference_nut;
        state.G_area += static_cast<long double>(f.area)*f.reference_G;
        if (!std::isfinite(state.out.total_area))
            throw std::overflow_error("V13 SST wall triplet: area overflow");
    }
    std::vector<SstWallCellTripletV13> results;
    results.reserve(by_cell.size());
    for (const auto& [cell, state]:by_cell) {
        auto out=state.out;
        const auto normalized=[&](long double sum) {
            const auto value=static_cast<double>(sum/out.total_area);
            if (!std::isfinite(value))
                throw std::overflow_error("V13 SST wall triplet: weighted reference overflow");
            return value;
        };
        out.reference_omega_area_weighted=normalized(state.omega_area);
        out.reference_nut_area_weighted=normalized(state.nut_area);
        out.reference_G_area_weighted=normalized(state.G_area);
        const auto gap=[](double reference, double owner) {
            const double norm=std::max({std::numeric_limits<double>::min(),reference,owner});
            return std::abs(reference-owner)/norm;
        };
        out.omega_owner_gap=gap(out.reference_omega_area_weighted,state.owner_omega);
        out.nut_owner_proxy_gap=gap(out.reference_nut_area_weighted,state.owner_nut);
        out.G_owner_proxy_gap=gap(out.reference_G_area_weighted,state.owner_G);
        results.push_back(out);
    }
    return results;
}

} // namespace cfdx::physics

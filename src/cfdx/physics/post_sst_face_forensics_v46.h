#pragma once
// V46: opt-in, read-only face-by-face audit of the POST-SST momentum operator.
// The authoritative Rhie-Chow flux, updated SST viscosity, cell state and
// boundary conditions must come from the SAME candidate. This reports evidence
// only: arithmetic closure is NOT a physical convergence criterion.
#include "cfdx/physics/momentum_face_budget.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <sstream>
#include <limits>
#include <string>
#include <vector>

namespace cfdx::physics {

struct PostSSTFaceTermV46 {
    std::size_t face = 0;
    std::size_t owner = 0;
    long neighbour = -1;
    int sign = 0;  // +1 for owner, -1 for neighbour
    std::string patch;
    double mass = 0.0;
    double convection = 0.0;
    double viscosity = 0.0;
    double transpose_stress = 0.0;
    // This difference is NOT pure Rhie-Chow: it includes pressure splitting.
    double collocated_convection_difference = 0.0;
};

struct PostSSTCellForensicsV46 {
    std::size_t cell = 0;
    MomentumFaceCellBudget balance;
    std::vector<PostSSTFaceTermV46> faces;
    bool face_sum_consistent = false;
    bool matrix_face_consistent = false;
};

struct PostSSTFaceForensicsV46 {
    std::vector<PostSSTCellForensicsV46> cells;
    bool all_cells_matrix_face_consistent = false;
    double global_matrix_face_linf = 0.0;
    double global_matrix_face_relative_linf = 0.0;
};

// Explicit opt-in parsing. The default 58..65 interval covers the first
// observed V44 phase orbit; callers may select other bounded windows.
struct PostSSTFaceSelectionV46 {
    std::size_t first = 58;
    std::size_t last = 65;
    std::vector<std::size_t> cells;
    bool contains(std::size_t iteration) const noexcept {
        return iteration >= first && iteration <= last;
    }
};
inline std::size_t post_sst_v46_unsigned(const std::string& token) {
    if (token.empty() || token.size()>20 ||
        !std::all_of(token.begin(),token.end(),[](char c){return c>='0' && c<='9';}))
        throw std::invalid_argument("V46 post-SST face trace: invalid unsigned integer");
    std::size_t used=0;
    const auto value=std::stoull(token,&used,10);
    if (used!=token.size() || value>std::numeric_limits<std::size_t>::max())
        throw std::invalid_argument("V46 post-SST face trace: integer overflow");
    return static_cast<std::size_t>(value);
}
inline PostSSTFaceSelectionV46 post_sst_face_selection_v46(
    const char* window, const char* cells, std::size_t nc) {
    if (nc==0) throw std::invalid_argument("V46 face trace: empty mesh");
    PostSSTFaceSelectionV46 result;
    if (window && *window) {
        const std::string text(window);
        const auto colon=text.find(':');
        if (colon==std::string::npos || text.find(':',colon+1)!=std::string::npos)
            throw std::invalid_argument("V46 face trace: window must be FIRST:LAST");
        result.first=post_sst_v46_unsigned(text.substr(0,colon));
        result.last=post_sst_v46_unsigned(text.substr(colon+1));
        if (!result.first || result.last<result.first ||
            result.last-result.first>255)
            throw std::invalid_argument("V46 face trace: window must cover <=256 iterations");
    }
    if (cells && *cells) {
        std::istringstream in(cells);
        std::string token;
        while (std::getline(in,token,',')) {
            const auto c=post_sst_v46_unsigned(token);
            if (c>=nc || std::find(result.cells.begin(),result.cells.end(),c)!=result.cells.end())
                throw std::invalid_argument("V46 face trace: duplicate/out-of-range cell");
            result.cells.push_back(c);
            if (result.cells.size()>32)
                throw std::invalid_argument("V46 face trace: at most 32 cells");
        }
        if (!cells[0] || cells[std::char_traits<char>::length(cells)-1]==',')
            throw std::invalid_argument("V46 face trace: empty cell token");
    } else {
        for (const auto c : {std::size_t(1),std::size_t(2),std::size_t(12),std::size_t(30)})
            if (c<nc) result.cells.push_back(c);
        if (result.cells.empty()) result.cells.push_back(0);
    }
    return result;
}

// Deliberately select a bounded set of cells. Never write an O(nFaces)
// diagnostic log by default. `physical_body` EXCLUDES the transpose-stress
// divergence, which is passed as an oriented face flux and integrated once.
inline PostSSTFaceForensicsV46 diagnose_post_sst_faces_v46(
    const cfdx::core::Mesh& mesh, const FvGeometry& geometry,
    const cfdx::core::Field<double,cfdx::core::Location::FACE>& mass_flux,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& component,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& pressure_gradient,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& physical_body,
    const std::vector<double>& updated_dynamic_viscosity,
    const ScalarBoundaryConditions& boundary_conditions,
    ConvectionScheme scheme, double high_order_fraction, bool bounded_convection,
    const ScalarEquation& post_sst_equation,
    const std::vector<std::size_t>& requested_cells,
    const cfdx::core::Field<double,cfdx::core::Location::FACE>* transpose_stress_flux=nullptr,
    const std::vector<double>* wall_face_viscosity=nullptr,
    const ScalarBoundaryFaceConditions* face_conditions=nullptr,
    const cfdx::core::Field<double,cfdx::core::Location::FACE>* collocated_flux=nullptr)
{
    using cfdx::core::Location;
    const auto nc=mesh.n_cells(), nf=mesh.n_faces();
    if (nc == 0 || post_sst_equation.matrix.n_rows() != nc ||
        !post_sst_equation.matrix.is_consistent() ||
        requested_cells.size() > 32 ||
        (transpose_stress_flux && (transpose_stress_flux->size()!=nf ||
                                   transpose_stress_flux->dimension()!=1)))
        throw std::invalid_argument("V46 post-SST face trace: invalid operator/selection");
    for (const auto c : requested_cells)
        if (c>=nc)
            throw std::invalid_argument("V46 post-SST face trace: cell out of range");

    PostSSTFaceForensicsV46 result;
    if (requested_cells.empty()) return result;
    auto budget=diagnose_momentum_face_budget(
        mesh,geometry,mass_flux,component,pressure_gradient,physical_body,
        updated_dynamic_viscosity,boundary_conditions,scheme,
        high_order_fraction,bounded_convection,collocated_flux,
        transpose_stress_flux,wall_face_viscosity,nullptr,face_conditions);
    reconcile_momentum_budget_with_matrix(budget,post_sst_equation,component);
    result.all_cells_matrix_face_consistent=momentum_matrix_face_closure_consistent(budget);
    result.global_matrix_face_linf=budget.matrix_closure_linf;
    result.global_matrix_face_relative_linf=budget.matrix_closure_relative_linf;

    // Use the FV flux reconstructor (not CSR multiplication) to attribute
    // convection and non-orthogonal viscous contributions to each face.
    const auto conv=reconstruct_scalar_transport_flux(
        mesh,geometry,mass_flux,component,0.0,boundary_conditions,
        bounded_convection,scheme,face_conditions,nullptr,nullptr,nullptr,
        high_order_fraction);
    const auto full=reconstruct_scalar_transport_flux(
        mesh,geometry,mass_flux,component,0.0,boundary_conditions,
        bounded_convection,scheme,face_conditions,&updated_dynamic_viscosity,
        nullptr,wall_face_viscosity,high_order_fraction);
    cfdx::core::Field<double,Location::FACE> colconv;
    if (collocated_flux)
        colconv=reconstruct_scalar_transport_flux(
            mesh,geometry,*collocated_flux,component,0.0,boundary_conditions,
            bounded_convection,scheme,face_conditions,nullptr,nullptr,nullptr,
            high_order_fraction);

    const auto* offsets=mesh.cells().offsets_data();
    const auto* face_ids=mesh.cells().faces_data();
    for (const auto cell : requested_cells) {
        if (std::any_of(result.cells.begin(),result.cells.end(),
            [cell](const auto& prior){return prior.cell==cell;}))
            throw std::invalid_argument("V46 post-SST face trace: repeated cell");
        PostSSTCellForensicsV46 selected;
        selected.cell=cell;
        selected.balance=budget.cells.at(cell);
        double sum_conv=0.,sum_visc=0.,sum_trans=0.,sum_mass=0.,sum_delta=0.;
        for (std::size_t j=offsets[cell]; j<offsets[cell+1]; ++j) {
            const auto f=static_cast<std::size_t>(face_ids[j]);
            if (f>=nf)
                throw std::invalid_argument("V46 post-SST face trace: bad cell face");
            const auto owner=mesh.ownership().owner(f);
            const auto neighbour=mesh.ownership().neighbour(f);
            const int sign=(owner==cell) ? 1 :
                (neighbour>=0 && static_cast<std::size_t>(neighbour)==cell) ? -1 : 0;
            if (sign==0 || owner>=nc ||
                (neighbour>=0 && static_cast<std::size_t>(neighbour)>=nc) ||
                std::any_of(selected.faces.begin(),selected.faces.end(),
                    [f](const auto& previous){return previous.face==f;}))
                throw std::invalid_argument("V46 post-SST face trace: inconsistent connectivity");
            PostSSTFaceTermV46 t;
            t.face=f; t.owner=owner; t.neighbour=neighbour; t.sign=sign;
            if (neighbour < 0) {
                const auto patch=geometry.face_patch.at(f);
                if (patch >= mesh.boundary().n_patches())
                    throw std::invalid_argument("V46 post-SST face trace: boundary patch missing");
                t.patch=mesh.boundary().patch(patch).name;
            } else t.patch="internal";
            t.mass=sign*mass_flux(f);
            t.convection=sign*conv(f);
            t.viscosity=sign*(full(f)-conv(f));
            t.transpose_stress=transpose_stress_flux
                ? -sign*(*transpose_stress_flux)(f) : 0.;
            t.collocated_convection_difference=collocated_flux
                ? sign*(conv(f)-colconv(f)) : 0.;
            for (const auto x : {t.mass,t.convection,t.viscosity,t.transpose_stress,
                                 t.collocated_convection_difference})
                if (!std::isfinite(x))
                    throw std::invalid_argument("V46 post-SST face trace: nonfinite flux");
            sum_mass+=t.mass;sum_conv+=t.convection;sum_visc+=t.viscosity;
            sum_trans+=t.transpose_stress;
            sum_delta+=t.collocated_convection_difference;
            selected.faces.push_back(std::move(t));
        }
        const auto& b=selected.balance;
        const double scale=std::max({1e-15,std::abs(b.mass_divergence),
            std::abs(b.convection),std::abs(b.viscous),
            std::abs(b.transpose_viscous),std::abs(b.collocated_difference)});
        selected.face_sum_consistent=
            std::abs(sum_mass-b.mass_divergence)<=2e-10*scale &&
            std::abs(sum_conv-b.convection)<=2e-10*scale &&
            std::abs(sum_visc-b.viscous)<=2e-10*scale &&
            std::abs(sum_trans-b.transpose_viscous)<=2e-10*scale &&
            std::abs(sum_delta-b.collocated_difference)<=2e-10*scale;
        selected.matrix_face_consistent=std::isfinite(b.matrix_face_difference) &&
            std::isfinite(b.matrix_face_scale) &&
            std::abs(b.matrix_face_difference)<=2e-10*std::max(1e-15,b.matrix_face_scale);
        result.cells.push_back(std::move(selected));
    }
    return result;
}
} // namespace cfdx::physics

// ORIGINAL instrumentation for OpenFOAM Foundation v11 (GPL-compatible).
// Hook in the native solver immediately after omegaEqn.boundaryManipulate(...)
// and BEFORE omegaEqn.solve(). This is NOT built by CFDX; compile it with
// the native OpenFOAM v11 toolchain, never substitute CFDX-computed fields.
#pragma once
#include "fvCFD.H"
#include "omegaWallFunctionFvPatchScalarField.H"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

inline void cfdxOmegaWallSnapshotV31(
    const Foam::fvMesh& mesh,
    const Foam::volScalarField& omega,
    const Foam::volScalarField& G,
    const Foam::volScalarField& nut,
    const Foam::fvScalarMatrix& omegaEqn,
    const std::string& filename)
{
    const auto n=static_cast<std::size_t>(mesh.nCells());
    if (!n || omega.size()!=mesh.nCells() || G.size()!=mesh.nCells() ||
        nut.size()!=mesh.nCells())
        throw std::runtime_error("V31 OpenFOAM snapshot: field sizes");
    std::vector<int> corner(n,0);
    std::set<Foam::label> cell_ids,face_ids;
    forAll(omega.boundaryField(), patchi) {
        const auto& patch_field=omega.boundaryField()[patchi];
        if (!Foam::isA<Foam::omegaWallFunctionFvPatchScalarField>(patch_field))
            continue;
        const auto& patch=mesh.boundary()[patchi];
        forAll(patch.faceCells(), f) {
            const auto owner=patch.faceCells()[f];
            if (owner<0 || static_cast<std::size_t>(owner)>=n)
                throw std::runtime_error("V31 OpenFOAM snapshot: invalid owner");
            ++corner[owner];
            cell_ids.insert(owner);
            face_ids.insert(patch.start()+f);
        }
    }
    if (cell_ids.empty() || face_ids.empty())
        throw std::runtime_error("V31 OpenFOAM snapshot: no omegaWallFunction face");
    // Assemble actual algebraic LDU rows; this must be called after boundary
    // matrix manipulations. Do not infer coefficients from wall formulas.
    std::vector<std::map<Foam::label,double>> rows(n);
    const auto& diag=omegaEqn.diag();
    const auto& rhs=omegaEqn.source();
    const auto& a=omegaEqn.lduAddr();
    const auto& lowerAddr=a.lowerAddr();
    const auto& upperAddr=a.upperAddr();
    const auto& lower=omegaEqn.lower();
    const auto& upper=omegaEqn.upper();
    for (std::size_t i=0;i<n;++i) rows[i][static_cast<Foam::label>(i)]+=diag[i];
    forAll(lowerAddr, i) {
        const auto lo=lowerAddr[i],hi=upperAddr[i];
        if (lo<0 || hi<0 || static_cast<std::size_t>(hi)>=n)
            throw std::runtime_error("V31 OpenFOAM snapshot: invalid LDU address");
        rows[lo][hi]+=upper[i];
        rows[hi][lo]+=lower[i];
    }
    std::ofstream out(filename);
    if (!out.good()) throw std::runtime_error("V31 OpenFOAM snapshot: cannot open file");
    out<<std::setprecision(17);
    out<<"{\"schema\":\"CFDX_OF11_NATIVE_SNAPSHOT_V31\",\"stage\":\"post_manipulate_pre_solve\",\"cells\":[";
    bool first=true;
    for (const auto cell:cell_ids) {
        if (!first) out<<','; first=false;
        out<<"{\"cell\":"<<cell<<",\"omega\":"<<omega[cell]
           <<",\"G\":"<<G[cell]<<",\"nut\":"<<nut[cell]
           <<",\"corner_faces\":"<<corner[cell]<<'}';
    }
    out<<"],\"faces\":[";
    first=true;
    forAll(omega.boundaryField(), patchi) {
        const auto& pf=omega.boundaryField()[patchi];
        if (!Foam::isA<Foam::omegaWallFunctionFvPatchScalarField>(pf)) continue;
        const auto& patch=mesh.boundary()[patchi];
        forAll(patch.faceCells(), f) {
            const auto owner=patch.faceCells()[f];
            const double coverage=patch.magSf()[f]/patch.patch().magFaceAreas()[f];
            const double weight=coverage<=.1?0.:(coverage-.1)/.9;
            if (!std::isfinite(weight) || weight<0. || weight>1.+1e-12)
                throw std::runtime_error("V31 OpenFOAM snapshot: invalid patch weight");
            if (!first) out<<',';first=false;
            out<<"{\"face\":"<<(patch.start()+f)<<",\"owner\":"<<owner
               <<",\"coverage\":"<<coverage<<",\"weight\":"<<weight
               <<",\"omega\":"<<omega[owner]<<",\"G\":"<<G[owner]
               <<",\"nut\":"<<nut[owner]<<'}';
        }
    }
    out<<"],\"matrix_rows\":[";
    first=true;
    // Include ALL rows affected by a constrained cell (its own row plus any
    // neighbour row carrying its off-diagonal). Matrix scope is declared in
    // the export and compared strictly, never silently truncated.
    std::set<Foam::label> affected=cell_ids;
    for (std::size_t i=0;i<n;++i)
        for (const auto& entry:rows[i])
            if (cell_ids.count(entry.first)) affected.insert(i);
    for (const auto row:affected) {
        if (!first) out<<',';first=false;
        out<<"{\"row\":"<<row<<",\"diagonal\":"<<diag[row]
           <<",\"rhs\":"<<rhs[row]<<",\"entries\":[";
        bool first_entry=true;
        for (const auto& [col,coefficient]:rows[row]) {
            if (!first_entry) out<<',';first_entry=false;
            out<<'['<<col<<','<<coefficient<<']';
        }
        out<<"]}";
    }
    out<<"]}\n";
    if (!out.good()) throw std::runtime_error("V31 OpenFOAM snapshot: write failed");
}

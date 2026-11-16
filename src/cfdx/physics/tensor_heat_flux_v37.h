#pragma once
// CFDX V37 experimental tensor heat-flux operator: face-based, SPD-controlled
// and conservative by construction. This is a face primitive for future
// production tensor-energy assembly, NOT a completed tensor energy solver.
#include "cfdx/physics/finite_volume_transport.h"
#include <array>
#include <cstdint>
#include <limits>
#include <numeric>

namespace cfdx::physics::v37 {
using cfdx::core::Vec3;
enum class ThermalBackend { Serial, CpuOpenMP, GpuRequested };
struct TensorK {
    double xx=1.,yy=1.,zz=1.,xy=0.,xz=0.,yz=0.;
};
inline Vec3 apply(TensorK k,Vec3 g) {
    return {k.xx*g.x+k.xy*g.y+k.xz*g.z,
            k.xy*g.x+k.yy*g.y+k.yz*g.z,
            k.xz*g.x+k.yz*g.y+k.zz*g.z};
}
inline bool tensor_is_spd(TensorK a,double minimum_relative_pivot=1e-14) {
    const double s=std::max({std::abs(a.xx),std::abs(a.yy),std::abs(a.zz),
        std::abs(a.xy),std::abs(a.xz),std::abs(a.yz)});
    if(!(s>0.)||!std::isfinite(s)||!(minimum_relative_pivot>=0.)||
       !std::isfinite(minimum_relative_pivot)) return false;
    const long double x=a.xx/s,y=a.yy/s,z=a.zz/s,
        b=a.xy/s,c=a.xz/s,d=a.yz/s;
    if(!(x>minimum_relative_pivot))return false;
    const long double l11=std::sqrt(x),l21=b/l11,l31=c/l11;
    const long double d22=y-l21*l21;
    if(!(d22>minimum_relative_pivot))return false;
    const long double l22=std::sqrt(d22),l32=(d-l31*l21)/l22;
    const long double d33=z-l31*l31-l32*l32;
    return d33>minimum_relative_pivot;
}
struct TensorFaceData {
    Vec3 centre_to_centre,area_owner_to_neighbour,gradient_at_face;
    double temperature_owner=0.,temperature_neighbour=0.;
};
struct TensorFaceFlux {
    double owner_outward_watts=0.,neighbour_outward_watts=0.;
    double implicit_conductance=0., deferred_correction=0.;
};
inline TensorFaceFlux tensor_diffusive_face(TensorFaceData f,TensorK k) {
    if(!tensor_is_spd(k))throw std::invalid_argument("thermal tensor must be SPD");
    const auto d=f.centre_to_centre,S=f.area_owner_to_neighbour,g=f.gradient_at_face;
    if(!std::isfinite(d.x)||!std::isfinite(d.y)||!std::isfinite(d.z)||
       !std::isfinite(S.x)||!std::isfinite(S.y)||!std::isfinite(S.z)||
       !std::isfinite(g.x)||!std::isfinite(g.y)||!std::isfinite(g.z)||
       !std::isfinite(f.temperature_owner)||!std::isfinite(f.temperature_neighbour))
        throw std::invalid_argument("thermal tensor face requires finite geometry and state");
    const double projection=S.dot(d);
    if(!(projection>0.)||!(S.mag()>0.)||!(d.mag()>0.))
        throw std::invalid_argument("thermal tensor face invalid orientation");
    const double conductance=S.dot(apply(k,S))/projection;
    const double ortho=-conductance*(f.temperature_neighbour-f.temperature_owner);
    const double full=-S.dot(apply(k,g));
    const double correction=full-ortho;
    if(!(conductance>0.)||!std::isfinite(conductance)||!std::isfinite(full)||
       !std::isfinite(correction))throw std::runtime_error("thermal tensor face overflow");
    return {full,-full,conductance,correction};
}
inline std::vector<TensorFaceFlux> evaluate_tensor_faces(
    const std::vector<TensorFaceData>& faces,TensorK tensor,ThermalBackend backend)
{
    if(backend==ThermalBackend::GpuRequested)
        throw std::runtime_error("tensor GPU backend not implemented; CPU fallback forbidden");
#ifndef _OPENMP
    if(backend==ThermalBackend::CpuOpenMP)
        throw std::runtime_error("tensor OpenMP backend not compiled; CPU fallback forbidden");
#endif
    if(!tensor_is_spd(tensor))throw std::invalid_argument("thermal tensor must be SPD");
    // Preflight outside parallel regions. Exceptions in OpenMP workers can
    // terminate the process, so the kernel must not throw inside the team.
    for(const auto& f:faces) {
        const auto d=f.centre_to_centre,S=f.area_owner_to_neighbour,g=f.gradient_at_face;
        if(!std::isfinite(d.x)||!std::isfinite(d.y)||!std::isfinite(d.z)||
           !std::isfinite(S.x)||!std::isfinite(S.y)||!std::isfinite(S.z)||
           !std::isfinite(g.x)||!std::isfinite(g.y)||!std::isfinite(g.z)||
           !std::isfinite(f.temperature_owner)||!std::isfinite(f.temperature_neighbour)||
           !(S.dot(d)>0.0))
            throw std::invalid_argument("thermal tensor invalid face preflight");
    }
    std::vector<TensorFaceFlux> result(faces.size());
    const auto unchecked=[&](std::size_t i) {
        const auto& f=faces[i];
        const auto& S=f.area_owner_to_neighbour;
        const double D=S.dot(apply(tensor,S))/S.dot(f.centre_to_centre);
        const double orth=-D*(f.temperature_neighbour-f.temperature_owner);
        const double full=-S.dot(apply(tensor,f.gradient_at_face));
        result[i]={full,-full,D,full-orth};
    };
#ifdef _OPENMP
    if(backend==ThermalBackend::CpuOpenMP) {
#pragma omp parallel for schedule(static)
        for(std::int64_t i=0;i<static_cast<std::int64_t>(faces.size());++i)
            unchecked(static_cast<std::size_t>(i));
    } else
#endif
        for(std::size_t i=0;i<faces.size();++i)unchecked(i);
    for(const auto& flux:result)
        if(!(flux.implicit_conductance>0.)||
           !std::isfinite(flux.implicit_conductance)||
           !std::isfinite(flux.owner_outward_watts)||
           !std::isfinite(flux.deferred_correction))
            throw std::runtime_error("thermal tensor kernel produced invalid flux");
    return result;
}
// Deterministic cell-owned CSR gathering of the ONE oriented flux per face.
// Each cell accumulates faces in monotonically increasing face ID order.
// OpenMP and serial execute each cell's reduction in the SAME arithmetic order.
struct FaceGatherGraph {
    std::vector<std::size_t> offsets,faces;
    std::vector<std::int8_t> signs;
};
inline FaceGatherGraph make_face_gather(const cfdx::core::Mesh& mesh) {
    const std::size_t nc=mesh.n_cells(),nf=mesh.n_faces();
    FaceGatherGraph g;g.offsets.assign(nc+1,0);
    for(std::size_t f=0;f<nf;++f) {
        const auto o=mesh.ownership().owner(f);
        if(o>=nc)throw std::invalid_argument("tensor gather invalid owner");
        ++g.offsets[o+1];
        const auto n=mesh.ownership().neighbour(f);
        if(n>=0) {
            if(static_cast<std::size_t>(n)>=nc)
                throw std::invalid_argument("tensor gather invalid neighbour");
            ++g.offsets[static_cast<std::size_t>(n)+1];
        }
    }
    std::partial_sum(g.offsets.begin(),g.offsets.end(),g.offsets.begin());
    g.faces.resize(g.offsets.back());g.signs.resize(g.offsets.back());
    auto cursor=g.offsets;
    for(std::size_t f=0;f<nf;++f) {
        const auto o=mesh.ownership().owner(f);
        const auto n=mesh.ownership().neighbour(f);
        auto j=cursor[o]++;g.faces[j]=f;g.signs[j]=+1;
        if(n>=0){j=cursor[static_cast<std::size_t>(n)]++;
            g.faces[j]=f;g.signs[j]=-1;}
    }
    return g;
}
inline std::vector<double> gather_face_watts(
    const FaceGatherGraph& graph,const std::vector<double>& oriented_face_watts,
    ThermalBackend backend)
{
    if(backend==ThermalBackend::GpuRequested)
        throw std::runtime_error("thermal GPU reduction not implemented; CPU fallback forbidden");
#ifndef _OPENMP
    if(backend==ThermalBackend::CpuOpenMP)
        throw std::runtime_error("thermal OpenMP reduction not compiled; CPU fallback forbidden");
#endif
    if(graph.offsets.empty()||graph.faces.size()!=graph.signs.size()||
       graph.offsets.back()!=graph.faces.size())
        throw std::invalid_argument("thermal face CSR dimensions invalid");
    for(std::size_t c=0;c+1<graph.offsets.size();++c)
        if(graph.offsets[c]>graph.offsets[c+1])
            throw std::invalid_argument("thermal face CSR nonmonotone offsets");
    for(const auto f:graph.faces)
        if(f>=oriented_face_watts.size()||!std::isfinite(oriented_face_watts[f]))
            throw std::invalid_argument("thermal face CSR invalid watts");
    for(const auto sign:graph.signs)
        if(sign!=1&&sign!=-1)
            throw std::invalid_argument("thermal face CSR invalid sign");
    std::vector<double> out(graph.offsets.size()-1,0.);
    const auto run=[&](std::size_t cell){
        double watts=0.0;
        for(std::size_t j=graph.offsets[cell];j<graph.offsets[cell+1];++j)
            watts+=graph.signs[j]*oriented_face_watts[graph.faces[j]];
        out[cell]=watts;
    };
#ifdef _OPENMP
    if(backend==ThermalBackend::CpuOpenMP) {
#pragma omp parallel for schedule(static)
        for(std::int64_t i=0;i<static_cast<std::int64_t>(out.size());++i)
            run(static_cast<std::size_t>(i));
        return out;
    }
#endif
    for(std::size_t i=0;i<out.size();++i)run(i);
    return out;
}
} // namespace cfdx::physics::v37

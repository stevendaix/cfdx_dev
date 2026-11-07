// CFDX V5: production-captured coupled preconditioner comparison.
// Only the evidence callback copies A, b, x0, and pressure operators.
// Linear phase replays the same captured system; physical phase solves
// the same initial-boundary-value problem with an unchanged nonlinear gate.
#include "cfdx/physics/steady_incompressible_solver.h"
#include "cfdx/core/linalg/coupled_amg_schur.h"
#include "cfdx/core/linalg/cell_overlapping_schwarz.h"
#include "cfdx/core/linalg/lsc_bfbt_schur.h"
#include "cfdx/core/linalg/pcd_schur.h"
#include "cfdx/core/linalg/linear_solver_dispatch.h"
#include "cfdx/core/linalg/gmres_solver.h"
#include "cfdx/core/linalg/coupled_block_extraction.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <exception>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <fstream>
#if defined(__GLIBC__)
#include <malloc.h>
#endif
#if defined(__linux__) || defined(__APPLE__)
#include <sys/resource.h>
#endif

using namespace cfdx::core;
using namespace cfdx::physics;
namespace {
using Clock = std::chrono::steady_clock;
using Tick = Clock::time_point;
double ms(Tick since) { return std::chrono::duration<double,std::milli>(Clock::now()-since).count(); }
Mesh make_cavity_mesh(std::size_t nx, std::size_t ny)
{
    if (nx < 4 || ny < 4) throw std::invalid_argument("cavity mesh is too small");
    Mesh mesh;
    mesh.points().resize((nx + 1) * (ny + 1) * 2);
    const auto point_id = [nx](std::size_t i, std::size_t j, std::size_t k) {
        return (j * (nx + 1) + i) * 2 + k;
    };
    for (std::size_t j=0;j<=ny;++j) for (std::size_t i=0;i<=nx;++i) {
        const double x=static_cast<double>(i)/static_cast<double>(nx);
        const double y=static_cast<double>(j)/static_cast<double>(ny);
        mesh.points().set(point_id(i,j,0),x,y,0.0);
        mesh.points().set(point_id(i,j,1),x,y,1.0);
    }

    std::map<std::vector<std::size_t>,std::size_t> face_map;
    std::vector<std::vector<std::size_t>> cell_faces(nx*ny);
    auto add_or_get_face=[&](std::initializer_list<std::size_t> vertices,std::size_t cell) {
        std::vector<std::size_t> sorted(vertices);
        std::sort(sorted.begin(),sorted.end());
        const auto found=face_map.find(sorted);
        if(found!=face_map.end()) {
            mesh.ownership().set_neighbour(found->second,static_cast<int>(cell));
            return found->second;
        }
        const std::size_t face=mesh.faces().n_faces();
        mesh.faces().push_face(vertices);
        mesh.ownership().resize(mesh.faces().n_faces());
        face_map.emplace(std::move(sorted),face);
        mesh.ownership().set_owner(face,cell);
        mesh.ownership().set_neighbour(face,FaceOwnership::BOUNDARY);
        return face;
    };

    for(std::size_t j=0;j<ny;++j) for(std::size_t i=0;i<nx;++i) {
        const std::size_t cell=j*nx+i;
        const std::size_t a=point_id(i,j,0), b=point_id(i+1,j,0);
        const std::size_t c=point_id(i+1,j+1,0), d=point_id(i,j+1,0);
        const std::size_t e=point_id(i,j,1), f=point_id(i+1,j,1);
        const std::size_t g=point_id(i+1,j+1,1), h=point_id(i,j+1,1);
        cell_faces[cell]={
            add_or_get_face({a,d,c,b},cell), add_or_get_face({e,f,g,h},cell),
            add_or_get_face({a,b,f,e},cell), add_or_get_face({d,h,g,c},cell),
            add_or_get_face({a,e,h,d},cell), add_or_get_face({b,c,g,f},cell)};
    }
    for(const auto& faces:cell_faces) mesh.cells().push_cell(faces);

    Patch bottom; bottom.name="bottom"; bottom.type=PatchType::WALL;
    Patch top; top.name="top"; top.type=PatchType::WALL;
    Patch left; left.name="left"; left.type=PatchType::WALL;
    Patch right; right.name="right"; right.type=PatchType::WALL;
    Patch front; front.name="front"; front.type=PatchType::WALL;
    Patch back; back.name="back"; back.type=PatchType::WALL;

    for(std::size_t f=0;f<mesh.n_faces();++f) {
        if(mesh.ownership().neighbour(f)>=0) continue;
        const auto& all_vertices=mesh.faces().vertices();
        const auto begin=all_vertices.begin()+static_cast<std::ptrdiff_t>(mesh.faces().face_offset(f));
        const auto end=begin+static_cast<std::ptrdiff_t>(mesh.faces().face_size(f));
        double x=0,y=0,z=0;
        for(auto it=begin;it!=end;++it) { x+=mesh.points().x(*it); y+=mesh.points().y(*it); z+=mesh.points().z(*it); }
        const double nverts=static_cast<double>(mesh.faces().face_size(f));
        x/=nverts; y/=nverts; z/=nverts;
        constexpr double tol=1e-12;
        if(std::abs(y)<tol) bottom.face_ids.push_back(f);
        else if(std::abs(y-1)<tol) top.face_ids.push_back(f);
        else if(std::abs(x)<tol) left.face_ids.push_back(f);
        else if(std::abs(x-1)<tol) right.face_ids.push_back(f);
        else if(std::abs(z)<tol) front.face_ids.push_back(f);
        else if(std::abs(z-1)<tol) back.face_ids.push_back(f);
        else throw std::runtime_error("cavity: unclassified boundary face");
    }
    mesh.boundary().add_patch(bottom); mesh.boundary().add_patch(top);
    mesh.boundary().add_patch(left); mesh.boundary().add_patch(right);
    mesh.boundary().add_patch(front); mesh.boundary().add_patch(back);
    return mesh;
}


struct Snapshot {
    SparseMatrix coupled, mass, kp, fp;
    Vector b, x0;
    std::size_t gauge = 0;
    bool fixed_pressure = false;
    bool valid = false;
    double max_old_velocity = 0.0;
};

// FNV1a on IEEE754 bytes: replay fidelity, not a cryptographic signature.
void hash_bytes(std::uint64_t& h, const void* ptr, std::size_t n) {
    const auto* p = static_cast<const unsigned char*>(ptr);
    for(std::size_t i=0;i<n;++i) { h ^= p[i]; h *= UINT64_C(1099511628211); }
}
template<typename T> void hash_scalar(std::uint64_t& h, T v){ hash_bytes(h,&v,sizeof(v)); }
void hash_matrix(std::uint64_t& h, const SparseMatrix& a){
    hash_scalar(h,a.n_rows()); hash_scalar(h,a.n_cols()); hash_scalar(h,a.nnz());
    for(std::size_t i=0;i<a.n_rows()+1;++i) hash_scalar(h,a.row_offsets_data()[i]);
    for(std::size_t i=0;i<a.nnz();++i){hash_scalar(h,a.columns_data()[i]);hash_scalar(h,a.values_data()[i]);}
}
void hash_vector(std::uint64_t& h, const Vector& v){
    hash_scalar(h,v.size()); for(std::size_t i=0;i<v.size();++i)hash_scalar(h,v(i));
}
std::string fingerprint(const Snapshot& s){
    std::uint64_t h=UINT64_C(14695981039346656037);
    hash_matrix(h,s.coupled); hash_vector(h,s.b); hash_vector(h,s.x0);
    hash_matrix(h,s.mass); hash_matrix(h,s.kp); hash_matrix(h,s.fp);
    hash_scalar(h,s.gauge); hash_scalar(h,s.fixed_pressure); hash_scalar(h,s.max_old_velocity);
    std::ostringstream os; os<<std::hex<<std::setw(16)<<std::setfill('0')<<h; return os.str();
}

// Versioned binary format for replay across isolated processes. The host
// representation and endianness are checked; the payload fingerprint covers
// matrix values, all pressure operators, RHS and initial iterate.
constexpr std::uint64_t kEndianTag=UINT64_C(0x0102030405060708);
template<class T> void put(std::ofstream& f,const T& x){
    f.write(reinterpret_cast<const char*>(&x),sizeof(x));
    if(!f)throw std::runtime_error("snapshot write failure");
}
template<class T> T get(std::ifstream& f){
    T x{};f.read(reinterpret_cast<char*>(&x),sizeof(x));
    if(!f)throw std::runtime_error("truncated capture");return x;
}
void put_matrix(std::ofstream& f,const SparseMatrix& a){
    put(f,static_cast<std::uint64_t>(a.n_rows()));
    put(f,static_cast<uint64_t>(a.n_cols()));
    put(f,static_cast<std::uint64_t>(a.nnz()));
    f.write(reinterpret_cast<const char*>(a.row_offsets_data()),
            static_cast<std::streamsize>((a.n_rows()+1)*sizeof(SparseMatrix::Index)));
    f.write(reinterpret_cast<const char*>(a.columns_data()),
            static_cast<std::streamsize>(a.nnz()*sizeof(SparseMatrix::Index)));
    f.write(reinterpret_cast<const char*>(a.values_data()),
            static_cast<std::streamsize>(a.nnz()*sizeof(double)));
    if(!f)throw std::runtime_error("snapshot matrix write failure");
}
SparseMatrix get_matrix(std::ifstream& f){
    const auto rows=get<std::uint64_t>(f),cols=get<std::uint64_t>(f),nnz=get<std::uint64_t>(f);
    if(rows>50000 || cols>50000 || nnz>8000000 ||
       nnz>std::numeric_limits<SparseMatrix::Index>::max() ||
       (rows!=0&&cols==0))throw std::runtime_error("snapshot dimensions out of bounds");
    std::vector<SparseMatrix::Index> offsets(static_cast<std::size_t>(rows+1));
    std::vector<SparseMatrix::Index> columns(static_cast<std::size_t>(nnz));
    std::vector<double> values(static_cast<std::size_t>(nnz));
    f.read(reinterpret_cast<char*>(offsets.data()),
           static_cast<std::streamsize>(offsets.size()*sizeof(SparseMatrix::Index)));
    f.read(reinterpret_cast<char*>(columns.data()),
           static_cast<std::streamsize>(columns.size()*sizeof(SparseMatrix::Index)));
    f.read(reinterpret_cast<char*>(values.data()),
           static_cast<std::streamsize>(values.size()*sizeof(double)));
    if(!f)throw std::runtime_error("truncated snapshot matrix");
    // Production CSR may contain repeated column entries. Their contributions
    // are summed by matvec but their *original storage* must be retained for
    // a byte-exact production replay. from_sorted_csr deliberately rejects
    // duplicates; COO finalize uses stable sorting and preserves them.
    if(offsets.empty()||offsets.front()!=0||offsets.back()!=nnz)
        throw std::runtime_error("invalid CSR offsets");
    SparseMatrix matrix(static_cast<std::size_t>(rows),static_cast<std::size_t>(cols));
    for(std::size_t row=0;row<rows;++row){
        if(offsets[row]>offsets[row+1]||offsets[row+1]>nnz)
            throw std::runtime_error("invalid CSR row offsets");
        for(std::size_t k=offsets[row];k<offsets[row+1];++k){
            if(columns[k]>=cols || (k>offsets[row]&&columns[k]<columns[k-1]))
                throw std::runtime_error("invalid CSR columns");
            matrix.push_back(row,columns[k],values[k]);
        }
    }
    matrix.finalize();
    if(!matrix.is_consistent())throw std::runtime_error("inconsistent reconstructed CSR");
    return matrix;
}
void put_vector(std::ofstream& f,const Vector& v){
    put(f,static_cast<std::uint64_t>(v.size()));
    for(std::size_t i=0;i<v.size();++i)put(f,v(i));
}
Vector get_vector(std::ifstream& f){
    const auto n=get<std::uint64_t>(f);
    if(n>50000)throw std::runtime_error("snapshot vector too large");
    Vector v(static_cast<std::size_t>(n));
    for(std::size_t i=0;i<v.size();++i)v(i)=get<double>(f);
    if(!v.is_valid())throw std::runtime_error("nonfinite snapshot vector");
    return v;
}
void save_snapshot(const Snapshot& s,const std::string& path){
    if(!s.valid)throw std::runtime_error("cannot save incomplete capture");
    std::ofstream f(path,std::ios::binary|std::ios::trunc);
    if(!f)throw std::runtime_error("unable to create snapshot file");
    f.write("CFDXV5A1",8);put(f,kEndianTag);
    put_matrix(f,s.coupled);put_vector(f,s.b);put_vector(f,s.x0);
    put_matrix(f,s.mass);put_matrix(f,s.kp);put_matrix(f,s.fp);
    put(f,static_cast<std::uint64_t>(s.gauge));
    put(f,static_cast<std::uint8_t>(s.fixed_pressure?1:0));
    put(f,s.max_old_velocity);
    const std::string digest=fingerprint(s);
    f.write(digest.data(),static_cast<std::streamsize>(digest.size()));
    if(!f)throw std::runtime_error("snapshot digest write failure");
}
Snapshot load_snapshot(const std::string& path){
    std::ifstream f(path,std::ios::binary);
    if(!f)throw std::runtime_error("unable to open snapshot");
    char magic[8]={};f.read(magic,8);
    if(!f||std::memcmp(magic,"CFDXV5A1",8)!=0 || get<std::uint64_t>(f)!=kEndianTag)
        throw std::runtime_error("wrong capture format or endianness");
    Snapshot s;
    s.coupled=get_matrix(f);s.b=get_vector(f);s.x0=get_vector(f);
    s.mass=get_matrix(f);s.kp=get_matrix(f);s.fp=get_matrix(f);
    s.gauge=static_cast<std::size_t>(get<std::uint64_t>(f));
    const auto flag=get<std::uint8_t>(f);
    if(flag>1)throw std::runtime_error("invalid pressure boundary flag");
    s.fixed_pressure=flag!=0;s.max_old_velocity=get<double>(f);
    char digest[16]={};f.read(digest,16);
    if(!f||f.peek()!=std::char_traits<char>::eof())
        throw std::runtime_error("truncated or appended capture file");
    s.valid=s.coupled.n_rows()!=0 && s.coupled.n_rows()%4==0 &&
       s.coupled.n_rows()==s.coupled.n_cols() &&
       s.b.size()==s.coupled.n_rows() && s.x0.size()==s.b.size() &&
       s.mass.n_rows()==s.coupled.n_rows()/4 &&
       s.kp.n_rows()==s.mass.n_rows() && s.fp.n_rows()==s.mass.n_rows() &&
       s.gauge<s.mass.n_rows() && s.max_old_velocity>1e-8 &&
       std::isfinite(s.max_old_velocity);
    if(!s.valid||fingerprint(s)!=std::string(digest,16))
        throw std::runtime_error("capture invariant or checksum failed");
    return s;
}

double true_relative_residual(const SparseMatrix& A,const Vector& b,const Vector& x){
    if(A.n_rows()!=b.size() || x.size()!=A.n_cols())return std::numeric_limits<double>::infinity();
    long double rr=0,bb=0;
    for(std::size_t i=0;i<A.n_rows();++i){
        long double sum=0;
        for(std::size_t k=A.row_offsets_data()[i]; k<A.row_offsets_data()[i+1];++k)
            sum += static_cast<long double>(A.values_data()[k])*x(A.columns_data()[k]);
        const long double r=static_cast<long double>(b(i))-sum;
        rr += r*r; bb += static_cast<long double>(b(i))*b(i);
    }
    if(bb==0) return rr==0 ? 0 : std::numeric_limits<double>::infinity();
    return static_cast<double>(std::sqrt(rr/bb));
}
std::size_t heap_live_kib(){
#if defined(__GLIBC__)
    const auto info=mallinfo2();
    return static_cast<std::size_t>(info.uordblks/1024);
#else
    return 0; // unavailable on this platform; never impute a memory value
#endif
}
std::size_t duplicate_column_entries(const SparseMatrix& a){
    std::size_t count=0;
    for(std::size_t row=0;row<a.n_rows();++row)
        for(std::size_t k=a.row_offsets_data()[row]+1;
            k<a.row_offsets_data()[row+1];++k)
            if(a.columns_data()[k]==a.columns_data()[k-1])++count;
    return count;
}
std::size_t rss_kib(){
#if defined(__linux__)
    std::ifstream f("/proc/self/status"); std::string name; std::size_t value=0;
    while(f>>name){if(name=="VmRSS:"){f>>value;return value;}std::string rest;std::getline(f,rest);}
#endif
    return 0;
}
std::size_t peak_kib(){
#if defined(__linux__) || defined(__APPLE__)
    struct rusage usage{};
    if(getrusage(RUSAGE_SELF,&usage)==0){
#if defined(__APPLE__)
        return static_cast<std::size_t>(usage.ru_maxrss/1024);
#else
        return static_cast<std::size_t>(usage.ru_maxrss);
#endif
    }
#endif
    return 0;
}
IncompressibleSolverControls common_controls(std::size_t nx,std::size_t ny){
    IncompressibleSolverControls c;
    c.algorithm=PressureVelocityAlgorithm::COUPLED;
    c.density=1.0; c.kinematic_viscosity=0.01; // Re=100, reference Ghia setup
    c.pressure_reference_cell=(ny/2)*nx+(nx/2); c.pressure_reference_value=0.0;
    c.coupling.alpha_u=0.7; c.coupling.alpha_p=0.3;
    c.coupling.coupled_max_iterations=400;
    c.coupling.coupled_linear_tolerance=1e-10;
    c.coupling.linear_solver.krylov=KrylovModel::FGMRES;
    c.coupling.linear_solver.preconditioner=PreconditionerModel::CoupledBlockSchur;
    c.coupling.linear_solver.gmres_restart=static_cast<int>(4*nx*ny);
    c.pcd_convection_diffusion_max_iterations=200;
    c.pcd_laplacian_max_iterations=100;
    c.pcd_convection_diffusion_linear_solver.gmres_restart=static_cast<int>(nx*ny);
    c.use_bounded_convection=true; c.convection_scheme=ConvectionScheme::UPWIND;
    c.convergence.max_iterations=350;
    c.convergence.relative_tolerance=1e-8;
    c.convergence.continuity_tolerance=1e-8;
    return c;
}
void configure_method(IncompressibleSolverControls& c, const std::string& model){
    if(model=="pcd")c.coupling.schur_model=CoupledSchurModel::PCD;
    else if(model=="lsc")c.coupling.schur_model=CoupledSchurModel::LSC;
    else if(model=="bfbt")c.coupling.schur_model=CoupledSchurModel::BFBT;
    else if(model=="amg")c.coupling.schur_model=CoupledSchurModel::BlockLocal;
    else if(model=="schwarz"){
        c.coupling.schur_model=CoupledSchurModel::BlockLocal;
        c.coupling.coupled_schwarz_experimental=true;
        c.coupling.coupled_schwarz_cells_per_subdomain=2;
        c.coupling.coupled_schwarz_overlap_hops=1;
        c.coupling.coupled_schwarz_max_local_unknowns=128;
        c.coupling.coupled_schwarz_weighted=false;
    }else throw std::invalid_argument("unknown preconditioner: "+model);
}

struct RunData {
    IncompressibleSolveResult result;
    Snapshot snapshot;
    double wall_ms=0;
    double rms_u=std::numeric_limits<double>::quiet_NaN();
    double rms_v=std::numeric_limits<double>::quiet_NaN();
    double pressure_range=std::numeric_limits<double>::quiet_NaN();
    double net_boundary_flux=std::numeric_limits<double>::quiet_NaN();
};
RunData run_physics(std::size_t nx, std::size_t ny,
                    const std::string& variant, std::size_t outer_budget, bool capture){
    Mesh mesh=make_cavity_mesh(nx,ny);
    Field<double,Location::CELL> U(mesh.n_cells(),"U","m/s",3);
    Field<double,Location::CELL> p(mesh.n_cells(),"p","Pa",1);
    U.fill(0.0);p.fill(0.0);
    VelocityBoundaryConditions ubc;
    for(const char* name:{"bottom","left","right","front","back"})
        ubc[name]={VelocityBoundaryCondition::Type::FIXED_VALUE,{0,0,0}};
    ubc["top"]={VelocityBoundaryCondition::Type::FIXED_VALUE,{1,0,0}};
    ScalarBoundaryConditions pbc;
    for(const char* name:{"bottom","top","left","right","front","back"})
        pbc[name]={ScalarBoundaryType::ZERO_GRADIENT,0,0};
    auto c=common_controls(nx,ny);
    configure_method(c,variant);
    c.convergence.max_iterations=outer_budget;
    RunData data;
    if(capture){
        std::size_t assembly_number=0;
        c.coupled_assembly_callback=[&,assembly_number](const CoupledAssemblySnapshot& s) mutable {
            ++assembly_number;
            if(assembly_number!=3)return; // stage 3: convection from nonzero interior velocity
            if(!s.matrix||!s.rhs||!s.initial_guess||!s.pcd_pressure_mass||
               !s.pcd_pressure_laplacian||!s.pcd_pressure_convection_diffusion)
                throw std::runtime_error("incomplete production PCD capture");
            auto& v=data.snapshot;
            v.coupled=*s.matrix; v.b=*s.rhs; v.x0=*s.initial_guess;
            v.mass=*s.pcd_pressure_mass;v.kp=*s.pcd_pressure_laplacian;
            v.fp=*s.pcd_pressure_convection_diffusion;
            v.gauge=s.reference_cell;v.fixed_pressure=s.has_fixed_pressure;
            if(!s.velocity_old) throw std::runtime_error("missing linearization velocity");
            for(std::size_t field=0;field<3;++field)for(std::size_t i=0;i<nx*ny;++i)
                v.max_old_velocity=std::max(v.max_old_velocity,
                      std::abs(s.velocity_old->component_data(field)[i]));
            v.valid=true;
        };
    }
    const auto t=Clock::now();
    data.result=solve_steady_incompressible(mesh,U,p,ubc,pbc,c);
    data.wall_ms=ms(t);
    long double uu=0,vv=0;
    double pmin=std::numeric_limits<double>::infinity();
    double pmax=-std::numeric_limits<double>::infinity();
    for(std::size_t c_id=0;c_id<mesh.n_cells();++c_id){
        const double ux=U.component_data(0)[c_id];
        const double uy=U.component_data(1)[c_id];
        uu+=static_cast<long double>(ux)*ux;
        vv+=static_cast<long double>(uy)*uy;
        pmin=std::min(pmin,p(c_id));pmax=std::max(pmax,p(c_id));
    }
    data.rms_u=std::sqrt(static_cast<double>(uu/mesh.n_cells()));
    data.rms_v=std::sqrt(static_cast<double>(vv/mesh.n_cells()));
    data.pressure_range=pmax-pmin;
    if(data.result.authoritative_mass_flux.size()==mesh.n_faces()){
        double net=0;
        for(std::size_t f=0;f<mesh.n_faces();++f)
            if(mesh.ownership().neighbour(f)<0)net+=data.result.authoritative_mass_flux(f);
        data.net_boundary_flux=net;
    }
    return data;
}

struct PreconditionerBundle {
    std::unique_ptr<Preconditioner> pc;
    // Callbacks for nested systems intentionally own pressure objects.
    std::shared_ptr<const SparseMatrix> operator_keepalive;
    std::shared_ptr<Preconditioner> nested_keepalive;
    std::size_t nested_iterations=0; // callbacks increment this
};

// Same V4 production algorithms and gauges, not a generic AMG on the
// entire indefinite saddle-point matrix.
std::unique_ptr<Preconditioner> prepare_preconditioner(
    const std::string& model, const Snapshot& s, std::size_t& inner_iterations){
    const std::size_t nc=s.coupled.n_rows()/4;
    const auto reference = s.fixed_pressure ? std::numeric_limits<std::size_t>::max():s.gauge;
    if(model=="schwarz"){
        CellSchwarzConfig cfg;
        cfg.cells_per_subdomain=2; cfg.overlap_hops=1;
        cfg.max_local_unknowns=128; cfg.restriction=CellSchwarzRestriction::Restricted;
        auto pre=std::make_unique<CellOverlappingSchwarzPreconditioner>(nc,cfg);
        if(!pre->setup(s.coupled)) throw std::runtime_error("Schwarz setup rejected");
        return pre;
    }
    CoupledBlockSchurOptions opts;
    opts.factorization=CoupledSchurFactorization::Full;
    opts.velocity_approximation=CoupledSchurVelocityApproximation::Block;
    auto pre=std::make_unique<CoupledBlockSchurAMGPreconditioner>(nc,opts,reference);
    if(model=="amg"){
        if(!pre->setup(s.coupled))throw std::runtime_error("block-local AMG setup failed: "+pre->last_error());
        return pre;
    }
    if(model=="pcd"){
        opts.schur_approximation=CoupledSchurApproximationModel::PCD;
        pre=std::make_unique<CoupledBlockSchurAMGPreconditioner>(nc,opts,reference);
        // Like production, the named PCD baseline is the V4 -Fp^-1 action
        // with nested FGMRES/ILU0, NOT the experimental PCD commutator.
        const auto nested=std::make_shared<SparseMatrix>(s.fp);
        auto nested_pc=std::shared_ptr<Preconditioner>(
            make_scalar_preconditioner(PreconditionerModel::ILU0).release());
        if(!nested_pc||!nested_pc->setup(*nested))throw std::runtime_error("PCD Fp ILU0 setup failed");
        auto solve_fp=[nested,nested_pc,&inner_iterations](const Vector& rhs,Vector& x){
            const auto result=solve_fgmres(*nested,rhs,x,
                static_cast<int>(std::min<std::size_t>(512,nested->n_rows())),
                200,1e-8,nested_pc.get(),{}, {},false);
            inner_iterations+=result.iterations;
            return nested_schur_pressure_solve_usable(result,x,nested->n_rows());
        };
        auto forbidden=[](const Vector&,Vector&)->bool{
            throw std::logic_error("legacy PCD unexpectedly called Kp solver");
        };
        auto pcd=std::make_unique<PcdSchurApproximation>(
            s.mass,s.kp,s.fp,forbidden,solve_fp,std::nullopt,
            s.fixed_pressure ? std::optional<std::size_t>() : std::optional<std::size_t>(s.gauge));
        pre->set_pcd_schur(std::move(pcd));
        if(!pre->setup(s.coupled))throw std::runtime_error("PCD coupled setup failed: "+pre->last_error());
        return pre;
    }
    if(model=="lsc"||model=="bfbt"){
        const auto Auu=extract_coupled_block(s.coupled,0,0,3*nc,3*nc);
        const auto G=extract_coupled_block(s.coupled,0,1,3*nc,nc);
        const auto D=extract_coupled_block(s.coupled,1,0,nc,3*nc);
        const auto C=extract_coupled_block(s.coupled,1,1,nc,nc);
        const BlockOperator blocks(Auu,G,D,C);
        const auto mode=model=="lsc"?LscBfbtSchurApproximation::Mode::LSC:
            LscBfbtSchurApproximation::Mode::BFBT;
        std::vector<double> q(3*nc,1.0);
        if(model=="lsc"){
            for(std::size_t i=0;i<3*nc;++i){
                double diagonal=0;
                for(std::size_t k=Auu.row_offsets_data()[i];k<Auu.row_offsets_data()[i+1];++k)
                    if(Auu.columns_data()[k]==i)diagonal+=Auu.values_data()[k];
                if(!(diagonal>0)&&std::isfinite(diagonal))throw std::runtime_error("invalid LSC diagonal");
                if(!std::isfinite(diagonal))throw std::runtime_error("non-finite LSC diagonal");
                q[i]=diagonal;
            }
        }
        auto H=std::make_shared<SparseMatrix>(
            LscBfbtSchurApproximation::assemble_stabilized_pressure_operator(blocks,mode,q));
        apply_nested_schur_pressure_gauge(*H,reference);
        auto nested_pc=std::shared_ptr<Preconditioner>(
            make_scalar_preconditioner(PreconditionerModel::NativeAMG,false,false).release());
        if(!nested_pc || !nested_pc->setup(*H)) throw std::runtime_error("LSC/BFBT H AMG setup failed");
        auto solve_h=[H,nested_pc,&inner_iterations](const Vector& rhs,Vector& x){
            const auto result=solve_gmres(*H,rhs,x,
                static_cast<int>(std::min<std::size_t>(64,H->n_rows())),
                400,1e-10,nested_pc.get(),{}, {},false);
            inner_iterations+=result.iterations;
            return nested_schur_pressure_solve_usable(result,x,H->n_rows());
        };
        auto approx=std::make_unique<LscBfbtSchurApproximation>(mode,solve_h,q);
        opts.schur_approximation=model=="lsc"?
            CoupledSchurApproximationModel::LSC:CoupledSchurApproximationModel::BFBT;
        pre=std::make_unique<CoupledBlockSchurAMGPreconditioner>(nc,opts,reference);
        pre->set_algebraic_schur(std::move(approx));
        if(!pre->setup(s.coupled)) throw std::runtime_error("LSC/BFBT coupled setup failed: "+pre->last_error());
        return pre;
    }
    throw std::invalid_argument("unknown method");
}

void write_line(const std::string& label, const std::string& text) {
    std::cout<<label<<" "<<text<<"\n";std::cout.flush();
}
void replay_snapshot(const std::string& variant,const Snapshot& s){
    const std::size_t nc=s.coupled.n_rows()/4;
    const std::size_t nx=static_cast<std::size_t>(std::sqrt(nc));
    if(!s.valid || s.coupled.n_rows()!=4*nx*nx || s.b.size()!=s.coupled.n_rows() ||
       s.x0.size()!=s.b.size() || s.mass.n_rows()!=nx*nx ||
       !(s.max_old_velocity>1e-8))
        throw std::runtime_error("capture was incomplete");
    const std::string digest=fingerprint(s);
    const std::size_t baseline_rss=rss_kib(), baseline_peak=peak_kib();
    const std::size_t baseline_heap=heap_live_kib();
    std::size_t inner=0;
    const auto setup_start=Clock::now();
    auto pre=prepare_preconditioner(variant,s,inner);
    const double setup_ms=ms(setup_start);
    const std::size_t rss_after_setup=rss_kib();
    const std::size_t heap_after_setup=heap_live_kib();
    Vector x=s.x0;
    KrylovControls opts;
    const int restart=static_cast<int>(std::min<std::size_t>(s.coupled.n_rows(),400));
    opts.restart_min=restart;opts.restart_max=restart;opts.adaptive_restart=false;
    const auto solve_start=Clock::now();
    const auto result=solve_fgmres(s.coupled,s.b,x,restart,400,1e-10,pre.get(),opts);
    const double solve_ms=ms(solve_start);
    const std::size_t rss_after_solve=rss_kib();
    const std::size_t heap_after_solve=heap_live_kib();
    const double rr=true_relative_residual(s.coupled,s.b,x);
    const bool valid=std::isfinite(rr)&&rr<=1e-10&&result.status==SolverStatus::CONVERGED;
    std::ostringstream os;
    os<<std::setprecision(15)<<"method="<<variant<<" grid="<<nx
      <<" fingerprint="<<digest<<" stage=3 max_old_velocity="<<s.max_old_velocity
      <<" n="<<s.coupled.n_rows()
      <<" nnz="<<s.coupled.nnz()
      <<" csr_duplicate_entries="<<duplicate_column_entries(s.coupled)
      <<" status="<<static_cast<int>(result.status)<<" converged="<<valid
      <<" outer_iterations="<<result.iterations<<" nested_iterations="<<inner
      <<" true_relative_residual="<<rr<<" setup_ms="<<setup_ms
      <<" solve_ms="<<solve_ms<<" total_linear_ms="<<setup_ms+solve_ms
      <<" baseline_capture_rss_kib="<<baseline_rss
      <<" baseline_capture_peak_kib="<<baseline_peak
      <<" rss_after_setup_kib="<<rss_after_setup
      <<" rss_after_solve_kib="<<rss_after_solve
      <<" heap_before_setup_kib="<<baseline_heap
      <<" heap_after_setup_kib="<<heap_after_setup
      <<" heap_after_solve_kib="<<heap_after_solve
      <<" heap_setup_delta_kib="<<(heap_after_setup>baseline_heap?heap_after_setup-baseline_heap:0)
      <<" setup_rss_delta_kib="<<(rss_after_setup>baseline_rss?rss_after_setup-baseline_rss:0)
      <<" process_peak_rss_kib="<<peak_kib();
    write_line("CFDX_LINEAR_V5",os.str());
    if(!valid)throw std::runtime_error("linear true-residual gate not satisfied");
}
void physical(const std::string& variant,std::size_t nx,std::size_t max_outer){
    auto data=run_physics(nx,nx,variant,max_outer,false);
    const bool has_final=!data.result.history.empty();
    const auto& final=data.result.history;
    std::ostringstream os;
    os<<std::setprecision(15)<<"method="<<variant<<" grid="<<nx
      <<" converged="<<data.result.converged
      <<" iterations="<<data.result.iterations
      <<" wall_ms="<<data.wall_ms
      <<" process_peak_rss_kib="<<peak_kib()
      <<" rss_end_kib="<<rss_kib()
      <<" heap_end_kib="<<heap_live_kib()
      <<" continuity_linf="<<(has_final?final.back().continuity_linf:std::numeric_limits<double>::infinity())
      <<" continuity_normalized="<<(has_final?final.back().continuity_normalized:std::numeric_limits<double>::infinity())
      <<" momentum_relative="<<(has_final?final.back().momentum_equation_residual_relative:std::numeric_limits<double>::infinity())
      <<" velocity_rms_u="<<data.rms_u
      <<" velocity_rms_v="<<data.rms_v
      <<" pressure_range="<<data.pressure_range
      <<" boundary_mass_flux="<<data.net_boundary_flux
      <<" convergence_reason="<<data.result.convergence_reason;
    write_line("CFDX_PHYSICAL_V5",os.str());
    if(!data.result.converged)throw std::runtime_error("physical nonlinear gates unmet");
}
} // namespace
int main(int argc,char** argv){
    try {
        if(argc<4)throw std::invalid_argument(
           "usage: linear METHOD GRID | capture GRID FILE | replay METHOD FILE | physical METHOD GRID [outer_budget]");
        const auto mode=std::string(argv[1]);
        const auto first=std::string(argv[2]);
        const auto last=std::string(argv[3]);
        if(mode=="capture"){
            const int nx=std::stoi(first);
            if(nx<4||nx>64)throw std::invalid_argument("grid out of bounds");
            auto s=run_physics(static_cast<std::size_t>(nx),static_cast<std::size_t>(nx),"pcd",3,true).snapshot;
            if(!s.valid||!(s.max_old_velocity>1e-8))throw std::runtime_error("stage3 snapshot missing convection");
            save_snapshot(s,last);
            const Snapshot check=load_snapshot(last); // byte exact, fail-closed
            std::cout<<"CFDX_CAPTURE_V5 grid="<<nx<<" n="<<check.coupled.n_rows()
                <<" nnz="<<check.coupled.nnz()
                <<" csr_duplicate_entries="<<duplicate_column_entries(check.coupled)
                <<" stage=3 max_old_velocity="
                <<check.max_old_velocity<<" fingerprint="<<fingerprint(check)<<"\n";
        }else if(mode=="replay"){
            replay_snapshot(first,load_snapshot(last));
        }else if(mode=="linear"||mode=="physical"){
            const int grid=std::stoi(last);
            if(grid<4||grid>64)throw std::invalid_argument("grid out of bounds");
            if(mode=="linear"){
                auto s=run_physics(static_cast<std::size_t>(grid),static_cast<std::size_t>(grid),"pcd",3,true).snapshot;
                replay_snapshot(first,s);
            }else{
                physical(first,static_cast<std::size_t>(grid),
                    argc>4?static_cast<std::size_t>(std::stoul(argv[4])):350);
            }
        }else throw std::invalid_argument("unknown mode");
        return 0;
    }catch(const std::exception& e){std::cerr<<"CFDX_V5_FAILURE "<<e.what()<<"\n";return 2;}
}

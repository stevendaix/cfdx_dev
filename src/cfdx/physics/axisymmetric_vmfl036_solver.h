#pragma once
#include <algorithm>
#include "cfdx/physics/forces.h"
#include <array>
#include <iostream>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace cfdx::physics::axisymmetric {

struct AxisymNode { double x{}, r{}; };
struct AxisymFace {
    std::size_t owner{};
    std::size_t neighbour{static_cast<std::size_t>(-1)};
    std::size_t n0{}, n1{};
    double cx{}, cr{}, ds{}, nx{}, nr{}, area{};
    bool sphere{}, outer{}, axis{};
};
struct AxisymCell {
    std::size_t i{}, j{};
    std::size_t nodes[4]{};
    double cx{}, cr{}, area{}, volume{};
    std::vector<std::size_t> faces;
};

struct AxisymMesh {
    std::size_t nr{}, nt{};
    double diameter{1.0}, outer_radius{50.0};
    std::vector<AxisymNode> nodes;
    std::vector<AxisymCell> cells;
    std::vector<AxisymFace> faces;
    std::vector<std::array<std::size_t,4>> cell_faces;

    static AxisymMesh read(const std::string& path)
    {
        AxisymMesh m;
        std::ifstream in(path);
        if (!in) throw std::runtime_error("cannot open axisymmetric mesh: " + path);
        std::string magic;
        in >> magic;
        if (magic != "CFDX_AXISYMMETRIC_POLAR_V2")
            throw std::runtime_error("invalid axisymmetric mesh format");
        in >> m.nr >> m.nt >> m.diameter >> m.outer_radius;
        if (m.nr < 2 || m.nt < 4) throw std::runtime_error("invalid axisymmetric mesh dimensions");
        m.nodes.resize((m.nr+1)*(m.nt+1));
        for (auto& p : m.nodes) if (!(in >> p.x >> p.r)) throw std::runtime_error("truncated axisymmetric mesh");
        m.cells.resize(m.nr*m.nt);
        m.cell_faces.resize(m.cells.size());

        auto nid=[&](std::size_t i,std::size_t j){return i*(m.nt+1)+j;};
        auto cid=[&](std::size_t i,std::size_t j){return i*m.nt+j;};

        // Build each quadrilateral once; face ownership is deterministic.
        for(std::size_t i=0;i<m.nr;++i) for(std::size_t j=0;j<m.nt;++j){
            auto& c=m.cells[cid(i,j)];
            c.i=i;c.j=j;
            // Counter-clockwise in the meridional (x,r) plane.
            c.nodes[0]=nid(i,j); c.nodes[1]=nid(i+1,j);
            c.nodes[2]=nid(i+1,j+1); c.nodes[3]=nid(i,j+1);
            const auto &a=m.nodes[c.nodes[0]],&b=m.nodes[c.nodes[1]],
                       &cc=m.nodes[c.nodes[2]],&d=m.nodes[c.nodes[3]];
            double cross=0.0,mx=0.0,mr=0.0;
            const AxisymNode q[4]={a,b,cc,d};
            for(int k=0;k<4;++k){
                const auto& p=q[k]; const auto& z=q[(k+1)%4];
                const double cr=p.x*z.r-z.x*p.r;
                cross+=cr; mx+=(p.x+z.x)*cr; mr+=(p.r+z.r)*cr;
            }
            c.area=0.5*cross;
            if(!(c.area>0.0)) throw std::runtime_error("axisymmetric cell has non-positive meridional area");
            c.cx=mx/(6.0*c.area); c.cr=mr/(6.0*c.area);
            const double first_moment=mr/6.0;
            c.volume=2.0*pi()*first_moment;
            if(!(c.volume>0.0) || !std::isfinite(c.volume))
                throw std::runtime_error("axisymmetric cell has non-positive revolution volume");
        }

        auto add_face=[&](std::size_t c,std::size_t n0,std::size_t n1,bool boundary){
            AxisymFace f;
            f.owner=c; f.n0=n0; f.n1=n1;
            const auto& p=m.nodes[n0]; const auto& q=m.nodes[n1];
            f.cx=0.5*(p.x+q.x); f.cr=0.5*(p.r+q.r);
            const double dx=q.x-p.x, dr=q.r-p.r;
            f.ds=std::hypot(dx,dr);
            if(!(f.ds>0.0)) throw std::runtime_error("degenerate axisymmetric face");
            // CCW cell boundary: right-hand normal is outward.
            f.nx=dr/f.ds; f.nr=-dx/f.ds;
            f.axis=boundary && std::abs(p.r)<1e-14 && std::abs(q.r)<1e-14;
            f.area=f.axis ? 0.0 : 2.0*pi()*std::max(0.0,f.cr)*f.ds;
            // Axis endpoints at the sphere/farfield corners belong only to the
            // zero-radius axis patch, never simultaneously to sphere/outer.
            f.sphere=boundary && m.cells[c].i==0 && !f.axis;
            f.outer=boundary && m.cells[c].i+1==m.nr && !f.axis;
            const std::size_t id=m.faces.size(); m.faces.push_back(f); m.cells[c].faces.push_back(id);
            return id;
        };
        // Build edges and pair internal faces by canonical endpoint key.
        struct Pending { std::size_t face; std::size_t cell; };
        std::vector<std::pair<std::pair<std::size_t,std::size_t>,Pending>> pending;
        for(std::size_t c=0;c<m.cells.size();++c){
            const auto &v=m.cells[c].nodes;
            const std::size_t edge[4][2]={{v[0],v[1]},{v[1],v[2]},{v[2],v[3]},{v[3],v[0]}};
            for(int e=0;e<4;++e){
                const std::pair<std::size_t,std::size_t> key{
                    std::min(edge[e][0],edge[e][1]), std::max(edge[e][0],edge[e][1])};
                auto it=std::find_if(pending.begin(),pending.end(),[&](const auto& z){return z.first==key;});
                if(it==pending.end()){
                    const auto fid=add_face(c,edge[e][0],edge[e][1],true);
                    pending.push_back({key,{fid,c}});
                } else {
                    const auto fid=it->second.face;
                    // Existing face belongs to the previous owner; current cell is neighbour.
                    m.faces[fid].neighbour=c;
                    // The face is internal now: clear the speculative boundary
                    // classification assigned when the edge was first created.
                    // A stale sphere flag on a ring-0 internal face (a radial or
                    // outer-arc edge first seen by the ring-0 cell) would
                    // otherwise pollute the drag integral, which must only
                    // sample genuine sphere boundary faces.
                    m.faces[fid].sphere=false;
                    m.faces[fid].outer=false;
                    m.faces[fid].axis=false;
                    m.cells[c].faces.push_back(fid);
                    pending.erase(it);
                }
            }
        }
        // Every unpaired edge is a boundary face. Validate the expected polar
        // topology (sphere=nt, outer=nt, axis=2*nr) fail-closed instead of
        // silently discarding the leftover pairing work.
        std::size_t bnd_sphere=0,bnd_outer=0,bnd_axis=0,bnd_other=0,other_cell=0;
        for(const auto& z:pending){
            const auto& f=m.faces[z.second.face];
            if(f.sphere) ++bnd_sphere;
            else if(f.outer) ++bnd_outer;
            else if(f.axis) ++bnd_axis;
            else { ++bnd_other; other_cell=z.second.cell; }
        }
        if(bnd_other!=0)
            throw std::runtime_error("axisymmetric mesh has an unclassified boundary face at cell "+
                std::to_string(other_cell));
        if(bnd_sphere!=m.nt || bnd_outer!=m.nt || bnd_axis!=2*m.nr)
            throw std::runtime_error("axisymmetric mesh boundary topology mismatch: sphere="+
                std::to_string(bnd_sphere)+" expected "+std::to_string(m.nt)+
                ", outer="+std::to_string(bnd_outer)+" expected "+std::to_string(m.nt)+
                ", axis="+std::to_string(bnd_axis)+" expected "+std::to_string(2*m.nr));
        // For paired faces, ensure the face orientation points from owner to neighbour.
        for(auto& f:m.faces){
            if(f.neighbour!=static_cast<std::size_t>(-1)){
                const auto& o=m.cells[f.owner]; const auto& n=m.cells[f.neighbour];
                const double dx=n.cx-o.cx, dr=n.cr-o.cr;
                if(f.nx*dx+f.nr*dr<0.0){f.nx=-f.nx;f.nr=-f.nr;std::swap(f.n0,f.n1);}
            }
        }
        for(std::size_t ci=0;ci<m.cells.size();++ci){
            const auto& cc=m.cells[ci];
            if(cc.faces.size()!=4) throw std::runtime_error("axisymmetric cell does not have four faces");
            m.cell_faces[ci]={cc.faces[0],cc.faces[1],cc.faces[2],cc.faces[3]};
        }
        return m;
    }

    static constexpr double pi(){return 3.141592653589793238462643383279502884;}
};

struct AxisymControls {
    std::size_t max_outer_iterations{1200};
    // momentum_sweeps counts true Gauss-Seidel passes over the whole domain
    // per outer iteration (not per cell).
    std::size_t momentum_sweeps{8};
    std::size_t pressure_sweeps{400};
    double momentum_tolerance{1e-8};
    double pressure_tolerance{1e-9};
    double continuity_tolerance{1e-7};
    double alpha_u{0.65};
    double alpha_p{0.25};
};

struct AxisymResult {
    bool converged{};
    std::size_t iterations{};
    double continuity{};
    double momentum_residual{};
    double pressure_correction{};
    double drag_pressure{};
    double drag_viscous{};
    double drag_total{};
    double cd_pressure{};
    double cd_viscous{};
    double cd_total{};
};

class VMFL036AxisymmetricSolver {
public:
    VMFL036AxisymmetricSolver(const AxisymMesh& mesh, double rho, double mu)
        : m_(mesh),rho_(rho),mu_(mu),n_(mesh.cells.size()),
          ux_(n_,1.0),ur_(n_,0.0),p_(n_,0.0),phi_(mesh.faces.size(),0.0),
          apx_(n_,1.0),apr_(n_,1.0),pc_(n_,0.0)
    {
        // Ring-0 sphere extrapolation stencil, shared by gradients() and the
        // pressure-correction gradient below: for every ring-0 cell, its
        // radially-outward neighbour (ring 1, same angular index) and the
        // geometric factor d(centroid -> sphere face) / d(centroid ->
        // neighbour centroid) used by the pressure-like sphere-wall sample.
        sphere_nb_.assign(n_,static_cast<std::size_t>(-1));
        sphere_w_.assign(n_,0.0);
        for(std::size_t c=0;c<n_;++c){
            if(m_.cells[c].i!=0) continue;
            for(const auto fid:m_.cells[c].faces){
                const auto& f=m_.faces[fid];
                if(f.neighbour==static_cast<std::size_t>(-1)) continue;
                const std::size_t other=(f.owner==c)?f.neighbour:f.owner;
                if(m_.cells[other].i==1 && m_.cells[other].j==m_.cells[c].j){
                    sphere_nb_[c]=other;
                    for(const auto fid2:m_.cells[c].faces){
                        const auto& f2=m_.faces[fid2];
                        if(f2.neighbour==static_cast<std::size_t>(-1) && f2.sphere)
                            sphere_w_[c]=std::hypot(f2.cx-m_.cells[c].cx,f2.cr-m_.cells[c].cr)
                                /std::max(1e-12,dist(m_.cells[c],m_.cells[other]));
                    }
                }
            }
        }
    }

    AxisymResult solve(const AxisymControls& c = {})
    {
        if(!(rho_>0&&mu_>0)) throw std::invalid_argument("invalid rho/mu");
        std::vector<double> gp_x(n_),gp_r(n_);
        double continuity=std::numeric_limits<double>::infinity();
        double mom=std::numeric_limits<double>::infinity();
        double pcorr=std::numeric_limits<double>::infinity();

        // Face fluxes are recomputed at the top of each outer iteration from
        // the freshly solved velocities plus the Rhie-Chow pressure-velocity
        // coupling term (see compute_fluxes). Carrying plain interpolated
        // face velocities re-injected the interpolation imbalance every
        // iteration and left the momentum residual (the per-iteration
        // velocity change) in a limit cycle ~1e-5 above the 2e-6 gate: the
        // campaign budget was exhausted at mom=2.1e-5 while continuity and
        // pcorr passed. The Rhie-Chow term restores the pressure-velocity
        // coupling at the face so the recomputed flux field is consistent
        // with the discrete momentum operator; the TS twin campaign
        // (diag-tools/port4_rhie_chow.ts, branch ci/long-campaign-diagnostics)
        // then converges in 75-180 outer iterations on coarse/medium/fine
        // with drag unchanged.
        compute_fluxes(gp_x,gp_r);
        for(std::size_t it=1;it<=c.max_outer_iterations;++it){
            gradients(p_,gp_x,gp_r,false,0.0);
            assemble_momentum(ux_,ur_,gp_x,gp_r,apx_,apr_,c);
            mom=momentum_residual(ux_,ur_,gp_x,gp_r,apx_,apr_);
            compute_fluxes(gp_x,gp_r);
            pcorr=pressure_correction(c);
            continuity=max_continuity();
            if(it%50==0 || it==1)
                std::cout<<"VMFL036_AXISYM ITER="<<it<<" continuity="<<continuity
                         <<" momentum="<<mom<<" pcorr="<<pcorr<<"\n";
            if(continuity<c.continuity_tolerance && mom<c.momentum_tolerance && pcorr<c.pressure_tolerance)
                return finish(true,it,continuity,mom,pcorr);
        }
        // Report the achieved metrics against their gates so a failed
        // campaign names the missing quantity instead of only "no convergence".
        std::cout<<"VMFL036_AXISYM NOT_CONVERGED budget="<<c.max_outer_iterations
                 <<" continuity="<<continuity<<" (gate "<<c.continuity_tolerance<<")"
                 <<" momentum="<<mom<<" (gate "<<c.momentum_tolerance<<")"
                 <<" pcorr="<<pcorr<<" (gate "<<c.pressure_tolerance<<")\n";
        return finish(false,c.max_outer_iterations,continuity,mom,pcorr);
    }

private:
    AxisymResult finish(bool converged,std::size_t it,double continuity,double mom,double pcorr)
    {
        std::vector<double> gux_x(n_),gux_r(n_),gur_x(n_),gur_r(n_);
        gradients(ux_,gux_x,gux_r,true,1.0); gradients(ur_,gur_x,gur_r,true,0.0);
        require_finite_field("ux",ux_);
        require_finite_field("ur",ur_);
        require_finite_field("p",p_);
        require_finite_field("grad_ux_x",gux_x);
        require_finite_field("grad_ux_r",gux_r);
        require_finite_field("grad_ur_x",gur_x);
        require_finite_field("grad_ur_r",gur_r);
        const auto d=drag(gux_x,gux_r,gur_x,gur_r);
        return {converged,it,continuity,mom,pcorr,d[0],d[1],d[2],d[3],d[4],d[5]};
    }

    static void require_finite_field(const char* name, const std::vector<double>& a)
    {
        for(std::size_t i=0;i<a.size();++i)
            if(!std::isfinite(a[i]))
                throw std::runtime_error(std::string("VMFL036 non-finite ")+name+
                    " at cell "+std::to_string(i)+" value="+std::to_string(a[i]));
    }

    const AxisymMesh& m_;
    double rho_,mu_;
    std::size_t n_;
    std::vector<double> ux_,ur_,p_,phi_,apx_,apr_,pc_;
    std::vector<std::size_t> sphere_nb_;
    std::vector<double> sphere_w_;

    static double dist(const AxisymCell&a,const AxisymCell&b){return std::hypot(a.cx-b.cx,a.cr-b.cr);}
    double boundary_ux(const AxisymFace& f,std::size_t o) const {
        if(f.sphere) return 0.0;
        if(f.outer && f.cx < 0.0) return 1.0; // documented inlet arc
        return ux_[o];
    }
    double boundary_ur(const AxisymFace& f,std::size_t o) const {
        if(f.sphere || f.axis) return 0.0;
        if(f.outer && f.cx < 0.0) return 0.0; // inlet arc
        return ur_[o]; // pressure outlet: zero-gradient velocity
    }

    void compute_fluxes(const std::vector<double>& gp_x,const std::vector<double>& gp_r){
        for(std::size_t k=0;k<m_.faces.size();++k){
            const auto& f=m_.faces[k];
            double u=0.0,v=0.0,rc=0.0;
            if(f.neighbour!=static_cast<std::size_t>(-1)){
                u=0.5*(ux_[f.owner]+ux_[f.neighbour]);
                v=0.5*(ur_[f.owner]+ur_[f.neighbour]);
                // Rhie-Chow and pressure correction use the same response
                // derived from the relaxed momentum equation:
                //   a_P^r = a_P / alpha_u,  U' = -(V/a_P^r) grad(p').
                // apx_/apr_ store a_P^r, so the response is V/apx_ or V/apr_.
                // Dividing by alpha_u again would undo relaxation a second time
                // and make the pressure matrix differ from the face-flux operator.
                const double d=dist(m_.cells[f.owner],m_.cells[f.neighbour]);
                const double Dox=m_.cells[f.owner].volume/apx_[f.owner];
                const double Dor=m_.cells[f.owner].volume/apr_[f.owner];
                const double Dnx=m_.cells[f.neighbour].volume/apx_[f.neighbour];
                const double Dnr=m_.cells[f.neighbour].volume/apr_[f.neighbour];
                const double ex=(m_.cells[f.neighbour].cx-m_.cells[f.owner].cx)/d;
                const double er=(m_.cells[f.neighbour].cr-m_.cells[f.owner].cr)/d;
                const double Dxf=0.5*(Dox+Dnx);
                const double Drf=0.5*(Dor+Dnr);
                // The implicit coefficient is (D_f e) dot n. The full
                // interpolated D grad(p) dot n remains explicit; using
                // D_f projected twice on n is inconsistent on skewed faces.
                const double response_s=Dxf*ex*f.nx+Drf*er*f.nr;
                rc=0.5*(Dox*gp_x[f.owner]+Dnx*gp_x[f.neighbour])*f.nx
                  +0.5*(Dor*gp_r[f.owner]+Dnr*gp_r[f.neighbour])*f.nr
                  -response_s*(p_[f.neighbour]-p_[f.owner])/d;
            }else{
                u=boundary_ux(f,f.owner); v=boundary_ur(f,f.owner);
            }
            if(!std::isfinite(f.area) || f.area<0.0 || !std::isfinite(u) || !std::isfinite(v))
                throw std::runtime_error("VMFL036 non-finite face flux input at face "+std::to_string(k));
            phi_[k]=rho_*f.area*(u*f.nx+v*f.nr+rc);
            if(!std::isfinite(phi_[k]))
                throw std::runtime_error("VMFL036 non-finite face flux at face "+std::to_string(k));
        }
    }

    void gradients(const std::vector<double>& q,std::vector<double>& gx,std::vector<double>& gr,bool velocity,double outer_value) const{
        gx.assign(n_,0.0); gr.assign(n_,0.0);
        for(std::size_t c=0;c<n_;++c){
            for(const auto fid:m_.cells[c].faces){
                const auto& f=m_.faces[fid];
                double qf=q[c];
                if(f.neighbour!=static_cast<std::size_t>(-1)){
                    // Face normals are stored owner->neighbour and the owner
                    // is always the lower-index cell, so on the low-theta
                    // radial face of every cell this cell IS f.neighbour and
                    // the previous arithmetic mean self-sampled qf=q[c] on
                    // ~25% of the internal faces, domain-wide (the gradient
                    // operator error on analytic linear fields was ~50% on
                    // gx and ~100% on gr, growing with r). Sample the actual
                    // opposite cell. Inverse-distance centroid weights are
                    // exact for linear fields on this geometrically stretched
                    // polar mesh, where the two centroids are not
                    // equidistant from the face centre (a plain 0.5/0.5 mean
                    // is only exact for equidistant stencils).
                    const std::size_t other=(f.owner==c)?f.neighbour:f.owner;
                    const double dc=std::max(1e-12,std::hypot(f.cx-m_.cells[c].cx,f.cr-m_.cells[c].cr));
                    const double dn=std::max(1e-12,std::hypot(f.cx-m_.cells[other].cx,f.cr-m_.cells[other].cr));
                    const double wc=1.0/dc, wn=1.0/dn;
                    qf=(q[c]*wc+q[other]*wn)/(wc+wn);
                }
                else if(velocity){
                    if(f.sphere) qf=0.0;
                    else if(f.outer && f.cx < 0.0) qf=outer_value;
                    else if(f.axis) qf=(outer_value==0.0?0.0:q[c]);
                } else if(f.outer && f.cx >= 0.0){
                    qf=0.0; // pressure outlet
                } else if(f.sphere && sphere_nb_[c]!=static_cast<std::size_t>(-1)){
                    // Pressure-like field on the sphere wall. The ring-0
                    // centroid sits deep inside a thick first ring, so the
                    // zero-gradient sample qf=q[c] is only first-order and
                    // biased the near-wall pressure force (Re=50 medium Cd
                    // -4.9% vs the Tabata-Itakura reference with the
                    // opposite-cell fix alone at 0.5/0.5). Linear
                    // extrapolation along the ring-0/ring-1 centroid
                    // direction is exact for linear fields (analytic
                    // gradient tests: operator error back to ~3%).
                    qf=q[c]+(q[c]-q[sphere_nb_[c]])*sphere_w_[c];
                }
                // All remaining pressure boundaries use zero normal gradient.
                // Face normals are stored oriented owner->neighbour. On faces
                // where this cell is the neighbour, the outward normal is the
                // negation of the stored one; without the sign flip the
                // Green-Gauss sum mixes inward-oriented contributions on half
                // of every stencil and the reconstructed gradient is wrong in
                // sign (every second ring/sector direction).
                const double s=(f.owner==c)?1.0:-1.0;
                gx[c]+=s*qf*f.nx*f.area;
                gr[c]+=s*qf*f.nr*f.area;
            }
            // Axisymmetric radial closure: the revolved face areas of a cell
            // do not sum to zero radially. For straight meridional edges each
            // face contributes n_r*2*pi*r_mid*ds = 2*pi*r_mid*dr exactly, so
            // sum_f n_r*A_f = 2*pi*A_meridional (divergence theorem on the
            // revolution volume: div(e_r)=1/r integrates to 2*pi*A_m). The
            // x-component closes to zero as usual. Without the subtraction
            // below, even a constant field produces a spurious radial
            // gradient q*2*pi*A_m/V ~ q/r_bar, which feeds the pressure force
            // of the momentum equation, the velocity correction and the wall
            // viscous tractions, and was the amplification mechanism behind
            // the divergent campaign runs.
            gr[c]-=q[c]*2.0*AxisymMesh::pi()*m_.cells[c].area;
            gx[c]/=m_.cells[c].volume; gr[c]/=m_.cells[c].volume;
        }
    }

    void assemble_momentum(std::vector<double>& ux,std::vector<double>& ur,
                           const std::vector<double>& gpx,const std::vector<double>& gpr,
                           std::vector<double>& apx,std::vector<double>& apr,
                           const AxisymControls& c){
        // Phase 1: assemble the under-relaxed momentum rows into flat
        // Gauss-Seidel stencils. Neighbour values are never read here, so
        // the sweep phase below sees a frozen, consistent operator.
        std::vector<std::array<std::size_t,4>> nb_id(n_);
        std::vector<std::array<double,4>> nb_an(n_);
        std::vector<std::size_t> nb_cnt(n_,0);
        std::vector<double> rhs_x(n_,0.0),rhs_r(n_,0.0);

        for(std::size_t cell=0;cell<n_;++cell){
            const auto& ce=m_.cells[cell];
            double ax=0.0, ar=0.0;
            double bx=-gpx[cell]*ce.volume;
            double br=-gpr[cell]*ce.volume;
            for(const auto fid:ce.faces){
                const auto& f=m_.faces[fid];
                const bool internal=f.neighbour!=static_cast<std::size_t>(-1);
                const double F=(f.owner==cell?phi_[fid]:-phi_[fid]);
                if(internal){
                    const std::size_t nb=(f.owner==cell?f.neighbour:f.owner);
                    const double d=dist(ce,m_.cells[nb]);
                    const double D=mu_*f.area/d;
                    const double an=D+std::max(-F,0.0);
                    ax+=D+std::max(F,0.0);
                    ar+=D+std::max(F,0.0);
                    nb_id[cell][nb_cnt[cell]]=nb;
                    nb_an[cell][nb_cnt[cell]]=an;
                    ++nb_cnt[cell];
                }else if(f.axis){
                    // Symmetry axis: ux has zero normal gradient, ur=0. The
                    // axis face has zero revolution area, so F and D vanish
                    // numerically; the terms are kept for residual consistency.
                    const double d=std::max(1e-12,std::abs((f.cx-ce.cx)*f.nx+(f.cr-ce.cr)*f.nr));
                    const double D=mu_*f.area/d;
                    ax+=F;
                    ar+=D+std::max(F,0.0);
                }else if(f.outer && f.cx >= 0.0){
                    // Pressure outlet: zero-gradient velocity, u_b = u_P. The
                    // exact upwind row contribution is F*u_P (signed): outflow
                    // max(F,0) and inflow max(-F,0) both cancel against u_b=u_P,
                    // and the diffusive term vanishes because u_b = u_P.
                    ax+=F;
                    ar+=F;
                }else{
                    // Sphere wall or inlet arc: Dirichlet velocity.
                    const double ubx=boundary_ux(f,cell), ubr=boundary_ur(f,cell);
                    const double d=std::max(1e-12,std::abs((f.cx-ce.cx)*f.nx+(f.cr-ce.cr)*f.nr));
                    const double D=mu_*f.area/d;
                    const double ab=D+std::max(-F,0.0);
                    ax+=D+std::max(F,0.0);
                    ar+=D+std::max(F,0.0);
                    bx+=ab*ubx;
                    br+=ab*ubr;
                }
            }
            // Cylindrical radial viscous term -nu*u_r/r^2, implicit.
            if(!std::isfinite(ce.cr) || !(ce.cr>0.0))
                throw std::runtime_error("VMFL036 invalid positive cell radius at cell "+std::to_string(cell));
            const double hoop=mu_*ce.volume/(ce.cr*ce.cr);
            if(!std::isfinite(hoop) || !(hoop>0.0))
                throw std::runtime_error("VMFL036 invalid radial geometric coefficient at cell "+std::to_string(cell));
            ar+=hoop;
            if(!std::isfinite(ax) || !std::isfinite(ar) || !(ax>0.0) || !(ar>0.0)) {
                std::string message="VMFL036 invalid momentum coefficients at cell "+
                    std::to_string(cell)+" ij=("+
                    std::to_string(ce.i)+","+std::to_string(ce.j)+") ax="+
                    std::to_string(ax)+" ar="+std::to_string(ar);
                // Keep the fail-closed error actionable: negative diagonals
                // can be caused by a backflowing zero-gradient pressure outlet.
                // Record each face's signed mass flux and classification so a
                // failure identifies the offending boundary without changing
                // the discretization or suppressing the negative coefficient.
                for(const auto fid:ce.faces) {
                    const auto& f=m_.faces[fid];
                    const double F=(f.owner==cell?phi_[fid]:-phi_[fid]);
                    message+=" [face="+std::to_string(fid)+
                        " F="+std::to_string(F)+
                        " internal="+std::to_string(
                            f.neighbour!=static_cast<std::size_t>(-1))+
                        " outer="+std::to_string(f.outer)+
                        " cx="+std::to_string(f.cx)+
                        " sphere="+std::to_string(f.sphere)+
                        " axis="+std::to_string(f.axis)+"]";
                }
                throw std::runtime_error(message);
            }
            // Under-relaxation: relaxed diagonal a_P/alpha_u on the left,
            // explicit source from the pre-solve values on the right.
            const double oldx=ux[cell], oldr=ur[cell];
            apx[cell]=ax/c.alpha_u;
            apr[cell]=ar/c.alpha_u;
            bx+=(1.0-c.alpha_u)/c.alpha_u*ax*oldx;
            br+=(1.0-c.alpha_u)/c.alpha_u*ar*oldr;
            if(!std::isfinite(bx) || !std::isfinite(br))
                throw std::runtime_error("VMFL036 non-finite momentum source at cell "+std::to_string(cell));
            rhs_x[cell]=bx;
            rhs_r[cell]=br;
        }

        // Phase 2: true Gauss-Seidel sweeps over the whole domain. Each sweep
        // updates every cell once with the latest available neighbour values.
        for(std::size_t s=0;s<c.momentum_sweeps;++s){
            for(std::size_t cell=0;cell<n_;++cell){
                double sumx=rhs_x[cell], sumr=rhs_r[cell];
                for(std::size_t k=0;k<nb_cnt[cell];++k){
                    const std::size_t nb=nb_id[cell][k];
                    const double an=nb_an[cell][k];
                    sumx+=an*ux[nb];
                    sumr+=an*ur[nb];
                }
                const double sx=sumx/apx[cell];
                const double sr=sumr/apr[cell];
                if(!std::isfinite(sx) || !std::isfinite(sr))
                    throw std::runtime_error(
                        "VMFL036 non-finite momentum update at cell "+std::to_string(cell)+
                        " sweep="+std::to_string(s)+
                        " ux="+std::to_string(sx)+" ur="+std::to_string(sr)+
                        " apx="+std::to_string(apx[cell])+" apr="+std::to_string(apr[cell]));
                ux[cell]=sx;
                ur[cell]=sr;
            }
        }
        // Boundary conditions are imposed through boundary-face fluxes and
        // diffusion stencils. Do not overwrite cell-centred unknowns here:
        // boundary cells still represent fluid volume and must remain solvable.
    }

    double momentum_residual(const std::vector<double>& ux,const std::vector<double>& ur,
                             const std::vector<double>& gpx,const std::vector<double>& gpr,
                             const std::vector<double>&,const std::vector<double>&) const{
        double r=0.0;
        for(std::size_t c=0;c<n_;++c){
            double rx=gpx[c]*m_.cells[c].volume, rr=gpr[c]*m_.cells[c].volume;
            for(const auto fid:m_.cells[c].faces){
                const auto& f=m_.faces[fid]; const double F=(f.owner==c?phi_[fid]:-phi_[fid]);
                const bool in=f.neighbour!=static_cast<std::size_t>(-1);
                const std::size_t nb=in?(f.owner==c?f.neighbour:f.owner):c;
                const double d=in?dist(m_.cells[c],m_.cells[nb]):std::max(1e-12,std::abs((f.cx-m_.cells[c].cx)*f.nx+(f.cr-m_.cells[c].cr)*f.nr));
                const double D=mu_*f.area/d;
                const double an=D+std::max(-F,0.0);
                const double ubx=in?ux[nb]:boundary_ux(f,c);
                const double ubr=in?ur[nb]:boundary_ur(f,c);
                rx += (D+std::max(F,0.0))*ux[c]-an*ubx;
                rr += (D+std::max(F,0.0))*ur[c]-an*ubr;
            }
            rr += mu_*m_.cells[c].volume*ur[c]/std::max(m_.cells[c].cr*m_.cells[c].cr,1e-20);
            r=std::max(r,std::max(std::abs(rx),std::abs(rr))/(1.0+std::abs(ux[c])+std::abs(ur[c])));
        }
        return r;
    }

    // Coefficient of the centre-to-centre pressure difference on an
    // internal face: (D_f e) dot n, using the relaxed diagonal response.
    double internal_face_response(
        std::size_t owner,std::size_t neighbour,double nx,double nr) const {
        const auto& o=m_.cells[owner];
        const auto& n=m_.cells[neighbour];
        const double dx=n.cx-o.cx, dr=n.cr-o.cr;
        const double d=std::hypot(dx,dr);
        if (!(d>0.0) || !std::isfinite(d))
            throw std::runtime_error("VMFL036 degenerate internal face centre distance");
        const double ex=dx/d, er=dr/d;
        const double Dox=o.volume/apx_[owner];
        const double Dor=o.volume/apr_[owner];
        const double Dnx=n.volume/apx_[neighbour];
        const double Dnr=n.volume/apr_[neighbour];
        return 0.5*(Dox+Dnx)*ex*nx + 0.5*(Dor+Dnr)*er*nr;
    }

    // Boundary counterpart: e points from the cell centroid to the face centre.
    double boundary_face_response(std::size_t cell,const AxisymFace& face) const {
        const auto& c=m_.cells[cell];
        const double dx=face.cx-c.cx, dr=face.cr-c.cr;
        const double d=std::hypot(dx,dr);
        if (!(d>0.0) || !std::isfinite(d))
            throw std::runtime_error("VMFL036 degenerate boundary face centre distance");
        const double ex=dx/d, er=dr/d;
        return (c.volume/apx_[cell])*ex*face.nx
             + (c.volume/apr_[cell])*er*face.nr;
    }

    double pressure_correction(const AxisymControls& c){
        pc_.assign(n_,0.0);
        std::vector<double> diag(n_,0.0),rhs(n_,0.0);
        // Flat off-diagonal stencils; every cell has at most four faces.
        std::vector<std::array<std::size_t,4>> prow_nb(n_);
        std::vector<std::array<double,4>> prow_a(n_);
        std::vector<std::size_t> prow_cnt(n_,0);
        for(std::size_t cell=0;cell<n_;++cell){
            double imbalance=0.0;
            for(const auto fid:m_.cells[cell].faces)
                imbalance += (m_.faces[fid].owner==cell?phi_[fid]:-phi_[fid]);
            rhs[cell]=-imbalance;
            for(const auto fid:m_.cells[cell].faces){
                const auto& f=m_.faces[fid];
                if(f.neighbour!=static_cast<std::size_t>(-1)){
                    const std::size_t nb=(f.owner==cell?f.neighbour:f.owner);
                    const double d=dist(m_.cells[cell],m_.cells[nb]);
                    // Use the same owner-to-neighbour centre-line response
                    // as compute_fluxes() and the conservative flux update.
                    // This keeps the implicit pressure coefficient identical
                    // to the pressure-jump part of the face-flux operator.
                    const double df=internal_face_response(
                        f.owner,f.neighbour,f.nx,f.nr);
                    const double a=rho_*f.area*df/d;
                    diag[cell]+=a;
                    prow_nb[cell][prow_cnt[cell]]=nb;
                    prow_a[cell][prow_cnt[cell]]=-a;
                    ++prow_cnt[cell];
                } else if(f.outer && f.cx >= 0.0){
                    const double d=std::hypot(
                        f.cx-m_.cells[cell].cx,f.cr-m_.cells[cell].cr);
                    const double response=boundary_face_response(cell,f);
                    if (!(response>0.0) || !std::isfinite(response))
                        throw std::runtime_error(
                            "VMFL036 invalid outlet pressure response at face "+
                            std::to_string(fid));
                    const double a=rho_*f.area*response/d;
                    // Fixed physical pressure at the outlet implies a
                    // homogeneous pressure-correction boundary condition:
                    // p'_b = 0. Do not inject -p_P into the correction RHS.
                    diag[cell]+=a;
                }
            }
        }
        for(std::size_t cell=0;cell<n_;++cell){
            if(!std::isfinite(diag[cell]) || !(diag[cell]>0.0))
                throw std::runtime_error("VMFL036 invalid pressure-correction diagonal at cell "+
                    std::to_string(cell)+" value="+std::to_string(diag[cell]));
            if(!std::isfinite(rhs[cell]))
                throw std::runtime_error("VMFL036 non-finite pressure-correction RHS at cell "+
                    std::to_string(cell)+" value="+std::to_string(rhs[cell]));
        }
        // The pressure-outlet arc supplies the Dirichlet pressure-correction
        // condition; no artificial reference cell is required.
        double maxcorr=std::numeric_limits<double>::infinity();
        for(std::size_t s=0;s<c.pressure_sweeps;++s){
            maxcorr=0.0;
            for(std::size_t cell=0;cell<n_;++cell){
                double sum=rhs[cell];
                // prow_a stores the assembled off-diagonal -a_N (see the
                // assembly above). Subtracting it adds +a_N*p'_N: the
                // Gauss-Seidel update of the assembled equation
                // a_P*p'_P - sum_N a_N*p'_N = b_P. Adding it instead would
                // solve the sign-flipped system while still satisfying maxcorr.
                for(std::size_t k=0;k<prow_cnt[cell];++k)
                    sum-=prow_a[cell][k]*pc_[prow_nb[cell][k]];
                // Successive over-relaxation, omega=1.8: plain Gauss-Seidel
                // contracts at q~0.998/sweep on the campaign meshes, leaving
                // relative residual 5.8e-3 after 400 sweeps on 36x54 - the
                // 2e-7 pcorr gate is unreachable and the momentum residual
                // inherits the floor mom ~ pres_tol*q/(1-q). SOR measured
                // stable on coarse/medium/fine: 36x54 relres 2.9e-5@400,
                // 1e-8@800; 48x72 6.5e-4@400, 8.2e-6@800.
                const double old=pc_[cell];
                const double v=old+1.8*(sum/diag[cell]-old);
                maxcorr=std::max(maxcorr,std::abs(v-old));
                pc_[cell]=v;
            }
            if(maxcorr<c.pressure_tolerance){
                // Combined stopping criterion. A small iterate increment
                // (maxcorr) only proves the sweep stopped moving; at a
                // maxcorr-only stop the assembled-equation row residual is
                // ~34*maxcorr. With the conservative flux correction above,
                // the post-correction cell imbalance IS that row residual, so
                // stopping on maxcorr alone would let the continuity metric of
                // the iteration sit far above the reported pcorr. Also require
                // the row residual at the continuity-like normalisation to be
                // below the tolerance before breaking.
                double res_rel=0.0;
                for(std::size_t cell=0;cell<n_;++cell){
                    double r=diag[cell]*pc_[cell];
                    for(std::size_t k=0;k<prow_cnt[cell];++k)
                        r+=prow_a[cell][k]*pc_[prow_nb[cell][k]];
                    double scale=0.0;
                    for(const auto fid:m_.cells[cell].faces)
                        scale+=std::abs(phi_[fid]);
                    r-=rhs[cell];
                    res_rel=std::max(res_rel,std::abs(r)/(1.0+scale));
                }
                if(res_rel<c.pressure_tolerance) break;
            }
        }
        // Algebraic consistency gate: a small iterate increment (maxcorr)
        // only proves the sweep stopped moving, not that it stopped on the
        // fixed point of the ASSEMBLED equation. When the sweep reports
        // convergence, verify the residual of diag*p'_P + sum(prow_a*p'_N)
        // - rhs directly: a converged solve leaves a row residual bounded
        // by the neighbour increments (<= sum a_N*maxcorr), while an
        // operator/sign regression converges per maxcorr onto the wrong
        // system and leaves a solution-sized residual. Fail closed on the
        // latter; this is the regression guard for the whole error class,
        // independent of the stopping criterion.
        if(maxcorr<c.pressure_tolerance){
            double res=0.0, diag_max=0.0;
            for(std::size_t cell=0;cell<n_;++cell){
                double r=diag[cell]*pc_[cell];
                for(std::size_t k=0;k<prow_cnt[cell];++k)
                    r+=prow_a[cell][k]*pc_[prow_nb[cell][k]];
                r-=rhs[cell];
                res=std::max(res,std::abs(r));
                diag_max=std::max(diag_max,diag[cell]);
            }
            if(res>100.0*diag_max*maxcorr)
                throw std::runtime_error("VMFL036 pressure-correction sweep converged per maxcorr="+
                    std::to_string(maxcorr)+" but violates the assembled pressure equation (residual "+
                    std::to_string(res)+"): operator/sign regression in the sweep");
        }
        // pressure correction gradient; boundary p'=0.
        for(std::size_t cell=0;cell<n_;++cell){
            double gx=0.0,gr=0.0;
            for(const auto fid:m_.cells[cell].faces){
                const auto& f=m_.faces[fid];
                double qf=pc_[cell];
                if(f.neighbour!=static_cast<std::size_t>(-1)){
                    // Same opposite-cell + inverse-distance face sampling as
                    // gradients(): the previous arithmetic mean
                    // 0.5*(pc_[cell]+pc_[f.neighbour]) self-sampled on every
                    // face where this cell is the face neighbour, corrupting
                    // the velocity correction on ~25% of the internal faces.
                    const std::size_t other=(f.owner==cell)?f.neighbour:f.owner;
                    const double dc=std::max(1e-12,std::hypot(f.cx-m_.cells[cell].cx,f.cr-m_.cells[cell].cr));
                    const double dn=std::max(1e-12,std::hypot(f.cx-m_.cells[other].cx,f.cr-m_.cells[other].cr));
                    const double wc=1.0/dc, wn=1.0/dn;
                    qf=(pc_[cell]*wc+pc_[other]*wn)/(wc+wn);
                }
                else if(f.outer && f.cx >= 0.0)
                    qf=0.0; // pressure correction is zero on a fixed-pressure outlet
                else if(f.sphere && sphere_nb_[cell]!=static_cast<std::size_t>(-1))
                    qf=pc_[cell]+(pc_[cell]-pc_[sphere_nb_[cell]])*sphere_w_[cell];
                // Owner->neighbour stored normal: flip the orientation for
                // faces where this cell is the neighbour, exactly like in
                // gradients().
                const double s=(f.owner==cell)?1.0:-1.0;
                gx+=s*qf*f.nx*f.area; gr+=s*qf*f.nr*f.area;
            }
            // Axisymmetric radial closure (see gradients()): subtract the exact
            // identity sum_f n_r*A_f = 2*pi*A_meridional applied to the local
            // p' value, so a uniform p' field produces no radial gradient of
            // the correction. Without it the velocity correction
            // u' = -rAU*grad p' carries a spurious p'/r_bar radial component
            // every outer iteration.
            gr-=pc_[cell]*2.0*AxisymMesh::pi()*m_.cells[cell].area;
            gx/=m_.cells[cell].volume; gr/=m_.cells[cell].volume;
            if(!std::isfinite(apx_[cell]) || !(apx_[cell]>0.0) ||
               !std::isfinite(apr_[cell]) || !(apr_[cell]>0.0))
                throw std::runtime_error("VMFL036 invalid momentum diagonal during pressure correction at cell "+
                    std::to_string(cell));
            ux_[cell]-=m_.cells[cell].volume/apx_[cell]*gx;
            ur_[cell]-=m_.cells[cell].volume/apr_[cell]*gr;
            if(!std::isfinite(ux_[cell]) || !std::isfinite(ur_[cell]))
                throw std::runtime_error("VMFL036 non-finite corrected velocity at cell "+
                    std::to_string(cell));
            p_[cell]+=c.alpha_p*pc_[cell];
            if(!std::isfinite(p_[cell]))
                throw std::runtime_error("VMFL036 non-finite corrected pressure at cell "+
                    std::to_string(cell));
        }
        // Do not impose inlet/outlet/wall values by overwriting cell centres.
        // The boundary-face treatment above is the actual finite-volume BC.
        // Conservative flux correction: instead of recomputing face fluxes
        // from the corrected velocities, apply the discrete continuity
        // statement of the assembled pressure equation directly to the
        // carried fluxes:
        //   internal faces: phi_f += rho*A*((D_f e) dot n)/d
        //                                  *(p'_owner - p'_neighbour)
        //   outlet faces   : phi_f += rho*A*((D_o e_b) dot n)/d*p'_owner
        //                   (p'_b = 0). These are the same coefficients used
        //                   in pressure_correction() above.
        //   sphere/axis    : no pressure coefficient was assembled, so no
        //                   flux correction (wall/symmetry velocities do not
        //                   respond to p' in this operator)
        // With these increments the post-correction face imbalance of every
        // cell equals the row residual of the pressure solve (exact
        // identity), so the continuity metric tracks the pressure solve
        // instead of stalling on the face-interpolation error of a recomputed
        // flux field (measured in the TS port: 5.5e-8..1e-12 with the
        // conservative correction vs a 1e-4 stall without it).
        for(std::size_t k=0;k<m_.faces.size();++k){
            const auto& f=m_.faces[k];
            if(!std::isfinite(phi_[k]))
                throw std::runtime_error("VMFL036 non-finite face flux before conservative correction at face "+
                    std::to_string(k));
            if(f.neighbour!=static_cast<std::size_t>(-1)){
                const double d=dist(m_.cells[f.owner],m_.cells[f.neighbour]);
                const double df=internal_face_response(
                    f.owner,f.neighbour,f.nx,f.nr);
                phi_[k]+=rho_*f.area*df/d*(pc_[f.owner]-pc_[f.neighbour]);
            } else if(f.outer && f.cx >= 0.0){
                const std::size_t o=f.owner;
                const double d=std::hypot(f.cx-m_.cells[o].cx,f.cr-m_.cells[o].cr);
                const double response=boundary_face_response(o,f);
                if (!(response>0.0) || !std::isfinite(response))
                    throw std::runtime_error(
                        "VMFL036 invalid outlet pressure response during flux correction at face "+
                        std::to_string(k));
                phi_[k]+=rho_*f.area*response/d*pc_[o];
            }
            if(!std::isfinite(phi_[k]))
                throw std::runtime_error("VMFL036 non-finite face flux after conservative correction at face "+
                    std::to_string(k));
        }
        return maxcorr;
    }

    double max_continuity() const{
        double r=0.0;
        for(std::size_t c=0;c<n_;++c){
            double b=0.0,scale=0.0;
            for(const auto fid:m_.cells[c].faces){b+=(m_.faces[fid].owner==c?phi_[fid]:-phi_[fid]);scale+=std::abs(phi_[fid]);}
            r=std::max(r,std::abs(b)/(1.0+scale));
        }
        return r;
    }

    std::array<double,6> drag(const std::vector<double>& gx,const std::vector<double>& gr,
                              const std::vector<double>& rx,const std::vector<double>& rr) const{
        std::vector<cfdx::physics::forces::AxisymmetricSample> samples;
        for(std::size_t c=0;c<n_;++c){
            if(m_.cells[c].i!=0) continue;
            for(const auto fid:m_.cells[c].faces){
                const auto& f=m_.faces[fid];
                // Only genuine sphere boundary faces carry body-surface drag:
                // internal faces must never reach the force integral even if a
                // stale boundary classification survived the mesh reader.
                if(f.neighbour!=static_cast<std::size_t>(-1) || !f.sphere) continue;
                // Mesh face normal points out of the fluid domain, i.e. into the
                // body. The force library uses the body->fluid normal.
                const double nx=-f.nx;
                const double nr=-f.nr;
                const double tau_xx=2.0*mu_*gx[c];
                const double tau_xr=mu_*(gr[c]+rx[c]);
                const double tau_rr=2.0*mu_*(rr[c]);
                samples.push_back({f.cx,f.cr,nx,nr,f.ds,p_[c],tau_xx,tau_xr,tau_rr});
            }
        }
        const auto result=cfdx::physics::forces::integrate_axisymmetric(
            samples,rho_,1.0,AxisymMesh::pi()*m_.diameter*m_.diameter/4.0);
        return {result.pressure_force,result.viscous_force,result.total_force,
                result.cd_pressure,result.cd_viscous,result.cd_total};
    }
};

} // namespace cfdx::physics::axisymmetric


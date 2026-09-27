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

        struct Key { std::size_t a,b; bool operator==(const Key& o) const{return a==o.a&&b==o.b;} };
        std::vector<std::pair<Key,std::size_t>> edge_map;
        edge_map.reserve(2*m.cells.size());
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
            f.area=2.0*pi()*std::max(0.0,f.cr)*f.ds;
            f.axis=boundary && std::abs(p.r)<1e-14 && std::abs(q.r)<1e-14;
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
                    m.cells[c].faces.push_back(fid);
                    pending.erase(it);
                }
            }
        }
        for(const auto& z:pending){
            (void)z;
        }
        // For paired faces, ensure the face orientation points from owner to neighbour.
        for(auto& f:m.faces){
            if(f.neighbour!=static_cast<std::size_t>(-1)){
                const auto& o=m.cells[f.owner]; const auto& n=m.cells[f.neighbour];
                const double dx=n.cx-o.cx, dr=n.cr-o.cr;
                if(f.nx*dx+f.nr*dr<0.0){f.nx=-f.nx;f.nr=-f.nr;std::swap(f.n0,f.n1);}
            }
        }
        for(auto& c:m.cells){
            if(c.faces.size()!=4) throw std::runtime_error("axisymmetric cell does not have four faces");
        }
        return m;
    }

    static constexpr double pi(){return 3.141592653589793238462643383279502884;}
};

struct AxisymControls {
    std::size_t max_outer_iterations{1200};
    std::size_t momentum_sweeps{80};
    std::size_t pressure_sweeps{250};
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
          apx_(n_,1.0),apr_(n_,1.0),pc_(n_,0.0) {}

    AxisymResult solve(const AxisymControls& c = {})
    {
        if(!(rho_>0&&mu_>0)) throw std::invalid_argument("invalid rho/mu");
        std::vector<double> gp_x(n_),gp_r(n_);
        std::vector<double> gux_x(n_),gux_r(n_),gur_x(n_),gur_r(n_);
        double continuity=std::numeric_limits<double>::infinity();
        double mom=std::numeric_limits<double>::infinity();
        double pcorr=std::numeric_limits<double>::infinity();

        for(std::size_t it=1;it<=c.max_outer_iterations;++it){
            compute_fluxes();
            gradients(p_,gp_x,gp_r,false,0.0);
            assemble_momentum(ux_,ur_,gp_x,gp_r,apx_,apr_,c);
            mom=momentum_residual(ux_,ur_,gp_x,gp_r,apx_,apr_);
            compute_fluxes();
            pcorr=pressure_correction(c);
            continuity=max_continuity();
            if(it%50==0 || it==1)
                std::cout<<"VMFL036_AXISYM ITER="<<it<<" continuity="<<continuity
                         <<" momentum="<<mom<<" pcorr="<<pcorr<<"\n";
            if(continuity<c.continuity_tolerance && mom<c.momentum_tolerance && pcorr<c.pressure_tolerance){
                gradients(ux_,gux_x,gux_r,true,1.0); gradients(ur_,gur_x,gur_r,true,0.0);
                require_finite_field("ux",ux_);
                require_finite_field("ur",ur_);
                require_finite_field("p",p_);
                require_finite_field("grad_ux_x",gux_x);
                require_finite_field("grad_ux_r",gux_r);
                require_finite_field("grad_ur_x",gur_x);
                require_finite_field("grad_ur_r",gur_r);
                const auto d=drag(gux_x,gux_r,gur_x,gur_r);
                return {true,it,continuity,mom,pcorr,d[0],d[1],d[2],d[3],d[4],d[5]};
            }
        }
        gradients(ux_,gux_x,gux_r,true,1.0); gradients(ur_,gur_x,gur_r,true,0.0);
        require_finite_field("ux",ux_);
        require_finite_field("ur",ur_);
        require_finite_field("p",p_);
        require_finite_field("grad_ux_x",gux_x);
        require_finite_field("grad_ux_r",gux_r);
        require_finite_field("grad_ur_x",gur_x);
        require_finite_field("grad_ur_r",gur_r);
        const auto d=drag(gux_x,gux_r,gur_x,gur_r);
        return {false,c.max_outer_iterations,continuity,mom,pcorr,d[0],d[1],d[2],d[3],d[4],d[5]};
    }

private:
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

    void compute_fluxes(){
        for(std::size_t k=0;k<m_.faces.size();++k){
            const auto& f=m_.faces[k];
            double u=0.0,v=0.0;
            if(f.neighbour!=static_cast<std::size_t>(-1)){
                u=0.5*(ux_[f.owner]+ux_[f.neighbour]);
                v=0.5*(ur_[f.owner]+ur_[f.neighbour]);
            }else{
                u=boundary_ux(f,f.owner); v=boundary_ur(f,f.owner);
            }
            if(!std::isfinite(f.area) || f.area<0.0 || !std::isfinite(u) || !std::isfinite(v))
                throw std::runtime_error("VMFL036 non-finite face flux input at face "+std::to_string(k));
            phi_[k]=rho_*f.area*(u*f.nx+v*f.nr);
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
                if(f.neighbour!=static_cast<std::size_t>(-1))
                    qf=0.5*(q[c]+q[f.neighbour]);
                else if(velocity){
                    if(f.sphere) qf=0.0;
                    else if(f.outer && f.cx < 0.0) qf=outer_value;
                    else if(f.axis) qf=(outer_value==0.0?0.0:q[c]);
                } else if(f.outer && f.cx >= 0.0){
                    qf=0.0; // pressure outlet
                }
                // All remaining pressure boundaries use zero normal gradient.
                gx[c]+=qf*f.nx*f.area;
                gr[c]+=qf*f.nr*f.area;
            }
            gx[c]/=m_.cells[c].volume; gr[c]/=m_.cells[c].volume;
        }
    }

    void assemble_momentum(std::vector<double>& ux,std::vector<double>& ur,
                           const std::vector<double>& gpx,const std::vector<double>& gpr,
                           std::vector<double>& apx,std::vector<double>& apr,
                           const AxisymControls& c){
        std::vector<double> nxux(n_,0.0),nxur(n_,0.0);
        for(std::size_t cell=0;cell<n_;++cell){
            double ax=0.0,ar=0.0,bx=-gpx[cell]*m_.cells[cell].volume;
            double br=-gpr[cell]*m_.cells[cell].volume;
            const auto& ce=m_.cells[cell];
            for(const auto fid:ce.faces){
                const auto& f=m_.faces[fid];
                const bool internal=f.neighbour!=static_cast<std::size_t>(-1);
                const std::size_t nb=internal?(f.owner==cell?f.neighbour:f.owner):cell;
                const double F=(f.owner==cell?phi_[fid]:-phi_[fid]);
                const double d=internal?dist(ce,m_.cells[nb]):std::max(1e-12, std::abs((f.cx-ce.cx)*f.nx+(f.cr-ce.cr)*f.nr));
                const double D=mu_*f.area/d;
                if(internal){
                    const double an=D+std::max(-F,0.0);
                    ax+=D+std::max(F,0.0);
                    ar+=D+std::max(F,0.0);
                    nxux[cell]+=an*ux[nb]; nxur[cell]+=an*ur[nb];
                }else if(f.axis){
                    // symmetry: ux has zero normal gradient, ur=0.
                    ax+=std::max(F,0.0);
                    ar+=D+std::max(F,0.0);
                }else if(f.outer && f.cx >= 0.0){
                    // Pressure outlet: zero normal velocity gradient.
                    ax+=std::max(F,0.0);
                    ar+=std::max(F,0.0);
                }else{
                    const double ubx=boundary_ux(f,cell), ubr=boundary_ur(f,cell);
                    const double ab=D+std::max(-F,0.0);
                    ax+=D+std::max(F,0.0); ar+=D+std::max(F,0.0);
                    bx+=ab*ubx; br+=ab*ubr;
                }
            }
            // Cylindrical radial viscous term -nu*u_r/r^2, implicit.
            if(!std::isfinite(ce.cr) || !(ce.cr>0.0))
                throw std::runtime_error("VMFL036 invalid positive cell radius at cell "+std::to_string(cell));
            const double rr=ce.cr;
            const double hoop=mu_*m_.cells[cell].volume/(rr*rr);
            if(!std::isfinite(hoop) || !(hoop>0.0))
                throw std::runtime_error("VMFL036 invalid radial geometric coefficient at cell "+std::to_string(cell));
            ar+=hoop;
            const double oldx=ux[cell], oldr=ur[cell];
            if(!std::isfinite(ax) || !std::isfinite(ar) || !(ax>0.0) || !(ar>0.0))
                throw std::runtime_error("VMFL036 invalid momentum coefficients at cell "+
                    std::to_string(cell)+" ax="+std::to_string(ax)+" ar="+std::to_string(ar));
            apx[cell]=ax/c.alpha_u; apr[cell]=ar/c.alpha_u;
            bx+=(1.0-c.alpha_u)/c.alpha_u*ax*oldx;
            br+=(1.0-c.alpha_u)/c.alpha_u*ar*oldr;
            // In-place Gauss-Seidel sweeps, diagonal uses relaxed aP.
            for(std::size_t s=0;s<c.momentum_sweeps;++s){
                const double nx=nxux[cell], nr=nxur[cell];
                const double sx=(nx+bx)/apx[cell], sr=(nr+br)/apr[cell];
                if(!std::isfinite(sx) || !std::isfinite(sr))
                    throw std::runtime_error(
                        "VMFL036 non-finite momentum sweep at cell "+std::to_string(cell)+
                        " sweep="+std::to_string(s)+
                        " ux="+std::to_string(sx)+" ur="+std::to_string(sr)+
                        " apx="+std::to_string(apx[cell])+" apr="+std::to_string(apr[cell]));
                ux[cell]=sx; ur[cell]=sr;
                // Remove the old neighbour contribution and refresh on the next sweep.
                // A full neighbour refresh is inexpensive for this deterministic mesh.
                double sumx=0.0,sumr=0.0;
                for(const auto fid:ce.faces){
                    const auto& ff=m_.faces[fid];
                    if(ff.neighbour!=static_cast<std::size_t>(-1)){
                        const auto nb2=(ff.owner==cell?ff.neighbour:ff.owner);
                        const double FF=(ff.owner==cell?phi_[fid]:-phi_[fid]);
                        const double dd=dist(ce,m_.cells[nb2]);
                        const double DD=mu_*ff.area/dd;
                        const double aa=DD+std::max(-FF,0.0);
                        sumx+=aa*ux[nb2]; sumr+=aa*ur[nb2];
                    }
                }
                const double ux_new=(sumx+bx)/apx[cell];
                const double ur_new=(sumr+br)/apr[cell];
                if(!std::isfinite(ux_new) || !std::isfinite(ur_new))
                    throw std::runtime_error(
                        "VMFL036 non-finite momentum update at cell "+std::to_string(cell)+
                        " ux="+std::to_string(ux_new)+" ur="+std::to_string(ur_new)+
                        " apx="+std::to_string(apx[cell])+" apr="+std::to_string(apr[cell])+
                        " bx="+std::to_string(bx)+" br="+std::to_string(br)+
                        " r="+std::to_string(ce.cr));
                ux[cell]=ux_new; ur[cell]=ur_new;
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

    double pressure_correction(const AxisymControls& c){
        pc_.assign(n_,0.0);
        std::vector<double> diag(n_,0.0),rhs(n_,0.0);
        std::vector<std::vector<std::pair<std::size_t,double>>> rows(n_);
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
                    const double df=0.5*(m_.cells[cell].volume/apx_[cell]+m_.cells[nb].volume/apx_[nb]);
                    const double a=rho_*f.area*df/d;
                    diag[cell]+=a; rows[cell].push_back({nb,-a});
                } else if(f.outer && f.cx >= 0.0){
                    const double d=std::max(1e-12,std::abs((f.cx-m_.cells[cell].cx)*f.nx+(f.cr-m_.cells[cell].cr)*f.nr));
                    const double a=rho_*f.area*(m_.cells[cell].volume/apx_[cell])/d;
                    diag[cell]+=a;
                    rhs[cell] += a*(-p_[cell]);
                }
            }
        }
        // The pressure-outlet arc supplies the Dirichlet pressure-correction
        // condition; no artificial reference cell is required.
        double maxcorr=std::numeric_limits<double>::infinity();
        for(std::size_t s=0;s<c.pressure_sweeps;++s){
            maxcorr=0.0;
            for(std::size_t cell=0;cell<n_;++cell){
                double sum=rhs[cell];
                for(const auto [nb,a]:rows[cell]) sum-=a*pc_[nb];
                const double v=sum/diag[cell];
                maxcorr=std::max(maxcorr,std::abs(v-pc_[cell]));
                pc_[cell]=v;
            }
            if(maxcorr<c.pressure_tolerance) break;
        }
        for(std::size_t cell=0;cell<n_;++cell){
            if(!std::isfinite(diag[cell]) || !(diag[cell]>0.0))
                throw std::runtime_error("VMFL036 invalid pressure-correction diagonal at cell "+
                    std::to_string(cell)+" value="+std::to_string(diag[cell]));
            if(!std::isfinite(rhs[cell]))
                throw std::runtime_error("VMFL036 non-finite pressure-correction RHS at cell "+
                    std::to_string(cell)+" value="+std::to_string(rhs[cell]));
        }
        // pressure correction gradient; boundary p'=0.
        for(std::size_t cell=0;cell<n_;++cell){
            double gx=0.0,gr=0.0;
            for(const auto fid:m_.cells[cell].faces){
                const auto& f=m_.faces[fid];
                double qf=pc_[cell];
                if(f.neighbour!=static_cast<std::size_t>(-1))
                    qf=0.5*(pc_[cell]+pc_[f.neighbour]);
                else if(f.outer && f.cx >= 0.0)
                    qf=-p_[cell]; // physical outlet pressure is fixed to zero
                gx+=qf*f.nx*f.area; gr+=qf*f.nr*f.area;
            }
            gx/=m_.cells[cell].volume; gr/=m_.cells[cell].volume;
            if(!std::isfinite(apx_[cell]) || !(apx_[cell]>0.0) ||
               !std::isfinite(apr_[cell]) || !(apr_[cell]>0.0))
                throw std::runtime_error("VMFL036 invalid momentum diagonal during pressure correction at cell "+
                    std::to_string(cell));
            ux_[cell]-=c.alpha_u*m_.cells[cell].volume/apx_[cell]*gx;
            ur_[cell]-=c.alpha_u*m_.cells[cell].volume/apr_[cell]*gr;
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
        compute_fluxes();
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
                if(!f.sphere) continue;
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

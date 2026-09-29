#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <queue>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace cfdx::physics {

struct WallDistanceVec3 {
    double x{0.0}, y{0.0}, z{0.0};
    WallDistanceVec3 operator+(const WallDistanceVec3& b) const { return {x+b.x,y+b.y,z+b.z}; }
    WallDistanceVec3 operator-(const WallDistanceVec3& b) const { return {x-b.x,y-b.y,z-b.z}; }
    WallDistanceVec3 operator*(double s) const { return {x*s,y*s,z*s}; }
};

inline double wd_dot(const WallDistanceVec3& a,const WallDistanceVec3& b) {
    return a.x*b.x+a.y*b.y+a.z*b.z;
}
inline double wd_norm2(const WallDistanceVec3& a) { return wd_dot(a,a); }
inline double wd_norm(const WallDistanceVec3& a) { return std::sqrt(wd_norm2(a)); }

struct WallTriangle { std::array<std::size_t,3> v{}; };

struct WallSurface {
    std::vector<WallDistanceVec3> points;
    std::vector<WallTriangle> triangles;
};

enum class WallDistanceConvergenceStatus {
    DIRECT,
    CONVERGED,
    MAX_ITER,
    RESIDUAL_TOO_HIGH,
    RESIDUAL_STAGNATION,
    NOT_REPORTED
};

struct WallDistanceResult {
    std::vector<double> distance;
    std::vector<unsigned char> valid;
    std::string method;
    std::size_t iterations{0};
    double residual_inf{0.0};
    std::size_t auxiliary_iterations{0};
    double auxiliary_residual_inf{0.0};
    std::size_t poisson_wall_ray_hits{0};
    std::size_t poisson_wall_ray_misses{0};
    std::size_t poisson_wall_fallbacks{0};
    std::size_t poisson_wall_bad_alignment{0};
    double poisson_wall_min_alignment{1.0};
    double poisson_phi_min{0.0};
    double poisson_phi_max{0.0};
    double poisson_grad_min{0.0};
    double poisson_grad_max{0.0};
    double poisson_distance_l2_error{0.0};
    double poisson_distance_linf_error{0.0};
    bool converged{false};
    std::string stopping_reason{"not_reported"};
    WallDistanceConvergenceStatus convergence_status{WallDistanceConvergenceStatus::NOT_REPORTED};
};

enum class WallDistanceMethod {
    EXACT_GEOMETRIC,
    SEARCH_BASED,
    MESH_WAVE,
    DIRECTIONAL_MESH_WAVE,
    POISSON,
    EIKONAL,
    HAMILTON_JACOBI,
    ADVECTION_DIFFUSION,
    HYBRID_POISSON_EIKONAL
};

inline const char* wall_distance_method_name(WallDistanceMethod m) {
    switch(m) {
        case WallDistanceMethod::EXACT_GEOMETRIC: return "exact_geometric";
        case WallDistanceMethod::SEARCH_BASED: return "search_based";
        case WallDistanceMethod::MESH_WAVE: return "mesh_wave";
        case WallDistanceMethod::DIRECTIONAL_MESH_WAVE: return "directional_mesh_wave";
        case WallDistanceMethod::POISSON: return "poisson";
        case WallDistanceMethod::EIKONAL: return "eikonal";
        case WallDistanceMethod::HAMILTON_JACOBI: return "hamilton_jacobi";
        case WallDistanceMethod::ADVECTION_DIFFUSION: return "advection_diffusion";
        case WallDistanceMethod::HYBRID_POISSON_EIKONAL: return "hybrid_poisson_eikonal";
    }
    return "unknown";
}

inline double point_triangle_distance2(const WallDistanceVec3& p,
                                       const WallDistanceVec3& a,
                                       const WallDistanceVec3& b,
                                       const WallDistanceVec3& c) {
    const WallDistanceVec3 ab=b-a, ac=c-a, ap=p-a;
    const double d1=wd_dot(ab,ap), d2=wd_dot(ac,ap);
    if(d1<=0.0 && d2<=0.0) return wd_norm2(ap);
    const WallDistanceVec3 bp=p-b;
    const double d3=wd_dot(ab,bp), d4=wd_dot(ac,bp);
    if(d3>=0.0 && d4<=d3) return wd_norm2(bp);
    const double vc=d1*d4-d3*d2;
    if(vc<=0.0 && d1>=0.0 && d3<=0.0) {
        const double v=d1/(d1-d3);
        return wd_norm2(p-(a+ab*v));
    }
    const WallDistanceVec3 cp=p-c;
    const double d5=wd_dot(ab,cp), d6=wd_dot(ac,cp);
    if(d6>=0.0 && d5<=d6) return wd_norm2(cp);
    const double vb=d5*d2-d1*d6;
    if(vb<=0.0 && d2>=0.0 && d6<=0.0) {
        const double w=d2/(d2-d6);
        return wd_norm2(p-(a+ac*w));
    }
    const double va=d3*d6-d5*d4;
    if(va<=0.0 && (d4-d3)>=0.0 && (d5-d6)>=0.0) {
        const double w=(d4-d3)/((d4-d3)+(d5-d6));
        return wd_norm2(p-(b+(c-b)*w));
    }
    const double denom=1.0/(va+vb+vc);
    const double v=vb*denom, w=vc*denom;
    return wd_norm2(p-(a+ab*v+ac*w));
}

struct WallDistanceAabb {
    WallDistanceVec3 lo{std::numeric_limits<double>::infinity(),
                         std::numeric_limits<double>::infinity(),
                         std::numeric_limits<double>::infinity()};
    WallDistanceVec3 hi{-std::numeric_limits<double>::infinity(),
                         -std::numeric_limits<double>::infinity(),
                         -std::numeric_limits<double>::infinity()};
};

inline WallDistanceAabb wall_distance_expand(WallDistanceAabb box,
                                              const WallDistanceVec3& p) {
    box.lo.x=std::min(box.lo.x,p.x); box.lo.y=std::min(box.lo.y,p.y); box.lo.z=std::min(box.lo.z,p.z);
    box.hi.x=std::max(box.hi.x,p.x); box.hi.y=std::max(box.hi.y,p.y); box.hi.z=std::max(box.hi.z,p.z);
    return box;
}

inline double wall_distance_aabb_distance2(const WallDistanceAabb& box,
                                            const WallDistanceVec3& p) {
    const auto axis=[](double x,double lo,double hi) {
        return x<lo ? lo-x : (x>hi ? x-hi : 0.0);
    };
    const double dx=axis(p.x,box.lo.x,box.hi.x);
    const double dy=axis(p.y,box.lo.y,box.hi.y);
    const double dz=axis(p.z,box.lo.z,box.hi.z);
    return dx*dx+dy*dy+dz*dz;
}

class WallDistanceBvh {
public:
    explicit WallDistanceBvh(const WallSurface& surface) : surface_(surface) {
        if(surface_.points.empty() || surface_.triangles.empty())
            throw std::invalid_argument("wall distance BVH surface is empty");
        indices_.resize(surface_.triangles.size());
        for(std::size_t i=0;i<indices_.size();++i) indices_[i]=i;
        nodes_.reserve(indices_.size()*2);
        build(0,indices_.size());
    }

    double nearest_distance(const WallDistanceVec3& p) const {
        double best2=std::numeric_limits<double>::infinity();
        nearest(0,p,best2);
        return std::sqrt(best2);
    }

    double nearest_distance_within(const WallDistanceVec3& p,double threshold) const {
        if(threshold<0.0) return std::numeric_limits<double>::infinity();
        double best2=threshold*threshold;
        bool found=false;
        nearest_within(0,p,best2,found);
        return found ? std::sqrt(best2) : std::numeric_limits<double>::infinity();
    }

    // Distance from p to the first wall intersection along a grid-face ray.
    // This is the correct geometric quantity for a cut Dirichlet face:
    // the boundary lies on the segment joining the fluid node to the solid
    // neighbour.  nearest_distance()/nearest_normal() are not sufficient
    // near edges/corners because the closest surface point need not lie on
    // that face-normal ray.
    double ray_distance(const WallDistanceVec3& p,
                        const WallDistanceVec3& direction,
                        double max_distance) const {
        const double dn=wd_norm(direction);
        if(!(dn>1e-30) || !(max_distance>0.0)) return std::numeric_limits<double>::infinity();
        const WallDistanceVec3 d=direction*(1.0/dn);
        const auto cross=[](const WallDistanceVec3& a,const WallDistanceVec3& b) {
            return WallDistanceVec3{
                a.y*b.z-a.z*b.y,
                a.z*b.x-a.x*b.z,
                a.x*b.y-a.y*b.x
            };
        };
        double best=std::numeric_limits<double>::infinity();
        const double eps=1e-12*std::max(1.0,max_distance);

        for(const auto& tri:surface_.triangles) {
            const auto& v0=surface_.points[tri.v[0]];
            const auto& v1=surface_.points[tri.v[1]];
            const auto& v2=surface_.points[tri.v[2]];
            const auto e1=v1-v0;
            const auto e2=v2-v0;
            const auto h=cross(d,e2);
            const double a=wd_dot(e1,h);
            if(std::abs(a)<1e-14) continue;
            const double inv_a=1.0/a;
            const auto s=p-v0;
            const double u=inv_a*wd_dot(s,h);
            if(u < -1e-12 || u > 1.0+1e-12) continue;
            const auto q=cross(s,e1);
            const double v=inv_a*wd_dot(d,q);
            if(v < -1e-12 || u+v > 1.0+1e-12) continue;
            const double t=inv_a*wd_dot(e2,q);
            if(t>=-eps && t<=max_distance+eps) best=std::min(best,std::max(0.0,t));
        }
        return best;
    }

    std::size_t triangle_count() const { return surface_.triangles.size(); }
    std::size_t node_count() const { return nodes_.size(); }

    WallDistanceVec3 nearest_normal(const WallDistanceVec3& p) const {
        double best2=std::numeric_limits<double>::infinity();
        std::size_t best=0;
        nearest_triangle(0,p,best2,best);
        const auto& t=surface_.triangles[best];
        const auto ab=surface_.points[t.v[1]]-surface_.points[t.v[0]];
        const auto ac=surface_.points[t.v[2]]-surface_.points[t.v[0]];
        const WallDistanceVec3 n{ab.y*ac.z-ab.z*ac.y,
                                 ab.z*ac.x-ab.x*ac.z,
                                 ab.x*ac.y-ab.y*ac.x};
        const double nn=wd_norm(n);
        return nn>1e-30 ? n*(1.0/nn) : WallDistanceVec3{1.0,0.0,0.0};
    }

private:
    struct Node {
        WallDistanceAabb box{};
        std::size_t begin{0}, end{0};
        std::size_t left{std::numeric_limits<std::size_t>::max()};
        std::size_t right{std::numeric_limits<std::size_t>::max()};
        bool leaf{false};
    };

    const WallSurface& surface_;
    std::vector<std::size_t> indices_;
    std::vector<Node> nodes_;
    static constexpr std::size_t leaf_size_=4;

    WallDistanceAabb triangle_box(std::size_t ti) const {
        const auto& t=surface_.triangles[ti];
        WallDistanceAabb box;
        box=wall_distance_expand(box,surface_.points[t.v[0]]);
        box=wall_distance_expand(box,surface_.points[t.v[1]]);
        box=wall_distance_expand(box,surface_.points[t.v[2]]);
        return box;
    }

    WallDistanceVec3 triangle_centroid(std::size_t ti) const {
        const auto& t=surface_.triangles[ti];
        return (surface_.points[t.v[0]]+surface_.points[t.v[1]]+surface_.points[t.v[2]])*(1.0/3.0);
    }

    std::size_t build(std::size_t begin,std::size_t end) {
        const std::size_t node_id=nodes_.size();
        nodes_.push_back({});
        Node& node=nodes_.back();
        node.begin=begin; node.end=end;
        for(std::size_t i=begin;i<end;++i) {
            const auto box=triangle_box(indices_[i]);
            node.box=wall_distance_expand(node.box,box.lo);
            node.box=wall_distance_expand(node.box,box.hi);
        }
        const std::size_t count=end-begin;
        if(count<=leaf_size_) {
            node.leaf=true;
            return node_id;
        }
        const auto extent=node.box.hi-node.box.lo;
        int axis=0;
        if(extent.y>extent.x && extent.y>=extent.z) axis=1;
        else if(extent.z>extent.x && extent.z>extent.y) axis=2;
        const std::size_t mid=begin+count/2;
        auto coord=[axis](const WallDistanceVec3& p) {
            return axis==0?p.x:(axis==1?p.y:p.z);
        };
        std::nth_element(indices_.begin()+static_cast<std::ptrdiff_t>(begin),
                         indices_.begin()+static_cast<std::ptrdiff_t>(mid),
                         indices_.begin()+static_cast<std::ptrdiff_t>(end),
                         [&](std::size_t a,std::size_t b) {
                             return coord(triangle_centroid(a))<coord(triangle_centroid(b));
                         });
        const std::size_t left=build(begin,mid);
        const std::size_t right=build(mid,end);
        nodes_[node_id].left=left;
        nodes_[node_id].right=right;
        return node_id;
    }

    void nearest_triangle(std::size_t node_id,const WallDistanceVec3& p,
                          double& best2,std::size_t& best) const {
        const Node& node=nodes_[node_id];
        const double box_d2=wall_distance_aabb_distance2(node.box,p);
        if(box_d2>=best2) return;
        if(node.leaf) {
            for(std::size_t i=node.begin;i<node.end;++i) {
                const std::size_t ti=indices_[i];
                const auto& t=surface_.triangles[ti];
                const double d2=point_triangle_distance2(
                    p,surface_.points[t.v[0]],surface_.points[t.v[1]],surface_.points[t.v[2]]);
                if(d2<best2) { best2=d2; best=ti; }
            }
            return;
        }
        const double dl=wall_distance_aabb_distance2(nodes_[node.left].box,p);
        const double dr=wall_distance_aabb_distance2(nodes_[node.right].box,p);
        if(dl<dr) { nearest_triangle(node.left,p,best2,best); nearest_triangle(node.right,p,best2,best); }
        else { nearest_triangle(node.right,p,best2,best); nearest_triangle(node.left,p,best2,best); }
    }

    void nearest_within(std::size_t node_id,const WallDistanceVec3& p,
                         double& best2,bool& found) const {
        const Node& node=nodes_[node_id];
        const double box_d2=wall_distance_aabb_distance2(node.box,p);
        if(box_d2>best2) return;
        if(node.leaf) {
            for(std::size_t i=node.begin;i<node.end;++i) {
                const auto& t=surface_.triangles[indices_[i]];
                const double d2=point_triangle_distance2(
                    p,surface_.points[t.v[0]],surface_.points[t.v[1]],surface_.points[t.v[2]]);
                if(d2<=best2) {
                    best2=d2;
                    found=true;
                }
            }
            return;
        }
        const Node& left=nodes_[node.left];
        const Node& right=nodes_[node.right];
        const double dl=wall_distance_aabb_distance2(left.box,p);
        const double dr=wall_distance_aabb_distance2(right.box,p);
        if(dl<dr) {
            nearest_within(node.left,p,best2,found);
            nearest_within(node.right,p,best2,found);
        } else {
            nearest_within(node.right,p,best2,found);
            nearest_within(node.left,p,best2,found);
        }
    }

    void nearest(std::size_t node_id,const WallDistanceVec3& p,double& best2) const {
        const Node& node=nodes_[node_id];
        if(wall_distance_aabb_distance2(node.box,p)>=best2) return;
        if(node.leaf) {
            for(std::size_t i=node.begin;i<node.end;++i) {
                const auto& t=surface_.triangles[indices_[i]];
                best2=std::min(best2,point_triangle_distance2(
                    p,surface_.points[t.v[0]],surface_.points[t.v[1]],surface_.points[t.v[2]]));
            }
            return;
        }
        const Node& left=nodes_[node.left];
        const Node& right=nodes_[node.right];
        const double dl=wall_distance_aabb_distance2(left.box,p);
        const double dr=wall_distance_aabb_distance2(right.box,p);
        if(dl<dr) {
            nearest(node.left,p,best2);
            nearest(node.right,p,best2);
        } else {
            nearest(node.right,p,best2);
            nearest(node.left,p,best2);
        }
    }
};

inline double exact_point_distance(const WallSurface& s,const WallDistanceVec3& p) {
    if(s.points.empty() || s.triangles.empty())
        throw std::invalid_argument("wall distance surface is empty");
    double d2=std::numeric_limits<double>::infinity();
    for(const auto& t:s.triangles) {
        d2=std::min(d2,point_triangle_distance2(p,s.points[t.v[0]],s.points[t.v[1]],s.points[t.v[2]]));
    }
    return std::sqrt(d2);
}

struct WallDistanceGrid {
    std::size_t nx{0}, ny{0}, nz{0};
    WallDistanceVec3 origin{};
    WallDistanceVec3 spacing{1.0,1.0,1.0};
    std::vector<WallDistanceVec3> points;
    std::vector<unsigned char> solid;
    std::vector<std::vector<std::size_t>> neighbours;

    std::size_t index(std::size_t i,std::size_t j,std::size_t k) const {
        return (k*ny+j)*nx+i;
    }
};

inline WallDistanceGrid make_wall_distance_grid(std::size_t nx,std::size_t ny,std::size_t nz,
                                                 const WallDistanceVec3& origin,
                                                 const WallDistanceVec3& spacing,
                                                 const std::function<bool(const WallDistanceVec3&)>& inside) {
    if(nx<3 || ny<3 || nz<3) throw std::invalid_argument("wall distance grid too small");
    WallDistanceGrid g; g.nx=nx; g.ny=ny; g.nz=nz; g.origin=origin; g.spacing=spacing;
    const std::size_t n=nx*ny*nz;
    g.points.resize(n); g.solid.resize(n,0); g.neighbours.resize(n);
    for(std::size_t k=0;k<nz;++k) for(std::size_t j=0;j<ny;++j) for(std::size_t i=0;i<nx;++i) {
        const std::size_t id=g.index(i,j,k);
        g.points[id]={origin.x+spacing.x*static_cast<double>(i),
                      origin.y+spacing.y*static_cast<double>(j),
                      origin.z+spacing.z*static_cast<double>(k)};
        g.solid[id]=inside(g.points[id])?1:0;
        constexpr std::array<std::array<int,3>,6> dirs={{{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1}}};
        for(const auto& d:dirs) {
            const long ii=static_cast<long>(i)+d[0], jj=static_cast<long>(j)+d[1], kk=static_cast<long>(k)+d[2];
            if(ii>=0 && jj>=0 && kk>=0 && ii<static_cast<long>(nx) && jj<static_cast<long>(ny) && kk<static_cast<long>(nz))
                g.neighbours[id].push_back(g.index(static_cast<std::size_t>(ii),static_cast<std::size_t>(jj),static_cast<std::size_t>(kk)));
        }
    }
    return g;
}

inline std::vector<double> exact_reference(const WallSurface& s,const WallDistanceGrid& g) {
    std::vector<double> d(g.points.size(),std::numeric_limits<double>::infinity());
    for(std::size_t i=0;i<g.points.size();++i) if(!g.solid[i]) d[i]=exact_point_distance(s,g.points[i]);
    return d;
}

inline std::vector<double> search_based_reference(const WallSurface& s,const WallDistanceGrid& g,double threshold) {
    // Exact surface search accelerated by a BVH. The threshold argument is
    // retained for API compatibility, but it must never switch the method to
    // a vertex-distance approximation: the wall-distance contract is exact
    // point-to-triangle distance for every fluid query point.
    (void)threshold;
    const WallDistanceBvh bvh(s);
    std::vector<double> d(g.points.size(),std::numeric_limits<double>::infinity());
    for(std::size_t i=0;i<g.points.size();++i)
        if(!g.solid[i]) d[i]=bvh.nearest_distance(g.points[i]);
    return d;
}

inline std::vector<std::size_t> wall_seed_nodes(const WallDistanceBvh& bvh,const WallDistanceGrid& g,double seed_width) {
    std::vector<std::size_t> seeds;
    for(std::size_t i=0;i<g.points.size();++i)
        if(!g.solid[i] && std::isfinite(bvh.nearest_distance_within(g.points[i],seed_width)))
            seeds.push_back(i);
    return seeds;
}

inline std::vector<std::size_t> wall_seed_nodes(const WallSurface& s,const WallDistanceGrid& g,double seed_width) {
    const WallDistanceBvh bvh(s);
    return wall_seed_nodes(bvh,g,seed_width);
}

inline WallDistanceResult graph_wave(const WallSurface& s,const WallDistanceGrid& g,bool directional) {
    const double h=std::min({g.spacing.x,g.spacing.y,g.spacing.z});
    const WallDistanceBvh bvh(s);
    const auto seeds=wall_seed_nodes(bvh,g,1.6*h);
    std::vector<double> d(g.points.size(),std::numeric_limits<double>::infinity());
    using Item=std::pair<double,std::size_t>;
    std::priority_queue<Item,std::vector<Item>,std::greater<Item>> q;
    for(const auto id:seeds) { d[id]=bvh.nearest_distance(g.points[id]); q.push({d[id],id}); }
    while(!q.empty()) {
        const auto [du,u]=q.top(); q.pop();
        if(du!=d[u]) continue;
        for(const auto v:g.neighbours[u]) {
            if(g.solid[v]) continue;
            const auto dv=g.points[v]-g.points[u];
            const double r=wd_norm(dv);
            double w=r;
            if(directional) {
                // Penalize propagation tangential to the local wall normal.
                // This is a real directional metric; the previous implementation
                // used dot(dv,-dv), which is identically -|dv|^2.
                const auto n=bvh.nearest_normal(g.points[u]);
                const double align=std::abs(wd_dot(dv,n))/std::max(r,1e-30);
                w*=1.0+0.15*(1.0-align);
            }
            if(du+w<d[v]) { d[v]=du+w; q.push({d[v],v}); }
        }
    }
    return {std::move(d),std::vector<unsigned char>(g.points.size(),1),
            directional?"directional_mesh_wave":"mesh_wave",0};
}

inline double poisson_wall_offset(const WallDistanceBvh& bvh,
                                  const WallDistanceVec3& p,
                                  const WallDistanceVec3& q,
                                  double h);

inline double eikonal_update(const std::array<double,3>& a,
                                  const std::array<double,3>& h,
                                  double rhs) {
    if(rhs<=0.0) return 0.0;
    std::array<std::pair<double,double>,3> v{};
    std::size_t n=0;
    for(std::size_t q=0;q<3;++q) if(std::isfinite(a[q]) && h[q]>0.0) v[n++]={a[q],h[q]};
    if(n==0) return std::numeric_limits<double>::infinity();
    // n<=3 by construction.  Use a bounded insertion sort instead of
    // std::sort on a fixed-size array; GCC otherwise emits a spurious
    // -Warray-bounds warning after aggressive inlining of the 3-element case.
    for(std::size_t q=1;q<n;++q) {
        auto key=v[q];
        std::size_t j=q;
        while(j>0 && key.first<v[j-1].first) {
            v[j]=v[j-1];
            --j;
        }
        v[j]=key;
    }
    double A=0.0,B=0.0,C=-rhs*rhs;
    double x=v[0].first+v[0].second*rhs;
    for(std::size_t m=0;m<n;++m) {
        const double ai=v[m].first, hi=v[m].second;
        A+=1.0/(hi*hi); B+=-2.0*ai/(hi*hi); C+=ai*ai/(hi*hi);
        const double disc=std::max(0.0,B*B-4.0*A*C);
        x=(-B+std::sqrt(disc))/(2.0*A);
        if(m+1==n || x<=v[m+1].first) return x;
    }
    return x;
}

inline std::vector<double> eikonal_fast_sweep(const WallDistanceBvh& bvh,const WallDistanceGrid& g,
                                               std::size_t max_iter,std::size_t* used_iter=nullptr) {
    const double h=std::min({g.spacing.x,g.spacing.y,g.spacing.z});
    const auto seeds=wall_seed_nodes(bvh,g,1.6*h);
    std::vector<unsigned char> fixed(g.points.size(),0);
    std::vector<double> d(g.points.size(),std::numeric_limits<double>::infinity());
    for(auto id:seeds) { d[id]=bvh.nearest_distance(g.points[id]); fixed[id]=1; }
    const std::array<int,2> signs={-1,1};
    std::size_t it_used=max_iter;
    for(std::size_t it=0;it<max_iter;++it) {
        double max_change=0.0;
        for(int sx:signs) for(int sy:signs) for(int sz:signs)
            for(std::size_t kk=0;kk<g.nz;++kk) for(std::size_t jj=0;jj<g.ny;++jj) for(std::size_t ii=0;ii<g.nx;++ii) {
                const std::size_t i=sx>0?ii:g.nx-1-ii, j=sy>0?jj:g.ny-1-jj, k=sz>0?kk:g.nz-1-kk;
                const std::size_t id=g.index(i,j,k);
                if(g.solid[id] || fixed[id]) continue;
                std::array<double,3> a={std::numeric_limits<double>::infinity(),
                                        std::numeric_limits<double>::infinity(),
                                        std::numeric_limits<double>::infinity()};
                if(i>0) a[0]=std::min(a[0],d[g.index(i-1,j,k)]);
                if(i+1<g.nx) a[0]=std::min(a[0],d[g.index(i+1,j,k)]);
                if(j>0) a[1]=std::min(a[1],d[g.index(i,j-1,k)]);
                if(j+1<g.ny) a[1]=std::min(a[1],d[g.index(i,j+1,k)]);
                if(k>0) a[2]=std::min(a[2],d[g.index(i,j,k-1)]);
                if(k+1<g.nz) a[2]=std::min(a[2],d[g.index(i,j,k+1)]);
                const double cand=eikonal_update(a,{g.spacing.x,g.spacing.y,g.spacing.z},1.0);
                if(std::isfinite(cand) && (!std::isfinite(d[id]) || cand<d[id])) {
                    const double old=d[id]; d[id]=cand;
                    if(std::isfinite(old)) max_change=std::max(max_change,std::abs(cand-old));
                    else max_change=std::max(max_change,cand);
                }
            }
        if(max_change<1e-10*h) { it_used=it+1; break; }
    }
    if(used_iter) *used_iter=it_used;
    return d;
}

inline std::vector<double> eikonal_fast_sweep(const WallSurface& s,const WallDistanceGrid& g,
                                               std::size_t max_iter,std::size_t* used_iter=nullptr) {
    const WallDistanceBvh bvh(s);
    return eikonal_fast_sweep(bvh,g,max_iter,used_iter);
}

inline double laplacian_at(const std::vector<double>& f,const WallDistanceGrid& g,
                           std::size_t id) {
    // Finite-volume-compatible seven-point operator on the active fluid graph.
    //
    // Every fluid-fluid or fluid-solid face contributes exactly once with
    // weight 1/h^2.  A missing outer neighbour is a zero-flux Neumann face and
    // therefore contributes nothing.  This is deliberately different from
    // the old mirrored-ghost stencil (2*(f_nb-f_i)): that stencil makes the
    // matrix non-symmetric at an outer boundary (boundary row -2/h^2 versus
    // interior row -1/h^2), so it is not a valid Euclidean-inner-product PCG
    // operator.
    const std::size_t k=id/(g.nx*g.ny), rem=id%(g.nx*g.ny), j=rem/g.nx, i=rem%g.nx;
    double l=0.0;
    auto add_face=[&](std::size_t q,bool exists,double h) {
        if(!exists) return; // homogeneous Neumann: zero flux
        const double w=1.0/(h*h);
        const double fq=g.solid[q] ? 0.0 : f[q];
        if(std::isfinite(fq)) l+=(fq-f[id])*w;
    };
    const std::size_t xm=i>0?g.index(i-1,j,k):0, xp=i+1<g.nx?g.index(i+1,j,k):0;
    const std::size_t ym=j>0?g.index(i,j-1,k):0, yp=j+1<g.ny?g.index(i,j+1,k):0;
    const std::size_t zm=k>0?g.index(i,j,k-1):0, zp=k+1<g.nz?g.index(i,j,k+1):0;
    add_face(xm,i>0,g.spacing.x); add_face(xp,i+1<g.nx,g.spacing.x);
    add_face(ym,j>0,g.spacing.y); add_face(yp,j+1<g.ny,g.spacing.y);
    add_face(zm,k>0,g.spacing.z); add_face(zp,k+1<g.nz,g.spacing.z);
    return l;
}

inline double laplacian_at(const WallDistanceBvh& bvh,
                           const std::vector<double>& f,
                           const WallDistanceGrid& g,
                           std::size_t id) {
    const std::size_t k=id/(g.nx*g.ny), rem=id%(g.nx*g.ny), j=rem/g.nx, i=rem%g.nx;
    double l=0.0;
    auto add_face=[&](std::size_t q,bool exists,double hh) {
        if(!exists) return;
        if(g.solid[q]) {
            const double delta=poisson_wall_offset(bvh,g.points[id],g.points[q],hh);
            if(delta>0.0 && std::isfinite(delta)) l-=f[id]/(hh*delta);
        } else if(std::isfinite(f[q])) {
            l+=(f[q]-f[id])/(hh*hh);
        }
    };
    const std::size_t xm=i>0?g.index(i-1,j,k):0, xp=i+1<g.nx?g.index(i+1,j,k):0;
    const std::size_t ym=j>0?g.index(i,j-1,k):0, yp=j+1<g.ny?g.index(i,j+1,k):0;
    const std::size_t zm=k>0?g.index(i,j,k-1):0, zp=k+1<g.nz?g.index(i,j,k+1):0;
    add_face(xm,i>0,g.spacing.x); add_face(xp,i+1<g.nx,g.spacing.x);
    add_face(ym,j>0,g.spacing.y); add_face(yp,j+1<g.ny,g.spacing.y);
    add_face(zm,k>0,g.spacing.z); add_face(zp,k+1<g.nz,g.spacing.z);
    return l;
}

inline double godunov_gradient_at(const WallDistanceBvh& bvh,
                                  const std::vector<double>& d,const WallDistanceGrid& g,
                                  std::size_t id) {
    const std::size_t k=id/(g.nx*g.ny), rem=id%(g.nx*g.ny), j=rem/g.nx, i=rem%g.nx;
    auto component=[&](int axis)->double {
        double dm=0.0, dp=0.0;
        auto face=[&](std::size_t q,bool minus_side,double hh) {
            if(g.solid[q]) {
                const double delta=poisson_wall_offset(bvh,g.points[id],g.points[q],hh);
                if(delta>0.0 && std::isfinite(delta)) {
                    if(minus_side) dm=d[id]/delta;
                    else dp=-d[id]/delta;
                }
            } else if(std::isfinite(d[q])) {
                if(minus_side) dm=(d[id]-d[q])/hh;
                else dp=(d[q]-d[id])/hh;
            }
        };
        if(axis==0) {
            if(i>0) face(g.index(i-1,j,k),true,g.spacing.x);
            if(i+1<g.nx) face(g.index(i+1,j,k),false,g.spacing.x);
        } else if(axis==1) {
            if(j>0) face(g.index(i,j-1,k),true,g.spacing.y);
            if(j+1<g.ny) face(g.index(i,j+1,k),false,g.spacing.y);
        } else {
            if(k>0) face(g.index(i,j,k-1),true,g.spacing.z);
            if(k+1<g.nz) face(g.index(i,j,k+1),false,g.spacing.z);
        }
        return std::sqrt(std::max(dm,0.0)*std::max(dm,0.0)
                       + std::min(dp,0.0)*std::min(dp,0.0));
    };
    const double gx=component(0), gy=component(1), gz=component(2);
    return std::sqrt(gx*gx+gy*gy+gz*gz);
}

inline double godunov_gradient_at(const std::vector<double>& d,const WallDistanceGrid& g,
                                  std::size_t id) {
    const std::size_t k=id/(g.nx*g.ny), rem=id%(g.nx*g.ny), j=rem/g.nx, i=rem%g.nx;
    auto component=[&](int axis)->double {
        double dm=0.0, dp=0.0;
        if(axis==0) {
            if(i>0) {
                const auto m=g.index(i-1,j,k);
                if(g.solid[m]) dm=d[id]/g.spacing.x;
                else if(std::isfinite(d[m])) dm=(d[id]-d[m])/g.spacing.x;
            }
            if(i+1<g.nx) {
                const auto p=g.index(i+1,j,k);
                if(g.solid[p]) dp=-d[id]/g.spacing.x;
                else if(std::isfinite(d[p])) dp=(d[p]-d[id])/g.spacing.x;
            }
        } else if(axis==1) {
            if(j>0) {
                const auto m=g.index(i,j-1,k);
                if(g.solid[m]) dm=d[id]/g.spacing.y;
                else if(std::isfinite(d[m])) dm=(d[id]-d[m])/g.spacing.y;
            }
            if(j+1<g.ny) {
                const auto p=g.index(i,j+1,k);
                if(g.solid[p]) dp=-d[id]/g.spacing.y;
                else if(std::isfinite(d[p])) dp=(d[p]-d[id])/g.spacing.y;
            }
        } else {
            if(k>0) {
                const auto m=g.index(i,j,k-1);
                if(g.solid[m]) dm=d[id]/g.spacing.z;
                else if(std::isfinite(d[m])) dm=(d[id]-d[m])/g.spacing.z;
            }
            if(k+1<g.nz) {
                const auto p=g.index(i,j,k+1);
                if(g.solid[p]) dp=-d[id]/g.spacing.z;
                else if(std::isfinite(d[p])) dp=(d[p]-d[id])/g.spacing.z;
            }
        }
        return std::sqrt(std::max(dm,0.0)*std::max(dm,0.0)
                       + std::min(dp,0.0)*std::min(dp,0.0));
    };
    const double gx=component(0), gy=component(1), gz=component(2);
    return std::sqrt(gx*gx+gy*gy+gz*gz);
}

inline std::vector<double> hamilton_jacobi_distance(const WallSurface& s,const WallDistanceGrid& g,
                                                     std::size_t max_iter,double epsilon=0.25,
                                                     double relaxation=0.25,
                                                     std::size_t* used_iter=nullptr,
                                                     double* residual_out=nullptr) {
    // Modified Hamilton-Jacobi equation:
    //   |grad d| = 1 + Gamma(d) laplacian(d), Gamma(d)=epsilon*d.
    // Solve it by pseudo-time relaxation with an upwind/Godunov gradient
    // using a CFL-scale explicit step for the advection/diffusion update.
    // and a central Laplacian. The local pseudo-time step is diffusion-limited.
    const double h=std::min({g.spacing.x,g.spacing.y,g.spacing.z});
    const WallDistanceBvh bvh(s);
    const auto seeds=wall_seed_nodes(bvh,g,1.6*h);
    std::vector<unsigned char> fixed(g.points.size(),0);
    std::vector<double> d=eikonal_fast_sweep(bvh,g,std::max<std::size_t>(8,max_iter));
    for(auto id:seeds) { d[id]=bvh.nearest_distance(g.points[id]); fixed[id]=1; }

    const std::size_t steps=std::max<std::size_t>(20,max_iter);
    std::size_t used=steps;
    double final_residual=std::numeric_limits<double>::infinity();
    const double inv_h2=1.0/(g.spacing.x*g.spacing.x)
                      +1.0/(g.spacing.y*g.spacing.y)
                      +1.0/(g.spacing.z*g.spacing.z);
    for(std::size_t it=0;it<steps;++it) {
        double max_change=0.0;
        for(std::size_t id=0;id<d.size();++id) {
            if(g.solid[id] || fixed[id] || !std::isfinite(d[id])) continue;
            const double grad=godunov_gradient_at(bvh,d,g,id);
            const double gamma=epsilon*std::max(d[id],0.0);
            const double lap=laplacian_at(bvh,d,g,id);
            const double residual=grad-1.0-gamma*lap;
            // Use the full explicit CFL limits of the upwind advection and
            // central diffusion terms. The previous 0.25 safety factor was
            // unnecessarily restrictive and made H-J convergence ~4x slower.
            // 3-D CFL: dt*h^{-1} must account for all active
            // upwind directions.  0.9*h is only a 1-D bound and can violate
            // the monotone explicit stability limit when several components
            // of the Godunov gradient are active.
            const double dt_adv=0.9/(1.0/g.spacing.x+1.0/g.spacing.y+1.0/g.spacing.z);
            const double dt_diff=gamma>0.0 ? 0.9/(2.0*gamma*inv_h2)
                                           : std::numeric_limits<double>::infinity();
            const double dt=std::clamp(relaxation,0.1,1.0)*std::min(dt_adv,dt_diff);
            const double nd=std::max(0.0,d[id]-dt*residual);
            max_change=std::max(max_change,std::abs(nd-d[id]));
            d[id]=nd;
        }
        // Pure Poisson mode deliberately does not overwrite the PDE solution with BVH seeds.
        if((it&3u)==3u || max_change<1e-10*h) {
            final_residual=0.0;
            for(std::size_t id=0;id<d.size();++id) {
                if(g.solid[id] || fixed[id] || !std::isfinite(d[id])) continue;
                final_residual=std::max(final_residual,std::abs(
                    godunov_gradient_at(bvh,d,g,id)-1.0-
                    epsilon*std::max(d[id],0.0)*laplacian_at(bvh,d,g,id)));
            }
            if(final_residual<1e-8 && max_change<1e-9*h) { used=it+1; break; }
        }
    }
    if(used_iter) *used_iter=used;
    if(residual_out) *residual_out=final_residual;
    return d;
}

inline double poisson_laplacian_at(const WallDistanceBvh& bvh,
                                   const std::vector<double>& f,
                                   const WallDistanceGrid& g,
                                   std::size_t id);

inline double poisson_residual_inf(const WallDistanceBvh& bvh,
                                      const std::vector<double>& phi,
                                      const WallDistanceGrid& g,
                                      const std::vector<unsigned char>& fixed) {
    // Public diagnostic: use exactly the same cut-face operator as PCG.
    double rmax=0.0;
    for(std::size_t id=0;id<phi.size();++id) {
        if(g.solid[id] || fixed[id] || !std::isfinite(phi[id])) continue;
        rmax=std::max(rmax,std::abs(-poisson_laplacian_at(bvh,phi,g,id)-1.0));
    }
    return rmax;
}

inline double poisson_wall_offset(const WallDistanceBvh& bvh,
                                      const WallDistanceVec3& p,
                                      const WallDistanceVec3& q,
                                      double h);
inline double poisson_laplacian_at(const WallDistanceBvh& bvh,
                                   const std::vector<double>& f,
                                   const WallDistanceGrid& g,
                                   std::size_t id);
inline double poisson_diagonal(const WallDistanceBvh& bvh,
                               const WallDistanceGrid& g,
                               std::size_t id);

struct WallDistancePoissonAudit {
    std::size_t fluid_nodes{0};
    std::size_t fluid_fluid_faces{0};
    std::size_t solid_faces{0};
    std::size_t outer_faces{0};
    double min_diagonal{std::numeric_limits<double>::infinity()};
    double max_diagonal{0.0};
    double min_diagonal_dominance{std::numeric_limits<double>::infinity()};
    double symmetry_error{0.0};
};

struct WallDistancePoissonOffsetAudit {
    double min_delta{std::numeric_limits<double>::infinity()};
    double max_delta{0.0};
    double min_delta_over_h{std::numeric_limits<double>::infinity()};
    double max_delta_over_h{0.0};
    double max_wall_coefficient{0.0};
    std::size_t max_wall_coefficient_cell{0};
    double max_wall_coefficient_delta{0.0};
    double max_wall_coefficient_h{0.0};
    WallDistanceVec3 max_wall_coefficient_point{};
    std::size_t degenerate_count{0};
    std::size_t solid_face_count{0};
    std::size_t ray_hit_count{0};
    std::size_t ray_miss_count{0};
    std::size_t fallback_count{0};
    std::size_t bad_alignment_count{0};
    double min_alignment{1.0};
};

struct WallDistancePoissonWallOffsetDiagnostic {
    double delta{0.0};
    bool ray_hit{false};
    bool fallback{false};
    bool bad_alignment{false};
    double alignment{1.0};
};

inline WallDistancePoissonWallOffsetDiagnostic poisson_wall_offset_diagnostic(
    const WallDistanceBvh& bvh,const WallDistanceVec3& p,
    const WallDistanceVec3& q,double h);

inline WallDistancePoissonOffsetAudit audit_poisson_wall_offsets(
    const WallDistanceBvh& bvh, const WallDistanceGrid& g) {
    WallDistancePoissonOffsetAudit a;
    for(std::size_t id=0; id<g.points.size(); ++id) {
        if(g.solid[id]) continue;
        const std::size_t k=id/(g.nx*g.ny), rem=id%(g.nx*g.ny), j=rem/g.nx, i=rem%g.nx;
        auto inspect=[&](std::size_t q, bool exists, double h) {
            if(!exists || !g.solid[q]) return;
            ++a.solid_face_count;
            const auto wd=poisson_wall_offset_diagnostic(bvh,g.points[id],g.points[q],h);
            const double delta=wd.delta;
            if(wd.ray_hit) ++a.ray_hit_count;
            else ++a.ray_miss_count;
            if(wd.fallback) ++a.fallback_count;
            if(wd.bad_alignment) ++a.bad_alignment_count;
            a.min_alignment=std::min(a.min_alignment,wd.alignment);
            if(!(std::isfinite(delta) && delta>0.0)) {
                ++a.degenerate_count;
                return;
            }
            const double ratio=delta/h;
            const double coeff=1.0/(h*delta);
            a.min_delta=std::min(a.min_delta,delta);
            a.max_delta=std::max(a.max_delta,delta);
            a.min_delta_over_h=std::min(a.min_delta_over_h,ratio);
            a.max_delta_over_h=std::max(a.max_delta_over_h,ratio);
            if(coeff>a.max_wall_coefficient) {
                a.max_wall_coefficient=coeff;
                a.max_wall_coefficient_cell=id;
                a.max_wall_coefficient_delta=delta;
                a.max_wall_coefficient_h=h;
                a.max_wall_coefficient_point=g.points[id];
            }
        };
        inspect(i>0?g.index(i-1,j,k):0,i>0,g.spacing.x);
        inspect(i+1<g.nx?g.index(i+1,j,k):0,i+1<g.nx,g.spacing.x);
        inspect(j>0?g.index(i,j-1,k):0,j>0,g.spacing.y);
        inspect(j+1<g.ny?g.index(i,j+1,k):0,j+1<g.ny,g.spacing.y);
        inspect(k>0?g.index(i,j,k-1):0,k>0,g.spacing.z);
        inspect(k+1<g.nz?g.index(i,j,k+1):0,k+1<g.nz,g.spacing.z);
    }
    if(!std::isfinite(a.min_delta)) a.min_delta=0.0;
    if(!std::isfinite(a.min_delta_over_h)) a.min_delta_over_h=0.0;
    return a;
}

inline WallDistancePoissonAudit audit_poisson_operator(const WallDistanceBvh& bvh,
                                                   const WallDistanceGrid& g) {
    WallDistancePoissonAudit a;
    const std::size_t n=g.points.size();

    auto offdiag=[&](std::size_t row,std::size_t col)->double {
        if(g.solid[row] || g.solid[col]) return 0.0;
        const std::size_t rk=row/(g.nx*g.ny), rr=row%(g.nx*g.ny), rj=rr/g.nx, ri=rr%g.nx;
        const std::size_t ck=col/(g.nx*g.ny), cr=col%(g.nx*g.ny), cj=cr/g.nx, ci=cr%g.nx;
        if(rj==cj && rk==ck && (ri+1==ci || ci+1==ri))
            return -1.0/(g.spacing.x*g.spacing.x);
        if(ri==ci && rk==ck && (rj+1==cj || cj+1==rj))
            return -1.0/(g.spacing.y*g.spacing.y);
        if(ri==ci && rj==cj && (rk+1==ck || ck+1==rk))
            return -1.0/(g.spacing.z*g.spacing.z);
        return 0.0;
    };

    for(std::size_t id=0;id<n;++id) {
        if(g.solid[id]) continue;
        ++a.fluid_nodes;
        const std::size_t k=id/(g.nx*g.ny), rem=id%(g.nx*g.ny), j=rem/g.nx, i=rem%g.nx;
        const double diag=poisson_diagonal(bvh,g,id);
        a.min_diagonal=std::min(a.min_diagonal,diag);
        a.max_diagonal=std::max(a.max_diagonal,diag);

        double abs_offdiag_sum=0.0;
        auto inspect=[&](std::size_t q,bool exists) {
            if(!exists) {
                ++a.outer_faces;
                return;
            }
            if(g.solid[q]) {
                ++a.solid_faces;
                return;
            }
            if(q>id) ++a.fluid_fluid_faces;
            abs_offdiag_sum += std::abs(offdiag(id,q));
            a.symmetry_error=std::max(
                a.symmetry_error,
                std::abs(offdiag(id,q)-offdiag(q,id)));
        };

        inspect(i>0?g.index(i-1,j,k):0,i>0);
        inspect(i+1<g.nx?g.index(i+1,j,k):0,i+1<g.nx);
        inspect(j>0?g.index(i,j-1,k):0,j>0);
        inspect(j+1<g.ny?g.index(i,j+1,k):0,j+1<g.ny);
        inspect(k>0?g.index(i,j,k-1):0,k>0);
        inspect(k+1<g.nz?g.index(i,j,k+1):0,k+1<g.nz);

        // For the actual operator A=-L:
        //   A_ii = sum(fluid-face weights) + sum(solid-face weights)
        //   A_ij = -fluid-face weight.
        // Therefore the diagonal dominance margin is exactly the total
        // Dirichlet contribution. This uses poisson_diagonal(), including
        // the BVH-based cut-face/wall-offset coefficient actually used by PCG.
        a.min_diagonal_dominance=std::min(
            a.min_diagonal_dominance, diag-abs_offdiag_sum);
    }

    if(!std::isfinite(a.min_diagonal)) a.min_diagonal=0.0;
    if(!std::isfinite(a.min_diagonal_dominance)) a.min_diagonal_dominance=0.0;
    return a;
}

inline WallDistancePoissonWallOffsetDiagnostic poisson_wall_offset_diagnostic(
    const WallDistanceBvh& bvh,const WallDistanceVec3& p,
    const WallDistanceVec3& q,double h) {
    WallDistancePoissonWallOffsetDiagnostic out;
    const auto segment=q-p;
    const double length=wd_norm(segment);
    const double ray=bvh.ray_distance(p,segment,length);
    if(std::isfinite(ray) && ray>1e-12) {
        out.delta=std::clamp(ray,1e-12,h);
        out.ray_hit=true;
        return out;
    }

    out.fallback=true;
    const double d=bvh.nearest_distance(p);
    if(!std::isfinite(d) || d<=0.0) {
        out.delta=std::max(1e-12,h);
        out.bad_alignment=true;
        out.alignment=0.0;
        return out;
    }

    const auto dir=segment*(1.0/std::max(length,1e-30));
    const auto n=bvh.nearest_normal(p);
    out.alignment=std::abs(wd_dot(n,dir));
    if(out.alignment<0.25) {
        out.bad_alignment=true;
        out.delta=std::min(d,h);
    } else {
        out.delta=std::clamp(d/out.alignment,1e-12,h);
    }
    return out;
}

inline double poisson_wall_offset(const WallDistanceBvh& bvh,
                                      const WallDistanceVec3& p,
                                      const WallDistanceVec3& q,
                                      double h) {
    return poisson_wall_offset_diagnostic(bvh,p,q,h).delta;
}

inline double poisson_laplacian_at(const WallDistanceBvh& bvh,
                                   const std::vector<double>& f,
                                   const WallDistanceGrid& g,
                                   std::size_t id) {
    const std::size_t k=id/(g.nx*g.ny), rem=id%(g.nx*g.ny), j=rem/g.nx, i=rem%g.nx;
    double l=0.0;
    auto add_face=[&](std::size_t q,bool exists,double h) {
        if(!exists) return;
        if(g.solid[q]) {
            const double delta=poisson_wall_offset(bvh,g.points[id],g.points[q],h);
            l-=f[id]/(h*delta);
        } else if(std::isfinite(f[q])) {
            l+=(f[q]-f[id])/(h*h);
        }
    };
    const std::size_t xm=i>0?g.index(i-1,j,k):0, xp=i+1<g.nx?g.index(i+1,j,k):0;
    const std::size_t ym=j>0?g.index(i,j-1,k):0, yp=j+1<g.ny?g.index(i,j+1,k):0;
    const std::size_t zm=k>0?g.index(i,j,k-1):0, zp=k+1<g.nz?g.index(i,j,k+1):0;
    add_face(xm,i>0,g.spacing.x); add_face(xp,i+1<g.nx,g.spacing.x);
    add_face(ym,j>0,g.spacing.y); add_face(yp,j+1<g.ny,g.spacing.y);
    add_face(zm,k>0,g.spacing.z); add_face(zp,k+1<g.nz,g.spacing.z);
    return l;
}

inline double poisson_diagonal(const WallDistanceBvh& bvh,
                               const WallDistanceGrid& g,
                               std::size_t id) {
    const std::size_t k=id/(g.nx*g.ny), rem=id%(g.nx*g.ny), j=rem/g.nx, i=rem%g.nx;
    double d=0.0;
    auto add_face=[&](std::size_t q,bool exists,double h) {
        if(!exists) return;
        d += g.solid[q] ? 1.0/(h*poisson_wall_offset(bvh,g.points[id],g.points[q],h))
                        : 1.0/(h*h);
    };
    const std::size_t xm=i>0?g.index(i-1,j,k):0, xp=i+1<g.nx?g.index(i+1,j,k):0;
    const std::size_t ym=j>0?g.index(i,j-1,k):0, yp=j+1<g.ny?g.index(i,j+1,k):0;
    const std::size_t zm=k>0?g.index(i,j,k-1):0, zp=k+1<g.nz?g.index(i,j,k+1):0;
    add_face(xm,i>0,g.spacing.x); add_face(xp,i+1<g.nx,g.spacing.x);
    add_face(ym,j>0,g.spacing.y); add_face(yp,j+1<g.ny,g.spacing.y);
    add_face(zm,k>0,g.spacing.z); add_face(zp,k+1<g.nz,g.spacing.z);
    return d;
}

inline std::vector<double> poisson_potential(const WallDistanceBvh& bvh,const WallDistanceGrid& g,
                                              std::size_t max_iter,double,
                                              std::size_t* used_iter=nullptr,
                                              double* residual_out=nullptr) {
    // Solve the actual elliptic boundary-value problem
    //
    //     -laplacian(phi) = 1,    phi = 0 on the wall,
    //
    // with homogeneous Neumann treatment at the outer computational boundary.
    //
    // Solid neighbours are Dirichlet values phi_b=0 in the discrete operator.
    // Outer boundaries are homogeneous Neumann zero-flux faces.  The stencil
    // is assembled as a symmetric face graph so PCG is mathematically valid.
    // No artificial zero-valued seed band is imposed on the Poisson solution.
    //
    // The operator is symmetric positive definite for a domain connected to
    // at least one Dirichlet wall. Use preconditioned conjugate gradients,
    // as OpenFOAM does for its Poisson wall-distance solve (PCG/DIC or GAMG).
    const std::size_t n=g.points.size();
    std::vector<double> phi(n,0.0);
    std::vector<double> r(n,0.0);
    std::vector<double> z(n,0.0);
    std::vector<double> p(n,0.0);
    std::vector<double> Ap(n,0.0);

    auto is_fluid=[](const WallDistanceGrid& grid,std::size_t id) {
        return !grid.solid[id];
    };

    auto diagonal=[&](std::size_t id) {
        return poisson_diagonal(bvh,g,id);
    };

    auto apply=[&](const std::vector<double>& x,std::vector<double>& y) {
        std::fill(y.begin(),y.end(),0.0);
        for(std::size_t id=0;id<n;++id) {
            if(!is_fluid(g,id)) continue;
            // A = -L, with the same solid Dirichlet and outer Neumann
            // treatment used by poisson_laplacian_at().
            y[id] = -poisson_laplacian_at(bvh,x,g,id);
        }
    };

    // b = 1 in every fluid degree of freedom.
    apply(phi,Ap);
    double residual_inf=0.0;
    std::size_t fluid_count=0;
    for(std::size_t id=0;id<n;++id) {
        if(!is_fluid(g,id)) continue;
        r[id]=1.0-Ap[id];
        residual_inf=std::max(residual_inf,std::abs(r[id]));
        ++fluid_count;
    }

    if(fluid_count==0) {
        if(used_iter) *used_iter=0;
        if(residual_out) *residual_out=0.0;
        return phi;
    }

    for(std::size_t id=0;id<n;++id) {
        if(is_fluid(g,id)) {
            z[id]=r[id]/diagonal(id);
            p[id]=z[id];
        }
    }

    double rz_old=0.0;
    for(std::size_t id=0;id<n;++id)
        if(is_fluid(g,id)) rz_old += r[id]*z[id];

    std::size_t used=0;
    constexpr double abs_tol=1e-10;
    constexpr double breakdown_tol=1e-30;

    if(residual_inf<=abs_tol) {
        if(used_iter) *used_iter=0;
        if(residual_out) *residual_out=residual_inf;
        return phi;
    }

    for(std::size_t it=0;it<max_iter;++it) {
        apply(p,Ap);

        double pAp=0.0;
        for(std::size_t id=0;id<n;++id)
            if(is_fluid(g,id)) pAp += p[id]*Ap[id];

        if(!(pAp>breakdown_tol) || !std::isfinite(pAp)) {
            used=it;
            break;
        }

        const double alpha=rz_old/pAp;
        if(!std::isfinite(alpha)) {
            used=it;
            break;
        }

        double max_change=0.0;
        for(std::size_t id=0;id<n;++id) {
            if(!is_fluid(g,id)) continue;
            const double old=phi[id];
            phi[id] += alpha*p[id];
            r[id] -= alpha*Ap[id];
            max_change=std::max(max_change,std::abs(phi[id]-old));
        }

        residual_inf=0.0;
        double rz_new=0.0;
        for(std::size_t id=0;id<n;++id) {
            if(!is_fluid(g,id)) continue;
            residual_inf=std::max(residual_inf,std::abs(r[id]));
            z[id]=r[id]/diagonal(id);
            rz_new += r[id]*z[id];
        }

        used=it+1;
        if(residual_inf<=abs_tol) break;

        if(!(std::isfinite(rz_new)) || rz_old<=breakdown_tol) break;

        const double beta=rz_new/rz_old;
        for(std::size_t id=0;id<n;++id)
            if(is_fluid(g,id)) p[id]=z[id]+beta*p[id];

        rz_old=rz_new;

        // A stagnating Krylov iteration is a numerical failure, not
        // convergence. Keep the residual visible to the caller.
        if(max_change==0.0 && residual_inf>abs_tol) break;
    }

    if(used_iter) *used_iter=used;
    if(residual_out) *residual_out=residual_inf;

    return phi;
}

struct WallDistancePoissonReconstructionAudit {
    double phi_min{std::numeric_limits<double>::infinity()};
    double phi_max{-std::numeric_limits<double>::infinity()};
    double grad_min{std::numeric_limits<double>::infinity()};
    double grad_max{0.0};
    double distance_l2_error{0.0};
    double distance_linf_error{0.0};
};

inline std::vector<double> poisson_distance(const WallDistanceBvh& bvh,const WallDistanceGrid& g,
                                             std::size_t max_iter,double smooth,
                                             std::size_t* used_iter=nullptr,
                                             double* residual_out=nullptr,
                                             WallDistancePoissonReconstructionAudit* audit_out=nullptr) {
    std::vector<double> phi=poisson_potential(bvh,g,max_iter,smooth,used_iter,residual_out);
    std::vector<double> d(g.points.size(),std::numeric_limits<double>::infinity());
    for(std::size_t id=0;id<g.points.size();++id) if(!g.solid[id]) {
        const std::size_t k=id/(g.nx*g.ny), rem=id%(g.nx*g.ny), j=rem/g.nx, i=rem%g.nx;
        auto deriv=[&](int axis)->double {
            if(axis==0) {
                if(i>0&&i+1<g.nx&&!g.solid[g.index(i-1,j,k)]&&!g.solid[g.index(i+1,j,k)]) return (phi[g.index(i+1,j,k)]-phi[g.index(i-1,j,k)])/(2*g.spacing.x);
                if(i+1<g.nx) {
                    const auto q=g.index(i+1,j,k);
                    if(g.solid[q]) return -phi[id]/poisson_wall_offset(bvh,g.points[id],g.points[q],g.spacing.x);
                    return (phi[q]-phi[id])/g.spacing.x;
                }
                if(i>0) {
                    const auto q=g.index(i-1,j,k);
                    if(g.solid[q]) return phi[id]/poisson_wall_offset(bvh,g.points[id],g.points[q],g.spacing.x);
                    return (phi[id]-phi[q])/g.spacing.x;
                }
            } else if(axis==1) {
                if(j>0&&j+1<g.ny&&!g.solid[g.index(i,j-1,k)]&&!g.solid[g.index(i,j+1,k)]) return (phi[g.index(i,j+1,k)]-phi[g.index(i,j-1,k)])/(2*g.spacing.y);
                if(j+1<g.ny) {
                    const auto q=g.index(i,j+1,k);
                    if(g.solid[q]) return -phi[id]/poisson_wall_offset(bvh,g.points[id],g.points[q],g.spacing.y);
                    return (phi[q]-phi[id])/g.spacing.y;
                }
                if(j>0) {
                    const auto q=g.index(i,j-1,k);
                    if(g.solid[q]) return phi[id]/poisson_wall_offset(bvh,g.points[id],g.points[q],g.spacing.y);
                    return (phi[id]-phi[q])/g.spacing.y;
                }
            } else {
                if(k>0&&k+1<g.nz&&!g.solid[g.index(i,j,k-1)]&&!g.solid[g.index(i,j,k+1)]) return (phi[g.index(i,j,k+1)]-phi[g.index(i,j,k-1)])/(2*g.spacing.z);
                if(k+1<g.nz) {
                    const auto q=g.index(i,j,k+1);
                    if(g.solid[q]) return -phi[id]/poisson_wall_offset(bvh,g.points[id],g.points[q],g.spacing.z);
                    return (phi[q]-phi[id])/g.spacing.z;
                }
                if(k>0) {
                    const auto q=g.index(i,j,k-1);
                    if(g.solid[q]) return phi[id]/poisson_wall_offset(bvh,g.points[id],g.points[q],g.spacing.z);
                    return (phi[id]-phi[q])/g.spacing.z;
                }
            }
            return 0.0;
        };
        const double gx=deriv(0),gy=deriv(1),gz=deriv(2),grad=std::sqrt(gx*gx+gy*gy+gz*gz);
        d[id]=std::max(0.0,std::sqrt(std::max(0.0,grad*grad+2*phi[id]))-grad);
    }
    if(audit_out) {
        WallDistancePoissonReconstructionAudit a;
        double e2=0.0, ref2=0.0, emax=0.0, refmax=0.0;
        for(std::size_t id=0;id<g.points.size();++id) {
            if(g.solid[id] || !std::isfinite(phi[id]) || !std::isfinite(d[id])) continue;
            const std::size_t k=id/(g.nx*g.ny), rem=id%(g.nx*g.ny), j=rem/g.nx, i=rem%g.nx;
            auto deriv_local=[&](int axis)->double {
                if(axis==0) {
                    if(i>0&&i+1<g.nx&&!g.solid[g.index(i-1,j,k)]&&!g.solid[g.index(i+1,j,k)])
                        return (phi[g.index(i+1,j,k)]-phi[g.index(i-1,j,k)])/(2*g.spacing.x);
                    if(i+1<g.nx) { const auto q=g.index(i+1,j,k); if(g.solid[q]) return -phi[id]/poisson_wall_offset(bvh,g.points[id],g.points[q],g.spacing.x); return (phi[q]-phi[id])/g.spacing.x; }
                    if(i>0) { const auto q=g.index(i-1,j,k); if(g.solid[q]) return phi[id]/poisson_wall_offset(bvh,g.points[id],g.points[q],g.spacing.x); return (phi[id]-phi[q])/g.spacing.x; }
                } else if(axis==1) {
                    if(j>0&&j+1<g.ny&&!g.solid[g.index(i,j-1,k)]&&!g.solid[g.index(i,j+1,k)])
                        return (phi[g.index(i,j+1,k)]-phi[g.index(i,j-1,k)])/(2*g.spacing.y);
                    if(j+1<g.ny) { const auto q=g.index(i,j+1,k); if(g.solid[q]) return -phi[id]/poisson_wall_offset(bvh,g.points[id],g.points[q],g.spacing.y); return (phi[q]-phi[id])/g.spacing.y; }
                    if(j>0) { const auto q=g.index(i,j-1,k); if(g.solid[q]) return phi[id]/poisson_wall_offset(bvh,g.points[id],g.points[q],g.spacing.y); return (phi[id]-phi[q])/g.spacing.y; }
                } else {
                    if(k>0&&k+1<g.nz&&!g.solid[g.index(i,j,k-1)]&&!g.solid[g.index(i,j,k+1)])
                        return (phi[g.index(i,j,k+1)]-phi[g.index(i,j,k-1)])/(2*g.spacing.z);
                    if(k+1<g.nz) { const auto q=g.index(i,j,k+1); if(g.solid[q]) return -phi[id]/poisson_wall_offset(bvh,g.points[id],g.points[q],g.spacing.z); return (phi[q]-phi[id])/g.spacing.z; }
                    if(k>0) { const auto q=g.index(i,j,k-1); if(g.solid[q]) return phi[id]/poisson_wall_offset(bvh,g.points[id],g.points[q],g.spacing.z); return (phi[id]-phi[q])/g.spacing.z; }
                }
                return 0.0;
            };
            const double gx=deriv_local(0), gy=deriv_local(1), gz=deriv_local(2);
            const double grad=std::sqrt(gx*gx+gy*gy+gz*gz);
            a.phi_min=std::min(a.phi_min,phi[id]); a.phi_max=std::max(a.phi_max,phi[id]);
            a.grad_min=std::min(a.grad_min,grad); a.grad_max=std::max(a.grad_max,grad);
            const double ref=bvh.nearest_distance(g.points[id]), e=d[id]-ref;
            e2+=e*e; ref2+=ref*ref; emax=std::max(emax,std::abs(e)); refmax=std::max(refmax,ref);
        }
        a.distance_l2_error=std::sqrt(e2/std::max(1e-30,ref2));
        a.distance_linf_error=emax/std::max(1e-30,refmax);
        *audit_out=a;
    }
    return d;
}

inline std::vector<double> poisson_distance(const WallSurface& s,const WallDistanceGrid& g,
                                             std::size_t max_iter,double smooth,
                                             std::size_t* used_iter=nullptr,
                                             double* residual_out=nullptr,
                                             WallDistancePoissonReconstructionAudit* audit_out=nullptr) {
    const WallDistanceBvh bvh(s);
    return poisson_distance(bvh,g,max_iter,smooth,used_iter,residual_out,audit_out);
}

inline double wall_distance_pde_residual_inf(WallDistanceMethod method,
                                                        const std::vector<double>& d,
                                                        const WallDistanceGrid& g,
                                                        const std::vector<unsigned char>* fixed=nullptr,
                                                        double advection_diffusion_gamma=0.05,
                                                        const WallDistanceBvh* bvh=nullptr);

inline std::vector<double> advection_diffusion_distance(const WallSurface& s,const WallDistanceGrid& g,
                                                        std::size_t max_iter,double gamma=0.05,
                                                        std::size_t* used_iter=nullptr,
                                                        double* residual_out=nullptr) {
    // NASA/Tucker transport form: U·∇d = 1 + Gamma ∇²d, with U derived
    // from the current distance field. Advection is first-order upwind and
    // diffusion is second-order central, matching the published formulation.
    const double h=std::min({g.spacing.x,g.spacing.y,g.spacing.z});
    const WallDistanceBvh bvh(s);
    const auto seeds=wall_seed_nodes(bvh,g,1.6*h);
    std::vector<unsigned char> fixed(g.points.size(),0);
    std::vector<double> d=eikonal_fast_sweep(bvh,g,std::max<std::size_t>(8,max_iter/2));
    for(std::size_t id=0;id<d.size();++id) if(g.solid[id]) d[id]=0.0;
    for(auto id:seeds) { d[id]=bvh.nearest_distance(g.points[id]); fixed[id]=1; }
    std::size_t used=max_iter;
    double final_residual=std::numeric_limits<double>::infinity();
    for(std::size_t it=0;it<max_iter;++it) {
        double max_change=0.0;
        for(std::size_t id=0;id<d.size();++id) {
            if(g.solid[id]||fixed[id]||!std::isfinite(d[id])) continue;
            const std::size_t k=id/(g.nx*g.ny), rem=id%(g.nx*g.ny), j=rem/g.nx, i=rem%g.nx;
            const auto one=[&](int axis)->double {
                double gm=0.0,gp=0.0;
                if(axis==0) { if(i>0) gm=(d[g.index(i,j,k)]-d[g.index(i-1,j,k)])/g.spacing.x; if(i+1<g.nx) gp=(d[g.index(i+1,j,k)]-d[g.index(i,j,k)])/g.spacing.x; }
                if(axis==1) { if(j>0) gm=(d[g.index(i,j,k)]-d[g.index(i,j-1,k)])/g.spacing.y; if(j+1<g.ny) gp=(d[g.index(i,j+1,k)]-d[g.index(i,j,k)])/g.spacing.y; }
                if(axis==2) { if(k>0) gm=(d[g.index(i,j,k)]-d[g.index(i,j,k-1)])/g.spacing.z; if(k+1<g.nz) gp=(d[g.index(i,j,k+1)]-d[g.index(i,j,k)])/g.spacing.z; }
                return std::max(gm,0.0)+std::min(gp,0.0);
            };
            const double gx=one(0),gy=one(1),gz=one(2),gn=std::sqrt(gx*gx+gy*gy+gz*gz);
            const double ux=gn>1e-14?gx/gn:0.0, uy=gn>1e-14?gy/gn:0.0, uz=gn>1e-14?gz/gn:0.0;
            double diag=0.0,rhs=1.0;
            auto add_axis=[&](std::size_t m,bool has_m,std::size_t p,bool has_p,
                                   double hh,double u) {
                auto face_delta=[&](std::size_t q)->double {
                    return g.solid[q] ? poisson_wall_offset(bvh,g.points[id],g.points[q],hh) : hh;
                };
                if(has_m) {
                    const double delta=face_delta(m);
                    const double w=gamma/(hh*delta);
                    const double dm=g.solid[m]?0.0:d[m];
                    diag+=w; rhs+=w*dm;
                    if(u>0.0) {
                        const double a=u/delta;
                        diag+=a; if(!g.solid[m]) rhs+=a*d[m];
                    }
                }
                if(has_p) {
                    const double delta=face_delta(p);
                    const double w=gamma/(hh*delta);
                    const double dp=g.solid[p]?0.0:d[p];
                    diag+=w; rhs+=w*dp;
                    if(u<0.0) {
                        const double a=-u/delta;
                        diag+=a; if(!g.solid[p]) rhs+=a*d[p];
                    }
                }
            };
            const std::size_t xm=i>0?g.index(i-1,j,k):0, xp=i+1<g.nx?g.index(i+1,j,k):0;
            const std::size_t ym=j>0?g.index(i,j-1,k):0, yp=j+1<g.ny?g.index(i,j+1,k):0;
            const std::size_t zm=k>0?g.index(i,j,k-1):0, zp=k+1<g.nz?g.index(i,j,k+1):0;
            add_axis(xm,i>0,xp,i+1<g.nx,g.spacing.x,ux);
            add_axis(ym,j>0,yp,j+1<g.ny,g.spacing.y,uy);
            add_axis(zm,k>0,zp,k+1<g.nz,g.spacing.z,uz);
            if(diag>0.0) { const double nd=rhs/diag; max_change=std::max(max_change,std::abs(nd-d[id])); d[id]=nd; }
        }
        for(auto id:seeds) d[id]=bvh.nearest_distance(g.points[id]);
        if((it&3u)==3u || max_change<1e-10*h) {
            final_residual=wall_distance_pde_residual_inf(
                WallDistanceMethod::ADVECTION_DIFFUSION,d,g,&fixed,gamma,&bvh);
            if(final_residual<1e-6 && max_change<1e-8*h) {
                used=it+1;
                break;
            }
        }
    }
    if(used_iter) *used_iter=used;
    if(residual_out) *residual_out=final_residual;
    return d;
}

inline std::vector<double> hybrid_poisson_hamilton_jacobi_distance(
    const WallSurface& s,const WallDistanceGrid& g,std::size_t max_iter,
    double epsilon=0.25,double relaxation=0.9,
    std::size_t* used_iter=nullptr,double* residual_out=nullptr,
    double poisson_weight=1.0,
    std::size_t* poisson_iter_out=nullptr,double* poisson_residual_out=nullptr) {
    // Tucker 2011: Poisson supplies an auxiliary front-propagation direction
    // (effectively a wall normal); H-J then propagates distance with that
    // direction, optionally blended with the evolving Eikonal direction.
    const double h=std::min({g.spacing.x,g.spacing.y,g.spacing.z});
    const WallDistanceBvh bvh(s);
    const auto seeds=wall_seed_nodes(bvh,g,1.6*h);
    std::vector<unsigned char> fixed(g.points.size(),0);
    for(auto id:seeds) fixed[id]=1;

    std::size_t poisson_iter=0; double poisson_residual=0.0;
    const std::size_t poisson_max_iter=std::max<std::size_t>(200,4*max_iter);
    const auto phi=poisson_potential(bvh,g,poisson_max_iter,1.5,
                                     &poisson_iter,&poisson_residual);
    if(poisson_iter_out) *poisson_iter_out=poisson_iter;
    if(poisson_residual_out) *poisson_residual_out=poisson_residual;

    std::vector<double> d(g.points.size(),std::numeric_limits<double>::infinity());
    auto grad_comp=[&](const std::vector<double>& f,std::size_t id,int a)->double {
        const std::size_t k=id/(g.nx*g.ny),rem=id%(g.nx*g.ny),j=rem/g.nx,i=rem%g.nx;
        if(a==0) {
            const bool hm=i>0, hp=i+1<g.nx;
            const auto im=hm?g.index(i-1,j,k):0, ip=hp?g.index(i+1,j,k):0;
            const bool fm=hm && !g.solid[im], fp=hp && !g.solid[ip];
            if(fm&&fp) return (f[ip]-f[im])/(2*g.spacing.x);
            if(fp) return (f[ip]-f[id])/g.spacing.x;
            if(fm) return (f[id]-f[im])/g.spacing.x;
        } else if(a==1) {
            const bool hm=j>0, hp=j+1<g.ny;
            const auto im=hm?g.index(i,j-1,k):0, ip=hp?g.index(i,j+1,k):0;
            const bool fm=hm && !g.solid[im], fp=hp && !g.solid[ip];
            if(fm&&fp) return (f[ip]-f[im])/(2*g.spacing.y);
            if(fp) return (f[ip]-f[id])/g.spacing.y;
            if(fm) return (f[id]-f[im])/g.spacing.y;
        } else {
            const bool hm=k>0, hp=k+1<g.nz;
            const auto im=hm?g.index(i,j,k-1):0, ip=hp?g.index(i,j,k+1):0;
            const bool fm=hm && !g.solid[im], fp=hp && !g.solid[ip];
            if(fm&&fp) return (f[ip]-f[im])/(2*g.spacing.z);
            if(fp) return (f[ip]-f[id])/g.spacing.z;
            if(fm) return (f[id]-f[im])/g.spacing.z;
        }
        return 0.0;
    };
    for(std::size_t id=0;id<d.size();++id) if(!g.solid[id]) {
        const double gx=grad_comp(phi,id,0),gy=grad_comp(phi,id,1),gz=grad_comp(phi,id,2);
        const double gn=std::sqrt(gx*gx+gy*gy+gz*gz);
        d[id]=gn>1e-14 ? std::max(0.0,std::sqrt(std::max(0.0,gn*gn+2*phi[id]))-gn) : 0.0;
    }
    for(auto id:seeds) d[id]=bvh.nearest_distance(g.points[id]);

    const double alpha=std::clamp(poisson_weight,0.0,1.0);
    const std::size_t steps=std::max<std::size_t>(20,max_iter);
    std::size_t used=steps; double final_residual=std::numeric_limits<double>::infinity();

    auto velocity=[&](std::size_t id,double& ux,double& uy,double& uz) {
        const double px=grad_comp(phi,id,0),py=grad_comp(phi,id,1),pz=grad_comp(phi,id,2);
        const double pn=std::sqrt(px*px+py*py+pz*pz);
        ux=pn>1e-14?px/pn:0.0; uy=pn>1e-14?py/pn:0.0; uz=pn>1e-14?pz/pn:0.0;
        if(alpha<1.0) {
            const double dx=grad_comp(d,id,0),dy=grad_comp(d,id,1),dz=grad_comp(d,id,2);
            const double dn=std::sqrt(dx*dx+dy*dy+dz*dz);
            if(dn>1e-14) {
                ux=alpha*ux+(1-alpha)*dx/dn; uy=alpha*uy+(1-alpha)*dy/dn; uz=alpha*uz+(1-alpha)*dz/dn;
                const double n=std::sqrt(ux*ux+uy*uy+uz*uz);
                if(n>1e-14){ux/=n;uy/=n;uz/=n;}
            }
        }
    };

    for(std::size_t it=0;it<steps;++it) {
        double max_change=0.0;
        for(std::size_t id=0;id<d.size();++id) {
            if(g.solid[id]||fixed[id]||!std::isfinite(d[id])) continue;
            const std::size_t k=id/(g.nx*g.ny),rem=id%(g.nx*g.ny),j=rem/g.nx,i=rem%g.nx;
            double ux,uy,uz; velocity(id,ux,uy,uz);
            const auto up=[&](double u,int a)->double {
                auto value=[&](std::size_t q)->double {
                    return g.solid[q] ? 0.0 : d[q];
                };
                if(a==0) {
                    if(u>=0) {
                        if(i==0) return 0.0;
                        const auto q=g.index(i-1,j,k);
                        const double hh=g.spacing.x;
                        const double delta=g.solid[q]?poisson_wall_offset(bvh,g.points[id],g.points[q],hh):hh;
                        return u*(d[id]-value(q))/delta;
                    }
                    if(i+1==g.nx) return 0.0;
                    const auto q=g.index(i+1,j,k);
                    const double hh=g.spacing.x;
                    const double delta=g.solid[q]?poisson_wall_offset(bvh,g.points[id],g.points[q],hh):hh;
                    return u*(value(q)-d[id])/delta;
                }
                if(a==1) {
                    if(u>=0) {
                        if(j==0) return 0.0;
                        const auto q=g.index(i,j-1,k);
                        const double hh=g.spacing.y;
                        const double delta=g.solid[q]?poisson_wall_offset(bvh,g.points[id],g.points[q],hh):hh;
                        return u*(d[id]-value(q))/delta;
                    }
                    if(j+1==g.ny) return 0.0;
                    const auto q=g.index(i,j+1,k);
                    const double hh=g.spacing.y;
                    const double delta=g.solid[q]?poisson_wall_offset(bvh,g.points[id],g.points[q],hh):hh;
                    return u*(value(q)-d[id])/delta;
                }
                if(u>=0) {
                    if(k==0) return 0.0;
                    const auto q=g.index(i,j,k-1);
                    const double hh=g.spacing.z;
                    const double delta=g.solid[q]?poisson_wall_offset(bvh,g.points[id],g.points[q],hh):hh;
                    return u*(d[id]-value(q))/delta;
                }
                if(k+1==g.nz) return 0.0;
                const auto q=g.index(i,j,k+1);
                const double hh=g.spacing.z;
                const double delta=g.solid[q]?poisson_wall_offset(bvh,g.points[id],g.points[q],hh):hh;
                return u*(value(q)-d[id])/delta;
            };
            const double adv=up(ux,0)+up(uy,1)+up(uz,2);
            const double gamma=epsilon*std::max(d[id],0.0);
            const double residual=adv-1.0-gamma*laplacian_at(bvh,d,g,id);
            const double inv_h2=1.0/(g.spacing.x*g.spacing.x)+1.0/(g.spacing.y*g.spacing.y)+1.0/(g.spacing.z*g.spacing.z);
            // The hybrid transport step is also a 3-D explicit upwind
            // update. Use its actual velocity components for the CFL bound.
            const double cfl_rate=std::abs(ux)/g.spacing.x+
                                  std::abs(uy)/g.spacing.y+
                                  std::abs(uz)/g.spacing.z;
            const double dt_adv=cfl_rate>1e-14 ? 0.9/cfl_rate :
                                                std::numeric_limits<double>::infinity();
            const double dt_diff=gamma>0.0?0.9/(2*gamma*inv_h2):std::numeric_limits<double>::infinity();
            const double dt=std::clamp(relaxation,0.1,1.0)*std::min(dt_adv,dt_diff);
            const double nd=std::max(0.0,d[id]-dt*residual);
            max_change=std::max(max_change,std::abs(nd-d[id])); d[id]=nd;
        }
        for(auto id:seeds) d[id]=bvh.nearest_distance(g.points[id]);
        if((it&3u)==3u||max_change<1e-10*h) {
            // Evaluate exactly the same upwind transport operator used by the
            // pseudo-time update. The previous diagnostic used a central
            // gradient here, so the reported residual could disagree with
            // the nonlinear operator actually being iterated.
            final_residual=0.0;
            for(std::size_t id=0;id<d.size();++id) if(!g.solid[id]&&!fixed[id]&&std::isfinite(d[id])) {
                double ux,uy,uz; velocity(id,ux,uy,uz);
                const std::size_t k=id/(g.nx*g.ny),rem=id%(g.nx*g.ny),j=rem/g.nx,i=rem%g.nx;
                auto upwind=[&](double u,int axis)->double {
                    auto value=[&](std::size_t q)->double { return g.solid[q] ? 0.0 : d[q]; };
                    if(axis==0) {
                        if(u>=0.0) { if(i==0) return 0.0; return u*(d[id]-value(g.index(i-1,j,k)))/g.spacing.x; }
                        if(i+1==g.nx) return 0.0;
                        return u*(value(g.index(i+1,j,k))-d[id])/g.spacing.x;
                    }
                    if(axis==1) {
                        if(u>=0.0) { if(j==0) return 0.0; return u*(d[id]-value(g.index(i,j-1,k)))/g.spacing.y; }
                        if(j+1==g.ny) return 0.0;
                        return u*(value(g.index(i,j+1,k))-d[id])/g.spacing.y;
                    }
                    if(u>=0.0) { if(k==0) return 0.0; return u*(d[id]-value(g.index(i,j,k-1)))/g.spacing.z; }
                    if(k+1==g.nz) return 0.0;
                    return u*(value(g.index(i,j,k+1))-d[id])/g.spacing.z;
                };
                const double adv=upwind(ux,0)+upwind(uy,1)+upwind(uz,2);
                final_residual=std::max(final_residual,
                    std::abs(adv-1.0-epsilon*std::max(d[id],0.0)*laplacian_at(bvh,d,g,id)));
            }
            if(final_residual<1e-8&&max_change<1e-9*h){used=it+1;break;}
        }
    }
    if(used_iter) *used_iter=used;
    if(residual_out) *residual_out=final_residual;
    return d;
}

inline double wall_distance_pde_residual_inf(WallDistanceMethod method,
                                                        const std::vector<double>& d,
                                                        const WallDistanceGrid& g,
                                                        const std::vector<unsigned char>* fixed,
                                                        double advection_diffusion_gamma,
                                                        const WallDistanceBvh* bvh) {
    double rmax=0.0;
    for(std::size_t id=0;id<d.size();++id) {
        if(g.solid[id] || !std::isfinite(d[id]) || (fixed && (*fixed)[id])) continue;
        if(method==WallDistanceMethod::EIKONAL) {
            const double grad=bvh ? godunov_gradient_at(*bvh,d,g,id)
                                  : godunov_gradient_at(d,g,id);
            rmax=std::max(rmax,std::abs(grad-1.0));
            continue;
        }
        if(method==WallDistanceMethod::HAMILTON_JACOBI ||
           method==WallDistanceMethod::HYBRID_POISSON_EIKONAL) {
            const double gamma=0.25*std::max(d[id],0.0);
            const double grad=bvh ? godunov_gradient_at(*bvh,d,g,id)
                                  : godunov_gradient_at(d,g,id);
            const double lap=bvh ? laplacian_at(*bvh,d,g,id)
                                 : laplacian_at(d,g,id);
            rmax=std::max(rmax,std::abs(grad-1.0-gamma*lap));
            continue;
        }
        if(method==WallDistanceMethod::ADVECTION_DIFFUSION) {
            const std::size_t k=id/(g.nx*g.ny), rem=id%(g.nx*g.ny), j=rem/g.nx, i=rem%g.nx;
            auto value=[&](std::size_t q)->double { return g.solid[q] ? 0.0 : d[q]; };
            auto one=[&](int axis)->double {
                double gm=0.0,gp=0.0;
                if(axis==0) {
                    if(i>0) gm=(d[id]-value(g.index(i-1,j,k)))/g.spacing.x;
                    if(i+1<g.nx) gp=(value(g.index(i+1,j,k))-d[id])/g.spacing.x;
                } else if(axis==1) {
                    if(j>0) gm=(d[id]-value(g.index(i,j-1,k)))/g.spacing.y;
                    if(j+1<g.ny) gp=(value(g.index(i,j+1,k))-d[id])/g.spacing.y;
                } else {
                    if(k>0) gm=(d[id]-value(g.index(i,j,k-1)))/g.spacing.z;
                    if(k+1<g.nz) gp=(value(g.index(i,j,k+1))-d[id])/g.spacing.z;
                }
                return std::max(gm,0.0)+std::min(gp,0.0);
            };
            const double gx=one(0), gy=one(1), gz=one(2);
            const double gn=std::sqrt(gx*gx+gy*gy+gz*gz);
            if(!(gn>1e-14)) continue;
            const double ux=gx/gn, uy=gy/gn, uz=gz/gn;
            auto upwind=[&](double u,int axis)->double {
                if(axis==0) {
                    if(u>0.0 && i>0) return u*(d[id]-value(g.index(i-1,j,k)))/g.spacing.x;
                    if(u<0.0 && i+1<g.nx) return u*(value(g.index(i+1,j,k))-d[id])/g.spacing.x;
                } else if(axis==1) {
                    if(u>0.0 && j>0) return u*(d[id]-value(g.index(i,j-1,k)))/g.spacing.y;
                    if(u<0.0 && j+1<g.ny) return u*(value(g.index(i,j+1,k))-d[id])/g.spacing.y;
                } else {
                    if(u>0.0 && k>0) return u*(d[id]-value(g.index(i,j,k-1)))/g.spacing.z;
                    if(u<0.0 && k+1<g.nz) return u*(value(g.index(i,j,k+1))-d[id])/g.spacing.z;
                }
                return 0.0;
            };
            const double transport=bvh ? godunov_gradient_at(*bvh,d,g,id)
                                        : upwind(ux,0)+upwind(uy,1)+upwind(uz,2);
            const double lap=bvh ? laplacian_at(*bvh,d,g,id)
                                 : laplacian_at(d,g,id);
            rmax=std::max(rmax,std::abs(
                transport-1.0-advection_diffusion_gamma*lap));
        }
    }
    return rmax;
}

inline WallDistanceResult compute_wall_distance(WallDistanceMethod method,const WallSurface& s,
                                                 const WallDistanceGrid& g,
                                                 std::size_t iterations=500) {
    WallDistanceResult r;
    r.method=wall_distance_method_name(method);
    r.valid.assign(g.points.size(),1);
    switch(method) {
        case WallDistanceMethod::EXACT_GEOMETRIC:
            r.distance=exact_reference(s,g);
            r.converged=true;
            r.stopping_reason="direct";
            r.convergence_status=WallDistanceConvergenceStatus::DIRECT;
            break;
        case WallDistanceMethod::SEARCH_BASED:
            r.distance=search_based_reference(s,g,4.0*std::min({g.spacing.x,g.spacing.y,g.spacing.z}));
            r.converged=true;
            r.stopping_reason="direct";
            break;
        case WallDistanceMethod::MESH_WAVE:
            r=graph_wave(s,g,false);
            r.converged=true;
            r.stopping_reason="graph_complete";
            r.convergence_status=WallDistanceConvergenceStatus::DIRECT;
            break;
        case WallDistanceMethod::DIRECTIONAL_MESH_WAVE:
            r=graph_wave(s,g,true);
            r.converged=true;
            r.stopping_reason="graph_complete";
            break;
        case WallDistanceMethod::POISSON: {
            WallDistancePoissonReconstructionAudit pa;
            r.distance=poisson_distance(s,g,iterations,1.5,&r.iterations,&r.residual_inf,&pa);
            const WallDistanceBvh pbvh(s);
            const auto wa=audit_poisson_wall_offsets(pbvh,g);
            r.poisson_wall_ray_hits=wa.ray_hit_count;
            r.poisson_wall_ray_misses=wa.ray_miss_count;
            r.poisson_wall_fallbacks=wa.fallback_count;
            r.poisson_wall_bad_alignment=wa.bad_alignment_count;
            r.poisson_wall_min_alignment=wa.min_alignment;
            r.poisson_phi_min=pa.phi_min;
            r.poisson_phi_max=pa.phi_max;
            r.poisson_grad_min=pa.grad_min;
            r.poisson_grad_max=pa.grad_max;
            r.poisson_distance_l2_error=pa.distance_l2_error;
            r.poisson_distance_linf_error=pa.distance_linf_error;
            r.converged = std::isfinite(r.residual_inf) && r.residual_inf < 1e-6;
            r.convergence_status = r.converged ? WallDistanceConvergenceStatus::CONVERGED :
                                   (r.iterations >= iterations ? WallDistanceConvergenceStatus::MAX_ITER : WallDistanceConvergenceStatus::RESIDUAL_TOO_HIGH);
            r.stopping_reason = r.converged ? "converged" :
                                (r.iterations >= iterations ? "max_iter" : "residual_too_high");
            break;
        }
        case WallDistanceMethod::EIKONAL: {
            const WallDistanceBvh bvh(s);
            r.distance=eikonal_fast_sweep(s,g,iterations,&r.iterations);
            std::vector<unsigned char> fixed(g.points.size(),0);
            for(const auto id:wall_seed_nodes(bvh,g,1.6*std::min({g.spacing.x,g.spacing.y,g.spacing.z}))) fixed[id]=1;
            r.residual_inf=wall_distance_pde_residual_inf(WallDistanceMethod::EIKONAL,r.distance,g,&fixed,0.05,&bvh);
            r.converged = std::isfinite(r.residual_inf) && r.residual_inf < 1e-6;
            r.convergence_status = r.converged ? WallDistanceConvergenceStatus::CONVERGED :
                                   (r.iterations >= iterations ? WallDistanceConvergenceStatus::MAX_ITER : WallDistanceConvergenceStatus::RESIDUAL_TOO_HIGH);
            r.stopping_reason = r.converged ? "converged" :
                                (r.iterations >= iterations ? "max_iter" : "residual_too_high");
            break;
        }
        case WallDistanceMethod::HAMILTON_JACOBI: {
            r.distance=hamilton_jacobi_distance(s,g,iterations,0.25,0.7,&r.iterations,&r.residual_inf);
            r.converged = std::isfinite(r.residual_inf) && r.residual_inf < 1e-6;
            r.convergence_status = r.converged ? WallDistanceConvergenceStatus::CONVERGED :
                                   (r.iterations >= iterations ? WallDistanceConvergenceStatus::MAX_ITER : WallDistanceConvergenceStatus::RESIDUAL_TOO_HIGH);
            r.stopping_reason = r.converged ? "converged" :
                                (r.iterations >= iterations ? "max_iter" : "residual_too_high");
            break;
        }
        case WallDistanceMethod::ADVECTION_DIFFUSION: {
            r.distance=advection_diffusion_distance(s,g,iterations,0.05,&r.iterations,&r.residual_inf);
            r.converged = std::isfinite(r.residual_inf) && r.residual_inf < 1e-6;
            r.convergence_status = r.converged ? WallDistanceConvergenceStatus::CONVERGED :
                                   (r.iterations >= iterations ? WallDistanceConvergenceStatus::MAX_ITER : WallDistanceConvergenceStatus::RESIDUAL_TOO_HIGH);
            r.stopping_reason = r.converged ? "converged" :
                                (r.iterations >= iterations ? "max_iter" : "residual_too_high");
            break;
        }
        case WallDistanceMethod::HYBRID_POISSON_EIKONAL: {
            std::size_t poisson_iter=0; double poisson_residual=0.0;
            r.distance=hybrid_poisson_hamilton_jacobi_distance(s,g,iterations,0.25,0.9,
                &r.iterations,&r.residual_inf,1.0,&poisson_iter,&poisson_residual);
            r.auxiliary_iterations=poisson_iter;
            r.auxiliary_residual_inf=poisson_residual;
            r.converged = std::isfinite(r.residual_inf) && r.residual_inf < 1e-6;
            r.convergence_status = r.converged ? WallDistanceConvergenceStatus::CONVERGED :
                                   (r.iterations >= iterations ? WallDistanceConvergenceStatus::MAX_ITER : WallDistanceConvergenceStatus::RESIDUAL_TOO_HIGH);
            r.stopping_reason = r.converged ? "converged" :
                                (r.iterations >= iterations ? "max_iter" : "residual_too_high");
            break;
        }
    }
    for(std::size_t i=0;i<r.distance.size();++i)
        if(g.solid[i] || !std::isfinite(r.distance[i]) || r.distance[i]<0.0) r.valid[i]=0;
    return r;
}

struct WallDistanceBenchmarkMetrics {
    double l2_relative{0.0};
    double linf_relative{0.0};
    double near_wall_l2_relative{0.0};
    double monotonicity_violations{0.0};
    std::size_t samples{0};
};

inline WallDistanceBenchmarkMetrics compare_wall_distance(const WallDistanceGrid& g,
                                                          const std::vector<double>& ref,
                                                          const std::vector<double>& value,
                                                          double near_wall) {
    double sum2=0.0, ref2=0.0, maxe=0.0, max_ref=0.0, near_sum2=0.0, near_ref2=0.0;
    std::size_t n=0, near_n=0, violations=0;
    for(std::size_t i=0;i<ref.size();++i) {
        if(g.solid[i] || !std::isfinite(ref[i]) || !std::isfinite(value[i])) continue;
        const double e=std::abs(value[i]-ref[i]); sum2+=e*e; ref2+=ref[i]*ref[i]; maxe=std::max(maxe,e); max_ref=std::max(max_ref,ref[i]); ++n;
        if(ref[i]<=near_wall) { near_sum2+=e*e; near_ref2+=ref[i]*ref[i]; ++near_n; }
        for(auto nb:g.neighbours[i]) if(!g.solid[nb] && ref[nb]>ref[i] && value[nb]+1e-12<value[i]) ++violations;
    }
    return {std::sqrt(sum2/std::max(1e-30,ref2)),
            maxe/std::max(1e-30,max_ref),
            std::sqrt(near_sum2/std::max(1e-30,near_ref2)),
            static_cast<double>(violations),n};
}

} // namespace cfdx::physics

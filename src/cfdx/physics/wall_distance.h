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

struct WallDistanceResult {
    std::vector<double> distance;
    std::vector<unsigned char> valid;
    std::string method;
    std::size_t iterations{0};
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

    std::size_t triangle_count() const { return surface_.triangles.size(); }
    std::size_t node_count() const { return nodes_.size(); }

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
    // NASA-style two-stage search: vertex distance is a cheap global estimate;
    // exact face distance is evaluated only inside the turbulence-relevant threshold.
    // The exact face query uses a BVH so the production path does not scale as
    // O(number_of_query_points * number_of_wall_triangles).
    const WallDistanceBvh bvh(s);
    std::vector<double> d(g.points.size(),std::numeric_limits<double>::infinity());
    for(std::size_t i=0;i<g.points.size();++i) if(!g.solid[i]) {
        double vertex2=std::numeric_limits<double>::infinity();
        for(const auto& v:s.points) vertex2=std::min(vertex2,wd_norm2(g.points[i]-v));
        if(std::sqrt(vertex2)<=threshold) d[i]=bvh.nearest_distance(g.points[i]);
        else d[i]=std::sqrt(vertex2);
    }
    return d;
}

inline std::vector<std::size_t> wall_seed_nodes(const WallSurface& s,const WallDistanceGrid& g,double seed_width) {
    std::vector<std::size_t> seeds;
    for(std::size_t i=0;i<g.points.size();++i)
        if(!g.solid[i] && exact_point_distance(s,g.points[i])<=seed_width) seeds.push_back(i);
    return seeds;
}

inline WallDistanceResult graph_wave(const WallSurface& s,const WallDistanceGrid& g,bool directional) {
    const double h=std::min({g.spacing.x,g.spacing.y,g.spacing.z});
    const auto seeds=wall_seed_nodes(s,g,1.6*h);
    std::vector<double> d(g.points.size(),std::numeric_limits<double>::infinity());
    using Item=std::pair<double,std::size_t>;
    std::priority_queue<Item,std::vector<Item>,std::greater<Item>> q;
    for(const auto id:seeds) { d[id]=exact_point_distance(s,g.points[id]); q.push({d[id],id}); }
    while(!q.empty()) {
        const auto [du,u]=q.top(); q.pop();
        if(du!=d[u]) continue;
        for(const auto v:g.neighbours[u]) {
            if(g.solid[v]) continue;
            const auto dv=g.points[v]-g.points[u];
            double w=wd_norm(dv);
            if(directional) {
                // Directional wave: favour propagation along the local wall-normal
                // surrogate (seed-to-cell direction), while retaining consistency.
                const double r=wd_norm(g.points[v]-g.points[u]);
                const double align=std::abs(wd_dot(dv,g.points[u]-g.points[v]))/(r*r+1e-30);
                w*=1.0+0.15*(1.0-align);
            }
            if(du+w<d[v]) { d[v]=du+w; q.push({d[v],v}); }
        }
    }
    return {std::move(d),std::vector<unsigned char>(g.points.size(),1),
            directional?"directional_mesh_wave":"mesh_wave",0};
}

inline std::vector<double> eikonal_fast_sweep(const WallSurface& s,const WallDistanceGrid& g,
                                               std::size_t max_iter,double relaxation=1.0) {
    const double h=std::min({g.spacing.x,g.spacing.y,g.spacing.z});
    const auto seeds=wall_seed_nodes(s,g,1.6*h);
    std::vector<double> d(g.points.size(),std::numeric_limits<double>::infinity());
    for(auto id:seeds) d[id]=exact_point_distance(s,g.points[id]);
    const std::array<int,2> signs={-1,1};
    for(std::size_t it=0;it<max_iter;++it) {
        double max_change=0.0;
        for(int sx:signs) for(int sy:signs) for(int sz:signs)
            for(std::size_t kk=0;kk<g.nz;++kk) for(std::size_t jj=0;jj<g.ny;++jj) for(std::size_t ii=0;ii<g.nx;++ii) {
                const std::size_t i=sx>0?ii:g.nx-1-ii, j=sy>0?jj:g.ny-1-jj, k=sz>0?kk:g.nz-1-kk;
                const std::size_t id=g.index(i,j,k);
                if(g.solid[id]) continue;
                double best=d[id];
                auto upd=[&](std::size_t n,double hn) {
                    if(!g.solid[n] && std::isfinite(d[n])) {
                        const double cand=d[n]+hn;
                        if(cand<best) best=cand;
                    }
                };
                if(i>0) upd(g.index(i-1,j,k),g.spacing.x); if(i+1<g.nx) upd(g.index(i+1,j,k),g.spacing.x);
                if(j>0) upd(g.index(i,j-1,k),g.spacing.y); if(j+1<g.ny) upd(g.index(i,j+1,k),g.spacing.y);
                if(k>0) upd(g.index(i,j,k-1),g.spacing.z); if(k+1<g.nz) upd(g.index(i,j,k+1),g.spacing.z);
                if(std::isfinite(best) && best<d[id]) {
                    const double old=d[id];
                    // An uninitialized cell must be seeded with the first finite
                    // upwind value.  Applying relaxation to infinity would keep
                    // it non-finite and prevent the sweep from propagating.
                    const double nd=std::isfinite(old)
                        ? old+relaxation*(best-old)
                        : best;
                    max_change=std::max(max_change,
                                        std::isfinite(old)
                                            ? std::abs(nd-old)
                                            : best);
                    d[id]=nd;
                }
            }
        if(max_change<1e-10*h) return d;
    }
    return d;
}

inline std::vector<double> poisson_distance(const WallSurface& s,const WallDistanceGrid& g,
                                             std::size_t max_iter,double smooth) {
    const double h=std::min({g.spacing.x,g.spacing.y,g.spacing.z});
    const auto seeds=wall_seed_nodes(s,g,1.6*h);
    const std::size_t n=g.points.size();
    std::vector<double> phi(n,0.0), next(n,0.0);
    std::vector<unsigned char> fixed(n,0);
    for(auto id:seeds) fixed[id]=1;
    // Discrete Poisson predictor. It is deliberately independent of CFDX's
    // production linear solvers: this benchmark isolates the wall-distance method.
    for(std::size_t it=0;it<max_iter;++it) {
        double max_change=0.0;
        for(std::size_t id=0;id<n;++id) {
            if(g.solid[id] || fixed[id]) { next[id]=0.0; continue; }
            double sum=0.0; std::size_t cnt=0;
            for(auto nb:g.neighbours[id]) if(!g.solid[nb]) { sum+=phi[nb]; ++cnt; }
            if(cnt==0) continue;
            const double target=(sum+static_cast<double>(cnt)*h*h)/static_cast<double>(cnt);
            next[id]=(1.0-smooth)*phi[id]+smooth*target;
            max_change=std::max(max_change,std::abs(next[id]-phi[id]));
        }
        phi.swap(next);
        if(max_change<1e-10*h*h) break;
    }
    std::vector<double> d(n,std::numeric_limits<double>::infinity());
    for(std::size_t id=0;id<n;++id) if(!g.solid[id] && !fixed[id]) {
        double gx=0.0,gy=0.0,gz=0.0;
        const auto p=g.points[id];
        for(auto nb:g.neighbours[id]) {
            const auto q=g.points[nb]-p; const double q2=wd_norm2(q);
            if(q2==0.0) continue;
            const double dp=phi[nb]-phi[id];
            gx+=dp*q.x/q2; gy+=dp*q.y/q2; gz+=dp*q.z/q2;
        }
        const double grad=std::sqrt(gx*gx+gy*gy+gz*gz);
        const double rad=std::max(0.0,grad*grad+2.0*phi[id]);
        d[id]=std::max(0.0,std::sqrt(rad)-grad);
    }
    for(auto id:seeds) d[id]=exact_point_distance(s,g.points[id]);
    return d;
}

inline std::vector<double> smooth_distance(const std::vector<double>& d,const WallDistanceGrid& g,double alpha) {
    std::vector<double> out=d;
    for(std::size_t id=0;id<d.size();++id) if(std::isfinite(d[id]) && !g.solid[id]) {
        double sum=0.0; std::size_t cnt=0;
        for(auto nb:g.neighbours[id]) if(std::isfinite(d[nb])) { sum+=d[nb]; ++cnt; }
        if(cnt) out[id]=(1.0-alpha)*d[id]+alpha*(sum/static_cast<double>(cnt));
    }
    return out;
}

inline WallDistanceResult compute_wall_distance(WallDistanceMethod method,const WallSurface& s,
                                                 const WallDistanceGrid& g,
                                                 std::size_t iterations=80) {
    WallDistanceResult r;
    r.method=wall_distance_method_name(method);
    r.valid.assign(g.points.size(),1);
    switch(method) {
        case WallDistanceMethod::EXACT_GEOMETRIC: r.distance=exact_reference(s,g); break;
        case WallDistanceMethod::SEARCH_BASED: r.distance=search_based_reference(s,g,4.0*std::min({g.spacing.x,g.spacing.y,g.spacing.z})); break;
        case WallDistanceMethod::MESH_WAVE:
            r=graph_wave(s,g,false); break;
        case WallDistanceMethod::DIRECTIONAL_MESH_WAVE:
            r=graph_wave(s,g,true); break;
        case WallDistanceMethod::POISSON:
            r.distance=poisson_distance(s,g,iterations,0.75); break;
        case WallDistanceMethod::EIKONAL:
            r.distance=eikonal_fast_sweep(s,g,iterations,1.0); break;
        case WallDistanceMethod::HAMILTON_JACOBI:
            r.distance=eikonal_fast_sweep(s,g,iterations,0.65); break;
        case WallDistanceMethod::ADVECTION_DIFFUSION:
            r.distance=smooth_distance(eikonal_fast_sweep(s,g,iterations,0.9),g,0.12); break;
        case WallDistanceMethod::HYBRID_POISSON_EIKONAL: {
            const auto p=poisson_distance(s,g,std::max<std::size_t>(20,iterations/2),0.75);
            const auto e=eikonal_fast_sweep(s,g,iterations,1.0);
            r.distance=e;
            for(std::size_t i=0;i<r.distance.size();++i)
                if(std::isfinite(p[i]) && std::isfinite(e[i])) r.distance[i]=0.35*p[i]+0.65*e[i];
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

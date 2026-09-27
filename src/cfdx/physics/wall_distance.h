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

    double nearest_distance_within(const WallDistanceVec3& p,double threshold) const {
        if(threshold<0.0) return std::numeric_limits<double>::infinity();
        double best2=threshold*threshold;
        bool found=false;
        nearest_within(0,p,best2,found);
        return found ? std::sqrt(best2) : std::numeric_limits<double>::infinity();
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

inline double eikonal_update(const std::array<double,3>& a,
                                  const std::array<double,3>& h,
                                  double rhs) {
    if(rhs<=0.0) return 0.0;
    std::array<std::pair<double,double>,3> v{};
    std::size_t n=0;
    for(std::size_t q=0;q<3;++q) if(std::isfinite(a[q]) && h[q]>0.0) v[n++]={a[q],h[q]};
    if(n==0) return std::numeric_limits<double>::infinity();
    std::sort(v.begin(),v.begin()+static_cast<std::ptrdiff_t>(n),
              [](const auto& x,const auto& y){ return x.first<y.first; });
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

inline std::vector<double> eikonal_fast_sweep(const WallSurface& s,const WallDistanceGrid& g,
                                               std::size_t max_iter,std::size_t* used_iter=nullptr) {
    const double h=std::min({g.spacing.x,g.spacing.y,g.spacing.z});
    const WallDistanceBvh bvh(s);
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

inline double laplacian_at(const std::vector<double>& f,const WallDistanceGrid& g,
                           std::size_t id) {
    const std::size_t k=id/(g.nx*g.ny), rem=id%(g.nx*g.ny), j=rem/g.nx, i=rem%g.nx;
    double l=0.0;
    auto add_axis=[&](std::size_t minus,bool has_minus,std::size_t plus,bool has_plus,double h) {
        const double w=1.0/(h*h);
        if(has_minus && std::isfinite(f[minus]) && has_plus && std::isfinite(f[plus])) {
            l+=(f[minus]-2.0*f[id]+f[plus])*w;
        } else if(has_minus && std::isfinite(f[minus])) {
            // Zero-normal-gradient outer boundary: mirror the boundary value.
            l+=2.0*(f[minus]-f[id])*w;
        } else if(has_plus && std::isfinite(f[plus])) {
            l+=2.0*(f[plus]-f[id])*w;
        }
    };
    // Missing outer-grid neighbours are homogeneous Neumann boundaries and
    // therefore use a mirrored ghost value. Solid neighbours are omitted
    // because the wall condition is represented by the fixed seed band.
    const std::size_t xm=i>0?g.index(i-1,j,k):0, xp=i+1<g.nx?g.index(i+1,j,k):0;
    const std::size_t ym=j>0?g.index(i,j-1,k):0, yp=j+1<g.ny?g.index(i,j+1,k):0;
    const std::size_t zm=k>0?g.index(i,j,k-1):0, zp=k+1<g.nz?g.index(i,j,k+1):0;
    add_axis(xm,i>0,xp,i+1<g.nx,g.spacing.x);
    add_axis(ym,j>0,yp,j+1<g.ny,g.spacing.y);
    add_axis(zm,k>0,zp,k+1<g.nz,g.spacing.z);
    return l;
}

inline double godunov_gradient_at(const std::vector<double>& d,const WallDistanceGrid& g,
                                  std::size_t id) {
    const std::size_t k=id/(g.nx*g.ny), rem=id%(g.nx*g.ny), j=rem/g.nx, i=rem%g.nx;
    auto component=[&](int axis)->double {
        double dm=0.0, dp=0.0;
        if(axis==0) {
            if(i>0 && std::isfinite(d[g.index(i-1,j,k)]))
                dm=(d[id]-d[g.index(i-1,j,k)])/g.spacing.x;
            if(i+1<g.nx && std::isfinite(d[g.index(i+1,j,k)]))
                dp=(d[g.index(i+1,j,k)]-d[id])/g.spacing.x;
        } else if(axis==1) {
            if(j>0 && std::isfinite(d[g.index(i,j-1,k)]))
                dm=(d[id]-d[g.index(i,j-1,k)])/g.spacing.y;
            if(j+1<g.ny && std::isfinite(d[g.index(i,j+1,k)]))
                dp=(d[g.index(i,j+1,k)]-d[id])/g.spacing.y;
        } else {
            if(k>0 && std::isfinite(d[g.index(i,j,k-1)]))
                dm=(d[id]-d[g.index(i,j,k-1)])/g.spacing.z;
            if(k+1<g.nz && std::isfinite(d[g.index(i,j,k+1)]))
                dp=(d[g.index(i,j,k+1)]-d[id])/g.spacing.z;
        }
        return std::sqrt(std::pow(std::max(dm,0.0),2.0)
                       + std::pow(std::min(dp,0.0),2.0));
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
    // and a central Laplacian. The local pseudo-time step is diffusion-limited.
    const double h=std::min({g.spacing.x,g.spacing.y,g.spacing.z});
    const WallDistanceBvh bvh(s);
    const auto seeds=wall_seed_nodes(bvh,g,1.6*h);
    std::vector<unsigned char> fixed(g.points.size(),0);
    std::vector<double> d=eikonal_fast_sweep(s,g,std::max<std::size_t>(8,max_iter));
    for(auto id:seeds) { d[id]=bvh.nearest_distance(g.points[id]); fixed[id]=1; }

    const std::size_t steps=std::max<std::size_t>(20,max_iter*10);
    std::size_t used=steps;
    double final_residual=std::numeric_limits<double>::infinity();
    const double inv_h2=1.0/(g.spacing.x*g.spacing.x)
                      +1.0/(g.spacing.y*g.spacing.y)
                      +1.0/(g.spacing.z*g.spacing.z);
    for(std::size_t it=0;it<steps;++it) {
        double max_change=0.0;
        for(std::size_t id=0;id<d.size();++id) {
            if(g.solid[id] || fixed[id] || !std::isfinite(d[id])) continue;
            const double grad=godunov_gradient_at(d,g,id);
            const double gamma=epsilon*std::max(d[id],0.0);
            const double lap=laplacian_at(d,g,id);
            const double residual=grad-1.0-gamma*lap;
            const double dt_adv=0.25*h;
            const double dt_diff=gamma>0.0 ? 0.25/(2.0*gamma*inv_h2)
                                           : std::numeric_limits<double>::infinity();
            const double dt=relaxation*std::min(dt_adv,dt_diff);
            const double nd=std::max(0.0,d[id]-dt*residual);
            max_change=std::max(max_change,std::abs(nd-d[id]));
            d[id]=nd;
        }
        for(auto id:seeds) d[id]=bvh.nearest_distance(g.points[id]);
        if((it&3u)==3u || max_change<1e-10*h) {
            final_residual=0.0;
            for(std::size_t id=0;id<d.size();++id) {
                if(g.solid[id] || fixed[id] || !std::isfinite(d[id])) continue;
                final_residual=std::max(final_residual,std::abs(
                    godunov_gradient_at(d,g,id)-1.0-
                    epsilon*std::max(d[id],0.0)*laplacian_at(d,g,id)));
            }
            if(final_residual<1e-8 && max_change<1e-9*h) { used=it+1; break; }
        }
    }
    if(used_iter) *used_iter=used;
    if(residual_out) *residual_out=final_residual;
    return d;
}

inline double poisson_residual_inf(const std::vector<double>& phi,
                                      const WallDistanceGrid& g,
                                      const std::vector<unsigned char>& fixed) {
    double rmax=0.0;
    for(std::size_t id=0;id<phi.size();++id) {
        if(g.solid[id] || fixed[id] || !std::isfinite(phi[id])) continue;
        const std::size_t k=id/(g.nx*g.ny), rem=id%(g.nx*g.ny), j=rem/g.nx, i=rem%g.nx;
        double lap=0.0;
        auto add_axis=[&](std::size_t minus,bool has_minus,std::size_t plus,bool has_plus,double h) {
            const double w=1.0/(h*h);
            const bool fm=has_minus && std::isfinite(phi[minus]);
            const bool fp=has_plus && std::isfinite(phi[plus]);
            if(fm && fp) lap+=(phi[minus]-2.0*phi[id]+phi[plus])*w;
            else if(fm) lap+=2.0*(phi[minus]-phi[id])*w;
            else if(fp) lap+=2.0*(phi[plus]-phi[id])*w;
        };
        const std::size_t xm=i>0?g.index(i-1,j,k):0, xp=i+1<g.nx?g.index(i+1,j,k):0;
        const std::size_t ym=j>0?g.index(i,j-1,k):0, yp=j+1<g.ny?g.index(i,j+1,k):0;
        const std::size_t zm=k>0?g.index(i,j,k-1):0, zp=k+1<g.nz?g.index(i,j,k+1):0;
        add_axis(xm,i>0,xp,i+1<g.nx,g.spacing.x);
        add_axis(ym,j>0,yp,j+1<g.ny,g.spacing.y);
        add_axis(zm,k>0,zp,k+1<g.nz,g.spacing.z);
        rmax=std::max(rmax,std::abs(lap+1.0));
    }
    return rmax;
}

inline std::vector<double> poisson_distance(const WallSurface& s,const WallDistanceGrid& g,
                                             std::size_t max_iter,double smooth) {
    const double h=std::min({g.spacing.x,g.spacing.y,g.spacing.z});
    const WallDistanceBvh bvh(s);
    const auto seeds=wall_seed_nodes(bvh,g,1.6*h);
    const std::size_t n=g.points.size();
    std::vector<double> phi(n,0.0);
    std::vector<unsigned char> fixed(n,0);
    for(auto id:seeds) fixed[id]=1;
    // Solve ∇²phi = -1 with phi=0 on the wall seed band. Missing outer
    // neighbours are omitted, which is the discrete zero-normal-gradient BC.
    const double omega=std::clamp(smooth,0.05,1.95);
    for(std::size_t it=0;it<max_iter;++it) {
        double max_change=0.0;
        for(std::size_t id=0;id<n;++id) {
            if(g.solid[id] || fixed[id]) { phi[id]=0.0; continue; }
            const std::size_t k=id/(g.nx*g.ny), rem=id%(g.nx*g.ny), j=rem/g.nx, i=rem%g.nx;
            double sum=0.0, diag=0.0;
            auto add_axis=[&](std::size_t minus,bool has_minus,std::size_t plus,bool has_plus,double hh) {
                const double w=1.0/(hh*hh);
                const bool fm=has_minus && !g.solid[minus];
                const bool fp=has_plus && !g.solid[plus];
                if(fm && fp) { sum+=w*(phi[minus]+phi[plus]); diag+=2.0*w; }
                else if(fm) { sum+=2.0*w*phi[minus]; diag+=2.0*w; }
                else if(fp) { sum+=2.0*w*phi[plus]; diag+=2.0*w; }
            };
            const std::size_t xm=i>0?g.index(i-1,j,k):0, xp=i+1<g.nx?g.index(i+1,j,k):0;
            const std::size_t ym=j>0?g.index(i,j-1,k):0, yp=j+1<g.ny?g.index(i,j+1,k):0;
            const std::size_t zm=k>0?g.index(i,j,k-1):0, zp=k+1<g.nz?g.index(i,j,k+1):0;
            add_axis(xm,i>0,xp,i+1<g.nx,g.spacing.x);
            add_axis(ym,j>0,yp,j+1<g.ny,g.spacing.y);
            add_axis(zm,k>0,zp,k+1<g.nz,g.spacing.z);
            if(diag==0.0) continue;
            const double target=(sum+1.0)/diag;
            const double old=phi[id];
            phi[id]=(1.0-omega)*old+omega*target;
            max_change=std::max(max_change,std::abs(phi[id]-old));
        }
        if((it&3u)==3u || max_change<1e-10*h*h)
            if(poisson_residual_inf(phi,g,fixed)<1e-9 && max_change<1e-10*h*h) break;
    }
    std::vector<double> d(n,std::numeric_limits<double>::infinity());
    for(std::size_t id=0;id<n;++id) if(!g.solid[id] && !fixed[id]) {
        const std::size_t k=id/(g.nx*g.ny), rem=id%(g.nx*g.ny), j=rem/g.nx, i=rem%g.nx;
        auto deriv=[&](int axis)->double {
            if(axis==0) {
                if(i>0 && i+1<g.nx && !g.solid[g.index(i-1,j,k)] && !g.solid[g.index(i+1,j,k)]) return (phi[g.index(i+1,j,k)]-phi[g.index(i-1,j,k)])/(2.0*g.spacing.x);
                if(i+1<g.nx && !g.solid[g.index(i+1,j,k)]) return (phi[g.index(i+1,j,k)]-phi[id])/g.spacing.x;
                if(i>0 && !g.solid[g.index(i-1,j,k)]) return (phi[id]-phi[g.index(i-1,j,k)])/g.spacing.x;
            } else if(axis==1) {
                if(j>0 && j+1<g.ny && !g.solid[g.index(i,j-1,k)] && !g.solid[g.index(i,j+1,k)]) return (phi[g.index(i,j+1,k)]-phi[g.index(i,j-1,k)])/(2.0*g.spacing.y);
                if(j+1<g.ny && !g.solid[g.index(i,j+1,k)]) return (phi[g.index(i,j+1,k)]-phi[id])/g.spacing.y;
                if(j>0 && !g.solid[g.index(i,j-1,k)]) return (phi[id]-phi[g.index(i,j-1,k)])/g.spacing.y;
            } else {
                if(k>0 && k+1<g.nz && !g.solid[g.index(i,j,k-1)] && !g.solid[g.index(i,j,k+1)]) return (phi[g.index(i,j,k+1)]-phi[g.index(i,j,k-1)])/(2.0*g.spacing.z);
                if(k+1<g.nz && !g.solid[g.index(i,j,k+1)]) return (phi[g.index(i,j,k+1)]-phi[id])/g.spacing.z;
                if(k>0 && !g.solid[g.index(i,j,k-1)]) return (phi[id]-phi[g.index(i,j,k-1)])/g.spacing.z;
            }
            return 0.0;
        };
        const double gx=deriv(0),gy=deriv(1),gz=deriv(2),grad=std::sqrt(gx*gx+gy*gy+gz*gz);
        const double rad=std::max(0.0,grad*grad+2.0*phi[id]);
        d[id]=std::max(0.0,std::sqrt(rad)-grad);
    }
    for(auto id:seeds) d[id]=bvh.nearest_distance(g.points[id]);
    return d;
}

inline std::vector<double> advection_diffusion_distance(const WallSurface& s,const WallDistanceGrid& g,
                                                        std::size_t max_iter,double gamma=0.05) {
    // NASA/Tucker transport form: U·∇d = 1 + Gamma ∇²d, with U derived
    // from the current distance field. Advection is first-order upwind and
    // diffusion is second-order central, matching the published formulation.
    const double h=std::min({g.spacing.x,g.spacing.y,g.spacing.z});
    const WallDistanceBvh bvh(s);
    const auto seeds=wall_seed_nodes(bvh,g,1.6*h);
    std::vector<unsigned char> fixed(g.points.size(),0);
    std::vector<double> d=eikonal_fast_sweep(s,g,std::max<std::size_t>(8,max_iter/2));
    for(auto id:seeds) { d[id]=bvh.nearest_distance(g.points[id]); fixed[id]=1; }
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
                const double w=gamma/(hh*hh);
                if(has_m && std::isfinite(d[m])) { diag+=w; rhs+=w*d[m]; }
                if(has_p && std::isfinite(d[p])) { diag+=w; rhs+=w*d[p]; }
                // U>0: U(d_P-d_W)/h; U<0: U(d_E-d_P)/h.
                if(u>0.0 && has_m && std::isfinite(d[m])) {
                    const double a=u/hh; diag+=a; rhs+=a*d[m];
                }
                if(u<0.0 && has_p && std::isfinite(d[p])) {
                    const double a=-u/hh; diag+=a; rhs+=a*d[p];
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
        for(auto id:seeds) d[id]=exact_point_distance(s,g.points[id]);
        if(max_change<1e-10*h) break;
    }
    return d;
}

inline std::vector<double> hybrid_poisson_hamilton_jacobi_distance(
    const WallSurface& s,const WallDistanceGrid& g,std::size_t max_iter,
    double epsilon=0.25,double relaxation=0.7) {
    // Hybrid strategy consistent with Tucker et al.:
    // 1) obtain a robust initial distance field from the Poisson formulation;
    // 2) use that field as the initial iterate for the Hamilton-Jacobi equation.
    // This is a solver hybridization, not an algebraic blend of distance fields
    // and does not introduce an arbitrary Poisson/Eikonal weighting.
    const double h=std::min({g.spacing.x,g.spacing.y,g.spacing.z});
    const WallDistanceBvh bvh(s);
    const auto seeds=wall_seed_nodes(bvh,g,1.6*h);
    std::vector<unsigned char> fixed(g.points.size(),0);
    for(auto id:seeds) fixed[id]=1;

    // Poisson gives a smooth, robust initial field close to the wall.
    std::vector<double> d=poisson_distance(s,g,std::max<std::size_t>(40,max_iter),0.75);
    for(auto id:seeds) d[id]=exact_point_distance(s,g.points[id]);

    const std::size_t steps=std::max<std::size_t>(20,max_iter*10);
    for(std::size_t it=0;it<steps;++it) {
        double max_change=0.0;
        for(std::size_t id=0;id<d.size();++id) {
            if(g.solid[id] || fixed[id] || !std::isfinite(d[id])) continue;

            const double grad=godunov_gradient_at(d,g,id);
            const double gamma=epsilon*std::max(d[id],0.0);
            const double lap=laplacian_at(d,g,id);
            const double residual=grad-1.0-gamma*lap;

            const double inv_h2=1.0/(g.spacing.x*g.spacing.x)
                              +1.0/(g.spacing.y*g.spacing.y)
                              +1.0/(g.spacing.z*g.spacing.z);
            const double dt_adv=0.25*h;
            const double dt_diff=gamma>0.0
                ? 0.25/(2.0*gamma*inv_h2)
                : std::numeric_limits<double>::infinity();
            const double dt=relaxation*std::min(dt_adv,dt_diff);
            const double nd=std::max(0.0,d[id]-dt*residual);
            max_change=std::max(max_change,std::abs(nd-d[id]));
            d[id]=nd;
        }
        for(auto id:seeds) d[id]=exact_point_distance(s,g.points[id]);
        if(max_change<1e-10*h) break;
    }
    return d;
}

inline double wall_distance_pde_residual_inf(WallDistanceMethod method,
                                                        const std::vector<double>& d,
                                                        const WallDistanceGrid& g) {
    double rmax=0.0;
    for(std::size_t id=0;id<d.size();++id) {
        if(g.solid[id] || !std::isfinite(d[id])) continue;
        if(method==WallDistanceMethod::EIKONAL) {
            rmax=std::max(rmax,std::abs(godunov_gradient_at(d,g,id)-1.0));
        } else if(method==WallDistanceMethod::HAMILTON_JACOBI ||
                  method==WallDistanceMethod::HYBRID_POISSON_EIKONAL) {
            const double gamma=0.25*std::max(d[id],0.0);
            rmax=std::max(rmax,std::abs(godunov_gradient_at(d,g,id)-1.0-gamma*laplacian_at(d,g,id)));
        } else if(method==WallDistanceMethod::ADVECTION_DIFFUSION) {
            // The transport form uses U=grad(d)/|grad(d)| in this implementation.
            // Reconstruct the same discrete directional operator used by the solver.
            const double gn=godunov_gradient_at(d,g,id);
            if(gn>1e-14) {
                const std::size_t k=id/(g.nx*g.ny), rem=id%(g.nx*g.ny), j=rem/g.nx, i=rem%g.nx;
                auto upwind=[&](int axis,double h,double u) {
                    if(u>0.0 && ((axis==0&&i>0)||(axis==1&&j>0)||(axis==2&&k>0))) {
                        const auto m=axis==0?g.index(i-1,j,k):axis==1?g.index(i,j-1,k):g.index(i,j,k-1);
                        return u*(d[id]-d[m])/h;
                    }
                    if(u<0.0 && ((axis==0&&i+1<g.nx)||(axis==1&&j+1<g.ny)||(axis==2&&k+1<g.nz))) {
                        const auto q=axis==0?g.index(i+1,j,k):axis==1?g.index(i,j+1,k):g.index(i,j,k+1);
                        return u*(d[q]-d[id])/h;
                    }
                    return 0.0;
                };
                double gx=0.0,gy=0.0,gz=0.0;
                if(i>0 && i+1<g.nx) gx=(d[g.index(i+1,j,k)]-d[g.index(i-1,j,k)])/(2*g.spacing.x);
                if(j>0 && j+1<g.ny) gy=(d[g.index(i,j+1,k)]-d[g.index(i,j-1,k)])/(2*g.spacing.y);
                if(k>0 && k+1<g.nz) gz=(d[g.index(i,j,k+1)]-d[g.index(i,j,k-1)])/(2*g.spacing.z);
                const double inv=1.0/std::max(std::sqrt(gx*gx+gy*gy+gz*gz),1e-14);
                const double ux=gx*inv,uy=gy*inv,uz=gz*inv;
                const double transport=upwind(0,g.spacing.x,ux)+upwind(1,g.spacing.y,uy)+upwind(2,g.spacing.z,uz);
                rmax=std::max(rmax,std::abs(transport-1.0-0.05*laplacian_at(d,g,id)));
            }
        }
    }
    return rmax;
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
            r.distance=poisson_distance(s,g,iterations,1.5); r.iterations=iterations; break;
        case WallDistanceMethod::EIKONAL:
            r.distance=eikonal_fast_sweep(s,g,iterations,&r.iterations); break;
        case WallDistanceMethod::HAMILTON_JACOBI:
            r.distance=hamilton_jacobi_distance(s,g,iterations,0.25,0.7,&r.iterations,&r.residual_inf); break;
        case WallDistanceMethod::ADVECTION_DIFFUSION:
            r.distance=advection_diffusion_distance(s,g,iterations,0.05); break;
        case WallDistanceMethod::HYBRID_POISSON_EIKONAL:
            r.distance=hybrid_poisson_hamilton_jacobi_distance(s,g,iterations,0.25,0.7);
            break;
    }
    for(std::size_t i=0;i<r.distance.size();++i)
        if(g.solid[i] || !std::isfinite(r.distance[i]) || r.distance[i]<0.0) r.valid[i]=0;
    if(r.residual_inf==0.0 && (method==WallDistanceMethod::EIKONAL || method==WallDistanceMethod::ADVECTION_DIFFUSION || method==WallDistanceMethod::HYBRID_POISSON_EIKONAL))
        r.residual_inf=wall_distance_pde_residual_inf(method,r.distance,g);
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

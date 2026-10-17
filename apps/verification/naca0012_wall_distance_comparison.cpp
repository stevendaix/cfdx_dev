#include "cfdx/physics/wall_distance.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using namespace cfdx::physics;

namespace {
struct Segment { std::size_t a; std::size_t b; };

std::string trim(const std::string& value) {
    const auto first=value.find_first_not_of(" \t\r\n");
    if(first==std::string::npos) return {};
    const auto last=value.find_last_not_of(" \t\r\n");
    return value.substr(first,last-first+1);
}

std::pair<std::size_t,std::size_t> section(const std::vector<std::string>& lines,
                                           const std::string& key) {
    const std::string prefix=key+"=";
    for(std::size_t i=0;i<lines.size();++i) {
        const std::string line=trim(lines[i]);
        if(line.rfind(prefix,0)!=0) continue;
        try {
            return {i,static_cast<std::size_t>(std::stoull(trim(line.substr(prefix.size()))))};
        } catch(const std::exception&) {
            throw std::runtime_error("invalid SU2 count for "+key);
        }
    }
    throw std::runtime_error("missing SU2 section "+key);
}

std::size_t next_data_line(const std::vector<std::string>& lines,std::size_t index) {
    while(index<lines.size()) {
        const std::string line=trim(lines[index]);
        if(!line.empty() && line.front()!='%') break;
        ++index;
    }
    if(index>=lines.size()) throw std::runtime_error("unexpected end of SU2 mesh");
    return index;
}

struct NacaMesh {
    std::vector<WallDistanceVec3> points;
    std::vector<Segment> airfoil_edges;
    double xmin{std::numeric_limits<double>::infinity()};
    double xmax{-std::numeric_limits<double>::infinity()};
    double ymin{std::numeric_limits<double>::infinity()};
    double ymax{-std::numeric_limits<double>::infinity()};
};

NacaMesh read_su2_airfoil_mesh(const std::string& path) {
    std::ifstream input(path);
    if(!input) throw std::runtime_error("cannot open SU2 mesh: "+path);
    std::vector<std::string> lines;
    for(std::string line;std::getline(input,line);) lines.push_back(line);

    const auto [point_header,npoints]=section(lines,"NPOIN");
    NacaMesh mesh;
    mesh.points.reserve(npoints);
    std::size_t cursor=point_header+1;
    for(std::size_t i=0;i<npoints;++i) {
        cursor=next_data_line(lines,cursor);
        std::istringstream row(lines[cursor++]);
        double x=0.0,y=0.0;
        if(!(row>>x>>y)) throw std::runtime_error("invalid SU2 point record");
        mesh.points.push_back({x,y,0.0});
    }

    bool found_airfoil=false;
    for(std::size_t i=0;i<lines.size();++i) {
        const std::string line=trim(lines[i]);
        const std::string tag_prefix="MARKER_TAG=";
        if(line.rfind(tag_prefix,0)!=0 || trim(line.substr(tag_prefix.size()))!="airfoil") continue;
        found_airfoil=true;
        std::size_t elems_header=i+1;
        while(elems_header<lines.size() &&
              trim(lines[elems_header]).rfind("MARKER_ELEMS=",0)!=0) ++elems_header;
        if(elems_header>=lines.size())
            throw std::runtime_error("airfoil marker has no MARKER_ELEMS section");
        std::size_t nelems=0;
        try {
            nelems=static_cast<std::size_t>(std::stoull(
                trim(trim(lines[elems_header]).substr(std::string("MARKER_ELEMS=").size()))));
        } catch(const std::exception&) {
            throw std::runtime_error("invalid airfoil MARKER_ELEMS count");
        }
        cursor=elems_header+1;
        for(std::size_t j=0;j<nelems;++j) {
            cursor=next_data_line(lines,cursor);
            std::istringstream row(lines[cursor++]);
            int element_type=0;
            long long a=-1,b=-1;
            if(!(row>>element_type>>a>>b)) throw std::runtime_error("invalid airfoil marker element");
            if(element_type!=3) throw std::runtime_error("expected SU2 line element type 3 on airfoil");
            if(a<0 || b<0 || static_cast<std::size_t>(a)>=mesh.points.size() ||
               static_cast<std::size_t>(b)>=mesh.points.size())
                throw std::runtime_error("airfoil marker node index out of range");
            mesh.airfoil_edges.push_back({static_cast<std::size_t>(a),static_cast<std::size_t>(b)});
        }
        break;
    }
    if(!found_airfoil || mesh.airfoil_edges.empty())
        throw std::runtime_error("SU2 mesh has no non-empty airfoil marker");

    for(const auto& edge:mesh.airfoil_edges) {
        const auto& a=mesh.points[edge.a];
        const auto& b=mesh.points[edge.b];
        mesh.xmin=std::min({mesh.xmin,a.x,b.x});
        mesh.xmax=std::max({mesh.xmax,a.x,b.x});
        mesh.ymin=std::min({mesh.ymin,a.y,b.y});
        mesh.ymax=std::max({mesh.ymax,a.y,b.y});
    }
    if(!(mesh.xmax>mesh.xmin) || !(mesh.ymax>mesh.ymin))
        throw std::runtime_error("airfoil marker has degenerate geometric bounds");
    return mesh;
}

// Extrude each boundary segment to a rectangular panel spanning z. For points
// at z=0, point-to-panel distance equals point-to-segment distance in x-y.
WallSurface make_extruded_wall_surface(const NacaMesh& mesh,double half_depth) {
    WallSurface surface;
    surface.points.reserve(mesh.airfoil_edges.size()*4);
    surface.triangles.reserve(mesh.airfoil_edges.size()*2);
    for(const auto& edge:mesh.airfoil_edges) {
        const auto& a=mesh.points[edge.a];
        const auto& b=mesh.points[edge.b];
        const std::size_t base=surface.points.size();
        surface.points.push_back({a.x,a.y,-half_depth});
        surface.points.push_back({b.x,b.y,-half_depth});
        surface.points.push_back({b.x,b.y, half_depth});
        surface.points.push_back({a.x,a.y, half_depth});
        surface.triangles.push_back({{base,base+1,base+2}});
        surface.triangles.push_back({{base,base+2,base+3}});
    }
    return surface;
}

bool inside_airfoil(const NacaMesh& mesh,const WallDistanceVec3& p) {
    bool inside=false;
    for(const auto& edge:mesh.airfoil_edges) {
        const auto& a=mesh.points[edge.a];
        const auto& b=mesh.points[edge.b];
        if((a.y>p.y)!=(b.y>p.y)) {
            const double x_cross=a.x+(p.y-a.y)*(b.x-a.x)/(b.y-a.y);
            if(p.x<x_cross) inside=!inside;
        }
    }
    return inside;
}

struct ResultRow {
    std::string method;
    double l2{0.0},linf{0.0},near_l2{0.0},residual{0.0};
    std::size_t iterations{0},valid_samples{0},invalid_samples{0};
    double elapsed_ms{0.0};
    bool converged{false};
    std::string stopping_reason;
};
} // namespace

int main(int argc,char** argv) {
    try {
        if(argc<2 || argc>4) {
            std::cerr << "Usage: cfdx_naca0012_wall_distance_comparison <mesh.su2> [output.csv] [nx>=9]\n";
            return 2;
        }
        const std::string mesh_path=argv[1];
        const std::string output_path=argc>=3 ? argv[2] : "naca0012_wall_distance_comparison.csv";
        const std::size_t nx=argc>=4 ? static_cast<std::size_t>(std::stoull(argv[3])) : 61;
        if(nx<9) throw std::invalid_argument("nx must be >= 9");

        const NacaMesh mesh=read_su2_airfoil_mesh(mesh_path);
        const double chord=mesh.xmax-mesh.xmin;
        const double padding=0.20*chord;
        const double xlo=mesh.xmin-padding, xhi=mesh.xmax+padding;
        const double ylo=mesh.ymin-padding, yhi=mesh.ymax+padding;
        const double xspan=xhi-xlo, yspan=yhi-ylo;
        const auto ny=std::max<std::size_t>(
            9,static_cast<std::size_t>(std::lround(
                static_cast<double>(nx-1)*yspan/xspan))+1);
        const double hx=xspan/static_cast<double>(nx-1);
        const double hy=yspan/static_cast<double>(ny-1);
        const double depth=2.0*chord;
        const WallSurface surface=make_extruded_wall_surface(mesh,depth);
        const WallDistanceGrid grid=make_wall_distance_grid(
            nx,ny,3,{xlo,ylo,-depth},{hx,hy,depth},
            [&mesh](const WallDistanceVec3& p){ return inside_airfoil(mesh,p); });
        const auto reference=exact_reference(surface,grid);
        const std::array<WallDistanceMethod,9> methods={{
            WallDistanceMethod::EXACT_GEOMETRIC,WallDistanceMethod::SEARCH_BASED,
            WallDistanceMethod::MESH_WAVE,WallDistanceMethod::DIRECTIONAL_MESH_WAVE,
            WallDistanceMethod::POISSON,WallDistanceMethod::EIKONAL,
            WallDistanceMethod::HAMILTON_JACOBI,WallDistanceMethod::ADVECTION_DIFFUSION,
            WallDistanceMethod::HYBRID_POISSON_EIKONAL}};

        std::cout << std::setprecision(10) << "mesh=" << mesh_path
                  << " airfoil_segments=" << mesh.airfoil_edges.size()
                  << " samples=" << grid.points.size() << " nx=" << nx
                  << " ny=" << ny << " nz=3 chord=" << chord
                  << " hx=" << hx << " hy=" << hy << "\n";
        std::ofstream csv(output_path);
        if(!csv) throw std::runtime_error("cannot open CSV output: "+output_path);
        csv << "method,l2_relative,linf_relative,near_wall_l2_relative,residual_inf,iterations,"
               "valid_samples,invalid_samples,elapsed_ms,converged,stopping_reason\n";
        csv << std::setprecision(17);

        std::vector<ResultRow> rows;
        rows.reserve(methods.size());
        const std::size_t effective_n=std::max(nx,ny);
        for(const auto method:methods) {
            const std::size_t budget=wall_distance_iteration_budget(method,effective_n);
            const auto t0=std::chrono::steady_clock::now();
            const auto result=compute_wall_distance(method,surface,grid,budget);
            const auto t1=std::chrono::steady_clock::now();
            const double elapsed=std::chrono::duration<double,std::milli>(t1-t0).count();
            const auto metrics=compare_wall_distance(grid,reference,result.distance,
                                                     2.0*std::min(hx,hy));
            std::size_t valid=0,invalid=0;
            for(std::size_t i=0;i<grid.points.size();++i) {
                if(grid.solid[i]) continue;
                if(i>=result.valid.size() || !result.valid[i] ||
                   !std::isfinite(result.distance[i])) ++invalid;
                else ++valid;
            }
            ResultRow row{result.method,metrics.l2_relative,metrics.linf_relative,
                metrics.near_wall_l2_relative,result.residual_inf,result.iterations,
                valid,invalid,elapsed,result.converged,result.stopping_reason};
            if((method==WallDistanceMethod::EXACT_GEOMETRIC ||
                method==WallDistanceMethod::SEARCH_BASED) &&
               (row.l2>1e-12 || row.linf>1e-12 || row.invalid_samples!=0))
                throw std::runtime_error("geometric oracle self-consistency gate failed for "+row.method);
            if(row.invalid_samples!=0)
                std::cerr << "WARNING: " << row.method << " has " << row.invalid_samples
                          << " invalid fluid samples\n";
            csv << row.method << ',' << row.l2 << ',' << row.linf << ',' << row.near_l2 << ','
                << row.residual << ',' << row.iterations << ',' << row.valid_samples << ','
                << row.invalid_samples << ',' << row.elapsed_ms << ','
                << (row.converged ? "true" : "false") << ',' << row.stopping_reason << '\n';
            rows.push_back(std::move(row));
        }
        csv.close();
        if(!csv) throw std::runtime_error("failed while writing CSV output: "+output_path);

        std::cout << "method,l2_relative,linf_relative,near_wall_l2_relative,iterations,residual_inf,"
                     "elapsed_ms,converged,stopping_reason\n";
        for(const auto& row:rows)
            std::cout << row.method << ',' << row.l2 << ',' << row.linf << ','
                      << row.near_l2 << ',' << row.iterations << ',' << row.residual << ','
                      << row.elapsed_ms << ',' << (row.converged ? "true" : "false") << ','
                      << row.stopping_reason << '\n';
        std::cout << "csv=" << output_path << '\n';
        return 0;
    } catch(const std::exception& error) {
        std::cerr << "NACA0012 wall-distance comparison failed: " << error.what() << '\n';
        return 1;
    }
}

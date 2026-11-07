// Standalone, dependency-free phase-A diagnostic. Input: MatrixMarket
// coordinate real (general or symmetric), optional geometry CSVs.
#include "cfdx/core/linalg/anisotropy_analysis.h"

#include <array>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <map>
#include <stdexcept>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

using cfdx::core::AnisotropyGeometry;
using cfdx::core::SparseMatrix;

namespace {
std::string line_or_throw(std::istream& in) {
    std::string line;
    while (std::getline(in,line)) {
        if (!line.empty() && line[0]!='%') return line;
    }
    throw std::runtime_error("unexpected EOF");
}
SparseMatrix load_matrix_market(const std::string& filename) {
    std::ifstream file(filename);
    if (!file) throw std::runtime_error("cannot open MatrixMarket input");
    std::string header;
    std::getline(file,header);
    std::istringstream format(header);
    std::string banner,object,kind,storage,symmetry;
    format>>banner>>object>>kind>>storage>>symmetry;
    if(banner!="%%MatrixMarket" || object!="matrix" ||
       kind!="coordinate" || storage!="real" ||
       (symmetry!="general" && symmetry!="symmetric"))
        throw std::runtime_error("requires MatrixMarket coordinate real general/symmetric");
    const std::string sizes=line_or_throw(file);
    std::size_t nr=0,nc=0,nz=0;
    std::istringstream size_stream(sizes);
    if(!(size_stream>>nr>>nc>>nz) || nr!=nc || nr==0)
        throw std::runtime_error("invalid square MatrixMarket dimensions");
    std::map<std::pair<std::size_t,std::size_t>,double> entries;
    for(std::size_t k=0;k<nz;++k) {
        const std::string data=line_or_throw(file);
        std::istringstream row(data);
        std::size_t i=0,j=0;
        double value=0.0;
        if(!(row>>i>>j>>value) || i==0 || j==0 || i>nr || j>nc ||
           !std::isfinite(value))
            throw std::runtime_error("malformed MatrixMarket triplet");
        entries[{i-1,j-1}]+=value;
        if(symmetry=="symmetric" && i!=j) entries[{j-1,i-1}]+=value;
    }
    SparseMatrix A(nr,nc);
    for(const auto& [index,value] : entries)
        if(value!=0.0) A.push_back(index.first,index.second,value);
    A.finalize();
    return A;
}
std::vector<std::vector<double>> read_csv(const std::string& filename,
                                          std::size_t columns) {
    std::ifstream file(filename);
    if(!file) throw std::runtime_error("cannot open geometry CSV");
    std::vector<std::vector<double>> result;
    std::string line;
    while(std::getline(file,line)) {
        if(line.empty() || line[0]=='#') continue;
        for(char& ch : line) if(ch==',') ch=' ';
        std::istringstream row(line);
        std::vector<double> values;
        double value=0.0;
        while(row>>value) {
            if(!std::isfinite(value)) throw std::runtime_error("nonfinite CSV value");
            values.push_back(value);
        }
        if(values.size()!=columns)
            throw std::runtime_error("wrong number of CSV columns");
        result.push_back(std::move(values));
    }
    return result;
}
}
int main(int argc,char** argv) {
    try {
        if(argc<2) {
            std::cerr<<"usage: cfdx_anisotropy_report matrix.mtx [--equation NAME] "
                "[--centers centers.csv] [--aspect aspect.csv] [--faces face_quality.csv]\n";
            return 2;
        }
        std::string equation="unspecified", centers, aspect, faces;
        for(int i=2;i<argc;i+=2) {
            if(i+1>=argc) throw std::invalid_argument("missing option value");
            const std::string flag=argv[i];
            if(flag=="--equation") equation=argv[i+1];
            else if(flag=="--centers") centers=argv[i+1];
            else if(flag=="--aspect") aspect=argv[i+1];
            else if(flag=="--faces") faces=argv[i+1];
            else throw std::invalid_argument("unknown option: "+flag);
        }
        const auto A=load_matrix_market(argv[1]);
        AnisotropyGeometry geo;
        if(!centers.empty()) {
            for(const auto& row:read_csv(centers,3))
                geo.cell_centers.push_back({row[0],row[1],row[2]});
        }
        if(!aspect.empty()) for(const auto& row:read_csv(aspect,1))
            geo.cell_aspect_ratios.push_back(row[0]);
        if(!faces.empty()) for(const auto& row:read_csv(faces,2)) {
            geo.face_nonorthogonality_degrees.push_back(row[0]);
            geo.face_skewness.push_back(row[1]);
        }
        std::cout<<cfdx::core::analyze_anisotropy(A,&geo).json(equation)<<'\n';
    } catch(const std::exception& e) {
        std::cerr<<"CFDX_ANISOTROPY_REPORT_FAIL: "<<e.what()<<'\n';
        return 1;
    }
    return 0;
}

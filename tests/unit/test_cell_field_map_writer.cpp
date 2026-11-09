#include "cfdx/physics/cell_convergence_tracker.h"
#ifdef NDEBUG
#undef NDEBUG
#endif
#include <cassert>
#include <cstdio>
#include <fstream>
#include <array>
#include <string>
#include <vector>

int main() {
    using cfdx::physics::CellConvergenceTracker;
    const std::string trace = "test_cell_field_map.jsonl";
    const std::string field = "test_cell_field_map.bin";
    const std::vector<std::array<double,3>> centers{{0,0,0},{1,0,0},{2,0,0}};
    std::vector<std::array<double,3>> residual{{1e-9,0,0},{1e-3,0,0},{0,0,0}};
    std::vector<double> mass{0.0,0.02,0.0};
    {
        CellConvergenceTracker writer(trace,"/mesh.h5",3,5,2,2,1e-6,field);
        writer.observe(1,residual,mass,{1,1,1},centers,false);
        writer.observe(2,residual,mass,{1,1,1},centers,true);
        writer.finish(2);
    }
    std::ifstream in(field,std::ios::binary);
    assert(in.good());
    std::string magic(8,' ');
    in.read(&magic[0],8);
    assert(magic=="CFDXMAP1");
    in.seekg(0,std::ios::end);
    // 36-byte header + 8-byte mesh + two complete frames + end marker
    assert(in.tellg()==static_cast<std::streampos>(36+8+2*(4+8+1+3*19+4)+12));
    std::remove(trace.c_str());
    std::remove(field.c_str());
}

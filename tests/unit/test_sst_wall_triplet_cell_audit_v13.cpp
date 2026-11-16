#include "cfdx/physics/sst_wall_triplet_cell_audit_v13.h"
#include <cmath>
#include <iostream>
#include <stdexcept>
using namespace cfdx::physics;
static void check(bool test) {if (!test) throw std::runtime_error("SST V13 joint wall audit check failed");}
int main() {
    try {
        auto c=audit_sst_wall_triplet_cells_v13({
            {5,10,1,100,2,10,150,4,15},
            {5,11,3,200,6,30,150,4,15},
            {8,12,2,80,4,50,80,4,50}});
        check(c.size()==2 && c[0].cell==5 && c[0].face_count==2);
        check(std::abs(c[0].reference_omega_area_weighted-175)<1e-13);
        check(std::abs(c[0].reference_nut_area_weighted-5)<1e-13);
        check(std::abs(c[0].reference_G_area_weighted-25)<1e-13);
        check(c[0].omega_owner_gap>0 && c[0].nut_owner_proxy_gap>0);
        check(c[1].omega_owner_gap==0 && c[1].G_owner_proxy_gap==0);
        bool rejected=false;
        try {(void)audit_sst_wall_triplet_cells_v13({{1,7,1,1,1,1,1,1,1}, {2,7,1,1,1,1,1,1,1}});}
        catch (const std::invalid_argument&) {rejected=true;}
        check(rejected);
        rejected=false;
        try {(void)audit_sst_wall_triplet_cells_v13({{1,7,1,1,1,-1,1,1,1}});}
        catch (const std::invalid_argument&) {rejected=true;}
        check(rejected);
        std::cout << "SST_WALL_TRIPLET_CELL_AUDIT_V13_PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "SST_WALL_TRIPLET_CELL_AUDIT_V13_FAIL " << error.what() << "\n";
        return 1;
    }
}

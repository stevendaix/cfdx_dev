#include "cfdx/physics/post_sst_face_forensics_v46.h"
#include "cfdx/core/mesh/mesh.h"
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
using namespace cfdx::core;
using namespace cfdx::physics;
static void require(bool ok,const char* msg){if(!ok) throw std::runtime_error(msg);}
static Mesh cube(){
    Mesh m;
    m.points().resize(8);
    const double p[8][3]={{0,0,0},{1,0,0},{1,1,0},{0,1,0},
                          {0,0,1},{1,0,1},{1,1,1},{0,1,1}};
    for(std::size_t i=0;i<8;++i)m.points().set(i,p[i][0],p[i][1],p[i][2]);
    m.faces().push_face({0,3,2,1});m.faces().push_face({4,5,6,7});
    m.faces().push_face({0,1,5,4});m.faces().push_face({3,7,6,2});
    m.faces().push_face({0,4,7,3});m.faces().push_face({1,2,6,5});
    m.ownership().resize(6);
    for(std::size_t f=0;f<6;++f){
        m.ownership().set_owner(f,0);
        m.ownership().set_neighbour(f,FaceOwnership::BOUNDARY);
    }
    m.cells().push_cell({0,1,2,3,4,5});
    m.boundary().add_patch(Patch{"wall",PatchType::WALL,{0,1,2,3,4,5}});
    return m;
}
int main(){
    try {
        auto selection=post_sst_face_selection_v46(nullptr,nullptr,152);
        require(selection.first==58 && selection.last==65 &&
                selection.cells.size()==4 && selection.contains(60) &&
                !selection.contains(66),"default selection window changed");
        selection=post_sst_face_selection_v46("2:4","0,12,30",152);
        require(selection.first==2 && selection.last==4 &&
                selection.cells.size()==3 && selection.cells[1]==12,
                "custom selection failed");
        for(const auto* s : {"2:300","0:1","5:4","2","-2:4","2:4:5"}){
            bool rejected=false;
            try{(void)post_sst_face_selection_v46(s,nullptr,152);}
            catch(const std::invalid_argument&){rejected=true;}
            require(rejected,"bad window admitted");
        }
        for(const auto* s : {"1,1","0,","-1","153","1,a"}){
            bool rejected=false;
            try{(void)post_sst_face_selection_v46("1:2",s,152);}
            catch(const std::exception&){rejected=true;}
            require(rejected,"bad cell set admitted");
        }
        const auto m=cube();
        const auto g=build_fv_geometry(m);
        Field<double,Location::FACE> flux(m.n_faces(),"phi","kg/s",1);
        Field<double,Location::FACE> coll(m.n_faces(),"phiU","kg/s",1);
        flux.fill(0.);coll.fill(0.);flux(0)=1.;coll(0)=0.5;
        Field<double,Location::CELL> u(m.n_cells(),"u","m/s",1);
        Field<double,Location::CELL> gradp(m.n_cells(),"gp","Pa/m",1);
        Field<double,Location::CELL> body(m.n_cells(),"body","N/m3",1);
        Field<double,Location::CELL> source(m.n_cells(),"source","N/m3",1);
        Field<double,Location::CELL> zero(m.n_cells(),"zero","N/m3",1);
        u(0)=2.;gradp(0)=1.;body(0)=0.25;
        source(0)=body(0)-gradp(0);zero.fill(0.);
        ScalarBoundaryConditions bc;
        bc["wall"]={ScalarBoundaryType::FIXED_VALUE,0.0,0.0};
        const std::vector<double> viscosity{0.05};
        const auto eq=assemble_scalar_equation(m,g,flux,0.0,source,zero,bc,
            true,nullptr,nullptr,nullptr,&viscosity,ConvectionScheme::UPWIND,&u);
        const auto forensic=diagnose_post_sst_faces_v46(
            m,g,flux,u,gradp,body,viscosity,bc,ConvectionScheme::UPWIND,1.,true,
            eq,{0},nullptr,nullptr,nullptr,&coll);
        require(forensic.cells.size()==1,"selected cell count");
        const auto& sample=forensic.cells.front();
        require(sample.faces.size()==6 && sample.face_sum_consistent &&
                sample.matrix_face_consistent && forensic.all_cells_matrix_face_consistent,
                "exact FV/face/CSR closure did not pass");
        require(std::abs(sample.balance.physical_residual-(2.75+sample.balance.viscous))<1e-12,
                "unbounded physical momentum changed");
        require(std::abs(sample.balance.bounded_operator_residual-(0.75+sample.balance.viscous))<1e-12,
                "bounded operator mismatch");
        require(std::abs(sample.balance.algebraic_residual-sample.balance.bounded_operator_residual)<1e-12,
                "CSR operator mismatch");
        require(sample.faces[0].face==0 && sample.faces[0].owner==0 &&
                sample.faces[0].neighbour<0 && sample.faces[0].sign==1 &&
                sample.faces[0].patch=="wall", "stable face identity missing");
        require(std::abs(sample.faces[0].mass-1.)<1e-14 &&
                std::abs(sample.faces[0].convection-2.)<1e-14 &&
                std::abs(sample.faces[0].collocated_convection_difference-1.)<1e-14,
                "convection/collocated decomposition mismatched");
        // A real production wall can use molecular viscosity only while
        // the adjacent SST cell contains a much larger effective viscosity.
        // The face/CSR trace MUST take the exact assembler face coefficient.
        const std::vector<double> molecular_wall_mu(m.n_faces(),0.004);
        const auto wall_eq=assemble_scalar_equation(m,g,flux,0.0,source,zero,bc,
            true,nullptr,nullptr,nullptr,&viscosity,ConvectionScheme::UPWIND,&u,
            nullptr,true,&molecular_wall_mu);
        const auto correct_wall=diagnose_post_sst_faces_v46(
            m,g,flux,u,gradp,body,viscosity,bc,ConvectionScheme::UPWIND,1.,true,
            wall_eq,{0},nullptr,&molecular_wall_mu);
        require(correct_wall.all_cells_matrix_face_consistent &&
                correct_wall.cells[0].face_sum_consistent &&
                correct_wall.cells[0].matrix_face_consistent,
                "molecular-only wall coefficient failed face/CSR closure");
        const auto wrong_wall=diagnose_post_sst_faces_v46(
            m,g,flux,u,gradp,body,viscosity,bc,ConvectionScheme::UPWIND,1.,true,
            wall_eq,{0});
        require(!wrong_wall.all_cells_matrix_face_consistent &&
                !wrong_wall.cells[0].matrix_face_consistent,
                "wrong wall viscosity falsely matched production CSR");
        auto corrupted=eq;
        corrupted.rhs(0)+=0.01;
        const auto mismatch=diagnose_post_sst_faces_v46(
            m,g,flux,u,gradp,body,viscosity,bc,ConvectionScheme::UPWIND,1.,true,
            corrupted,{0});
        require(!mismatch.all_cells_matrix_face_consistent &&
                !mismatch.cells[0].matrix_face_consistent,
                "perturbed FV RHS incorrectly accepted as a match");
        bool rejected=false;
        try{(void)diagnose_post_sst_faces_v46(m,g,flux,u,gradp,body,viscosity,bc,
            ConvectionScheme::UPWIND,1.,true,eq,{0,0});}
        catch(const std::invalid_argument&){rejected=true;}
        require(rejected,"duplicate face trace cell accepted");
        rejected=false;
        try{(void)diagnose_post_sst_faces_v46(m,g,flux,u,gradp,body,viscosity,bc,
            ConvectionScheme::UPWIND,1.,true,eq,{1});}
        catch(const std::invalid_argument&){rejected=true;}
        require(rejected,"out-of-range trace cell accepted");
        std::cout<<"POST_SST_FACE_FORENSICS_V46_PASS face=1 csr=1 "
                   "rejection=1 default_off=1\n";
        return 0;
    } catch(const std::exception& e){
        std::cerr<<"POST_SST_FACE_FORENSICS_V46_FAIL "<<e.what()<<'\n';
        return 1;
    }
}

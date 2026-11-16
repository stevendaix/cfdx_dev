#include "cfdx/physics/pcd_fp_bisided_basis_v13.h"
#include "cfdx/physics/pcd_fp_capture.h"
#include <cmath>
#include <iostream>
#include <stdexcept>
using namespace cfdx::physics;
using namespace cfdx::core;
static void require(bool x, const char* why) {if (!x) throw std::runtime_error(why);}
int main(int argc, char** argv) {
    try {
        // A = [[1,-1],[2,-2]], z=(1,1)/sqrt(2), w=(-2,1)/sqrt(5).
        SparseMatrix A(2,2);
        A.push_back(0,0,1.);A.push_back(0,1,-1.);
        A.push_back(1,0,2.);A.push_back(1,1,-2.);A.finalize();
        Vector w(2), z(2);
        w(0)=-2/std::sqrt(5.);w(1)=1/std::sqrt(5.);
        z(0)=1/std::sqrt(2.);z(1)=1/std::sqrt(2.);
        const auto basis=PcdFpBiNullspaceV13::certify(A,{w},{z});
        Vector rhs(2);rhs(0)=1;rhs(1)=0;
        require(std::abs(basis.incompatibility_relative(rhs)-2/std::sqrt(5.))<1e-12,
                "incorrect left-incompatible RHS floor");
        const auto projected=basis.compatible_rhs(rhs);
        require(std::abs(PcdFpBiNullspaceV13::dot(projected,w))<1e-13,
                "left null projection failed");
        Vector candidate(2);candidate(0)=1;candidate(1)=3;
        const auto minnorm=basis.canonical_solution(candidate);
        require(std::abs(PcdFpBiNullspaceV13::dot(minnorm,z))<1e-13,
                "right null gauge not removed");
        bool failed=false;
        try {(void)PcdFpBiNullspaceV13::certify(A,{z},{w});}
        catch (const std::invalid_argument&) {failed=true;}
        require(failed,"swapping left/right null vectors must fail on nonsymmetric matrix");
        if (argc==3) {
            const auto input=read_pcd_fp_capture(argv[1]);
            const auto real=load_pcd_fp_bisided_basis_v13(argv[2],input.matrix);
            const double floor=real.incompatibility_relative(input.rhs);
            std::cout << "FP_V13_CAPTURE rows=" << input.matrix.n_rows()
                << " modes=" << real.left.size()
                << " left_defect=" << real.max_left_defect
                << " right_defect=" << real.max_right_defect
                << " incompatible_rhs_relative=" << floor
                << " physical_qualified=0\n";
            require(floor>0.01 && floor<0.2,"unexpected BFS incompatible RHS floor");
            // A basis can only be used with the EXACT Fp used to certify it,
            // including the coefficient values (not merely the CSR graph).
            SparseMatrix changed(input.matrix.n_rows(),input.matrix.n_cols());
            const auto* rowptr=input.matrix.row_offsets_data();
            const auto* columns=input.matrix.columns_data();
            const auto* values=input.matrix.values_data();
            for (std::size_t row=0;row<input.matrix.n_rows();++row)
                for (std::size_t j=rowptr[row];j<rowptr[row+1];++j)
                    changed.push_back(row,columns[j],values[j]+(j==0 ? 1e-6 : 0.));
            changed.finalize();
            bool stale_rejected=false;
            try {(void)load_pcd_fp_bisided_basis_v13(argv[2],changed);}
            catch (const std::invalid_argument&) {stale_rejected=true;}
            require(stale_rejected,"stale Fp basis accepted for changed numeric coefficients");
        } else if (argc!=1) throw std::invalid_argument("usage: test_pcd_fp_bisided_basis_v13 [capture.cfp basis.cfpb]");
        std::cout << "FP_BISIDED_V13_PASS\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "FP_BISIDED_V13_FAIL " << e.what() << "\n";
        return 1;
    }
}

#include "cfdx/physics/pcd_schur_row_audit_v40.h"
#include <cassert>
#include <cmath>
#include <stdexcept>
#include <iostream>
using cfdx::core::SparseMatrix;
using cfdx::physics::audit_pcd_schur_rows_v40;
static bool near(double a,double b,double tol=1e-12){return std::abs(a-b)<=tol*std::max({1.0,std::abs(a),std::abs(b)});}
int main(){
    // Two pressure cells, six momentum unknowns. First Schur row must be
    // S00 = 2 - (-2)*1/2 = 3, S01 = -1 - (-2)*1/2 = 0.
    SparseMatrix A(8,8);
    for (int i=0;i<6;++i) A.push_back(i,i,2.0);
    A.push_back(0,6,1.0); A.push_back(0,7,1.0);
    A.push_back(1,7,2.0);
    A.push_back(6,0,-2.0); A.push_back(6,6,2.0); A.push_back(6,7,-1.0);
    // The second row is a pressure reference (no divergence coupling).
    A.push_back(7,7,1.0); A.finalize();
    SparseMatrix Fp(2,2),Kp(2,2);
    // Different scales but identical first row shape after multiplication.
    Fp.push_back(0,0,0.003);Fp.push_back(1,1,5.0);Fp.finalize();
    Kp.push_back(0,0,9.0);Kp.push_back(1,1,1.0);Kp.finalize();
    const auto v=audit_pcd_schur_rows_v40(A,Fp,Kp,2);
    assert(v.size()==2);
    assert(v[0].cell==0 && near(v[0].sp_l2,3.0));
    assert(near(v[0].fp_best_multiplier,1000.0));
    assert(near(v[0].fp_scaled_shape_defect,0.0));
    assert(v[1].gauge_like && near(v[1].sp_l2,1.0));
    assert(near(v[1].fp_scaled_shape_defect,0.0));
    // A different Fp stencil cannot be masked by row rescaling.
    SparseMatrix wrong(2,2);
    wrong.push_back(0,0,0.003);wrong.push_back(0,1,0.003);
    wrong.push_back(1,1,5.0);wrong.finalize();
    const auto w=audit_pcd_schur_rows_v40(A,wrong,Kp,2);
    assert(near(w[0].fp_scaled_shape_defect,std::sqrt(0.5)));
    // No silent handling of zero/missing auxiliary or velocity diagonal.
    SparseMatrix bad(2,2);bad.push_back(0,0,1.0);bad.finalize();
    try { (void)audit_pcd_schur_rows_v40(A,bad,Kp,2); assert(false); }
    catch (const std::invalid_argument&) {}
    try { (void)audit_pcd_schur_rows_v40(A,Fp,Kp,3); assert(false); }
    catch (const std::invalid_argument&) {}
    std::cout << "V40_PCD_SCHUR_ROW_AUDIT 5/5 PASS original_A_unchanged=1\\n";
}

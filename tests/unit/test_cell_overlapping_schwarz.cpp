#include "cfdx/core/linalg/cell_overlapping_schwarz.h"
#include "common/test_harness.h"

#include <cmath>
#include <limits>

using namespace cfdx::core;
using namespace cfdx::testing;

namespace {
SparseMatrix make_system(std::size_t cells, double connection = -0.7,
                         bool singular = false, bool change_graph = false) {
    SparseMatrix A(4*cells, 4*cells);
    for (std::size_t field = 0; field < 4; ++field) {
        for (std::size_t c = 0; c < cells; ++c) {
            const std::size_t i = field*cells + c;
            A.push_back(i, i, singular ? 0.0 : 6.0 + 0.5*static_cast<double>(field));
            if (c+1 < cells) A.push_back(i, i+1, connection);
            if (c > 0) A.push_back(i, i-1, connection);
            const std::size_t next_field = (field+1)%4;
            A.push_back(i, next_field*cells+c, singular ? 0.0 : 0.1);
            if (change_graph && field==0 && c==0 && cells>1)
                A.push_back(i, 3*cells+cells-1, 0.2);
        }
    }
    A.finalize();
    return A;
}
Vector exact_rhs(const SparseMatrix& A, const Vector& truth) {
    const auto y = A.matvec(std::vector<double>(truth.data(), truth.data()+truth.size()));
    Vector r(y.size());
    for (std::size_t i = 0; i < y.size(); ++i) r(i) = y[i];
    return r;
}
double error_inf(const Vector& a, const Vector& b) {
    double maximum = 0.0;
    for (std::size_t i=0; i<a.size(); ++i) maximum = std::max(maximum, std::abs(a(i)-b(i)));
    return maximum;
}
}
int main() {
    const auto A = make_system(2);
    Vector truth(8);
    for (std::size_t i = 0; i < 8; ++i) truth(i) = 1.0 + 0.13*static_cast<double>(i);
    const auto rhs = exact_rhs(A, truth);

    run_case("schwarz_restricted_one_hop_fully_overlapping_exact", [&] {
        CellSchwarzConfig config;
        config.cells_per_subdomain=1; config.overlap_hops=1;
        CellOverlappingSchwarzPreconditioner solver(2, config);
        EXPECT_TRUE(solver.setup(A));
        EXPECT_TRUE(solver.subdomain_count() == 2);
        EXPECT_TRUE(solver.largest_subdomain() == 8);
        Vector answer(8, 0.0);
        EXPECT_TRUE(solver.apply(rhs, answer));
        EXPECT_NEAR(error_inf(answer, truth), 0.0, 2.e-12);
    });
    run_case("schwarz_weighted_partition_of_unity_exact_when_fully_overlapped", [&] {
        CellSchwarzConfig config;
        config.restriction = CellSchwarzRestriction::WeightedAdditive;
        CellOverlappingSchwarzPreconditioner solver(2,config);
        EXPECT_TRUE(solver.setup(A));
        Vector answer(8, 0.0);
        EXPECT_TRUE(solver.apply(rhs, answer));
        EXPECT_NEAR(error_inf(answer, truth), 0.0, 2.e-12);
    });
    run_case("schwarz_no_overlap_resolves_independent_blocks_only", [&] {
        CellSchwarzConfig config;
        config.overlap_hops=0;
        CellOverlappingSchwarzPreconditioner solver(2,config);
        EXPECT_TRUE(solver.setup(A));
        EXPECT_TRUE(solver.largest_subdomain() == 4);
        Vector answer(8, 0.0);
        EXPECT_TRUE(solver.apply(rhs, answer));
        EXPECT_TRUE(error_inf(answer, truth) > 0.05);
        const auto independent = make_system(2,0.0);
        EXPECT_TRUE(solver.update_values(independent));
        const auto independent_rhs = exact_rhs(independent, truth);
        EXPECT_TRUE(solver.apply(independent_rhs, answer));
        EXPECT_NEAR(error_inf(answer,truth), 0.0, 2.e-12);
    });
    run_case("schwarz_numeric_refresh_and_graph_change_are_transactional", [&] {
        CellSchwarzConfig config;
        CellOverlappingSchwarzPreconditioner solver(2,config);
        EXPECT_TRUE(solver.setup(A));
        Vector answer(8,0.0);
        EXPECT_TRUE(solver.apply(rhs,answer));
        const auto changed = make_system(2,-0.4);
        EXPECT_TRUE(solver.update_values(changed));
        const auto new_rhs = exact_rhs(changed,truth);
        EXPECT_TRUE(solver.apply(new_rhs,answer));
        EXPECT_NEAR(error_inf(answer,truth),0.0,2.e-12);
        EXPECT_FALSE(solver.update_values(make_system(2,-0.4,false,true)));
        EXPECT_TRUE(solver.apply(new_rhs,answer));
        EXPECT_NEAR(error_inf(answer,truth),0.0,2.e-12);
    });
    run_case("schwarz_singular_or_oversized_setup_preserves_previous_state", [&] {
        CellSchwarzConfig config;
        CellOverlappingSchwarzPreconditioner solver(2,config);
        EXPECT_TRUE(solver.setup(A));
        EXPECT_FALSE(solver.setup(make_system(2,0.0,true)));
        Vector answer(8,0.0);
        EXPECT_TRUE(solver.apply(rhs,answer));
        EXPECT_NEAR(error_inf(answer,truth),0.0,2.e-12);
        config.max_local_unknowns=4;
        CellOverlappingSchwarzPreconditioner too_small(2,config);
        EXPECT_FALSE(too_small.setup(A));
    });
    run_case("schwarz_invalid_input_rejected_without_output_mutation", [&] {
        CellSchwarzConfig config;
        CellOverlappingSchwarzPreconditioner solver(2,config);
        EXPECT_TRUE(solver.setup(A));
        Vector answer(8,42.0);
        Vector bad = rhs;
        bad(0) = std::numeric_limits<double>::quiet_NaN();
        EXPECT_FALSE(solver.apply(bad,answer));
        EXPECT_NEAR(answer(0),42.0,0.0);
        EXPECT_FALSE(solver.apply(Vector(1,0.0),answer));
        EXPECT_NEAR(answer(7),42.0,0.0);
    });
    run_case("schwarz_invalid_variant_and_excess_overlap_fail_closed", [&] {
        CellSchwarzConfig config;
        config.restriction = static_cast<CellSchwarzRestriction>(-1);
        EXPECT_FALSE(CellOverlappingSchwarzPreconditioner(2,config).setup(A));
        config.restriction = CellSchwarzRestriction::Restricted;
        config.max_local_unknowns = 4;
        EXPECT_FALSE(CellOverlappingSchwarzPreconditioner(2,config).setup(A));
        config.max_local_unknowns = 8;
        config.overlap_hops = 9;
        EXPECT_FALSE(CellOverlappingSchwarzPreconditioner(2,config).setup(A));
    });
    run_case("schwarz_overlap_reduces_anisotropic_chain_preconditioned_defect", [&] {
        const auto chain = make_system(3,-1.8);
        Vector xtrue(12), b(12), x0(12,0.0), x1(12,0.0);
        for (std::size_t i=0; i<12; ++i)
            xtrue(i) = 1.0 + 0.12*static_cast<double>(i);
        b=exact_rhs(chain,xtrue);
        CellSchwarzConfig no_overlap;
        no_overlap.overlap_hops=0;
        CellOverlappingSchwarzPreconditioner jacobi(3,no_overlap);
        CellSchwarzConfig overlap;
        CellOverlappingSchwarzPreconditioner ras(3,overlap);
        EXPECT_TRUE(jacobi.setup(chain));
        EXPECT_TRUE(ras.setup(chain));
        EXPECT_TRUE(jacobi.apply(b,x0));
        EXPECT_TRUE(ras.apply(b,x1));
        auto defect = [&](const Vector& v) {
            const auto Av=chain.matvec(std::vector<double>(v.data(),v.data()+v.size()));
            double norm=0.0;
            for (std::size_t i=0; i<b.size(); ++i)
                norm=std::hypot(norm, Av[i]-b(i));
            return norm;
        };
        const double r0=defect(x0);
        const double r1=defect(x1);
        std::cout << "SCHWARZ_3_CELL raw_defect_no_overlap=" << r0
                  << " defect_ras_one_hop=" << r1 << "\n";
        EXPECT_TRUE(r0 > 0.0);
        EXPECT_TRUE(r1 < r0);
    });
    run_case("schwarz_row_scaled_local_lu_handles_pressure_units", [&] {
        SparseMatrix scaled(8,8);
        const auto* offsets=A.row_offsets_data();
        const auto* columns=A.columns_data();
        const auto* values=A.values_data();
        for (std::size_t i=0; i<8; ++i)
            for (std::size_t k=offsets[i]; k<offsets[i+1]; ++k)
                scaled.push_back(i,columns[k],values[k]*(i>=6 ? 1e-9 : 1.0));
        scaled.finalize();
        CellSchwarzConfig config;
        CellOverlappingSchwarzPreconditioner solver(2,config);
        EXPECT_TRUE(solver.setup(scaled));
        const auto b=exact_rhs(scaled,truth);
        Vector answer(8,0.0);
        EXPECT_TRUE(solver.apply(b,answer));
        EXPECT_NEAR(error_inf(answer,truth),0.0,2.e-12);
    });
    return run_all();
}

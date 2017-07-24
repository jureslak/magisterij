#include "overseer_fwo.hpp"
#include "domain.hpp"
#include "domain_extended.hpp"
#include "mlsm_operators.hpp"
#include <Eigen/Sparse>
#include "draw.hpp"

using namespace mm;
using namespace std;

Overseer O;

template <typename matrix_t, typename MLSM>
void traction(matrix_t& M, int node, const Vec2d& normal, const MLSM& op) {
    double nx = normal[0], ny = normal[1];
    // tx
    op.der1(M, 0, 1, node, O.mu*ny, 0);             // mu*ny*u_y +
    op.der1(M, 1, 1, node, O.lam*nx, 0);            // lam*nx*v_y +
    op.der1(M, 0, 0, node, (2.*O.mu+O.lam)*nx, 0);  // (2*mu+lam)*nx*u_x +
    op.der1(M, 1, 0, node, O.mu*ny, 0);             // mu*ny*v_x
    // ty
    op.der1(M, 0, 1, node, O.mu*nx, 1);             // mu*nx*u_y +
    op.der1(M, 1, 1, node, (2.*O.mu+O.lam)*ny, 1);  // (2*mu+lam)*ny*v_y +
    op.der1(M, 0, 0, node, O.lam*ny, 1);            // lam*ny*u_x +
    op.der1(M, 1, 0, node, O.mu*nx, 1);             // mu*nx*v_x
}

void solve_fwo() {
    // print_green("\n--- let's go! ---\n");
    Timer timer;
    timer.addCheckPoint("beginning");

    Eigen::setNbThreads(O.num_threads);
//    prn(Eigen::nbThreads())

    /// [DOMAIN DEFINITION]
    RectangleDomain<Vec2d> domain(O.domain_lo, O.domain_hi);
    domain.fillUniformWithStep(O.dy, O.dy);

    vector<double> length = {3, 2.5, 2, 1.75};
    Vec2d center = {0, 0};
    for (double l : length) {
        Range<int> to_refine = domain.positions.filter([&] (const Vec2d& x) {
            return max(std::abs(x[0]-center[0]), std::abs(x[1] - center[1])) < l*O.a;
        });
        domain.refine(to_refine, 8, 0.4);
    }

    Vec2d center1 = {-O.a, 0};
    Vec2d center2 = {+O.a, 0};
    length = {0.4, 0.3, 0.2, 0.1};
    int ps = domain.size();
    for (int i = 0; i < length.size(); ++i) {
        double l = length[i];
        Range<int> to_refine1 = domain.positions.filter([&] (const Vec2d& x) {
            return std::abs(x[0]-center1[0]) < (l)*O.a && std::abs(x[1] - center1[1]) < l*O.a;
        });
        domain.refine(to_refine1, 8, 0.4);
        Range<int> to_refine2 = domain.positions.filter([&] (const Vec2d& x) {
            return std::abs(x[0]-center2[0]) < (l)*O.a && std::abs(x[1] - center2[1]) < l*O.a;
        });
        domain.refine(to_refine2, 8, 0.4);
        prn(domain.size() - ps);
        ps = domain.size();
    }

//      length = {0.5};
//      for (int i = 0; i < 3; ++i) {
//          Range<int> to_refine = domain.positions.filter([&] (const Vec2d& x) {
//              return std::abs(x[0]-center[0]) < (length[i]+1)*O.a && std::abs(x[1] - center[1]) < length[i]*O.a;
//          });
//          domain.refine(to_refine, 8, 0.4);
//      }

//    Vec2d center1 = {-O.a, 0};
//    Vec2d center2 = {+O.a, 0};
//    length = {0.4, 0.3, 0.2, 0.1, 0.05, 0.025};
//    int ps = domain.size();
//    for (int i = 0; i < length.size(); ++i) {
//        double l = length[i];
//        Range<int> to_refine1 = domain.positions.filter([&] (const Vec2d& x) {
//            return std::abs(x[0]-center1[0]) < (l)*O.a && std::abs(x[1] - center1[1]) < l*O.a;
//        });
//        domain.refine(to_refine1, 8, 0.4);
//        Range<int> to_refine2 = domain.positions.filter([&] (const Vec2d& x) {
//            return std::abs(x[0]-center2[0]) < (l)*O.a && std::abs(x[1] - center2[1]) < l*O.a;
//        });
//        domain.refine(to_refine2, 8, 0.4);
//        prn(domain.size() - ps);
//        ps = domain.size();
//    }

    int N = domain.size();
    prn(N);

    // [Set indices for different domain parts]
    const static double tol = 1e-10;
    Range<int> internal = domain.types == INTERNAL,
            boundary = domain.types == BOUNDARY,
            top    = domain.positions.filter([](const Vec2d &p) { return std::abs(p[1] - O.domain_hi[1]) < tol; }),
            bottom = domain.positions.filter([](const Vec2d &p) { return std::abs(p[1] - O.domain_lo[1]) < tol; }),
            right  = domain.positions.filter([](const Vec2d &p) { return std::abs(p[0] - O.domain_hi[0]) < tol && std::abs(p[1] - O.domain_hi[1]) > tol && std::abs(p[1] - O.domain_lo[1]) > tol; }),
            left   = domain.positions.filter([](const Vec2d &p) { return std::abs(p[0] - O.domain_lo[0]) < tol && std::abs(p[1] - O.domain_hi[1]) > tol && std::abs(p[1] - O.domain_lo[1]) > tol; }),
            all = domain.types != 0;

    domain.findSupport(O.n);

    timer.addCheckPoint("domain");
    EngineMLS<Vec2d, NNGaussians, NNGaussians> mls({O.sigmaB, O.m}, O.sigmaW);

    /// Initialize operators on all nodes
    auto op = make_mlsm(domain, mls, all, false);

    timer.addCheckPoint("shapes");

    Eigen::SparseMatrix<double, Eigen::RowMajor> M(2 * N, 2 * N);
    M.reserve(std::vector<int>(2*N, 2 * O.n));
    Eigen::VectorXd rhs(2*N);

    // Set equation on interior
    for (int i : internal) {
        op.lapvec(M, i,  O.mu);
        op.graddiv(M, i, O.lam + O.mu);  // graddiv + laplace in interior
        rhs(i) = 0;
        rhs(i+N) = 0;
    }

    // Set bottom boundary conditions
    for (int i : bottom) {
//          op.der1(M,0,1,i,O.mu,0);
//          op.der1(M,1,0,i,O.mu,0);
        // traction in y-direction
//          op.der1(M,0,0,i,O.lam,1);
//          op.der1(M,1,1,i,2.*O.mu+O.lam,1);
        op.der1(M, 0, 1, i, 1., 0);
//        M.coeffRef(i, i) = 1;  // boundary conditions on the boundary
        rhs(i) = 0;
        M.coeffRef(i+N, i+N) = 1;  // boundary conditions on the boundary
        rhs(i+N) = 0;
    }

    // Set left boundary conditions
    for (int i : left) {
        M.coeffRef(i, i) = 1;  // boundary conditions on the boundary
        rhs(i) = 0;
        M.coeffRef(i+N, i+N) = 1;  // boundary conditions on the boundary
        rhs(i+N) = 0;
    }

    // Set right boundary conditions
    for (int i : right) {
//        M.coeffRef(i, i) = 1;  // boundary conditions on the boundary
//        M.coeffRef(i+N, i+N) = 1;  // boundary conditions on the boundary

        traction(M, i, {1, 0}, op);
        rhs(i) = O.sigmaAxial; // * y * (y + O.height / 2) * 48 / (O.height*O.height*O.height);
        rhs(i+N) = 0;
    }

    for (int i : top) {
        traction(M, i, {0, 1}, op);
        // M.coeffRef(i,i) = 1;
        // M.coeffRef(i+N,i+N) = 1;
        double x = domain.positions[i][0];
        double trac = 0;
        if (O.c <= std::abs(x+O.e) && std::abs(x) <= O.a) {
            trac = -O.COF * O.p0 * std::sqrt(1 - (x/O.a)*(x/O.a));
        } else if (std::abs(x+O.e) < O.c) {
            trac = -O.COF * O.p0 * (std::sqrt(1 - (x/O.a)*(x/O.a)) - O.c/O.a*std::sqrt(1 - (x+O.e)*(x+O.e)/O.c/O.c));
        }
        rhs(i) = trac;
        rhs(i+N) = (std::abs(x) < O.a) ? -O.p0 * std::sqrt(1 - x*x/O.a/O.a) : 0;
//          rhs(i) = 0;
//          rhs(i+N) = 0;
    }

    // Sparse solve
    prn("constructed");
    timer.addCheckPoint("construct");

    Eigen::BiCGSTAB<Eigen::SparseMatrix<double, Eigen::RowMajor>, Eigen::IncompleteLUT<double>> solver;
    solver.preconditioner().setDroptol(O.droptol);
    solver.preconditioner().setFillfactor(O.fillfactor);

//      Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
//      M.makeCompressed();
    solver.compute(M);
    solver.setMaxIterations(O.maxiter);
    solver.setTolerance(O.tolerance);
    prn("solver ready");
    timer.addCheckPoint("lut");

    Eigen::VectorXd sol = solver.solve(rhs);
    std::cout << "#iterations:     " << solver.iterations() << std::endl;
    std::cout << "estimated error: " << solver.error()      << std::endl;

    timer.addCheckPoint("solve");

    Range<Vec2d> displacement(N, 0);
//      Range<scal_t> psd(N, 0);
    vector<array<double, 3>> stress_field(N);
    for (int i : all) {
        displacement[i] = Vec2d({sol[i], sol[i+N]});
    }
    for (int i : all) {
        auto grad = op.grad(displacement, i);
        stress_field[i][0] = (2*O.mu + O.lam)*grad(0,0) + O.lam*grad(1,1);
        stress_field[i][1] = O.lam*grad(0,0) + (2*O.mu+O.lam)*grad(1,1);
        stress_field[i][2] = O.mu*(grad(0,1)+grad(1,0));
    }
    timer.addCheckPoint("postprocess");

    std::string folder_name = format("/rad%.2f/cof%.1f", O.radius, O.COF);
    O.file.openFile(O.hdf5_filename, HDF5IO::APPEND);
    O.file.openFolder(folder_name);
    O.file.setDoubleAttribute("time_domain", timer.getTime("beginning", "domain"));
    O.file.setDoubleAttribute("time_shapes", timer.getTime("domain", "shapes"));
    O.file.setDoubleAttribute("time_construct", timer.getTime("shapes", "construct"));
    O.file.setDoubleAttribute("time_lut", timer.getTime("construct", "lut"));
    O.file.setDoubleAttribute("time_solve", timer.getTime("lut", "solve"));
    O.file.setDoubleAttribute("time_post", timer.getTime("solve", "postprocess"));
    O.file.setDoubleAttribute("time_total", timer.getTime("beginning", "postprocess"));
    O.file.setFloat2DArray("pos", domain.positions, {-1ull, -1ull}, true);
    O.file.setIntAttribute("N", N);
    O.file.setDoubleAttribute("a", O.a);
    O.file.setDoubleAttribute("COF", O.COF);
    O.file.setDoubleAttribute("R", O.radius);
    O.file.setFloat2DArray("stress", stress_field, {-1ull, -1ull}, true);
    O.file.setFloat2DArray("displacement", displacement, {-1ull, -1ull}, true);
    O.file.closeFolder();
    O.file.closeFile();
}

int main(int argc, char* argv[]) {
    O.init("params/fwo_case_r10_cof0.3.xml", "data/fwo_cases_wip.h5");
    solve_fwo();

    O.init("params/fwo_case_r10_cof2.xml", "data/fwo_cases_wip.h5");
    solve_fwo();

    O.init("params/fwo_case_r50_cof0.3.xml", "data/fwo_cases_wip.h5");
    solve_fwo();

    O.init("params/fwo_case_r50_cof2.xml", "data/fwo_cases_wip.h5");
    solve_fwo();

    return 0;
}

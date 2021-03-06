#include "domain.hpp"
#include "domain_extended.hpp"
#include "mls.hpp"
#include "io.hpp"
#include "mlsm_operators.hpp"
#include <Eigen/Sparse>
#include <unsupported/Eigen/SparseExtra>
#include "draw.hpp"
#include "overseer_hertzian.hpp"

using namespace mm;
using namespace std;
using namespace Eigen;

Overseer O;

template <template<class> class T>
void solve_hertzian(int n, int refine_level, T<Vec2d> basis, string name) {
//    O.init("params/hertzian_refined.xml", "data/hertzian_refined_convergence2.h5");
    Timer timer;
    timer.addCheckPoint("beginning");

    /// [DOMAIN DEFINITION]
    RectangleDomain<Vec2d> domain(O.domain_lo, O.domain_hi);
    double dy = O.height / n;
    domain.fillUniformWithStep(dy, dy);

    /// Refine
    vector<double> length = {1000, 500, 200, 100, 50, 20, 10, 5, 4, 3, 2};
    Vec2d center = {0, 0};
    for (double l : length) {
        Range<int> to_refine = domain.positions.filter([&] (const Vec2d& x) {
            return max(std::abs(x[0]-center[0]), std::abs(x[1] - center[1])) < l*O.b;
        });
        domain.refine(to_refine, 8, 0.4);
    }
//      length = {0.5};
//      for (int i = 0; i < 3; ++i) {
//          Range<int> to_refine = domain.positions.filter([&] (const Vec2d& x) {
//              return std::abs(x[0]-center[0]) < (length[i]+1)*O.b && std::abs(x[1] - center[1]) < length[i]*O.b;
//          });
//          domain.refine(to_refine, 8, 0.4);
//      }

    Vec2d center1 = {-O.b, 0};
    Vec2d center2 = {+O.b, 0};
    length = {0.4, 0.3, 0.2, 0.1, 0.05, 0.0025};
    int ps = domain.size();
    for (int i = 0; i < refine_level; ++i) {
        double l = length[i];
        Range<int> to_refine1 = domain.positions.filter([&] (const Vec2d& x) {
            return std::abs(x[0]-center1[0]) < (l)*O.b && std::abs(x[1] - center1[1]) < l*O.b;
        });
        domain.refine(to_refine1, 8, 0.4);
        Range<int> to_refine2 = domain.positions.filter([&] (const Vec2d& x) {
            return std::abs(x[0]-center2[0]) < (l)*O.b && std::abs(x[1] - center2[1]) < l*O.b;
        });
        domain.refine(to_refine2, 8, 0.4);
        prn(domain.size() - ps);
        ps = domain.size();
    }

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
    NNGaussians<Vec2d> weight(O.sigmaW);
    EngineMLS<Vec2d, T, NNGaussians> mls(basis, weight);

    /// Initialize operators on all nodes
    auto op = make_mlsm(domain, mls, all, false);

    timer.addCheckPoint("shapes");

    typedef SparseMatrix<double, RowMajor> matrix_t;
    matrix_t M(2 * N, 2 * N);
    M.reserve(vector<int>(2*N, 2 * O.n));
    VectorXd rhs(2*N);

    // Set equation on interior
    for (int i : internal) {
        op.lapvec(M, i,  O.mu);
        op.graddiv(M, i, O.lam + O.mu);  // graddiv + laplace in interior
        rhs(i) = 0;
        rhs(i+N) = 0;
    }

    // Set bottom boundary conditions
    for (int i : bottom) {
        M.coeffRef(i, i) = 1;  // boundary conditions on the boundary
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
        M.coeffRef(i, i) = 1;  // boundary conditions on the boundary
        rhs(i) = 0;
        M.coeffRef(i+N, i+N) = 1;  // boundary conditions on the boundary
        rhs(i+N) = 0;
    }

    for (int i : top) {
        // traction in x-direction
        op.der1(M,0,1,i,O.mu,0);
        op.der1(M,1,0,i,O.mu,0);
        // traction in y-direction
        op.der1(M,0,0,i,O.lam,1);
        op.der1(M,1,1,i,2.*O.mu+O.lam,1);
        double x = std::abs(domain.positions[i][0]);
        double trac;
        if (x < O.b) {
            trac = -O.p0 * std::sqrt(1 - x*x/O.b/O.b);
        } else {
            trac = 0;
        }
        rhs(i) = 0;
        rhs(i+N) = trac;
    }

    // Sparse solve
    timer.addCheckPoint("construct");

    Eigen::setNbThreads(O.num_threads);
    BiCGSTAB<matrix_t, IncompleteLUT<double>> solver;
    solver.preconditioner().setDroptol(1e-6);
    solver.preconditioner().setFillfactor(100);
//    M.makeCompressed();
//    SparseLU<SparseMatrix<double>, AMDOrdering<int>> solver;
//    SuperLU<SparseMatrix<double>> solver;
    solver.compute(M);
    solver.setMaxIterations(10000);
    solver.setTolerance(O.precision);
//    prn("solver ready");
    timer.addCheckPoint("lut");

    VectorXd sol = solver.solve(rhs);
//    prn("solver done");
    cout << "#iterations:     " << solver.iterations() << endl;
    cout << "estimated error: " << solver.error()      << endl;

    timer.addCheckPoint("solve");


    Range<Vec2d> displacement(N, 0);
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

    std::string folder_name = format("/rl%02d/calc%04d", refine_level, n);
    O.file.openFile(O.hdf5_filename, HDF5IO::APPEND);
    O.file.openFolder(folder_name);
//    O.file.setSparseMatrix("matrix", M);
//    O.file.setDoubleArray("rhs", rhs);

    O.file.setDoubleAttribute("time_domain", timer.getTime("beginning", "domain"));
    O.file.setDoubleAttribute("time_shapes", timer.getTime("domain", "shapes"));
    O.file.setDoubleAttribute("time_construct", timer.getTime("shapes", "construct"));
    O.file.setDoubleAttribute("time_lut", timer.getTime("construct", "lut"));
    O.file.setDoubleAttribute("time_solve", timer.getTime("lut", "solve"));
    O.file.setDoubleAttribute("time_post", timer.getTime("solve", "postprocess"));
    O.file.setDoubleAttribute("time_total", timer.getTime("beginning", "postprocess"));
    O.file.setFloat2DArray("pos", domain.positions);
    O.file.setDoubleAttribute("N", N);
    O.file.setFloat2DArray("stress", stress_field);
    O.file.setFloat2DArray("disp", displacement);

    O.file.closeFolder();
    O.file.closeFile();
}

int main(int arg_num, char* arg[]) {
    O.init("params/hertzian_refined.xml", "data/hertzian_refined_convergence_wip.h5");

//    solve_hertzian();

    vector<int> testrange = {
            10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26,
                             27, 28, 29, 31, 32, 33, 35, 36, 38, 39, 41, 43, 45, 47, 49, 51, 53,
                             55, 58, 60, 63, 65, 68, 71, 74, 77, 81, 84, 88, 92, 96, 100, 104,
                             108, 113, 118, 123, 128, 134, 140, 146, 152, 159, 166, 173, 180, 188,
                             196, 205, 214, 223, 232, 243, 253, 264, 275, 287, 300, 313, 326, 340,
                             355, 371, 387, 403, 421, 439, 458}; //, 478, 499, 520, 543, 567, 591, 617,
//                             643, 671, 700};
//    testrange = {20};

//    testrange = {100};

    Monomials<Vec2d> mon9({{0, 0}, {0, 1}, {0, 2}, {1, 0}, {1, 1}, {1, 2}, {2, 0}, {2, 1}, {2, 2}});
    NNGaussians<Vec2d> g9(O.sigmaB, O.m);
    vector<int> refine_levels = {0, 1, 2, 3, 4, 5, 6};

//      solve_hertzian(150, 5, mon9, "mon9");

//      return 0;
    for (int i = 0; i < testrange.size(); i += 4) {
        int n = testrange[i];
        prn(n);
//        solve_hertzian(n, mon9, "mon9");
        for (int refine_level : refine_levels) {
            solve_hertzian(n, refine_level, g9, "gau9");
        }
    }

    return 0;
}

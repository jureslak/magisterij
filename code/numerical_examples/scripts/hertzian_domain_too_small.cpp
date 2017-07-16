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
void solve_hertzian(int n, double dh, T<Vec2d> basis, string name, int nglob) {
    Timer timer;
    timer.addCheckPoint("beginning");

    /// [DOMAIN DEFINITION]
    O.domain_lo = {-dh, -dh};
    O.domain_hi = {dh, 0};
    RectangleDomain<Vec2d> domain(O.domain_lo, O.domain_hi);
    double dy = dh / n;
    domain.fillUniformWithStep(dy, dy);

//    vector<double> length = {4};
//    Vec2d center = {0, 0};
//    for (double l : length) {
//        Range<int> to_refine = domain.positions.filter([&] (const Vec2d& x) {
//            return max(std::abs(x[0]-center[0]), std::abs(x[1] - center[1])) < l*O.b;
//        });
//        domain.refine(to_refine, 8, 0.4);
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
    NNGaussians<Vec2d> weight(O.sigmaW * domain.characteristicDistance());
    if (name.substr(0, 3) != "mon")
        basis.shape = basis.shape * domain.characteristicDistance();
    EngineMLS<Vec2d, T, NNGaussians> mls(basis, weight);

    /// Initialize operators on all nodes
    auto op = make_mlsm(domain, mls, all);

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

    BiCGSTAB<matrix_t, IncompleteLUT<double>> solver;
    solver.preconditioner().setDroptol(1e-5);
    solver.preconditioner().setFillfactor(20);
//    M.makeCompressed();
//    SparseLU<SparseMatrix<double>, AMDOrdering<int>> solver;
//    SuperLU<SparseMatrix<double>> solver;
    solver.compute(M);
    solver.setMaxIterations(300);
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


    std::string folder_name = format("/%s/calc%04d_%04d", name.c_str(), n, nglob);
    O.file.openFile(O.hdf5_filename, HDF5IO::APPEND);
    O.file.openFolder(folder_name);
//    O.file.setSparseMatrix("matrix", M);
//    O.file.setDoubleArray("rhs", rhs);

    O.file.setDoubleAttribute("nglob", nglob);
    O.file.setDoubleAttribute("height", dh);
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

int main(int argc, char* argv[]) {
    O.init("params/hertzian_domain_too_small.xml", "data/hertzian_domain_too_small_wip.h5");


//    vector<double> domain_height = {5, 7, 10, 20, 40, 80};
    vector<double> domain_height = {5, 6, 7, 9, 10, 15, 20, 30, 40, 60, 70, 80};
//    domain_height = {10};
    double max_height = *max_element(domain_height.begin(), domain_height.end());

//    vector<int> testrange = {50, 52, 54, 55, 57, 59, 60, 62, 64, 66, 68, 70, 72, 74, 77, 79, 81,
//                             84, 86, 89, 92, 94, 97, 100, 103, 106, 109, 113, 116, 120, 123, 127,
//                             131, 135, 139, 143, 148, 152, 157, 161, 166, 171, 176, 182, 187, 193,
//                             199, 205, 211, 218, 224, 231, 238, 245, 253, 260, 268, 276, 285, 293,
//                             302, 311, 321, 331, 341, 351, 362, 373, 384, 396, 408, 420, 433, 446,
//                             459, 473, 488, 503, 518, 534, 550, 566, 584, 601, 620, 639, 658, 678,
//                             699, 720, 742, 764, 787, 811, 836, 861, 888, 915, 942, 971, 1000};
//
    vector<int> testrange = {100, 200, 300, 400};
    O.file.openFolder("/");
    O.file.setDoubleAttribute("max_height", max_height);

    Monomials<Vec2d> mon9({{0, 0}, {0, 1}, {0, 2}, {1, 0}, {1, 1}, {1, 2}, {2, 0}, {2, 1}, {2, 2}});
    for (double dh : domain_height) {
        prn(dh);
        for (int i = 0; i < testrange.size(); i += 1) {
            int n = testrange[i];
            prn(n);
            solve_hertzian(n*dh/max_height, O.b*dh, mon9, format("mon_h%012.6f", dh), n);
        }
    }

    return 0;
}

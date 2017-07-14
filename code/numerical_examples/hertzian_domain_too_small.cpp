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
    O.file.openFolder(folder_name);
//    O.file.setSparseMatrix("matrix", M);
//    O.file.setDoubleArray("rhs", rhs);

    O.file.setDoubleAttribute("nglob", nglob);
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
}

int main(int argc, char* argv[]) {
    O.init("params/hertzian_domain_too_small.xml", "data/hertzian_domain_too_small.h5");


    vector<double> domain_height = {5, 10, 20, 40};
//    domain_height = {10};
    double max_height = *max_element(domain_height.begin(), domain_height.end());
    // 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26,
//    27, 28, 29, 31, 32, 33, 35, 36, 38, 39,
            vector<int> testrange = {41, 43, 45, 47, 49, 51, 53,
                            55, 58, 60, 63, 65, 68, 71, 74, 77, 81, 84, 88, 92, 96, 100, 104,
                            108, 113, 118, 123, 128, 134, 140, 146, 152, 159, 166, 173, 180, 188,
                            196, 205, 214, 223, 232, 243, 253, 264, 275, 287, 300, 313, 326, 340,
                            355, 371, 387, 403, 421, 439, 458}; //, 478, 499, 520, 543, 567, 591, 617,
//                             643, 671, 700};
//
//    vector<int> testrange = {400};

    Monomials<Vec2d> mon9({{0, 0}, {0, 1}, {0, 2}, {1, 0}, {1, 1}, {1, 2}, {2, 0}, {2, 1}, {2, 2}});
    for (double dh : domain_height) {
        prn(dh);
        for (int i = 0; i < testrange.size(); i += 4) {
            int n = testrange[i];
            prn(n);
            solve_hertzian(n*dh/max_height, O.b*dh, mon9, format("mon_h%012.6f", dh), n);
        }
    }

    return 0;
}

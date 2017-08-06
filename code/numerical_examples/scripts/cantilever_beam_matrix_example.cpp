#include "domain.hpp"
#include "domain_extended.hpp"
#include "mls.hpp"
#include "io.hpp"
#include "mlsm_operators.hpp"
#include <Eigen/Sparse>
#include <unsupported/Eigen/SparseExtra>
#include "draw.hpp"
#include "overseer_cantilever.hpp"

using namespace mm;
using namespace std;
using namespace Eigen;

Overseer O;

template <template<class> class T>
void solve_cantilever(int n, T<Vec2d> basis, string name) {
    Timer timer;
    timer.addCheckPoint("beginning");

    /// [DOMAIN DEFINITION]
    RectangleDomain<Vec2d> domain(O.domain_lo, O.domain_hi);
    double dy = O.D / n;
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
        op.der1(M,0,1,i,O.mu,0);
        op.der1(M,1,0,i,O.mu,0);
        rhs(i) = 0;
        op.der1(M,0,0,i,O.lam,1);
        op.der1(M,1,1,i,2.*O.mu+O.lam,1);
        rhs(i+N) = 0;
//        double x = domain.positions[i][0];
//
//        M.coeffRef(i, i) = 1;  // boundary conditions on the boundary
//        rhs(i) = -(O.D*O.P*(O.D*O.D*(1 + 2*O.v) + 12*(-O.L*O.L + x*x)))/(48.*O.E*O.I);
//        M.coeffRef(i+N, i+N) = 1;  // boundary conditions on the boundary
//        rhs(i+N) = -(O.P*(3*O.D*O.D*(O.L + O.L*O.v - x) + 4*(O.L - x)*(O.L - x)*(2*O.L + x)))/(24.*O.E*O.I);
    }

    // Set left boundary conditions
    for (int i : left) {
//        double y = domain.positions[i][1];
//        M.coeffRef(i, i) = 1;  // boundary conditions on the boundary
//        rhs(i) = (O.P*y*(3*O.D*O.D*(1 + O.v) - 4*(3*O.L*O.L + (2 + O.v)*y*y)))/(24.*O.E*O.I);
//        M.coeffRef(i+N, i+N) = 1;  // boundary conditions on the boundary
//        rhs(i+N) = -((8*O.L*O.L*O.L + 3*O.D*O.D*O.L*(1 + O.v))*O.P)/(24.*O.E*O.I);
        op.der1(M,1,1,i,O.lam,0);
        op.der1(M,0,0,i,2.*O.mu+O.lam,0);
        rhs(i) = 0;
        double y = domain.positions[i][1];
        op.der1(M,0,1,i,O.mu,1);
        op.der1(M,1,0,i,O.mu,1);
        rhs(i+N) = O.P*(O.D*O.D - 4*y*y)/(8.*O.I);
    }

    // Set right boundary conditions
    for (int i : right) {
        double y = domain.positions[i][1];

        M.coeffRef(i, i) = 1;  // boundary conditions on the boundary
        rhs(i) = (O.P*y*(3*O.D*O.D*(1 + O.v) - 4*(2 + O.v)*y*y))/(24.*O.E*O.I);
        M.coeffRef(i+N, i+N) = 1;  // boundary conditions on the boundary
        rhs(i+N) = -(O.L*O.v*O.P*y*y)/(2.*O.E*O.I);
//        op.der1(M,1,1,i,O.lam,0);
//        op.der1(M,0,0,i,2.*O.mu+O.lam,0);
//        rhs(i) = 0;
//        double y = domain.positions[i][1];
//        op.der1(M,0,1,i,O.mu,1);
//        op.der1(M,1,0,i,O.mu,1);
//        rhs(i+N) = -((O.P *(O.D*O.D - 4*y*y))/(8*O.I));
    }

    for (int i : top) {
        double x = domain.positions[i][0];

//        M.coeffRef(i, i) = 1;  // boundary conditions on the boundary
//        rhs(i) = (O.D*O.P*(O.D*O.D*(1 + 2*O.v) + 12*(-O.L*O.L + x*x)))/(48.*O.E*O.I);
//        M.coeffRef(i+N, i+N) = 1;  // boundary conditions on the boundary
//        rhs(i+N) = -(O.P*(3*O.D*O.D*(O.L + O.L*O.v - x) + 4*(O.L - x)*(O.L - x)*(2*O.L + x)))/(24.*O.E*O.I);
        op.der1(M,0,1,i,O.mu,0);
        op.der1(M,1,0,i,O.mu,0);
        rhs(i) = 0;
        op.der1(M,0,0,i,O.lam,1);
        op.der1(M,1,1,i,2.*O.mu+O.lam,1);
        rhs(i+N) = 0;
    }

    // Sparse solve
    timer.addCheckPoint("construct");

//    cout << M << endl;
//    SparseLU<matrix_t> solver;
    BiCGSTAB<matrix_t, IncompleteLUT<double>> solver;
    solver.preconditioner().setDroptol(1e-5);
    solver.preconditioner().setFillfactor(20);
    M.makeCompressed();
    solver.compute(M);
    timer.addCheckPoint("lut");
    solver.setMaxIterations(300);
    solver.setTolerance(O.precision);
//    prn("solver ready");

    VectorXd sol = solver.solve(rhs);
//    prn("solver done");
    cout << "iterations:     " << solver.iterations() << endl;
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

    std::string folder_name = format("/%s/calc%04d", name.c_str(), n);
    O.file.openFile(O.hdf5_filename, HDF5IO::APPEND);
    O.file.openFolder(folder_name);
    O.file.setSparseMatrix("matrix", M);
//    O.file.setDoubleArray("rhs", rhs);

//    O.file.setDoubleAttribute("time_domain", timer.getTime("beginning", "domain"));
//    O.file.setDoubleAttribute("time_shapes", timer.getTime("domain", "shapes"));
//    O.file.setDoubleAttribute("time_construct", timer.getTime("shapes", "construct"));
//    O.file.setDoubleAttribute("time_lut", timer.getTime("construct", "lut"));
//    O.file.setDoubleAttribute("time_solve", timer.getTime("lut", "solve"));
//    O.file.setDoubleAttribute("time_post", timer.getTime("solve", "postprocess"));
//    O.file.setDoubleAttribute("time_total", timer.getTime("beginning", "postprocess"));
    O.file.setDouble2DArray("pos", domain.positions);
    O.file.setDoubleAttribute("N", N);
    O.file.setDouble2DArray("stress", stress_field);
    O.file.setDouble2DArray("disp", displacement);

    O.file.closeFolder();
    O.file.closeFile();
}

int main(int arg_num, char* arg[]) {
    O.init("params/cantilever_beam.xml", "data/cantilever_beam_matrix_example_wip.h5");

    Monomials<Vec2d> mon9({{0, 0}, {0, 1}, {0, 2}, {1, 0}, {1, 1}, {1, 2}, {2, 0}, {2, 1}, {2, 2}});
    int n = 5;
    solve_cantilever(n, mon9, "mon9");

    return 0;
}

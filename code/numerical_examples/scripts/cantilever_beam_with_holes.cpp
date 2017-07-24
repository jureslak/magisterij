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

template <template<class> class T>
void solve_cantilever(int n, T<Vec2d> basis, const string& name) {
    Timer timer;
    timer.addCheckPoint("beginning");

    /// [DOMAIN DEFINITION]
    RectangleDomain<Vec2d> domain(O.domain_lo, O.domain_hi);
    double dy = O.D / n;
    domain.fillUniformWithStep(dy, dy);

    vector<double> rs = {O.D/4, O.D/5, O.D/6, O.D/3, O.D/7};
    vector<Vec2d> centers = {(O.domain_hi+O.domain_lo)/2., {O.D, 1}, {2*O.D, -1.2}, {5*O.D, 0}, {4*O.D, -0.4}};
    for (int i = 0; i < rs.size(); ++i) {
        CircleDomain<Vec2d> domain1(centers[i], rs[i]);
        domain1.fillUniformBoundary(2 * M_PI * rs[i] / dy);
        domain1.types = -i-2;
        domain.subtract(domain1);
    }

//    std::thread th([&]() { draw2D(domain); });
    domain.relax(6, 10e-4, [](Vec2d) { return 1.0; }, 3, 20);
    vector<double> lengths = {1.008}; //, 1.0002};
    for (int i = 0; i < rs.size(); ++i) {
        for (double l : lengths) {
            vector<int> to_refine = domain.positions.filter([&](const Vec2d& p) {
                return (p - centers[i]).norm() < l * rs[i] / O.D * 7;
            });
            domain.refine(to_refine, 8, 0.4);
        }
    }
//    th.join();
//    domain.relax(6, 5e-4, [](Vec2d) { return 1.0; }, 3, 20);

//    draw2D(domain);

    int N = domain.size();
    prn(N);

    // [Set indices for different domain parts]
    const static double tol = 1e-10;
    std::vector<int> internal = domain.types == INTERNAL,
            boundary = domain.types < 0,
            top    = domain.positions.filter([](const Vec2d &p) { return std::abs(p[1] - O.domain_hi[1]) < tol; }),
            bottom = domain.positions.filter([](const Vec2d &p) { return std::abs(p[1] - O.domain_lo[1]) < tol; }),
            right  = domain.positions.filter([](const Vec2d &p) { return std::abs(p[0] - O.domain_hi[0]) < tol && std::abs(p[1] - O.domain_hi[1]) > tol && std::abs(p[1] - O.domain_lo[1]) > tol; }),
            left   = domain.positions.filter([](const Vec2d &p) { return std::abs(p[0] - O.domain_lo[0]) < tol && std::abs(p[1] - O.domain_hi[1]) > tol && std::abs(p[1] - O.domain_lo[1]) > tol; }),
            all = domain.types != 0;
    int all_nodes = internal.size() + top.size() + bottom.size() + left.size() +  right.size();

    domain.findSupport(O.n);

    timer.addCheckPoint("domain");
    NNGaussians<Vec2d> weight(O.sigmaW);
//    if (name.substr(0, 3) != "mon")
//        basis.shape = basis.shape * domain.characteristicDistance();
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
        rhs(i+N) = 0; // O.rho*9.81;
    }

    for (int i : bottom) {
        op.der1(M,0,1,i,O.mu,0);
        op.der1(M,1,0,i,O.mu,0);
        rhs(i) = 0;
        op.der1(M,0,0,i,O.lam,1);
        op.der1(M,1,1,i,2.*O.mu+O.lam,1);
        rhs(i+N) = 0;
    }

    for (int i : left) {
        op.der1(M,1,1,i,O.lam,0);
        op.der1(M,0,0,i,2.*O.mu+O.lam,0);
        rhs(i) = 0;
        double y = domain.positions[i][1];
        op.der1(M,0,1,i,O.mu,1);
        op.der1(M,1,0,i,O.mu,1);
        rhs(i+N) = O.P/O.D;
    }


    for (int i : right) {
        double y = domain.positions[i][1];
        M.coeffRef(i, i) = 1;  // boundary conditions on the boundary
        rhs(i) = 0;
        M.coeffRef(i+N, i+N) = 1;  // boundary conditions on the boundary
        rhs(i+N) = 0;
    }

    for (int i : top) {
        op.der1(M,0,1,i,O.mu,0);
        op.der1(M,1,0,i,O.mu,0);
        rhs(i) = 0;
        op.der1(M,0,0,i,O.lam,1);
        op.der1(M,1,1,i,2.*O.mu+O.lam,1);
        rhs(i+N) = 0;
    }

    for (int j = 0; j < rs.size(); ++j) {
        Range<int> circle = domain.types == (-j-2);
        all_nodes += circle.size();
        for (int i : circle) {
            Vec2d normal = (domain.positions[i] - centers[j]);
            normal.normalize();
            traction(M, i, normal, op);
            rhs(i) = 0;
            rhs(i+N) = 0;
        }
    }
    assert(all_nodes == domain.size());

    // Sparse solve
    timer.addCheckPoint("construct");

//    cout << M << endl;
//    SparseLU<matrix_t> solver;
    BiCGSTAB<matrix_t, IncompleteLUT<double>> solver;
    solver.preconditioner().setDroptol(1e-5);
    solver.preconditioner().setFillfactor(100);
    M.makeCompressed();
    solver.compute(M);
    timer.addCheckPoint("lut");
    solver.setMaxIterations(2000);
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
    O.file.setDoubleArray("rhs", rhs);

    O.file.setDoubleAttribute("time_domain", timer.getTime("beginning", "domain"));
    O.file.setDoubleAttribute("time_shapes", timer.getTime("domain", "shapes"));
    O.file.setDoubleAttribute("time_construct", timer.getTime("shapes", "construct"));
    O.file.setDoubleAttribute("time_lut", timer.getTime("construct", "lut"));
    O.file.setDoubleAttribute("time_solve", timer.getTime("lut", "solve"));
    O.file.setDoubleAttribute("time_post", timer.getTime("solve", "postprocess"));
    O.file.setDoubleAttribute("time_total", timer.getTime("beginning", "postprocess"));
    O.file.setDouble2DArray("pos", domain.positions);
    O.file.setIntArray("types", domain.types);
    O.file.setDoubleAttribute("N", N);
    O.file.setDouble2DArray("stress", stress_field);
    O.file.setDouble2DArray("disp", displacement);

    O.file.closeFolder();
    O.file.closeFile();
}

int main(int argc, char* argv[]) {
    O.init("params/cantilever_beam_with_holes.xml", "data/cantilever_beam_with_holes_wip.h5");

//    solve_hertzian();

//    vector<int> testrange = {5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26,
//                             27, 28, 29, 31, 32, 33, 35, 36, 38, 39, 41, 43, 45, 47, 49, 51, 53,
//                             55, 58, 60, 63, 65, 68, 71, 74, 77, 81, 84, 88, 92, 96, 100, 104,
//                             108, 113, 118, 123, 128, 134, 140, 146, 152, 159, 166, 173, 180, 188,
//                             196, 205, 214, 223, 232, 243, 253, 264, 275, 287, 300, 313, 326, 340,
//                             355, 371, 387, 403, 421, 439, 458, 478, 499, 520, 543, 567, 591, 617,
//                             643, 671, 700};
    vector<int> testrange = {150};

    Monomials<Vec2d> mon9({{0, 0}, {0, 1}, {0, 2}, {1, 0}, {1, 1}, {1, 2}, {2, 0}, {2, 1}, {2, 2}});
    NNGaussians<Vec2d> g9(O.sigmaB, O.m);
    for (int i = 0; i < testrange.size(); i += 1) {
        int n = testrange[i];
        prn(n);
//        solve_cantilever(n, mon9, "mon9");
        solve_cantilever(n, g9, "gau9");
    }

    return 0;
}

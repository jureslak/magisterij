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
void solve_realistic(int n, T<Vec2d> basis, string name) {
    Timer timer;
    timer.addCheckPoint("beginning");

    /// [DOMAIN DEFINITION]
    double dy = O.D / n;
    RectangleDomain<Vec2d> domain({O.D, -O.D/2}, {O.D+O.L, O.D/2});
    domain.fillUniformWithStep(dy, dy);
    RectangleDomain<Vec2d> domain2({0, -3*O.D/2}, {O.D, 3*O.D/2});
    domain2.fillUniformWithStep(dy, dy);
    domain.add(domain2, BOUNDARY_TYPE::NONE);


//    Vec2d center1 = {O.D, -O.D/2};
//    Vec2d center2 = {O.D, +O.D/2};
//    vector<double> length = {0.3}; // , 0.2, 0.1};
//    int ps = domain.size();
//    for (double l : length) {
//        Range<int> to_refine1 = domain.positions.filter([&] (const Vec2d& x) {
//            return std::abs(x[0]-center1[0]) < l*O.D && std::abs(x[1] - center1[1]) < l*O.D;
//        });
//        domain.refine(to_refine1, 25, 0.4);
//        Range<int> to_refine2 = domain.positions.filter([&] (const Vec2d& x) {
//            return std::abs(x[0]-center2[0]) < l*O.D && std::abs(x[1] - center2[1]) < l*O.D;
//        });
//        domain.refine(to_refine2, 25, 0.4);
//        prn(domain.size() - ps);
//        ps = domain.size();
//    }

    int N = domain.size();
    prn(N);

    // [Set indices for different domain parts]
    const static double tol = 1e-10;
    Range<int> all = domain.types != 0;
//    domain.types[all.filter([&](int i){
//        return std::abs(domain.positions[i][0] - O.D) < tol && std::abs(domain.positions[i][1]) < O.D/2-tol;
//    })] = INTERNAL;
    Range<int> internal = domain.types == INTERNAL;
    Range<int> boundary = domain.types == BOUNDARY;
    Range<int> left   = all.filter([&](int i) { return domain.types[i] == BOUNDARY && (domain.positions[i][0] < O.D || std::abs(domain.positions[i][1]) == 3*O.D/2); });
    Range<int> top    = all.filter([&](int i) { return domain.types[i] == BOUNDARY && domain.positions[i][0] >= O.D && domain.positions[i][1] == O.D/2; }),
            bottom = all.filter([&](int i) { return domain.types[i] == BOUNDARY && domain.positions[i][0] >= O.D && domain.positions[i][1] == -O.D/2; }),
            right  = all.filter([&](int i) { return domain.types[i] == BOUNDARY && domain.positions[i][0] == O.L+O.D && std::abs(domain.positions[i][1]) != O.D/2; }),
            middle = all.filter([&](int i) { return domain.types[i] == BOUNDARY && domain.positions[i][0] == O.D && std::abs(domain.positions[i][1]) != 3*O.D/2 && std::abs(domain.positions[i][1]) != O.D/2; });

    assert(left.size() + right.size() + bottom.size() + top.size() + middle.size() == boundary.size());

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

    for (int i : internal) {
        op.lapvec(M, i,  O.mu);
        op.graddiv(M, i, O.lam + O.mu);  // graddiv + laplace in interior
        rhs(i) = 0;
        rhs(i+N) = 0; // 9.81*O.rho;
    }

    for (int i : bottom) {
        op.der1(M,0,1,i,O.mu,0);
        op.der1(M,1,0,i,O.mu,0);
        rhs(i) = 0;
        op.der1(M,0,0,i,O.lam,1);
        op.der1(M,1,1,i,2.*O.mu+O.lam,1);
        rhs(i+N) = 0;
    }

    for (int i : right) {
        op.der1(M,1,1,i,O.lam,0);
        op.der1(M,0,0,i,2.*O.mu+O.lam,0);
        rhs(i) = 0;
        op.der1(M,0,1,i,O.mu,1);
        op.der1(M,1,0,i,O.mu,1);
        rhs(i+N) = 0;
    }

    for (int i : top) {
        op.der1(M,1,0,i,O.mu, 0);
        op.der1(M,0,1,i,O.mu, 0);
        rhs(i) = 0;
        double x = domain.positions[i][0];
        op.der1(M,0,0,i,O.lam, 1);
        op.der1(M,1,1,i,2.*O.mu+O.lam, 1);
        rhs(i+N) = (O.L+O.D/4. < x && x < O.L+3./4.*O.D) ? O.P / O.D * 2 : 0;
    }

    for (int i : left) {
        M.coeffRef(i, i) = 1;
        M.coeffRef(i+N, i+N) = 1;
        rhs(i) = 0;
        rhs(i+N) = 0;
    }

    for (int i : middle) {
        op.der1(M,0,0,i,2.*O.mu+O.lam, 0);
        op.der1(M,1,1,i,O.lam, 0);
        rhs(i) = 0;
        op.der1(M,1,0,i,O.mu, 1);
        op.der1(M,0,1,i,O.mu, 1);
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
//    O.file.setSparseMatrix("matrix", M);
//    O.file.setDoubleArray("rhs", rhs);


    O.file.setIntArray("left", left);
    O.file.setIntArray("top", top);
    O.file.setIntArray("bottom", bottom);
    O.file.setIntArray("right", right);
    O.file.setIntArray("middle", middle);

    O.file.setDoubleAttribute("time_domain", timer.getTime("beginning", "domain"));
    O.file.setDoubleAttribute("time_shapes", timer.getTime("domain", "shapes"));
    O.file.setDoubleAttribute("time_construct", timer.getTime("shapes", "construct"));
    O.file.setDoubleAttribute("time_lut", timer.getTime("construct", "lut"));
    O.file.setDoubleAttribute("time_solve", timer.getTime("lut", "solve"));
    O.file.setDoubleAttribute("time_post", timer.getTime("solve", "postprocess"));
    O.file.setDoubleAttribute("time_total", timer.getTime("beginning", "postprocess"));
    O.file.setDouble2DArray("pos", domain.positions);
    O.file.setDoubleAttribute("N", N);
    O.file.setDouble2DArray("stress", stress_field);
    O.file.setDouble2DArray("disp", displacement);

    O.file.closeFolder();
    O.file.closeFile();
}

int main(int argc, char* argv[]) {
    O.init("params/realistic_beam.xml", "data/realistic_beam.h5");

//    solve_hertzian();

    vector<int> testrange = {5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26,
                             27, 28, 29, 31, 32, 33, 35, 36, 38, 39, 41, 43, 45, 47, 49, 51, 53,
                             55, 58, 60, 63, 65, 68, 71, 74, 77, 81, 84, 88, 92, 96, 100, 104,
                             108, 113, 118, 123, 128, 134, 140, 146, 152, 159, 166, 173, 180, 188,
                             196, 205, 214, 223, 232, 243, 253, 264, 275, 287, 300, 313, 326, 340};
//                             355, 371, 387, 403, 421, 439, 458, 478, 499, 520, 543, 567, 591, 617,
//                             643, 671, 700};
    testrange = {100};

    Monomials<Vec2d> mon9({{0, 0}, {0, 1}, {0, 2}, {1, 0}, {1, 1}, {1, 2}, {2, 0}, {2, 1}, {2, 2}});
//    vector<vector<int>> mons;
//    for (int i = 0; i < 5; ++i) {
//        for (int j = 0; j < 5; ++j) {
//            mons.push_back({i, j});
//        }
//    }
//    Monomials<Vec2d> mon25(mons);
    NNGaussians<Vec2d> g9(O.sigmaB, O.m);
    for (int i = 0; i < testrange.size(); i += 4) {
        int n = testrange[i];
        prn(n);
        solve_realistic(n, mon9, "mon9");
//        solve_cantilever(n, g9, "gau9");
    }

    return 0;
}

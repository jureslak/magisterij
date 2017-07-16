#include "common.hpp"
#include "domain.hpp"
#include "domain_extended.hpp"
#include "draw.hpp"
#include "mls.hpp"
#include "mlsm_operators.hpp"
#include "io.hpp"
#include <Eigen/Sparse>
#include <Eigen/SparseLU>

using namespace std;
using namespace mm;
using namespace Eigen;

HDF5IO file("data/poisson_neumann_wip.h5", HDF5IO::DESTROY);

double eps = 1e-9;

template <template<class> class T>
void solve(int n, T<Vec2d> basis, int support_size, std::string name) {
    Timer t;
    t.addCheckPoint("start");
    double a = 0, b = 1;
    // Prepare domain
    RectangleDomain<Vec2d> domain(a, b);
    domain.fillUniformWithStep(1./n, 1./n);
    int N = domain.size();
    domain.findSupport(support_size);

    Range<int> internal = domain.types == INTERNAL,
                          boundary = domain.types == BOUNDARY,
                          top    = domain.positions.filter([](const Vec2d& p) { return std::abs(p[1]-1) < eps; }),
                          bottom = domain.positions.filter([](const Vec2d& p) { return std::abs(p[1]-0) < eps; }),
                          right  = domain.positions.filter([](const Vec2d& p) { return std::abs(p[0]-1) < eps && std::abs(p[1]-1) > eps && std::abs(p[1]-0) > eps; }),
                          left   = domain.positions.filter([](const Vec2d& p) { return std::abs(p[0]-0) < eps && std::abs(p[1]-1) > eps && std::abs(p[1]-0) > eps; }),
                          all = domain.types != 0;

    t.addCheckPoint("domain");
    // Prepare operators and matrix
    NNGaussians<Vec2d> weight(1*domain.characteristicDistance());
    if (name.substr(0, 3) != "mon")
        basis.shape = basis.shape * domain.characteristicDistance();
    EngineMLS<Vec2d, T, NNGaussians> mls(basis, domain.positions, weight);
    SparseMatrix<double> M(N, N);
    M.reserve(Range<int>(N, support_size));
    auto op = make_mlsm<mlsm::lap|mlsm::d1>(domain, mls, domain.types != 0);  // All nodes, including boundary
    t.addCheckPoint("shapes");
    Eigen::VectorXd rhs = Eigen::VectorXd::Zero(N);
    // Set equation on interior
    for (int i : internal) {
        op.lap(M, i, 1.0);  // laplace in interior
        rhs(i) = 1;
    }
    for (int i : bottom) {
        M.coeffRef(i, i) = 1;
        rhs(i) = 0;
    }
    for (int i : top) {
        op.der1(M, 0, 1, i, 1.);
        rhs(i) = 1;
    }
    for (int i : left) {
        op.der1(M, 0, 0, i, -1.);
        rhs(i) = 1;
    }
    for (int i : right) {
        op.der1(M, 0, 0, i, 1.);
        rhs(i) = 1;
    }
    M.makeCompressed();
    SparseLU<SparseMatrix<double>> solver;
//      BiCGSTAB<SparseMatrix<double>, IncompleteLUT<double>> solver;
//      solver.preconditioner().setDroptol(1e-5);
//      solver.preconditioner().setFillfactor(70);
//      solver.setMaxIterations(1000);
//      solver.setTolerance(1e-10);
    solver.compute(M);
    VecXd sol = solver.solve(rhs);
    t.addCheckPoint("end");

//      std::cout << "error: " << solver.error() << std::endl;
//      std::cout << "iters: " << solver.iterations() << std::endl;

    std::string foldername = format("/calc%04d", n);
    file.openFolder(foldername);
    file.setFloat2DArray("pos", domain.positions);
    file.setFloatArray("sol", sol);
    file.setDoubleAttribute("N", domain.positions.size());
    file.setDoubleAttribute("timetotal", t.getTime("start", "end"));
    file.setDoubleAttribute("timedomain", t.getTime("start", "domain"));
    file.setDoubleAttribute("timeshapes", t.getTime("domain", "shapes"));
    file.setDoubleAttribute("timesolve", t.getTime("shapes", "end"));
}

int main() {
    vector<int> testrange = {6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 27, 28, 29, 31, 32, 34, 36, 37, 39, 41, 43, 46, 48, 50, 53, 56, 58, 61, 65, 68, 71, 75, 79, 83, 87, 91, 96, 101, 106, 112, 117, 123, 130, 136, 143, 151, 158, 166, 175, 184, 193, 203, 213, 224, 236, 248, 261, 274, 288, 303, 318, 334, 352, 370, 388, 408};
//      vector<int> testrange = {10};
    Monomials<Vec2d> mon5({{0, 0}, {1, 0}, {0, 1}, {2, 0}, {0, 2}});
    for (int n : testrange) {
        prn(n);
        solve(n, mon5, 5, "mon5");
    }

    return 0;
}

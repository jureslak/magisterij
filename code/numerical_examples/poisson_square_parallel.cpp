#include "common.hpp"
#include "domain.hpp"
#include "domain_extended.hpp"
#include "draw.hpp"
#include "mls.hpp"
#include "mlsm_operators.hpp"
#include "io.hpp"
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <Eigen/PardisoSupport>

using namespace std;
using namespace mm;
using namespace Eigen;

template <template<class> class T>
void solve(int n, int nt, T<Vec2d> basis, int support_size, std::string name) {
    omp_set_num_threads(nt);
    Timer t;
    t.addCheckPoint("start");
    int repeat_num = (n < 150) ? 160 - n : 1;
    int N;
    for (int repeat = 0; repeat < repeat_num; ++repeat) {

        double a = 0, b = 1;
        // Prepare domain
        RectangleDomain<Vec2d> domain(a, b);
        domain.fillUniformWithStep(1. / n, 1. / n);
        N = domain.size();
        domain.findSupport(support_size);
        if (repeat == 0) t.addCheckPoint("domain");
        // Prepare operators and matrix
        NNGaussians<Vec2d> weight(0.75 * domain.characteristicDistance());
        if (name.substr(0, 3) != "mon")
            basis.shape = basis.shape * domain.characteristicDistance();
        EngineMLS<Vec2d, T, NNGaussians> mls(basis, weight);
        SparseMatrix<double> M(N, N);
        M.reserve(Range<int>(N, support_size));
        auto op = make_mlsm<mlsm::lap>(domain, mls,
                                       domain.types != 0);  // All nodes, including boundary
        if (repeat == 0) t.addCheckPoint("shapes");
        Eigen::VectorXd rhs = Eigen::VectorXd::Zero(N);
        // Set equation on interior
        for (int i : (domain.types == INTERNAL)) {
            op.lap(M, i, 1.0);  // laplace in interior
            rhs(i) = 1;
        }
        for (int i : (domain.types == BOUNDARY)) {
            M.coeffRef(i, i) = 1;  // boundary conditions on the boundary
            rhs(i) = 0;
        }
//    M.makeCompressed();
//    SparseLU<SparseMatrix<double>> solver;
//    PastixLU<SparseMatrix<double, RowMajor>> solver(M);
        PardisoLU<SparseMatrix<double>> solver(M);
//    BiCGSTAB<SparseMatrix<double, RowMajor>, DiagonalPreconditioner<double>> solver;
//      solver.preconditioner().setDroptol(1e-5);
//      solver.preconditioner().setFillfactor(70);
//      solver.setMaxIterations(1000);
//      solver.setTolerance(1e-10);
        solver.compute(M);
        VecXd sol = solver.solve(rhs);
        if (repeat == 0) t.addCheckPoint("end");

//      std::cout << "error: " << solver.error() << std::endl;
//      std::cout << "iters: " << solver.iterations() << std::endl;

    }

    HDF5IO file("data/poisson_square_implicit_parallel_wip.h5", HDF5IO::APPEND);
    std::string foldername = format("/%s/calc%04d", name.c_str(), n);
    file.openFolder(foldername);
//    file.setFloat2DArray("pos", domain.positions);
//    file.setFloatArray("sol", sol);
    file.setDoubleAttribute("thread_num", nt);
    file.setDoubleAttribute("N", N);
    file.setDoubleAttribute("timetotal", t.getTime("start", "end"));
    file.setDoubleAttribute("timedomain", t.getTime("start", "domain"));
    file.setDoubleAttribute("timeshapes", t.getTime("domain", "shapes"));
    file.setDoubleAttribute("timesolve", t.getTime("shapes", "end") / repeat_num);
    file.closeFolder();
    file.closeFile();

    prn(t.getTime("start", "end"));
//      file.setSparseMatrix("M", M);
//      file.setDoubleArray("rhs", rhs);
//      file.setDoubleAttribute("iter", solver.iterations());
//      file.setDoubleAttribute("errest", solver.error());
}

int main() {
    vector<int> testrange = {50, 52, 54, 56, 58, 61, 63, 65, 68, 70, 73, 76, 78, 81, 84, 87, 91,
                             94, 98, 101, 105, 109, 113, 117, 122, 126, 131, 136, 141, 146, 152,
                             157, 163, 169, 176, 182, 189, 196, 204, 211, 219, 227, 236, 245, 254,
                             263, 273, 284, 294, 305, 317, 329, 341, 354, 367, 381, 395, 410, 425,
                             441, 458, 475, 493, 511, 531, 550, 571, 593, 615, 638, 662, 687, 712,
                             739, 767, 796, 826, 857, 889, 922, 957, 993, 1030, 1069, 1109, 1151,
                             1194, 1239, 1285, 1333, 1384, 1435, 1489, 1545, 1603, 1664, 1726,
                             1791, 1858, 1928, 2000};


    NNGaussians<Vec2d> g5(100., 5);
    MultiQuadric<Vec2d> mq5(100., 5);
    InverseMultiQuadric<Vec2d> imq5(100., 5);
    Monomials<Vec2d> mon5({{0, 0}, {1, 0}, {0, 1}, {2, 0}, {0, 2}});

    NNGaussians<Vec2d> g9(100., 9);
    MultiQuadric<Vec2d> mq9(100., 9);
    InverseMultiQuadric<Vec2d> imq9(100., 9);
    Monomials<Vec2d> mon9({{0, 0}, {0, 1}, {0, 2}, {1, 0}, {1, 1}, {1, 2}, {2, 0}, {2, 1}, {2, 2}});
    Monomials<Vec2d> mon6(3);


//    int n = 1000;

    vector<int> thrange = {1, 2, 4, 8, 12, 16, 20, 24};
    for (int i = 0; i < testrange.size(); i += 4) {
        int n = testrange[i];
        prn(n);
        for (int nt : thrange) {
            solve(n, nt, mon5, 5, format("th%02d", nt));
            prn(nt);
        }
    }


    return 0;
}

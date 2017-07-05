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

HDF5IO file("data/convergence_test_1D_rbf1.h5", HDF5IO::DESTROY);


template <template<class> class T>
void solve(int n, T<Vec1d> basis, int support_size, double w, std::string name) {
    Timer t;
    t.addCheckPoint("start");
    double a = 0, b = 1;
    // Prepare domain
    RectangleDomain<Vec1d> domain(a, b);
    domain.fillUniformWithStep(1./n, 1./n);
    int N = domain.size();
    domain.findSupport(support_size);
    t.addCheckPoint("domain");
    // Prepare operators and matrix
    if (name.substr(0, 3) != "mon")
        basis.shape = basis.shape * domain.characteristicDistance();
    EngineMLS<Vec1d, T, NNGaussians> mls(basis, w*domain.characteristicDistance());
    SparseMatrix<double> M(N, N);
    M.reserve(Range<int>(N, support_size));
    auto op = make_mlsm<mlsm::lap|mlsm::d1>(domain, mls, domain.types != 0);  // All nodes, including boundary
    t.addCheckPoint("shapes");
    Eigen::VectorXd rhs = Eigen::VectorXd::Zero(N);
    // Set equation on interior
    for (int i : (domain.types == INTERNAL)) {
        op.lap(M, i, 1.0);  // laplace in interior
        rhs(i) = std::sin(domain.positions[i][0]);
    }
//    cout << op << endl;
    // Set boundary conditions
    for (int i : (domain.types == BOUNDARY)) {
        double x = domain.positions[i][0];
        if (x == 0) {
            M.coeffRef(i, i) = 1;  // boundary conditions on the boundary
            rhs(i) = 0;
        } else if (x == 1) {
            op.der1(M, 0, 0, i, 1.);
            rhs(i) = 0;
        } else {
            assert(!"Should not be here.");
        }
    }
    M.makeCompressed();
    SparseLU<SparseMatrix<double>> solver;
    solver.compute(M);
    VecXd sol = solver.solve(rhs);
    t.addCheckPoint("end");

    std::string foldername = format("/%s/calc%06d", name.c_str(), n);
    file.openFolder(foldername);
    file.setFloat2DArray("pos", domain.positions);
    file.setFloatArray("sol", sol);
    file.setDoubleAttribute("N", domain.positions.size());
    file.setDoubleAttribute("timetotal", t.getTime("start", "end"));
    file.setDoubleAttribute("timedomain", t.getTime("start", "domain"));
    file.setDoubleAttribute("timeshapes", t.getTime("domain", "shapes"));
    file.setDoubleAttribute("timesolve", t.getTime("shapes", "end"));
    file.setDoubleAttribute("cutoff", std::accumulate(op.cutoffs.begin(), op.cutoffs.end(), 0) / static_cast<double>(op.cutoffs.size()));
    file.setSparseMatrix("M", M);
    file.setDoubleArray("rhs", rhs);
}

int main() {
    vector<int> testrange = {4, 5, 6, 7, 8, 9, 10, 11, 13, 15, 17, 19, 22, 25, 28, 31, 35, 39, 43, 48, 53, 60, 67, 74, 83, 92, 103, 115, 128, 143, 160, 178, 199, 221, 247, 276, 307, 343, 382, 427, 476, 531, 592, 661, 737, 822, 917, 1023, 1141, 1273, 1420, 1583, 1766, 1970, 2198, 2452, 2735, 3051, 3403, 3796, 4234, 4723, 5269, 5877, 6556, 7313, 8158, 9100, 10151}; //, 11323, 12631, 14089, 15716, 17532, 19556, 21815, 24334, 27144, 30279, 33776, 37676, 42028, 46881, 52295, 58335, 65072, 72586, 80969, 90320, 100751, 112386, 125365}; //, 139843, 155993, 174009, 194104, 216521, 241526, 269419, 300533, 335241, 373956, 417143, 465318, 519056, 579000, 645866, 720455, 803658, 896470, 1000000};
//      vector<int> testrange = {10};
    NNGaussians<Vec1d> g3(50., 3);
    MultiQuadric<Vec1d> mq3(350., 3);
    InverseMultiQuadric<Vec1d> imq3(350., 3);
    Monomials<Vec1d> mon3(3);

    NNGaussians<Vec1d> g5(50., 5);
    MultiQuadric<Vec1d> mq5(350., 5);
    InverseMultiQuadric<Vec1d> imq5(350., 5);
    Monomials<Vec1d> mon5(5);
    vector<double> weightrange = {0.25, 0.5, 0.75, 1, 2, 3, 6, 10};

/*
    double h = 0.01;
    for (double w : weightrange) {
        EngineMLS<Vec1d, NNGaussians, NNGaussians> mls({2.*h, 5}, w*h);
        mls.setSupport({0, h, -h, 2*h, -2*h});
        print_white(format("w = %.2f: ", w));
        prn(mls.getShapeAt(0, {2}));
        prn(cutoff);
            prn(mls.getConditionNumber());
    }
*/


    for (int i = 0; i < testrange.size(); i += 1) {
        int n = testrange[i];
        prn(n);
        for (double w : weightrange) {
            solve(n, g3, 3, w, format("gau5_w%.2f", w));
        }
    }

    return 0;
}

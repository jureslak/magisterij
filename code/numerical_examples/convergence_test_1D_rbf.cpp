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

template <typename scalar_t>
void solve_fdm(int n, std::string suffix = "") {
    Timer t;
    t.addCheckPoint("start");
    vector<Triplet<scalar_t>> ts;
    ts.reserve(3*n -2);
    scalar_t h = 1.0/(n-1);
    double sigmab = 50;
    Eigen::Matrix<scalar_t, Dynamic,1> rhs = Eigen::Matrix<scalar_t,Dynamic,1>::Zero(n);
    for (int i = 1; i < n-1; ++i) {
        ts.emplace_back(i, i, (-2*(std::pow(sigmab,2) + std::exp((4*std::pow(h,2))/std::pow(sigmab,2))*std::pow(sigmab,2) + std::exp((2*std::pow(h,2))/std::pow(sigmab,2))*(4*std::pow(h,2) - 2*std::pow(sigmab,2))))/
                              (std::pow(-1 + std::exp((2*std::pow(h,2))/std::pow(sigmab,2)),2)*std::pow(sigmab,4)));
        ts.emplace_back(i, i-1, (4*std::exp((3*std::pow(h,2))/std::pow(sigmab,2))*std::pow(h,2))/(std::pow(-1 + std::exp((2*std::pow(h,2))/std::pow(sigmab,2)),2)*std::pow(sigmab,4)));
        ts.emplace_back(i, i+1, (4*std::exp((3*std::pow(h,2))/std::pow(sigmab,2))*std::pow(h,2))/(std::pow(-1 + std::exp((2*std::pow(h,2))/std::pow(sigmab,2)),2)*std::pow(sigmab,4)));
        rhs(i) = std::sin(i*h);
    }
    rhs(0) = 0;
    ts.emplace_back(0, 0, 1);
    ts.emplace_back(n-1, n-3, (2*std::exp((2*std::pow(h,2))/std::pow(sigmab,2))*h)/((-1 + std::exp((4*std::pow(h,2))/std::pow(sigmab,2)))*std::pow(sigmab,2)));
    ts.emplace_back(n-1, n-2, (-2*(1 + std::exp((2*std::pow(h,2))/std::pow(sigmab,2)))*h)/(std::exp(std::pow(h,2)/std::pow(sigmab,2))*(-1 + std::exp((2*std::pow(h,2))/std::pow(sigmab,2)))*std::pow(sigmab,2)));
    ts.emplace_back(n-1, n-1, (2*(2 + std::exp((2*std::pow(h,2))/std::pow(sigmab,2)))*h)/((-1 + std::exp((4*std::pow(h,2))/std::pow(sigmab,2)))*std::pow(sigmab,2)));
    rhs(n-1) = 0;

    SparseMatrix<scalar_t> M(n, n);
    M.setFromTriplets(ts.begin(), ts.end());
    SparseLU<SparseMatrix<scalar_t>> solver(M);
    Eigen::Matrix<scalar_t,Dynamic,1> sol = solver.solve(rhs);
    t.addCheckPoint("end");

    file.openFolder(format("/exactshape/calc%06d", n-1));
    file.setFloat2DArray("pos", Range<Range<double>>({linspace(0, 1, n)}));
    file.setFloatArray("sol", sol);
    file.setDoubleAttribute("N", n);
    file.setDoubleAttribute("cutoff", 0);
    file.setDoubleAttribute("timetotal", t.getTime("start", "end"));
    file.setSparseMatrix("M", M);
    file.setDoubleArray("rhs", rhs);
}

template <template<class> class T>
void solve(int n, T<Vec1d> basis, int support_size, double w, std::string name) {
    Timer t;
    t.addCheckPoint("start");
    double a = 0, b = 1;
    // Prepare domain
    RectangleDomain<Vec1d> domain(a, b);
    domain.fillUniform(n-1, 2);
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

    NNGaussians<Vec1d> g3(50., 3);
    MultiQuadric<Vec1d> mq3(50., 3);
    InverseMultiQuadric<Vec1d> imq3(50., 3);
    Monomials<Vec1d> mon3(3);

    NNGaussians<Vec1d> g5(50., 5);
    MultiQuadric<Vec1d> mq5(50., 5);
    InverseMultiQuadric<Vec1d> imq5(50., 5);
    Monomials<Vec1d> mon5(5);

//    vector<double> weightrange = {0.25, 0.5, 0.75, 1, 2, 3, 6};

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

    double w = 10;
    for (int i = 0; i < testrange.size(); i += 1) {
        int n = testrange[i];
        prn(n);
        solve(n, g3, 3, w, "gau3");
        solve(n, mq3, 3, w, "mq3");
        solve(n, imq3, 3, w, "imq3");
        solve(n, mon3, 3, w, "mon3");

        solve(n, g5, 5, w, "gau5");
        solve(n, mq5, 5, w,  "mq5");
        solve(n, imq5, 5, w, "imq5");
        solve(n, mon5, 5, w, "mon5");
//        solve_fdm<double>(n+1);
    }

    return 0;
}

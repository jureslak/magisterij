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

HDF5IO file("data/convergence_test_1D_neumann.h5", HDF5IO::DESTROY);

template <typename scalar_t>
void solve_fdm(int n, std::string suffix = "") {
    Timer t;
    t.addCheckPoint("start");
    vector<Triplet<scalar_t>> ts;
    ts.reserve(3*n -2);
    scalar_t h = 1.0/(n-1);
    Eigen::Matrix<scalar_t,Dynamic,1> rhs = Eigen::Matrix<scalar_t,Dynamic,1>::Zero(n);
    for (int i = 1; i < n-1; ++i) {
        ts.emplace_back(i, i, -2/h/h);
        ts.emplace_back(i, i-1, 1/h/h);
        ts.emplace_back(i, i+1, 1/h/h);
        rhs(i) = std::sin(i*h);
    }
    rhs(0) = 0;
    ts.emplace_back(0, 0, 1);
    ts.emplace_back(n-1, n-3, -1./2/h);
    ts.emplace_back(n-1, n-2, 2/h);
    ts.emplace_back(n-1, n-1, -3./2/h);
    rhs(n-1) = 0;

    SparseMatrix<scalar_t> M(n, n);
    M.setFromTriplets(ts.begin(), ts.end());
    SparseLU<SparseMatrix<scalar_t>> solver(M);
    Eigen::Matrix<scalar_t,Dynamic,1> sol = solver.solve(rhs);
    t.addCheckPoint("end");

    file.openFolder(format("/calc%06d", n));
    file.setFloatArray("solfdm"+suffix, sol);
    file.setDoubleAttribute("timefdm"+suffix, t.getTime("start", "end"));
}

template <typename scalar_t>
void solve_neu(int n, std::string suffix = "", int support_size = 3) {
    Timer t;
    t.addCheckPoint("start");
    scalar_t a = 0, b = 1;
    // Prepare domain
    RectangleDomain<Vec<scalar_t, 1>> domain(a, b);
    domain.fillUniform(n-1, 2);
    int N = domain.size();
    assert(N == n+1);
    domain.findSupport(support_size);
    // Prepare operators and matrix
    EngineMLS<Vec<scalar_t, 1>, Monomials, Monomials> mls(support_size, domain.positions, 1);
    SparseMatrix<scalar_t> M(N, N);
    M.reserve(Range<int>(N, support_size));
    auto op = make_mlsm<mlsm::d1|mlsm::lap>(domain, mls, domain.types != 0);  // All nodes, including boundary
    t.addCheckPoint("shapes");
    Eigen::Matrix<scalar_t,Dynamic,1> rhs = Eigen::Matrix<scalar_t,Dynamic,1>::Zero(N);
    // Set equation on interior
    for (int i : (domain.types == INTERNAL)) {
        op.lap(M, i, 1.0);  // laplace in interior
        rhs(i) = std::sin(domain.positions[i][0]);
    }
    // Set boundary conditions
    for (int i : (domain.types == BOUNDARY)) {
        scalar_t x = domain.positions[i][0];
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
    SparseLU<SparseMatrix<scalar_t>> solver;
    solver.compute(M);
    Eigen::Matrix<scalar_t,Dynamic,1> sol = solver.solve(rhs);
    t.addCheckPoint("end");

    file.openFolder(format("/calc%06d", N));
    file.setFloat2DArray("pos", domain.positions);
    file.setFloatArray("sol"+suffix, sol);
    file.setDoubleAttribute("time"+suffix, t.getTime("start", "end"));
    file.setDoubleAttribute("timepart"+suffix, t.getTime("start", "shapes"));
}

int main() {
    vector<int> testrange = {4, 5, 6, 7, 8, 9, 10, 11, 13, 15, 17, 19, 22, 25, 28, 31, 35, 39, 43, 48, 53, 60, 67, 74, 83, 92, 103, 115, 128, 143, 160, 178, 199, 221, 247, 276, 307, 343, 382, 427, 476, 531, 592, 661, 737, 822, 917, 1023, 1141, 1273, 1420, 1583, 1766, 1970, 2198, 2452, 2735, 3051, 3403, 3796, 4234, 4723, 5269, 5877, 6556, 7313, 8158, 9100, 10151, 11323, 12631, 14089, 15716, 17532, 19556, 21815, 24334, 27144, 30279, 33776, 37676, 42028, 46881, 52295, 58335, 65072, 72586, 80969, 90320, 100751, 112386, 125365}; //, 139843, 155993, 174009, 194104, 216521, 241526, 269419, 300533, 335241, 373956, 417143, 465318, 519056, 579000, 645866, 720455, 803658, 896470, 1000000};
    for (int n : testrange) {
        prn(n);
        solve_neu<double>(n);
        solve_fdm<double>(n+1);
        solve_neu<double>(n, "_s5", 5);
        if (n < 500) {
            solve_neu<float>(n, "_float");
            solve_fdm<float>(n+1, "_float");
        }
    }

    return 0;
}

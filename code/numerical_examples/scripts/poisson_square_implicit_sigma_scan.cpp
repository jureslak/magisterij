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

HDF5IO file("data/poisson_square_implicit_sigma_scan.h5", HDF5IO::DESTROY);

template <typename T>
double avg(const vector<T>& v) {
    double x = 0;
    for (T y : v) x += y;
    return x / v.size();
}

void solve(double sb, double sw) {
    Timer t;
    t.addCheckPoint("start");
    double a = 0, b = 1;
    int support_size = 9;
    // Prepare domain
    RectangleDomain<Vec2d> domain(a, b);
    int n = 50;
    domain.fillUniformWithStep(1./n, 1./n);
    int N = domain.size();
    domain.findSupport(support_size);
    t.addCheckPoint("domain");
    // Prepare operators and matrix
    NNGaussians<Vec2d> weight(sw*domain.characteristicDistance());
    NNGaussians<Vec2d> basis(sb*domain.characteristicDistance(), support_size);
    EngineMLS<Vec2d, NNGaussians, NNGaussians> mls(basis, weight);
    SparseMatrix<double> M(N, N);
    M.reserve(Range<int>(N, support_size));
    auto op = make_mlsm<mlsm::lap>(domain, mls, domain.types != 0);  // All nodes, including boundary
    t.addCheckPoint("shapes");
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
    M.makeCompressed();
    SparseLU<SparseMatrix<double>> solver;
    solver.compute(M);
    VecXd sol = solver.solve(rhs);
    t.addCheckPoint("end");

    std::string folder_name = format("/sb%013.6f/sw%013.6f", sb, sw);
    file.openFolder(folder_name);
    file.setFloat2DArray("pos", domain.positions);
    file.setFloatArray("sol", sol);
    file.setDoubleAttribute("N", domain.positions.size());
    file.setDoubleAttribute("cutoff", avg(op.cutoffs));
    file.setDoubleAttribute("sigmaW", sw);
    file.setDoubleAttribute("sigmaB", sb);
}

int main() {
    vector<double> sws = {0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3, 3.2, 3.4, 3.6, 3.8, 4, 4.2, 4.4, 4.6, 4.8, 5};
    vector<double> sbs = {0.5, 1, 6, 11, 16, 21, 26, 31, 36, 41, 46, 51, 56, 61, 66, 71, 76, 81, 86, 91, 96, 101, 106, 111, 116, 121, 126, 131, 136, 141, 146, 151, 156, 161, 166, 171, 176, 181, 186, 191, 196, 201, 206, 211, 216, 221, 226, 231, 236, 241, 246, 251, 256, 261, 266, 271, 276, 281, 286, 291, 296, 301, 306, 311, 316, 321, 326, 331, 336, 341, 346, 351, 356, 361, 366, 371, 376, 381, 386, 391, 396, 401, 406, 411, 416, 421, 426, 431, 436, 441, 446};
    for (int i = 0; i < sbs.size(); i += 1) {
        prn(sbs[i]);
        for (int j = 0; j < sws.size(); j += 1) {
            prn(sws[j]);
            solve(sbs[i], sws[j]);
        }
    }

    return 0;
}

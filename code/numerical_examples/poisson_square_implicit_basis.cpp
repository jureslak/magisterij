#include "common.hpp"
#include "domain.hpp"
#include "domain_extended.hpp"
#include "draw.hpp"
#include "mls.hpp"
#include "mlsm_operators.hpp"
#include "io.hpp"
#include <Eigen/Sparse>

using namespace std;
using namespace mm;
using namespace Eigen;

HDF5IO file("data/poisson_square_implicit_basis.h5", HDF5IO::DESTROY);


template <typename T>
double avg(const vector<T>& v) {
    double x = 0;
    for (T y : v) x += y;
    return x / v.size();
}

template <template<class> class T>
void solve(int n, T<Vec2d> basis, int support_size, double w, std::string name) {
    Timer t;
    t.addCheckPoint("start");
    double a = 0, b = 1;
    // Prepare domain
    RectangleDomain<Vec2d> domain(a, b);
    domain.fillUniformWithStep(1./n, 1./n);
    int N = domain.size();
    domain.findSupport(support_size);
    t.addCheckPoint("domain");
    // Prepare operators and matrix
    NNGaussians<Vec2d> weight(w*domain.characteristicDistance());
    if (name.substr(0, 3) != "mon")
        basis.shape = basis.shape * domain.characteristicDistance();
    EngineMLS<Vec2d, T, NNGaussians> mls(basis, weight);
    SparseMatrix<double> M(N, N);
    M.reserve(Range<int>(N, support_size));
    auto op = make_mlsm<mlsm::lap>(domain, mls, domain.types > 0);
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

    std::string foldername = format("/%s/calc%04d", name.c_str(), n);
    file.openFolder(foldername);
    file.setFloat2DArray("pos", domain.positions);
    file.setFloatArray("sol", sol);
    file.setDoubleAttribute("N", domain.positions.size());
    file.setDoubleAttribute("timetotal", t.getTime("start", "end"));
    file.setDoubleAttribute("timedomain", t.getTime("start", "domain"));
    file.setDoubleAttribute("timeshapes", t.getTime("domain", "shapes"));
    file.setDoubleAttribute("timesolve", t.getTime("shapes", "end"));
    file.setFloatArray("cutoff", op.cutoffs);
    file.setSparseMatrix("M", M);
    file.setDoubleArray("rhs", rhs);
//      file.setDoubleAttribute("iter", solver.iterations());
//      file.setDoubleAttribute("errest", solver.error());
}

int main() {
    vector<int> testrange = {6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 27, 28, 29, 31, 32, 34, 36, 37, 39, 41, 43, 46, 48, 50, 53, 56, 58, 61, 65, 68, 71, 75, 79, 83, 87, 91, 96, 101, 106, 112, 117, 123, 130, 136, 143, 151, 158, 166, 175, 184, 193, 203, 213, 224, 236, 248, 261, 274, 288, 303, 318, 334, 352, 370, 388, 408};
//      vector<int> testrange = {10};
    vector<double> sbs = {30, 90, 120, 150, 300};
    for (int i = 0; i < testrange.size(); i += 2) {
        int n = testrange[i];
        prn(n);
        for (double s : sbs) {
            solve(n, NNGaussians<Vec2d>(s, 9), 9, 0.75, format("w0.75/gau_s%08.2f", s));
            solve(n, NNGaussians<Vec2d>(s, 9), 9, 5, format("w5.00/gau_s%08.2f", s));
        }
    }

    return 0;
}

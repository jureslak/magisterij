#include "domain.hpp"
#include "domain_extended.hpp"
#include "mls.hpp"
#include "io.hpp"
#include "mlsm_operators.hpp"
#include <Eigen/Sparse>
#include <unsupported/Eigen/SparseExtra>
#include "draw.hpp"
#include "overseer_hertzian.hpp"

using namespace mm;
using namespace std;
using namespace Eigen;

Overseer O;

template <template<class> class T>
void solve_hertzian(int n, T<Vec2d> basis, string name) {
    /// [DOMAIN DEFINITION]
    RectangleDomain<Vec2d> domain(O.domain_lo, O.domain_hi);
    double dy = O.height / n;
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

    NNGaussians<Vec2d> weight(O.sigmaW * domain.characteristicDistance());
    if (name.substr(0, 3) != "mon")
        basis.shape = basis.shape * domain.characteristicDistance();
    EngineMLS<Vec2d, T, NNGaussians> mls(basis, weight);

    /// Initialize operators on all nodes
    auto op = make_mlsm(domain, mls, all);

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
        M.coeffRef(i, i) = 1;  // boundary conditions on the boundary
        rhs(i) = 0;
        M.coeffRef(i+N, i+N) = 1;  // boundary conditions on the boundary
        rhs(i+N) = 0;
    }

    // Set left boundary conditions
    for (int i : left) {
        M.coeffRef(i, i) = 1;  // boundary conditions on the boundary
        rhs(i) = 0;
        M.coeffRef(i+N, i+N) = 1;  // boundary conditions on the boundary
        rhs(i+N) = 0;
    }

    // Set right boundary conditions
    for (int i : right) {
        M.coeffRef(i, i) = 1;  // boundary conditions on the boundary
        rhs(i) = 0;
        M.coeffRef(i+N, i+N) = 1;  // boundary conditions on the boundary
        rhs(i+N) = 0;
    }



    for (int i : top) {
        // traction in x-direction
        op.der1(M,0,1,i,O.mu,0);
        op.der1(M,1,0,i,O.mu,0);
        // traction in y-direction
        op.der1(M,0,0,i,O.lam,1);
        op.der1(M,1,1,i,2.*O.mu+O.lam,1);
        double x = std::abs(domain.positions[i][0]);
        double trac;
        if (x < O.b) {
            trac = -O.p0 * std::sqrt(1 - x*x/O.b/O.b);
        } else {
            trac = 0;
        }
        rhs(i) = 0;
        rhs(i+N) = trac;
    }

    O.file.openFolder("/");
    O.file.setFloat2DArray("pos", domain.positions);
    O.file.setDoubleAttribute("N", N);
    O.file.setSparseMatrix("M", M);
}

int main(int arg_num, char* arg[]) {
    O.init("data/hertzian.xml", "data/hertzian_matrix_example.h5");

    Monomials<Vec2d> mon9({{0, 0}, {0, 1}, {0, 2}, {1, 0}, {1, 1}, {1, 2}, {2, 0}, {2, 1}, {2, 2}});
    int n = 5;
    solve_hertzian(n, mon9, "mon9");

    return 0;
}

#include "domain.hpp"
#include "io.hpp"
#include "draw.hpp"
#include "domain_extended.hpp"
#include <thread>

using namespace std;
using namespace mm;

double calculate_S(CircleDomain<Vec2d>& d) {
    d.findSupport(2);
    double dist = std::numeric_limits<double>::infinity();
    for (int i = 0; i < d.size(); ++i) {
        dist = min(dist, (d.positions[i] - d.positions[d.support[i][1]]).norm());
    }
    return dist;
}
double calculate_h(CircleDomain<Vec2d>& d) {
    KDTree<Vec2d> kdtree(d.positions);
    double dist = 0;
    int N = 500;
    for (int i = 0; i <= N; ++i) {
    for (int j = 0; j <= N; ++j) {
        double x = double(i) / N;
        double y = double(j) / N;
        if (x*x + y*y <= 1) {
            double d = kdtree.query({x, y}).second[0];
            dist = max(dist, d);
        }
    }
    }
    return dist;
}

int main() {
    HDF5IO file("../hs_data.h5", HDF5IO::DESTROY);
    file.openFolder("circle_refine");

    {
    vector<int> Ns;
    vector<double> hs, Ss;

    for (int N = 10; N < 2800; ++N) {

        CircleDomain<Vec2d> domain({0, 0}, 1);
        domain.fillUniform(N, int(3*4/3.14*std::sqrt(N)));
        domain.relax(6, 1e-2, [](Vec2d) { return 1.0; }, 3, 10);

        double S = calculate_S(domain);
        double h = calculate_h(domain);

        Ns.push_back(domain.size());
        hs.push_back(h);
        Ss.push_back(S);
        cout << N << " " << h << " " << S << endl;
    }

    file.setIntArray("N", Ns);
    file.setDoubleArray("S", Ss);
    file.setDoubleArray("h", hs);
    }

    file.openFolder("circle");

    {
    vector<int> Ns;
    vector<double> hs, Ss;

    for (int N = 10; N < 2800; ++N) {

        CircleDomain<Vec2d> domain({0, 0}, 1);
        domain.fillUniform(N, int(3*4/3.14*std::sqrt(N)));
        double S = calculate_S(domain);
        double h = calculate_h(domain);

        Ns.push_back(domain.size());
        hs.push_back(h);
        Ss.push_back(S);
        cout << N << " " << h << " " << S << endl;
    }

    file.setIntArray("N", Ns);
    file.setDoubleArray("S", Ss);
    file.setDoubleArray("h", hs);
    }
    return 0;
}


#include "domain.hpp"
#include "io.hpp"
#include "draw.hpp"
#include "domain_extended.hpp"
#include <thread>

using namespace std;
using namespace mm;

#include <X11/Xlib.h>

int main() {
    XInitThreads();
    int N = 300;

    CircleDomain<Vec2d> domain({0, 0}, 1);
    domain.fillRandomInterior(N);
    domain.fillUniformBoundary(int(3*4/3.14*std::sqrt(N)));

    auto d = Domain<Vec2d>::makeClone(domain);
    thread th ([&](){ draw2D(d); });
    domain.relax(6, 1e-2, [](Vec2d) { return 1.0; }, 3, 10);
    thread th2 ([&](){ draw2D(domain); });

    th.join();
    th2.join();

    return 0;
}



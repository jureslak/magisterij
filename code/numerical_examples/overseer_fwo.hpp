#ifndef __OVERSEER__
#define __OVERSEER__

#include "io.hpp"
#include "mls.hpp"
#include "types.hpp"
#include "util.hpp"

using namespace mm;

class Overseer {
  public:
    int num_threads;

    int n;                // support size;
    int m;                // number of basis functions or order of basis

    double sigmaW;        // squared sigma for weight
    double sigmaB;        // squared sigma for Gaussian basis
    double precision;     // time step, implicit solver precision
    int ny;               // num of discretization nodes in y direction
    double dy;            // discretization step in y direction

    // physics params
    double rho, mu, E, lam, v, Q, COF, sigmaAxial, a, p0, Estar, c, e; // material properties
    double time;          // time
    double width, height, thickness;         // half width of domain
    Vec2d domain_lo;      // domain lo point
    Vec2d domain_hi;      // domain hi point
    double radius;        // pressing cylinder radius

    // case dependant
    double force;

    int maxiter, fillfactor;
    double tolerance, droptol;

    std::string hdf5_filename;
    HDF5IO file;

    void init(const std::string& input_filename, const std::string& output_filename) {
        // [load other parameters from xml]
        hdf5_filename = output_filename;
        XMLloader xml(input_filename);

        width =      xml.getAttribute({"params", "case"}, "width");
        height =     xml.getAttribute({"params", "case"}, "height");
        thickness =     xml.getAttribute({"params", "case"}, "thickness");
        force =      xml.getAttribute({"params", "case"}, "force");
        radius =     xml.getAttribute({"params", "case"}, "radius");
        rho =        xml.getAttribute({"params", "phy"}, "rho");
        v =          xml.getAttribute({"params", "phy"},  "v");
        E =          xml.getAttribute({"params", "phy"}, "E");
        Q =          xml.getAttribute({"params", "phy"}, "Q");
        COF =        xml.getAttribute({"params", "phy"}, "COF");
        sigmaAxial = xml.getAttribute({"params", "phy"}, "sigmaAxial");

        m =          xml.getAttribute({"params", "mls"}, "m");
        n =          xml.getAttribute({"params", "mls"}, "n");
        sigmaW =     xml.getAttribute({"params", "mls"}, "sigmaW");
        sigmaB =     xml.getAttribute({"params", "mls"}, "sigmaB");

        ny =         xml.getAttribute({"params", "num"}, "ny");
        dy = height / ny;

        num_threads= xml.getAttribute({"params", "sys"}, "num_threads");

        droptol = xml.getAttribute({"params", "solver"}, "droptol");
        fillfactor = xml.getAttribute({"params", "solver"}, "fillfactor");
        tolerance = xml.getAttribute({"params", "solver"}, "tolerance");
        maxiter = xml.getAttribute({"params", "solver"}, "maxiter");

        // hdf5 output
        file.openFile(hdf5_filename, HDF5IO::APPEND);
        std::ifstream conf(input_filename);
        assert(conf && "Could not open conf file.");
        std::stringstream buffer;
        buffer << conf.rdbuf();

        // parameter logic
        mu = E / 2. / (1+v);
        lam = E * v / (1-2*v) / (1+v);
//        lam = 2 * mu * lam / (2 * mu + lam);  // plane stress
        Estar = E/(2*(1-v*v));
        a = 2*std::sqrt(std::abs(force*radius/(thickness*M_PI*Estar)));
        p0 = std::sqrt(std::abs(force*Estar/(thickness*M_PI*radius)));
        c = a * std::sqrt(1 - Q / COF / std::abs(force));
        e = a * sigmaAxial / 4 / COF / p0;

        file.openFolder("/");
        file.setStringAttribute("conf", buffer.str());
        file.setDoubleAttribute("width"      , width      );
        file.setDoubleAttribute("height"     , height     );
        file.setDoubleAttribute("force"      , force      );
        file.setDoubleAttribute("radius"     , radius     );
        file.setIntAttribute("m"             , m          );
        file.setIntAttribute("n"             , n          );
        file.setDoubleAttribute("sigmaW"     , sigmaW     );
        file.setDoubleAttribute("sigmaB"     , sigmaB     );
        file.setDoubleAttribute("precision"  , precision  );
        file.setIntAttribute("ny"            , ny         );
        file.setDoubleAttribute("v"          , v          );
        file.setDoubleAttribute("E"          , E          );
        file.setIntAttribute("num_threads"   , num_threads);
        file.setDoubleAttribute("mu"         , mu         );
        file.setDoubleAttribute("lam"        , lam        );
        file.setDoubleAttribute("COF"        , COF        );
        file.setDoubleAttribute("p0"         , p0         );
        file.setDoubleAttribute("sigmaAxial" , sigmaAxial );
        file.setDoubleAttribute("Q"          , Q          );
        file.setDoubleAttribute("a"          , a          );
        file.setDoubleAttribute("c"          , c          );
        file.setDoubleAttribute("e"          , e          );
        file.setDoubleAttribute("Estar"      , Estar      );

        file.closeFolder();
        file.closeFile();

        domain_lo = Vec2d(-width/2, -height/2);
        domain_hi = Vec2d(width/2, 0);

        //[System setup]
//        omp_set_num_threads(num_threads);
    }
};

#endif
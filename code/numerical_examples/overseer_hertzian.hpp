#ifndef __OVERSEER__
#define __OVERSEER__

#include "io.hpp"
#include "mls.hpp"
#include "types.hpp"
#include "util.hpp"

using namespace mm;

class Overseer {
private:
    XMLloader xml;         //param files
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
    double mu, E, lam, v, b, p0, E_; // material properties
    double time;          // time
    double width;         // half width of domain
    double height;        // height of domain
    Vec2d domain_lo;      // domain lo point
    Vec2d domain_hi;      // domain hi point
    double radius;        // pressing cylinder radius

    // case dependant
    double force;

    std::string hdf5_filename;
    HDF5IO file;

    void init(std::string input_filename, std::string output_filename) {
        // [load other parameters from xml]
        hdf5_filename = output_filename;
        xml(input_filename);

        width =      xml.getAttribute({"params", "case"}, "width");
        height =     xml.getAttribute({"params", "case"}, "height");
        force =      xml.getAttribute({"params", "case"}, "force");
        radius =     xml.getAttribute({"params", "case"}, "radius");

        m =          xml.getAttribute({"params", "mls"}, "m");
        n =          xml.getAttribute({"params", "mls"}, "n");
        sigmaW =     xml.getAttribute({"params", "mls"}, "sigmaW");
        sigmaB =     xml.getAttribute({"params", "mls"}, "sigmaB");

        precision =  xml.getAttribute({"params", "num"}, "precision");
        ny =         xml.getAttribute({"params", "num"}, "ny");
        dy = height / ny;

        v =          xml.getAttribute({"params", "phy"},  "v");
        E =          xml.getAttribute({"params", "phy"}, "E");

        num_threads= xml.getAttribute({"params", "sys"}, "num_threads");

        // hdf5 output
        file.openFile(hdf5_filename, HDF5IO::DESTROY);
        std::ifstream conf(input_filename);
        assert(conf && "Could not open conf file.");
        std::stringstream buffer;
        buffer << conf.rdbuf();

        E_ = E/(2*(1-v*v));
        b  = 2*std::sqrt(std::abs(force*radius/(M_PI*E_)));
        p0 = std::sqrt(std::abs(force*E_/(M_PI*radius)));

        file.openFolder("/");
        file.setStringAttribute("conf", buffer.str());
        file.setDoubleAttribute("b", b);
        file.setDoubleAttribute("p0", p0);
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

        //[params LOGIC]
        mu = E / 2. / (1+v);
        lam = E * v / (1-2*v) / (1+v);

        domain_lo = Vec2d(-width, -height);
        domain_hi = Vec2d(width, 0);

        //[System setup]
//        omp_set_num_threads(num_threads);
    }
};

#endif
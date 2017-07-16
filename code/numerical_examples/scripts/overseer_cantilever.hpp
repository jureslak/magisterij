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
    double mu, E, lam, v, I; // material properties
    double L;             // width of domain
    double D;             // height of domain
    Vec2d domain_lo;      // domain lo point
    Vec2d domain_hi;      // domain hi point

    // case dependant
    double P;

    std::string hdf5_filename;
    HDF5IO file;

    void init(std::string input_filename, std::string output_filename) {
        // [load other parameters from xml]
        hdf5_filename = output_filename;
        xml(input_filename);

        L =      xml.getAttribute({"params", "case"}, "L");
        D =      xml.getAttribute({"params", "case"}, "D");
        P =      xml.getAttribute({"params", "case"}, "P");
        I = D*D*D/12;

        m =          xml.getAttribute({"params", "mls"}, "m");
        n =          xml.getAttribute({"params", "mls"}, "n");
        sigmaW =     xml.getAttribute({"params", "mls"}, "sigmaW");
        sigmaB =     xml.getAttribute({"params", "mls"}, "sigmaB");

        precision =  xml.getAttribute({"params", "num"}, "precision");
        ny =         xml.getAttribute({"params", "num"}, "ny");
        dy = D / ny;

        v =          xml.getAttribute({"params", "phy"},  "v");
        E =          xml.getAttribute({"params", "phy"}, "E");

        num_threads= xml.getAttribute({"params", "sys"}, "num_threads");

        // hdf5 output
        file.openFile(hdf5_filename, HDF5IO::DESTROY);
        std::ifstream conf(input_filename);
        assert(conf && "Could not open conf file.");
        std::stringstream buffer;
        buffer << conf.rdbuf();

        file.openFolder("/");
        file.setStringAttribute("conf", buffer.str());
        file.setDoubleAttribute("L"      , L      );
        file.setDoubleAttribute("D"     , D     );
        file.setDoubleAttribute("P"      , P      );
        file.setDoubleAttribute("I"      , I      );
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
        // plane stress
        lam = 2 * mu * lam / (2 * mu + lam);

        file.setDoubleAttribute("mu"          , mu        );
        file.setDoubleAttribute("lambda"          , lam          );

        domain_lo = Vec2d(0, -D/2);
        domain_hi = Vec2d(L, D/2);

        //[System setup]
//        omp_set_num_threads(num_threads);
    }
};

#endif
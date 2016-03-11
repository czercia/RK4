#ifndef _RK4_H_
#define _RK4_H_
#include <iostream>
#include <cstdlib>
#include <vector>
#include <fstream>
#include <sstream>
#include <complex>
#include "Parameters.h"
//using namespace std;



namespace RK4
{
    typedef std::vector<double> VecDoub;
    typedef std::vector<std::complex<double>> VecComp;
    typedef std::complex<double> Complex;
    void saveToFile(Parameters p, VecDoub xx, VecComp yy);
    void solve (Parameters p);
    void createGnuplotScript(std::vector<std::string> names, Parameters p, int nx, std::string filetype, int nt);

}
#endif

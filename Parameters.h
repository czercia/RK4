#ifndef _PARAMETERS_H_
#define _PARAMETERS_H_

#include <iostream>
#include <cstdlib>
#include <math.h>
#include <complex>
#include <string>
#include <vector>

typedef std::complex<double> Complex;
typedef std::vector<std::complex<double>> VecComp;

class Parameters {

    double t0;
    double tMax;
    double dt;
    double L;
    int init;

    void initialPsi(Complex *y, VecComp &V, int y_length);

    Complex secondDerivative(int i, Complex *y, int y_length);

public:
    double getT0();

    double getTMax();

    double getDT();

    double getL();

    Parameters(double t_0, double t_max, double d_t, double l, int);

    Complex function(double x, double t, int j, Complex *y, int y_length);

    Complex psi0(double x, double a, double C);

    double DX();

};

#endif

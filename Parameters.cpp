#include "Parameters.h"

Parameters::Parameters(double t_0, double t_max, double d_t, double l, int initial) {
    t0 = t_0;
    tMax = t_max;
    dt = d_t;
    L = l;
    init = initial;
    Psi0 = 0;
}

double Parameters::getT0() {
    return this->t0;
}

double Parameters::getTMax() {
    return this->tMax;
}

double Parameters::getDT() {
    return this->dt;
}

double Parameters::getL() {
    return this->L;
}

Complex Parameters::function(double x, double t, int j, Complex *y, int y_length, Complex psi) {
    Complex a(0, 1);
    //Complex value =  a  * abs(y[j])*abs(y[j]) * y[j];
    Complex value = a * secondDerivative(j, y, y_length) - (a + (1.0, 0.0)) * abs(psi) * abs(psi) * psi + 0.1 * psi;

    return value;
}

Complex Parameters::psi0(double x, double a, double C) {
    if (this->Psi0 == 0)
        return exp(-x * x / 2);
    if (this->Psi0 == 1)
        return generateGaussianNoise(1);
}

double Parameters::DX() {
    return 1 / (400 * sqrt(this->dt));
}

void Parameters::initialPsi(Complex *y, VecComp &V, int y_length) {
    int in;
    in = this->init;
    switch (in) {
        case 0:    //warunki Dirichleta
        {
            Complex value(0, 0);
            V.push_back(value);
            V.push_back(value);
            break;
        }

        case 1: //periodic
        {
            int n = y_length - 1;
            V.push_back(y[n]);
            V.push_back(y[0]);
            break;
        }
        default: // dirichlet
        {
            Complex value(0, 0);
            V.push_back(value);
            V.push_back(value);
            break;
        }
    }
}

Complex Parameters::secondDerivative(int i, Complex *y, int y_length) {
    if (i > 0 && i < y_length - 2) {
        double dx2 = this->DX() * this->DX();
        Complex value = (y[i - 1] - 2.0 * y[i] + y[i + 1]) / dx2;
        return value;
    }
    else {

        VecComp cond;
        initialPsi(y, cond, y_length);
        if (i == 0)
            return (cond[0] - 2.0 * y[i] + y[i + 1]) / ((this->DX()) * (this->DX()));
        if (i == y_length - 1)
            return (y[i - 1] - 2.0 * y[i] + cond[1]) / ((this->DX()) * (this->DX()));
    }
}

double Parameters::generateGaussianNoise(double variance) {

    static bool haveSpare = false;
    static double rand1, rand2;

    if (haveSpare) {
        haveSpare = false;
        return sqrt(variance * rand1) * sin(rand2);
    }

    haveSpare = true;

    rand1 = rand() / ((double) RAND_MAX);
    if (rand1 < 1e-100) rand1 = 1e-100;
    rand1 = -2 * log(rand1);
    rand2 = (rand() / ((double) RAND_MAX)) * (2 * M_PI);

    return sqrt(variance * rand1) * cos(rand2);

}
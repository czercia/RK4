#include "Parameters.h"

Parameters::Parameters(double t_0, double t_max, double d_t, double l, int initial)
{
    t0 = t_0;
    tMax = t_max;
    dt = d_t;
    L = l;
    init = initial;
}

double Parameters::getT0()
{
    return this -> t0;
}
double Parameters::getTMax()
{
    return this ->tMax;
}
double Parameters::getDT()
{
    return this -> dt;
}
double Parameters::getL()
{
    return this->L;
}
Complex Parameters::function(double x, double t, int j, /*VecComp y*/Complex *y, int y_length)
{
    Complex a(0,1);
    //Complex value =  a  * abs(y[j])*abs(y[j]) * y[j];
    Complex value = a * secondDerivative(j, y, y_length)+a * abs(y[j])*abs(y[j]) * y[j];
    //std::cout << y[j+1] <<std::endl;
    return value;
}
Complex Parameters::psi0(double x, double a, double C)
{
    Complex k(0,0);
    return C * exp (- a * x * x/2);
     //return k;
}
double Parameters::DX() 
{
    return 1/(200 * sqrt(this -> dt));
    //return 0.1;
}

void Parameters::initialPsi(Complex *y, VecComp &V, int y_length)
{
    int in;
    in = this -> init;
    switch (in)
    {
        case 0:	//warunki Dirichleta
        {
            Complex value(0,0);
            V.push_back(value);
            V.push_back(value);
            break;
        }

        case 1: //periodyczne
        {
            //int n = y.size() - 1;
            int n = y_length-1;
            V.push_back(y[n]);
            V.push_back(y[0]);
            break;
        }
    }
}

Complex Parameters::secondDerivative(int i, /*VecComp y*/Complex *y, int y_length)
{
	if (i > 0 && i < y_length - 2)
	{
        double dx2 = this -> DX()*this -> DX();
	     Complex value =  (y[i - 1] - 2.0 * y[i] + y[i + 1])/dx2;
       // std::cout << (y[i+1] - 2.0*y[i] + y[i-1])/dx2 << std::endl;
        return value;
	}
else
    {

      VecComp cond;
    initialPsi(y, cond, y_length);
    //std::cout << cond[0] << " " <<cond[1] << std::endl;
        if (i == 0)
            return (cond[0] - 2.0 * y[i] + y[i + 1]) / ((this -> DX()) * (this -> DX()));
        if (i == y_length - 1)
            return (y[i - 1] - 2.0 * y[i] + cond[1]) / ((this -> DX()) * (this -> DX()));
    }
    //std::cout << y[i] <<std::endl;
   // return y[i];
}

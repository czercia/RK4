
#include "Parameters.h"
#include "RK4.h"


int main()
{
    std::srand(std::time(NULL));
    Parameters params1(0, 5, 0.00005, 60, 0); //t0 tmax dt L 0-dirichlet/1-periodic
    RK4::solve(params1);

    return 0;
}

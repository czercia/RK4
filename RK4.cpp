#include "RK4.h"

void RK4::saveToFile(Parameters p, VecDoub xx, VecComp yy) {
    std::ostringstream filename;
    filename << "//home//marta//rungekuttacomplex//results.dat";
    std::ofstream f(filename.str().c_str());

    f << "# Parameters: " << "\n";
    f << "# t0, t_max, dt: " << p.getT0() << " " << p.getTMax() << " " << p.getDT() << " \n";
    f << "# L, dx: " << p.getL() << " " << p.DX() << "\n";
    f << "# x	Re(psi(x))   Im(Psi(x))" << "\n";
    int i = 0;
    while (i < xx.size()) {
        f << xx[i] << "\t" << yy[i].real() << "\t" << yy[i].imag() << "\t" << abs(yy[i]) * abs(yy[i]) << "\n";
        i++;
    }
    f.close();

    std::cout << "Psi(x, t_Max) saved to \"" << filename.str().c_str() << "\"" << std::endl;


}


void RK4::solve(Parameters p) {
    Complex k1(0, 0);
    Complex k2(0, 0);
    Complex k3(0, 0);
    Complex k4(0, 0);
    double h = p.getDT();


    int nt = (p.getTMax() - p.getT0()) / h;
    int nx = 2 * p.getL() / p.DX();
    std::cout << "nx " << nx << std::endl;
    std::cout << "nt " << nt << std::endl;


    std::ostringstream nameA, nameR, nameI, nameP, nameRes;
    std::vector<std::string> names;
    nameA << "//home//marta//RK4//AbsPsi2.dat";
    nameR << "//home//marta//RK4//RePsi2.dat";
    nameI << "//home//marta//RK4//ImPsi2.dat";
    nameP << "//home//marta//RK4//ArgPsi2.dat";
    nameRes << "//home//marta//RK4//Results2.dat";
    std::ofstream ab(nameA.str().c_str());
    names.push_back(nameA.str().c_str());
    std::ofstream re(nameR.str().c_str());
    names.push_back(nameR.str().c_str());
    std::ofstream im(nameI.str().c_str());
    names.push_back(nameI.str().c_str());
    std::ofstream ar(nameP.str().c_str());
    names.push_back(nameP.str().c_str());
    std::ofstream res(nameRes.str().c_str());
    names.push_back(nameRes.str().c_str());

    ab << "# Parameters: " << "\n";
    ab << "# t0, t_max, dt: " << p.getT0() << " " << p.getTMax() << " " << p.getDT() << " \n";
    ab << "# L, dx: " << p.getL() << " " << p.DX() << "\n";
    re << "# Parameters: " << "\n";
    re << "# t0, t_max, dt: " << p.getT0() << " " << p.getTMax() << " " << p.getDT() << " \n";
    re << "# L, dx: " << p.getL() << " " << p.DX() << "\n";
    im << "# Parameters: " << "\n";
    im << "# t0, t_max, dt: " << p.getT0() << " " << p.getTMax() << " " << p.getDT() << " \n";
    im << "# L, dx: " << p.getL() << " " << p.DX() << "\n";
    ar << "# Parameters: " << "\n";
    ar << "# t0, t_max, dt: " << p.getT0() << " " << p.getTMax() << " " << p.getDT() << " \n";
    ar << "# L, dx: " << p.getL() << " " << p.DX() << "\n";
    res << "# Parameters: " << "\n";
    res << "# t0, t_max, dt: " << p.getT0() << " " << p.getTMax() << " " << p.getDT() << " \n";
    res << "# L, dx: " << p.getL() << " " << p.DX() << "\n";


    VecDoub T; //Time
    T.push_back(p.getT0());

    Complex *YBefore = NULL;
    YBefore = new Complex[nx];

    for (int i = 0; i < nt; i++) {

        double t = p.getT0() + i * h;
        T.push_back(t);

        double *X = NULL;
        X = new double[nx];

        Complex *Y = NULL;
        Y = new Complex[nx];
        if (i == 0) {
            for (int j = 0; j < nx; j++) {
                X[j] = -p.getL() + j * p.DX();
                Complex y;

                y = p.psi0(X[j], 1, 1);
                Y[j] = y;
                YBefore[j] = Y[j];
            }


        }
        else {
            int j;
            Complex y;
            //#pragma omp parallel for private (y)
            for (j = 0; j < nx; j++) {
                X[j] = -p.getL() + j * p.DX();
                if (i == 0) {
                    y = p.psi0(X[j], 1, 1);
                }
                else if (i > 0) {
                    k1 = h * p.function(X[j], T[i - 1], j, YBefore, nx, YBefore[j]);
                    k2 = h * p.function(X[j], T[i - 1] + h / 2, j, YBefore, nx, YBefore[j] + 0.5 * k1 * h);
                    //std::cout << YBefore[j] + 0.5 * k1 * h << std::endl;
                    k3 = h * p.function(X[j], T[i - 1] + h / 2, j, YBefore, nx, YBefore[j] + 0.5 * k2 * h);
                    k4 = h * p.function(X[j], T[i - 1] + h, j, YBefore, nx, YBefore[j] + k3 * h);
                    y = YBefore[j] +
                        (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
                }

                Y[j] = y;
                YBefore[j] = Y[j];
            }
        }
        if (i % 100 == 0) {
            for (int k = 0; k < nx; k++) {
                ab << abs(YBefore[k]) * abs(YBefore[k]) << "\t";
                re << Y[k].real() << "\t";
                im << Y[k].imag() << "\t";
                ar << std::arg(Y[k]) << "\t";

            }
        }
        if (i == nt - 1) {
            for (int l = 0; l < nx; l++) {
                res << X[l] << "\t";
                res << abs(Y[l]) * abs(Y[l]) << "\t";
                res << Y[l].real() << "\t";
                res << Y[l].imag() << "\t";
                res << std::arg(Y[l]) << "\n";
            }
        }

        delete Y;
        delete X;
        if (i % 100 == 0) {
            ab << "\n";
            re << "\n";
            im << "\n";
            ar << "\n";
        }
    }
    delete YBefore;


    ab.close();
    re.close();
    im.close();
    ar.close();
    res.close();


    std::cout << "calculated" << std::endl;
    std::cout << "[Abs(Psi(x, t))]^2 saved to " << nameA.str().c_str() << std::endl;
    std::cout << "Im(Psi(x, t)) saved to " << nameI.str().c_str() << std::endl;
    std::cout << "Re(Psi(x, t)) saved to " << nameR.str().c_str() << std::endl;


    createGnuplotScript(names, p, nx,
                        "png", nt / 100);

}

void ::RK4::createGnuplotScript(std::vector<std::string> names, Parameters p, int nx, std::string filetype, int nt) {
    std::ostringstream fn;
    fn << "//home//marta//RK4//script.txt";
    std::ofstream f(fn.str().c_str());

    f << "reset \n";
    f << "set pm3d map \n";
    f << "set xlabel \"x\" \n";
    f << "set ylabel \"t\" \n";
    f << "set terminal " << filetype << "\t";
    f << "size 1000, 1000 \n";
    f << "unset key\n";

    f << "set output \"AbsPsi." << filetype << "\"\n";
    f << "set title \" [Abs(Psi(x,t))]^2\" \n";
    f << "splot \"" << names[0] << "\" u (($1-" << nx / 2 << ")/" << nx / (2 * p.getL()) << "):(($2)/" <<
    (nt / p.getTMax()) << "):3 matrix w image\n";
    f << "unset title \n";
    f << "unset output \n";

    f << "set output \"RePsi." << filetype << "\" \n";
    f << "set title \" Re(Psi(x,t))\" \n";
    // f << "splot \""<< names[1] << "\" u (($2)/"<< p.getTMax() <<"):(($1-" << nx/2<<")/" << nx/ (2*p.getL())<<"):3 matrix w image \n";
    f << "splot \"" << names[1] << "\" u (($1-" << nx / 2 << ")/" << nx / (2 * p.getL()) << "):(($2)/" <<
    (nt / p.getTMax()) << "):3 matrix w image\n";
    f << "unset title \n";
    f << "unset output \n";

    f << "set output \"ImPsi." << filetype << "\" \n";
    f << "set title \" Im(Psi(x,t))\" \n";
    f << "splot \"" << names[2] << "\" u (($1-" << nx / 2 << ")/" << nx / (2 * p.getL()) << "):(($2)/" <<
    (nt / p.getTMax()) << "):3 matrix w image\n";
    f << "unset title \n";
    f << "unset output \n";

    f << "set output \"ArgPsi." << filetype << "\" \n";
    f << "set title \" Arg(Psi(x,t))\" \n";
    f << "splot \"" << names[3] << "\" u (($1-" << nx / 2 << ")/" << nx / (2 * p.getL()) << "):(($2)/" <<
    (nt / p.getTMax()) << "):3 matrix w image\n";
    f << "unset title \n";
    f << "unset output \n";


    f << "unset pm3d \n";
    f << "set grid \n";
    f << "set output \"Results." << filetype << "\" \n";
    f << "set xlabel \" x \" \n";
    f << "unset ylabel \n";
    f << "set multiplot layout 4,1 rowsfirst \n";
    f << "set title  \" [Abs(Psi(x,t = " << p.getTMax() << "))]^2\" \n";
    f << "plot \"" << names[4] << "\" using 1:2 with lines\n";
    f << "set title  \" Re(Psi(x,t = " << p.getTMax() << "))\" \n";
    f << "plot \"" << names[4] << "\" using 1:3 with lines\n";
    f << "set title  \" Im(Psi(x,t = " << p.getTMax() << "))\" \n";
    f << "plot \"" << names[4] << "\" using 1:4 with lines\n";
    f << "set title  \" Arg(Psi(x,t = " << p.getTMax() << "))\" \n";
    f << "plot \"" << names[4] << "\" using 1:5 with lines\n";
    f << "unset multiplot \n";
    f << "unset output \n";

    f << "unset terminal \n";


}

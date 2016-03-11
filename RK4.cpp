#include "RK4.h"

void RK4::saveToFile(Parameters p, VecDoub xx, VecComp yy)
{
    std::ostringstream filename;
    filename << "//home//marta//rungekuttacomplex//results.dat";
    std::ofstream f (filename.str().c_str());

    f << "# Parameters: " << "\n";
    f << "# t0, t_max, dt: " << p.getT0() << " "  << p.getTMax() << " " << p.getDT() << " \n";
    f << "# L, dx: " << p.getL() << " " << p.DX() <<"\n";
    f << "# x	Re(psi(x))   Im(Psi(x))" << "\n";
    int i = 0;
    while (i < xx.size())
    {
        f << xx[i] << "\t" << yy[i].real() << "\t" << yy[i].imag() << "\t" << abs(yy[i])*abs(yy[i]) <<"\n";
        i++;
    }
    f.close();

    std::cout <<"Psi(x, t_Max) saved to \"" << filename.str().c_str()<< "\"" << std::endl;


}


void RK4::solve(Parameters p)
{
    Complex k1(0, 0);
    Complex k2(0, 0);
    Complex k3(0, 0);
    Complex k4(0, 0);
    double h = p.getDT();

   // Vec tt;
    VecComp yy2;
    //VecDoub xx;



    int nt = (p.getTMax() - p.getT0() ) / h;
    int nx = 2 * p.getL()/p.DX();
	std::cout << "nx " << nx << std::endl;
    std::cout << "nt " << nt << std::endl;
    VecComp psi;

   // znalezc jakis madrzejszy sposob zapisnia tego
   //funkcja od pliku?
    std::ostringstream nameA, nameR, nameI, nameRes;
    std::vector<std::string> names;
    nameA << "//home//marta//RKponiedz//rungekuttacomplex//AbsPsi.dat";
    nameR << "//home//marta//RKponiedz//rungekuttacomplex//RePsi.dat";
    nameI << "//home//marta//RKponiedz//rungekuttacomplex//ImPsi.dat";
    nameRes << "//home//marta//RKponiedz//rungekuttacomplex//Results.dat";
    std::ofstream ab (nameA.str().c_str()); names.push_back(nameA.str().c_str());
    std::ofstream re (nameR.str().c_str()); names.push_back(nameR.str().c_str());
    std::ofstream im (nameI.str().c_str()); names.push_back(nameI.str().c_str());
    std::ofstream res (nameRes.str().c_str()); names.push_back(nameRes.str().c_str());

    ab << "# Parameters: " << "\n";
    ab << "# t0, t_max, dt: " << p.getT0() << " "  << p.getTMax() << " " << p.getDT() << " \n";
    ab << "# L, dx: " << p.getL() << " " << p.DX() <<"\n";
    re << "# Parameters: " << "\n";
    re << "# t0, t_max, dt: " << p.getT0() << " "  << p.getTMax() << " " << p.getDT() << " \n";
    re << "# L, dx: " << p.getL() << " " << p.DX() <<"\n";
    im << "# Parameters: " << "\n";
    im << "# t0, t_max, dt: " << p.getT0() << " "  << p.getTMax() << " " << p.getDT() << " \n";
    im << "# L, dx: " << p.getL() << " " << p.DX() <<"\n";
    res << "# Parameters: " << "\n";
    res << "# t0, t_max, dt: " << p.getT0() << " "  << p.getTMax() << " " << p.getDT() << " \n";
    res << "# L, dx: " << p.getL() << " " << p.DX() <<"\n";


    
//    for (int j = 0; j < nx; j++)
//    {
//        xx.push_back(-p.getL() + j * p.DX());
//
//        double t = p.getT0();
//        Complex y = p.psi0(xx[j], 1, 1);
//
//        VecDoub tt;
//	    VecComp yy;
//        yy.push_back(y);
//        tt.push_back(p.getT0());
    VecDoub T; //wektor wartości czasów
    T.push_back(p.getT0()); //pierwszy czas
    //VecComp YBefore; //poprzedni y
    Complex *YBefore = NULL;
    YBefore = new Complex [nx];

    //trzeba wypelnic wektor Y wartosciami poczatkowymi psi - dla kazdego x
    for (int i = 0; i < nt; i++) {

        double t = p.getT0() + i * h; //czas
        T.push_back(t); //dopisujemy czas[i] do wektora T
        //VecDoub X; //wektor wartosci x dla czasu T[i]
        double *X = NULL;
        X = new double[nx];
        //X.push_back(-p.getL()); //pierwszy x = -L
        //VecComp Y; //Psi(t) dla kolejnych x
        Complex *Y = NULL;
        Y = new Complex[nx];
        if (i==0)
        {
            for ( int j = 0; j < nx; j++)
            {
                //X.push_back(-p.getL() + j * p.DX()); //wartosci x
                X[j] = -p.getL() + j * p.DX();
                Complex y;

               //dla czasu poczatkowego zapelniamy wektor Y wartosciami poczatkowymi psi
                y = p.psi0(X[j], 1, 1);
                Y[j] = y;
                YBefore[j] = Y[j];
            }


        }
        else
        {
            int j;
            Complex y;
           //#pragma omp parallel for private (y)
            for ( j = 0; j < nx; j++)
            {
                X[j] = -p.getL() + j * p.DX();
                //X.push_back(-p.getL() + j * p.DX()); //wartosci x
//                Complex y;

                if (i == 0) {  //dla czasu poczatkowego zapelniamy wektor Y wartosciami poczatkowymi psi
                    y = p.psi0(X[j], 1, 1);
                }
                else if (i > 0) {
                    k1 = p.function(X[j], T[i - 1], j, YBefore, nx);
                    k2 = p.function(X[j], T[i - 1] + h / 2, j, YBefore, nx);
                    k3 = p.function(X[j], T[i - 1] + h / 2, j, YBefore, nx);
                    k4 = p.function(X[j], T[i - 1] + h, j, YBefore, nx);
                    y = YBefore[j] +
                        h * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
                }

                //Y.push_back(y);
                Y[j] = y;
                YBefore[j] = Y[j];
            }
        }
        if (i % 100 == 0) {
            for (int k = 0; k < nx; k++) {
                ab << abs(YBefore[k]) * abs(YBefore[k]) << "\t";
                re << Y[k].real() << "\t";
                im << Y[k].imag() << "\t";

            }
        }
        if (i == nt - 1) {
            for (int l = 0; l < nx; l++) {
                res << X[l] << "\t";
                res << abs(Y[l]) * abs(Y[l]) << "\t";
                res << Y[l].real() << "\t";
                res << Y[l].imag() << "\n";
            }
        }

        delete Y;
        delete X;
        //YBefore = Y;
        //YBefore[j] = Y[j];
        // if (YBefore.size() != Y.size()) std::cout << "ERROR" <<std::endl;//wektor Y zostanie stworzony nowy dla kolejnego czasu a sa potrzebne wartosci tego dla poprzedniego czasu, wiec
        //przypisanie ich do YBefore
        if (i % 100 == 0)
        {
            ab << "\n";
            re << "\n";
            im << "\n";
        }

    }
    delete YBefore;

      // psi.push_back(y);
        //yy2 = yy

    ab.close();
    re.close();
    im.close();
    res.close();


    std::cout << "calculated" << std::endl;
    std::cout << "[Abs(Psi(x, t))]^2 saved to " << nameA.str().c_str() <<std::endl;
    std::cout << "Im(Psi(x, t)) saved to " << nameI.str().c_str() <<std::endl;
    std::cout << "Re(Psi(x, t)) saved to " << nameR.str().c_str() <<std::endl;

    //std::cout << "rozmiar: " << xx.size() << " " << psi.size() << std::endl;
    //saveToFile(p, X, psi);
    createGnuplotScript(names, p, nx,
                        "png", nt/100);

}

void ::RK4::createGnuplotScript(std::vector<std::string> names, Parameters p, int nx, std::string filetype, int nt) {
    std::ostringstream fn;
    fn << "//home//marta//RKponiedz//rungekuttacomplex//script.txt";
    std::ofstream f(fn.str().c_str());

    f << "reset \n";
    f << "set pm3d map \n";
    f << "set xlabel \"x\" \n";
    f << "set ylabel \"t\" \n";
    f << "set terminal " << filetype << "\t" ;
    f << "size 1000, 1000 \n" ;
    f << "unset key\n";

    f << "set output \"AbsPsi." << filetype << "\"\n";
    f << "set title \" [Abs(Psi(x,t))]^2\" \n";
    f << "splot \"" << names[0] << "\" u (($1-" << nx/2 << ")/" << nx/ (2*p.getL()) << "):(($2)/"<< (nt/p.getTMax()) <<"):3 matrix w image\n";
    f << "unset title \n";
    f << "unset output \n";

    f << "set output \"RePsi." << filetype << "\" \n";
    f << "set title \" Re(Psi(x,t))\" \n";
   // f << "splot \""<< names[1] << "\" u (($2)/"<< p.getTMax() <<"):(($1-" << nx/2<<")/" << nx/ (2*p.getL())<<"):3 matrix w image \n";
    f << "splot \"" << names[1] << "\" u (($1-" << nx/2 << ")/" << nx/ (2*p.getL()) << "):(($2)/"<< (nt/p.getTMax()) <<"):3 matrix w image\n";
    f << "unset title \n";
    f << "unset output \n";

    f << "set output \"ImPsi." << filetype << "\" \n";
    f << "set title \" Im(Psi(x,t))\" \n";
    f << "splot \"" << names[2] << "\" u (($1-" << nx/2 << ")/" << nx/ (2*p.getL()) << "):(($2)/"<< (nt/p.getTMax()) <<"):3 matrix w image\n";
    f << "unset title \n";
    f << "unset output \n";


    f << "unset pm3d \n";
    f << "set grid \n";
    f << "set output \"Results." << filetype << "\" \n";
    f << "set xlabel \" x \" \n";
    f << "unset ylabel \n";
    f << "set multiplot layout 3,1 rowsfirst \n";
    f << "set title  \" [Abs(Psi(x,t = " << p.getTMax() <<"))]^2\" \n";
    f << "plot \"" << names[3] << "\" using 1:2 with lines\n";
    f << "set title  \" Re(Psi(x,t = " << p.getTMax() <<"))\" \n";
    f << "plot \"" << names[3] << "\" using 1:3 with lines\n";
    f << "set title  \" Im(Psi(x,t = " << p.getTMax() <<"))\" \n";
    f << "plot \"" << names[3] << "\" using 1:4 with lines\n";
    f << "unset multiplot \n";
    f << "unset output \n";

    f << "unset terminal \n";









}

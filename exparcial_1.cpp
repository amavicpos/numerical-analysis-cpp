#include <iostream>
#include <cmath>
#include <iomanip>

using namespace std;
double F (double x);
double Fprima (double x);
double bisection(double (*f) (double), double a, double b, double EP);
double newton(double (*f) (double), double (*fprima) (double), double x0, double tolerancia);
//double secant(double (*f) (double), double x1, double x2, double EP);

int main() {
	double a=1/sqrt(1e-4), b=1/sqrt(1.0), x0=1/sqrt(1e-4), EP=1e-8;
	cout << setprecision(7);
	cout << "Coeficiente de friccion segun el metodo de la biseccion: " << bisection(F, a, b, EP) << endl;
	cout << "Coeficiente de friccion segun el metodo de Newton: " << newton(F, Fprima, x0, EP) << endl;
	//cout << "Coeficiente de friccion segun el metodo de la secante: " << secant(F, a, b, EP) << endl;
	return 0;
}

double F (double x) {
	double K=0.25e-3, D=0.3, R_e=2e5, w, s;
	w=K/(D*3.71); s=2.51/R_e;
	return (-2.0*log10(w+s*x)-x);
}

double Fprima (double x) {
	double K=0.25e-3, D=0.3, R_e=2e5, w, s;
	w=K/D; s=2.51/R_e;
	return (-2.0/log(10)*s/(w+s*x)-1);
}

double bisection(double (*f) (double), double a, double b, double EP) {
 if (f(a) * f(b) >= 0) {
      cout << "No se pudo encontrar una raiz en el intervalo\n";
      return -1;
   }
   double c = a;
   int contador=0;
   while (fabs(f(c)) >= EP) {
      c = (a+b)/2;
      ++contador;
      if (f(c) == 0.0) break;
      else if (f(c)*f(a) < 0) b = c;
      else a = c;
   }
   cout << "Numero de iteraciones biseccion: " << contador << endl;
    return (1/pow(c,2));
}

double newton(double (*f) (double), double (*fprima) (double), double x0, double tolerancia) {
	double x1,xin;
	xin=x0;
	int contador=0;
	do { 	
	   if ((f(xin) !=0) && ( fprima(xin) !=0)) {
	     	x1 = xin -(f(xin)/fprima(xin));
	        xin = x1; 
	       	contador=contador+1;}
			else {
				cout << "No se pudo encontrar una raiz" << endl;
				return -1;
			}
	} while (fabs(f(x1)) > tolerancia);
	cout << endl << "Numero de iteraciones Newton: " << contador << endl;
    return (1/pow(x1,2));
}

/*
double secant(double (*f) (double), double x1, double x2, double EP) {
	float n = 0, xm, x0, c;
   if (f(x1) * f(x2) < 0) {
      do {
         x0 = (x1 * f(x2) - x2 * f(x1)) / (f(x2) - f(x1));
         c = f(x1) * f(x0);
         x1 = x2;
         x2 = x0;
         ++n;
         if (c == 0) {
         	cout << "c=0" << endl;
         	break;
		 }
         xm = (x1 * f(x2) - x2 * f(x1)) / (f(x2) - f(x1));
         cout << endl << fabs(xm - x0) << endl;
      } while (fabs(xm - x0) >= EP);
    return 1/pow(x1,2);
   } else {
   cout << "No se pudo encontrar una raiz en el intervalo\n";
   return -1; } } */

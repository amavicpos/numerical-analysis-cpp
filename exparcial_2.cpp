#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>

using namespace std;
const int n=10;
void lu(double A[n][n],double b[n],double x[n],int n);
double norma(double a[n], int n);
void jacobi(double a[n][n], double b[n], int n, double tolerancia);
void gauss_seidel(double a[n][n], double b[n], int n, double tolerancia);

int main() {
	double A[n][n], b[n]={0,1,2,3,4,5,6,7,8,9}, x[n];
	
	//Construyo la matriz A
	for (int i=0; i<=n-1; i++) {
	for (int j=0; j<=n-1; j++) {
		if (i==j)
		A[i][j]=5;
		else
		A[i][j]=0.5;
	}
	}
	
	cout << "Matriz A:" << endl;
	for (int i=0; i<=n-1; i++) {
	for (int j=0; j<=n-1; j++) {
		cout << A[i][j] << " ";
	}
	cout << endl;
	}
	cout << endl;
	
	//Solución por LU
	lu(A,b,x,n);
	cout << "x_LU = ";
	for (int i=0; i<n; ++i) {
			cout << setprecision(13) << x[i] << " ";
	}
	cout << endl << endl;
	
	ofstream fborrar;
	fborrar.open("resultados_parcial_1.txt", ofstream::out | ofstream::trunc);
	fborrar << "Tolerancia: Iteraciones (Jacobi):" << endl;
	fborrar.close();
	
	//Solución por Jacobi
	cout << "Solucion con Jacobi:" << endl;
	for (double tol = 1e-3; tol >= 1e-12; tol*=0.1) {
	cout << endl << "Para tolerancia = " << tol << endl;
	jacobi(A,b,n,tol);
	}
	
	ofstream fborrar2;
	fborrar2.open("resultados_parcial_2.txt", ofstream::out | ofstream::trunc);
	fborrar2 << "Tolerancia: Iteraciones (Gauss-Seidel):" << endl;
	fborrar2.close();
	
	//Solución por Gauss-Seidel
	cout << "Solucion con Gauss-Seidel:" << endl;
	for (double tol = 1e-3; tol >= 1e-12; tol=0.1*tol) {
	cout << endl << "Para tolerancia = " << tol << endl;
	gauss_seidel(A,b,n,tol);
	}
	
	
	return 0;
}

void lu(double A[n][n],double b[n],double x[n],int n)
{ // Resuelve sistema lineal por el metodo LU. Ax=b > LUx=b > Lz=b > Ux=z > x

double U[n][n], L[n][n], z[n];
double sum;
int i, j, k;
//Obtener L y U
for(i=0; i<=n-1; i++){ // i barre las filas
for(j=0; j<=n-1; j++){ // j barre las columnas
if(i<=j) { // elementos de la diagonal y por encima de la diagonal (i <= j)
sum=0.0;
for(k=0; k<=i-1; k++) {
sum=sum+L[i][k]*U[k][j];
}
U[i][j]=A[i][j]-sum;
if(i<j) {L[i][j]=0.0;}
else {L[i][j]=1.0;} }
else { // elementos por debajo de la diagonal (i > j)
sum=0.0;
for(k=0; k<=j-1; k++){
sum=sum+L[i][k]*U[k][j];
}
L[i][j]=(A[i][j]-sum)/U[j][j];
U[i][j]=0.0;
}}}
// Obtener z < Lz=b. SustituciÃ³n progresiva
for(i=0; i<=n-1; i++){
sum=0.0;
for(j=0; j<=i-1; j++){
sum=sum+L[i][j]*z[j];
}
z[i]=b[i]-sum;
}
// Obtener x < Ux=z. SustituciÃ³n regresiva
for(i=n-1 ;i>=0; i--){
sum=0.0;
for(j=i+1; j<=n-1; j++){
sum=sum+U[i][j]*x[j];
}
x[i]=(z[i]-sum)/U[i][i]; }}

double norma(double a[n], int n){
double aux=0;
for (int i = 0; i < n; i++)
aux += a[i] * a[i];

aux = sqrt(aux);
return aux; }

void jacobi(double a[n][n], double b[n], int n, double tolerancia){
	ofstream fsalida("resultados_parcial_1.txt", ios::app);
	double c[n], d[n]; // old x, x
   int i = 0, j = 0, l=0;
   for(int i=0; i<n; i++)
   {d[i]=0; } //initial guess
   do {
	for  (int i=0; i<n; i++){
	c[i]=d[i];}
	for (int i=0; i<n; i++){d[i]=b[i]/a[i][i];
	for (int j=0; j<n; j++){
	if (j!=i){d[i]=d[i]-(a[i][j]/a[i][i])*c[j];}}}
	l=l+1;
   } while ((fabs(norma(d, n)-norma(c,n)))>tolerancia);
   for (int i=0; i<n; i++)
   { cout<<"d"<<i+1<<"= "<< setprecision(13) << d[i]<<endl; }
   cout<<"Iteraciones realizadas con Jacobi: "<<l<<endl;
   fsalida << setprecision(1) << tolerancia << "      " << l << endl;
   fsalida.close();
}

void gauss_seidel(double a[n][n], double b[n], int n, double tolerancia){
	ofstream fsalida("resultados_parcial_2.txt", ios::app);
	double c[n], d[n];
	int i = 0, j = 0, l=0;
	for(int i=0; i<n; i++)
	{d[i]=0; }
   do {
	for  (int i=0; i<n; i++){
	c[i]=d[i];}
	for (int i=0; i<n; i++){d[i]=b[i]/a[i][i];
	for (int j=0; j<n; j++){
	if (j!=i){d[i]=d[i]-(a[i][j]/a[i][i])*d[j];}}}
	l=l+1;
   } while ((fabs(norma(d, n)-norma(c,n)))>tolerancia);
	for (int i=0; i<n; i++)
	{ cout<<"d"<<i+1<<"= "<< setprecision(13) << d[i]<<endl; }
	cout<<"Iteraciones realizadas con Gauss-Seidel: "<<l<<endl;
	fsalida << tolerancia << "      " << setprecision(1) << l << endl;
	fsalida.close();
}

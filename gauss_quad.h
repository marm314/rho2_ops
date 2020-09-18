#ifndef _GAUSS_HERM_QUADRATURE_H_
#define _GAUSS_HERM_QUADRATURE_H_

# include <cstdlib>
# include <cmath>
# include <iostream>
# include <fstream>
# include <iomanip>
# include <ctime>
# include <cstring>

using namespace std;

void gauss_hermite_rule(string filename,double &alpha, double &a, double &b,int &order);
void cdgqfgh ( int nt, int kind, double alpha, double beta, double t[],
  double wts[] );
void cgqfgh ( int nt, int kind, double alpha, double beta, double a, double b,
  double t[], double wts[] );
double class_matrixgh ( int kind, int m, double alpha, double beta, double aj[],
  double bj[] );
void imtqlxgh ( int n, double d[], double en[], double z[] );
void parchkgh ( int kind, int m, double alpha, double beta );
double r8_epsilongh ( );
double r8_hugegh( );
double r8_signgh ( double x );
void r8mat_writegh ( string output_filename, int m, int n, double table[] );
void rule_writegh ( int order, string filename, double x[], double w[],  double r[] );
void scqfgh ( int nt, double t[], int mlt[], double wts[], int nwts, int ndx[],double swts[], double st[], int kind, double alpha, double beta, double a,
  double b );
void sgqfgh ( int nt, double aj[], double bj[], double zemu, double t[],double wts[] );
#endif // _GAUSS_HERM_QUADRATURE_H_


#ifndef _MATHEMATICAL_FUNCTIONS_H_
#define _MATHEMATICAL_FUNCTIONS_H_

#include<iostream>
#include<stdlib.h>
#include<string>
#include<stdio.h>
#include<iomanip>
#include<cmath>
#include<complex.h>
#include<algorithm>
#include"Numbers.h"

using namespace std;
//////////////////////////
//Functions declaration //
//////////////////////////
double scalarP(double vec[3]);
double sign(double a,double b);
void second_order_polynom(double abc[3],double &root1,double &root2);//ax^2+bx+c=0
int factorial(int n);//factorial
int dfact(int n);//double factorial (2k-1)!!
void permutation(int *indexes,int order,int &permutation_term);//Permutes the indexes (Orders indexes first)
double norm3D(double vec[3]);//Norm for 3D vectors
double dot(int &n,double *a,double *b);//Dot product
void proyect(int &n,double *a,double *b,double *res);//Project vector a in b.
void print_mat(int &n,double **A);
void matmul(int n,double **A,double **B, double **RES);//RES=AxB
void mat_inverse(int n,double **A,double **Ainv);
void mat_transpose(int n,double **A, double **RES);//Transpose A->A^t==RES
void mat_equal(int n,double **A,double **RES);//RES=A
void jacobi(int n, double **m, double **v);//Diagonalization
double Gamma(double x);//Gamma function of x
double PochhammerS(double n,double x);//Rising Pochhammer symbol (x)n
double Y_lm(double l,double m,double teta, double phi);//Evaluates the spherical harmonic at teta and phi
double P_lm(double l, double m, double point);//Asociated legendre polynoms P(l,m,x)
double H_n(double n,double x);//Hermite polynomials
double L_n_alpha(double n,double alpha,double x);//Generalized Laguerre polynoms (alpha E Reals)
double T_n(double n,double x);//Chebyshev of first kind
double U_n(double n,double x);//Chebyshev of second kind
double J_v(double v, double x);//Bessel first kind function J_v(x)
double Y_v(double v, double x);//Bessel second kind function Y_v(x)
double j_n(double n,double x);//Spherical Bessel jn(x)
double y_n(double n,double x);//Spherical Bessel yn(x)
double Q_HARM_OSC_3D(int nx, int ny, int nz, double omega, double xyz[3]);//Sol for 3D Quantum Harmonic Oscillator
void nan(string in,bool &checked);//Check NAN
#endif // _MATHEMATICAL_FUNCTIONS_H_



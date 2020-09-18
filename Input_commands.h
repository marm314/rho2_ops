#ifndef _INPUT_COMMANDS_H_
#define _INPUT_COMMANDS_H_

#include<iostream>
#include<algorithm>
#include<fstream>
#include<stdlib.h>
#include<string>
#include<stdio.h>
#include<iomanip>
#include"String_ops.h"
#include"Numbers.h"

using namespace std;
//////////////////////////
//Functions declaration //
//////////////////////////
class Input
{
 public:
 int order_ang,order_ang2,nthreads;
 bool intracule,extracule,second_moments,reduce_terms,legendre,parallel,time_intra,rweight,time_extra;
 double radial_Init,radial_Step,radial_Last,radial_grid,threshold_in,tau,lambda_rs,lambda_scr;
 string name_dm2,name_basis;
 Input();
 Input(string);
 ~Input();
};
#endif // _INPUT_COMMANDS_H_

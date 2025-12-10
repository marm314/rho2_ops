#ifndef _INPUT_COMMANDS_DIRACINTFACE_H_
#define _INPUT_COMMANDS_DIRACINTFACE_H_

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
class Input_diracintface
{
 public:
 bool is_gaunt,large_mem,oneMOwfx,transf_cplx,print_coef_files,rotate_krammers;
 int OneMO_wfx;
 double threshold,maxmem;
 string name_dm2,name_out;
 Input_diracintface();
 Input_diracintface(string);
 ~Input_diracintface();
};
#endif // _INPUT_COMMANDS_DIRACINTFACE_H_

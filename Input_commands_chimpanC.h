#ifndef _INPUT_COMMANDS_CHIMPANC_H_
#define _INPUT_COMMANDS_CHIMPANC_H_

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
class Input_chimpanC
{
 public:
 bool two_dm2_mat,index_iiii,reduce,reduce_sym,all_dm2_are_given,donof,int8alldm2;
 double threshold;
 string name_dm2,name_fchk;
 Input_chimpanC();
 Input_chimpanC(string);
 ~Input_chimpanC();
};
#endif // _INPUT_COMMANDS_CHIMPANC_H_

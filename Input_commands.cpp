#include"Input_commands.h"
Input::Input(){cout<<"Not allowed default constructor"<<endl;}
Input::Input(string I_in)
{
 order_ang=0;order_ang2=0;nthreads=1;radial_Last=1e99;
 intracule=false;legendre=false;parallel=false;
 time_intra=false;second_moments=false;rweight=false;extracule=false;
 time_extra=false;sym_red=false;
 tau=ZERO;
 lambda_rs=1.0e99; 
 lambda_scr=0.0e0; 
 ifstream I_input_file;
 I_in.erase(std::remove_if(I_in.begin(),I_in.end(),::isspace),I_in.end());
 I_input_file.open((I_in).c_str());
 if(I_input_file.good()) //Check existence of file
 {
  while(getline(I_input_file,I_in))
  {
   I_in.erase(std::remove_if(I_in.begin(),I_in.end(),::isspace),I_in.end());
   lowercase(I_in);
   if(I_in=="$name_basis")
   {
    getline(I_input_file,I_in);
    I_in.erase(std::remove_if(I_in.begin(),I_in.end(),::isspace),I_in.end());
    name_basis=I_in;
   }
   else if(I_in=="$name_dm2")
   {
    getline(I_input_file,I_in);
    I_in.erase(std::remove_if(I_in.begin(),I_in.end(),::isspace),I_in.end());
    name_dm2=I_in;
   }
   else if(I_in=="$intracule"){intracule=true;second_moments=false;}
   else if(I_in=="$moments"){second_moments=true;intracule=false;extracule=false;}
   else if(I_in=="$extracule"){extracule=true;second_moments=false;}
   else if(I_in=="$symred"){sym_red=true;}
   else if(I_in=="$threshold"){I_input_file>>threshold_in;}
   else if(I_in=="$rangesep"){I_input_file>>lambda_rs;}
   else if(I_in=="$screening"){I_input_file>>lambda_scr;}
   else if(I_in=="$parallel")
   {
    parallel=true;
    I_input_file>>nthreads;
   }
   else if(I_in=="$tau"){I_input_file>>tau;}
   else if(I_in=="$radial"){I_input_file>>radial_Init>>radial_Step>>radial_Last;}
   else if(I_in=="$legendre")
   {
    legendre=true;
    I_input_file>>radial_grid>>radial_Last;
   }
   else if(I_in=="$angular"){I_input_file>>order_ang;}
   else if(I_in=="$angular2")
   {
    I_input_file>>order_ang;
    I_input_file>>order_ang2;
   }
   else if(I_in=="$rweights"){rweight=true;}
   else if(I_in=="$time"){time_intra=true;time_extra=true;}
   else{}
  }
  if(radial_Init<ZERO){radial_Init=ZERO;}
 }
 I_input_file.close();
}
Input::~Input()
{}

#include"Input_commands_diracintface.h"
Input_diracintface::Input_diracintface(){cout<<"Not allowed default constructor"<<endl;}
Input_diracintface::Input_diracintface(string I_in)
{
 is_gaunt=false;large_mem=false;oneMOwfx=false;transf_cplx=false;print_coef_files=false;
 rotate_krammers=false;only_wfx=false;
 threshold=pow(TEN,-TEN);
 maxmem=ONE;
 ifstream I_input_file;
 I_in.erase(std::remove_if(I_in.begin(),I_in.end(),::isspace),I_in.end());
 I_input_file.open((I_in).c_str());
 if(I_input_file.good()) //Check existence of file
 {
  while(getline(I_input_file,I_in))
  {
   I_in.erase(std::remove_if(I_in.begin(),I_in.end(),::isspace),I_in.end());
   lowercase(I_in);
   if(I_in=="$name_dm2")
   {
    getline(I_input_file,I_in);
    I_in.erase(std::remove_if(I_in.begin(),I_in.end(),::isspace),I_in.end());
    name_dm2=I_in;
   }
   else if(I_in=="$name_out")
   {
    getline(I_input_file,I_in);
    I_in.erase(std::remove_if(I_in.begin(),I_in.end(),::isspace),I_in.end());
    name_out=I_in;
   }
   else if(I_in=="$large_mem"){large_mem=true;}
   else if(I_in=="$gaunt"){is_gaunt=true;}
//   else if(I_in=="$transf_cplx"){transf_cplx=true;} // Currently switched off
   else if(I_in=="$print_coef_files"){print_coef_files=true;}
   else if(I_in=="$threshold")
   {
    I_input_file>>threshold;    
   }
   else if(I_in=="$onemowfx")
   {
    oneMOwfx=true;
    I_input_file>>OneMO_wfx;    
   }
   else if(I_in=="$maxmem")
   {
    I_input_file>>maxmem;    
   }
   else if(I_in=="$rotate_krammers"){rotate_krammers=true;}
   else if(I_in=="$only_wfx"){only_wfx=true;}
   else
   {
    // Nth
   }
  }
 }
 I_input_file.close();
}
Input_diracintface::~Input_diracintface()
{}

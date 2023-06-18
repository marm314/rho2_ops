#include"Input_commands_chimpanC.h"
Input_chimpanC::Input_chimpanC(){cout<<"Not allowed default constructor"<<endl;}
Input_chimpanC::Input_chimpanC(string I_in)
{
 large_mem=false,index_iiii=false,reduce=false,reduce_sym=false,all_dm2_are_given=false;
 int8alldm2=false,donof=false;method1=true,sl2rdm=false;
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
   else if(I_in=="$name_fchk")
   {
    getline(I_input_file,I_in);
    I_in.erase(std::remove_if(I_in.begin(),I_in.end(),::isspace),I_in.end());
    name_fchk=I_in;
   }
   else if(I_in=="$large_mem"){large_mem=true;}
   else if(I_in=="$index_iiii"){index_iiii=true;}
   else if(I_in=="$reduce"){reduce=true;}
   else if(I_in=="$reduce_sym"){reduce_sym=true;}
   else if(I_in=="$all_dm2_are_given"){all_dm2_are_given=true;}
   else if(I_in=="$int8alldm2"){int8alldm2=true;}
   else if(I_in=="$only_trace"){only_trace=true;}
   else if(I_in=="$donof"){donof=true;}
   else if(I_in=="$sl2rdm"){sl2rdm=true;}
   else if(I_in=="$method2"){method1=false;}
   else if(I_in=="$threshold")
   {
    I_input_file>>threshold;    
   }
   else if(I_in=="$maxmem")
   {
    I_input_file>>maxmem;    
   }
   else
   {
    // Nth
   }
  }
 }
 I_input_file.close();
}
Input_chimpanC::~Input_chimpanC()
{}

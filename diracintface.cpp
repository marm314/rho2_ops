#include<iostream>
#include<algorithm>
#include<fstream>
#include<stdlib.h>
#include<string>
#include<stdio.h>
#include<iomanip>
#include<vector>
#include<complex>
#include"String_ops.h"
#include"Numbers.h"
#include"Input_commands_diracintface.h"

using namespace std;

// Functions
void print_basis_file();
void read_dirac_out();
void clean_shell2aos();
void print_wfx();
void read_2rdm4cMO_and_transf();
void dm2_4c2LS(int &index_4cMO1,int &index_4cMO2,int &index_4cMO3,int &index_4cMO4,double Dijkl_4cMO); 
void transform_Dijkl2Dpqrs();
void reduce_getreal_print();
void reduce_print();
//void transform_Dijkl2Dpqrs_cplx(); // Currently switched off because we do not need the imaginary part of the 2-RDM in primitives 

// Global variables
complex<double>CZERO=(ZERO,ZERO);
bool transf_cplx=false,large_mem=false;//symmrr_prime=true;
int Ncenters,Nprimitives,Nbasis,Nbasis_L,Nbasis_S,Nshell,Nshell_L,Nshell_S,NMOs,NMOs_LS,NMOs_occ,Largest_Prim;
int OneMO_wfx=-1;
long int NMOs_LS_1,NMOs_LS_2,NMOs_LS_3,NMOs_LS_4,Nterms_printed=0;
long int Nprims1,Nprims2,Nprims3,Nprims4;       // Powers of Nprimitives
int RECORD_DELIMITER_LENGTH=4;
struct Shell2AOs
{
 int styp,atom,nprim,naos,paired2=-1;
 int **index_min2max; 
 double *Coef,*Expon,Coord[3];
};
Shell2AOs *shell2aos;
double threshold,maxmem,sumMEM=ZERO;
double Quaternion_coef[4];
double *OCCs,**Prim2AO_Coef;
double ****Dpqrs_ALL;
complex<double> **AO2MO_Coef,**Prim2MO_Coef;
complex<double> *Dpqrs_ALL_cplx;
double *Dijkl_MOsLS;
string dirac_output_name,dirac_output_file,dm2_file;
vector<int>Prim2Center_map;
vector<int>shell_types;
vector<double>Center_x;
vector<double>Center_y;
vector<double>Center_z;
vector<double>prim_exponents;
vector<double>MOsLS_occ;
vector<double>One_Prim2MO_Coef_RE;
vector<double>One_Prim2MO_Coef_IM;

// Main  
int main(int argc, char *argv[])
{
 cout<<"--------------------------------------------"<<endl;
 cout<<"--------------------------------------------"<<endl;
 cout<<"---  DIRAC+DM2 INTERFACE FOR RHO2_OPS    ---"<<endl;
 cout<<"--------------------------------------------"<<endl;
 cout<<"--------------------------------------------"<<endl;
 cout<<"-- Developed by: Dr. M. Rodriguez-Mayorga --"<<endl;
 cout<<"--     email: marm3.14@gmail.com          --"<<endl;
 cout<<"--------------------------------------------"<<endl;
 cout<<"--------------------------------------------"<<endl;
 cout<<endl;
 if(argc!=2)
 {
  cout<<endl;
  cout<<"Please, include the input file"<<endl;
  cout<<endl;
  cout<<endl;
  cout<<"--------------------------------------------"<<endl;
  cout<<"                       .-.                  "<<endl;  
  cout<<"                      |_:_|                 "<<endl; 
  cout<<"                     /(_Y_)\\               "<<endl; 
  cout<<".                   ( \\/M\\/ )             "<<endl;
  cout<<" '.               _.'-/'-'\\-'._            "<<endl;
  cout<<"   ':           _/.--'[[[[]'--.\\_          "<<endl;
  cout<<"     ':        /_'  : |::\"| :  '.\\        "<<endl;
  cout<<"       ':     //   ./ |oUU| \\.'  :\\       "<<endl;
  cout<<"         ':  _:'..' \\_|___|_/ :   :|       "<<endl;
  cout<<"           ':.  .'  |_[___]_|  :.':\\       "<<endl;
  cout<<"            [::\\ |  :  | |  :   ; : \\     "<<endl;
  cout<<"             '-'   \\/'.| |.' \\  .;.' |    "<<endl;
  cout<<"             |\\_    \\  '-'   :       |    "<<endl;
  cout<<"             |  \\    \\ .:    :   |   |    "<<endl;
  cout<<"             |   \\    | '.   :    \\  |    "<<endl;
  cout<<"             /       \\   :. .;       |     "<<endl;
  cout<<"            /     |   |  :__/     :  \\\\   "<<endl;
  cout<<"           |  |   |    \\:   | \\   |   ||  "<<endl;
  cout<<"          /    \\  : :  |:   /  |__|   /|   "<<endl;
  cout<<"      snd |     : : :_/_|  /'._\\  '--|_\\  "<<endl;
  cout<<"          /___.-/_|-'   \\  \\              "<<endl;  
  cout<<"                         '-'                "<<endl;
  cout<<"                                            "<<endl; 
  cout<<"  'I find your lack of faith disturbing.'   "<<endl; 
  cout<<"                             Dart Vader     "<<endl; 
  cout<<"                                            "<<endl; 
  cout<<"--------------------------------------------"<<endl;
  cout<<"--          Normal termination            --"<<endl;
  cout<<"--------------------------------------------"<<endl;
  return -1;
 } 
 bool repeated_prims;
 int ishell,ishell1,iprim,iprim1,iaos,iaos1,imos,imos1,imos2; // imos for Scalar MOs
 int naos;
 string aux(argv[1]);
 Input_diracintface Input_commands(aux);
 dirac_output_file=Input_commands.name_out;
 dm2_file=Input_commands.name_dm2;
 threshold=Input_commands.threshold;
 dirac_output_name=dirac_output_file.substr(0,dirac_output_file.length()-3);
 maxmem=Input_commands.maxmem;
 maxmem=maxmem*pow(TEN,NINE); // Gb to bytes
 if(Input_commands.oneMOwfx)
 {
  OneMO_wfx=Input_commands.OneMO_wfx;
  OneMO_wfx=(OneMO_wfx-1)*4; // even->unbar, odd->bar
  cout<<"Orbital selection is swittched on. Scalar (LS) orbitals to print in the WFX file: ";
  cout<<setw(5)<<OneMO_wfx+1<<" to "<<setw(5)<<OneMO_wfx+4<<endl;
 }
 if(Input_commands.transf_cplx)
 {
  cout<<"Warning! Currently the transf. of the Scalar (LS) 2-RDM to a complex 2-RDM in primitives is switched off!"<<endl;
 //transf_cplx=Input_commands.transf_cplx;
 }
 if(Input_commands.large_mem){large_mem=Input_commands.large_mem;}
 // TODO: 
 /*
  symmrr_prime=false;
 */
 // Read Dirac output
 read_dirac_out();
 // Find repeated shells ("pairing")
 for(ishell=0;ishell<Nshell-1;ishell++)
 {
  for(ishell1=ishell+1;ishell1<Nshell;ishell1++)
  {
   repeated_prims=false;
   if(shell2aos[ishell1].styp==shell2aos[ishell].styp && shell2aos[ishell1].nprim==shell2aos[ishell].nprim)
   {
    if(shell2aos[ishell1].Coord[0]==shell2aos[ishell].Coord[0] && 
       shell2aos[ishell1].Coord[1]==shell2aos[ishell].Coord[1] && 
       shell2aos[ishell1].Coord[2]==shell2aos[ishell].Coord[2])
    {
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      if(shell2aos[ishell1].Expon[iprim]==shell2aos[ishell].Expon[iprim])
      {
       repeated_prims=true;
      }
      else
      {
       repeated_prims=false;
       iprim=shell2aos[ishell].nprim;
      }
     }
    }
   }
   // Save pairing info
   if(repeated_prims && shell2aos[ishell1].paired2==-1)
   {
    shell2aos[ishell1].paired2=ishell;
   }
  }
 }
 // Print .basis file with unique primitives information
 print_basis_file();
 cout<<"Number of Primitives : "<<setw(12)<<Nprimitives<<endl;
 // Find the total number of AOs and allocate indexing arrays
 naos=0;
 for(ishell=0;ishell<Nshell;ishell++)
 {
  naos=naos+shell2aos[ishell].naos;
  shell2aos[ishell].index_min2max=new int*[shell2aos[ishell].naos];
  for(iaos=0;iaos<shell2aos[ishell].naos;iaos++)
  {
   shell2aos[ishell].index_min2max[iaos]=new int[shell2aos[ishell].nprim];
  }
 }
 cout<<"Number of AOs (calc) : "<<setw(12)<<naos<<endl;
 // Allocate the primitive coefs. matrix (used to build AOs from primitives)
 Prim2AO_Coef=new double*[naos];
 for(iaos=0;iaos<naos;iaos++)
 {
  Prim2AO_Coef[iaos]=new double[Nprimitives];
  for(iprim=0;iprim<Nprimitives;iprim++)
  {
   Prim2AO_Coef[iaos][iprim]=ZERO;
  }
 }
 // Initialize the indexing for the "unpaired" shells
 iprim1=0;
 for(ishell=0;ishell<Nshell;ishell++)
 {
  if(shell2aos[ishell].paired2==-1)
  {
   for(iaos=0;iaos<shell2aos[ishell].naos;iaos++)
   {
    for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
    {
     shell2aos[ishell].index_min2max[iaos][iprim]=iprim1;
     iprim1++;
    }
   }
  }
 }
 // Copy indexing from "unpaired" shells to the "paired" ones
 for(ishell=0;ishell<Nshell;ishell++)
 {
  if(shell2aos[ishell].paired2!=-1)
  {
   for(iaos=0;iaos<shell2aos[ishell].naos;iaos++)
   {
    for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
    {
     shell2aos[ishell].index_min2max[iaos][iprim]=shell2aos[shell2aos[ishell].paired2].index_min2max[iaos][iprim];
    }
   }
  } 
 }
 // Organize Prim 2 AO Coefs matrix using index_min2max
 cout<<endl;
 cout<<"AO to Primitives map (per shell):"<<endl;
 iaos1=1;
 for(ishell=0;ishell<Nshell;ishell++)
 {
  for(iaos=0;iaos<shell2aos[ishell].naos;iaos++)
  {
   cout<<"AO"<<setw(4)<<iaos1<<":";
   for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
   {
    cout<<"\t"<<shell2aos[ishell].index_min2max[iaos][iprim]+1;
   }
   cout<<endl;
   iaos1++;
  }
 }
 // Fill the Prim2AO_Coef matrix
 iaos1=0;
 for(ishell=0;ishell<Nshell;ishell++)
 {
  for(iaos=0;iaos<shell2aos[ishell].naos;iaos++)
  {
   for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
   {
    Prim2AO_Coef[iaos1][shell2aos[ishell].index_min2max[iaos][iprim]]=shell2aos[ishell].Coef[iprim];
   }
   iaos1++;
  }
 }
 // Deallocate arrays that are not used anymore
 clean_shell2aos();
 // Compute Primitives to MO coefficients
 cout<<endl;
 cout<<"Building the matrix Prim2MO_Coefs[MO][Primitives] ((2x4MO) x 2(L+S)) = AO2MO_Coef[MO][AO]*Prim2AO_Coef[AO][Primitives]"<<endl;
 Prim2MO_Coef=new complex<double>*[NMOs_LS]; 
 for(imos=0;imos<NMOs_LS;imos++)
 {
  Prim2MO_Coef[imos]=new complex<double>[Nprimitives];
  for(iprim=0;iprim<Nprimitives;iprim++)
  {
   Prim2MO_Coef[imos][iprim]=CZERO;
   for(iaos=0;iaos<Nbasis;iaos++)
   {
    Prim2MO_Coef[imos][iprim]=Prim2MO_Coef[imos][iprim]+AO2MO_Coef[imos][iaos]*Prim2AO_Coef[iaos][iprim];
   }
  }
 }
 cout<<"Matrix Prim2MO_Coefs[MO][Primitives] for bar and unbar orbs. built."<<endl;
 cout<<endl;
 // Deallocate arrays related to AOs
 for(imos=0;imos<NMOs_LS;imos++)
 {
  delete [] AO2MO_Coef[imos];AO2MO_Coef[imos]=NULL;
 }
 delete[] AO2MO_Coef;AO2MO_Coef=NULL;
 for(iaos=0;iaos<naos;iaos++)
 {
  delete[] Prim2AO_Coef[iaos];Prim2AO_Coef[iaos]=NULL;
 }
 delete[] Prim2AO_Coef;Prim2AO_Coef=NULL;
 // Print the Prim2MO_Coef matrix (coefficients are rows). 
 // WARNING! Below we may overwrite the positronic states (initial ones) with the occ. electronic states (i.e. put them on top of the list)!
 ofstream coefs_file;
 ofstream coefs_file_pos;
 if(Input_commands.print_coef_files)
 {
  coefs_file.open((dirac_output_name+"electronic.coef").c_str());
  coefs_file_pos.open((dirac_output_name+"positronic.coef").c_str());
  coefs_file<<setprecision(12)<<fixed<<scientific;
  coefs_file_pos<<setprecision(12)<<fixed<<scientific;
 }
 imos1=0;
 for(imos=0;imos<NMOs_LS;imos++)
 {
  if(OneMO_wfx!=-1 && OneMO_wfx==imos) // Check if we have to store this particular MO (plus the rest of components) for the WFX file
  {
   for(imos2=0;imos2<4;imos2++)
   {
    for(iprim=0;iprim<Nprimitives;iprim++)
    {
     One_Prim2MO_Coef_RE.push_back(Prim2MO_Coef[imos+imos2][iprim].real());
     One_Prim2MO_Coef_IM.push_back(Prim2MO_Coef[imos+imos2][iprim].imag());
    }
   }
   cout<<"Scalar (LS) orbitals selected for the WFX file from "<<setw(5)<<imos+1<<" to "<<setw(5)<<imos+4<<endl; 
  }
  if(Input_commands.print_coef_files)
  {
   if(OCCs[imos]==-TEN) // Write positronic MO coefs to file
   {
    for(iprim=0;iprim<Nprimitives;iprim++)
    {
     coefs_file_pos<<setw(20)<<Prim2MO_Coef[imos][iprim].real()<<setw(20)<<Prim2MO_Coef[imos][iprim].imag();
    }
   }
  }
  if(OCCs[imos]>pow(TEN,-EIGHT)) // Overwrite positronic MO coefs with electronic ones and print the electronic ones
  {
   for(iprim=0;iprim<Nprimitives;iprim++)
   {
    Prim2MO_Coef[imos1][iprim]=Prim2MO_Coef[imos][iprim];
    if(Input_commands.print_coef_files)
    {
     coefs_file<<setw(20)<<Prim2MO_Coef[imos1][iprim].real()<<setw(20)<<Prim2MO_Coef[imos1][iprim].imag();
    }
   }
   MOsLS_occ.push_back(OCCs[imos]);
   imos1++;
  }
 }
 NMOs_occ=imos1;
 if(Input_commands.print_coef_files)
 {
  coefs_file.close();
  coefs_file_pos.close();
 }
 cout<<"Num. of occ MO (Scalar): "<<setw(12)<<NMOs_occ<<endl;
 // Print a WFX file for RHO_OPS (currently only available for ATOMS) and delete OCCs array.
 if(OneMO_wfx==-1)
 {
  print_wfx();
 }
 else
 {
  if(OneMO_wfx+4<=NMOs_LS)
  {
   print_wfx();
  }
  else
  {
   cout<<"Unable to print the (un)bar MO orb: "<<setw(5)<<(OneMO_wfx/4)+1<<" is not present."<<endl;
  }
 }
 delete[] OCCs;OCCs=NULL;
 // Read 2-RDM (binary) and store it in the scalar (2(L+S)) MO basis (using only Scalar occupied MOs for the size of the 2-RDM at this stage!).
 ifstream check_dm2;
 check_dm2.open(dm2_file);
 if(check_dm2.good())
 {
  check_dm2.close();
  read_2rdm4cMO_and_transf();
 }
 else
 {
  check_dm2.close();
  cout<<endl;
  cout<<"Unable to find the file: "<<dm2_file<<endl;
  cout<<"Not transforming the 2-RDM matrix"<<endl;
  cout<<endl;
 }
 // Deallocate array Primitive to MO coefs (rows).
 for(imos=0;imos<NMOs_LS;imos++)
 {
  delete [] Prim2MO_Coef[imos];Prim2MO_Coef[imos]=NULL;
 }
 delete[] Prim2MO_Coef;Prim2MO_Coef=NULL;
 cout<<"--------------------------------------------"<<endl;
 cout<<"                       .-.                  "<<endl;
 cout<<"                      |_:_|                 "<<endl;
 cout<<"                     /(_Y_)\\               "<<endl;
 cout<<".                   ( \\/M\\/ )             "<<endl;
 cout<<" '.               _.'-/'-'\\-'._            "<<endl;
 cout<<"   ':           _/.--'[[[[]'--.\\_          "<<endl;
 cout<<"     ':        /_'  : |::\"| :  '.\\        "<<endl;
 cout<<"       ':     //   ./ |oUU| \\.'  :\\       "<<endl;
 cout<<"         ':  _:'..' \\_|___|_/ :   :|       "<<endl;
 cout<<"           ':.  .'  |_[___]_|  :.':\\       "<<endl;
 cout<<"            [::\\ |  :  | |  :   ; : \\     "<<endl;
 cout<<"             '-'   \\/'.| |.' \\  .;.' |    "<<endl;
 cout<<"             |\\_    \\  '-'   :       |    "<<endl;
 cout<<"             |  \\    \\ .:    :   |   |    "<<endl;
 cout<<"             |   \\    | '.   :    \\  |    "<<endl;
 cout<<"             /       \\   :. .;       |     "<<endl;
 cout<<"            /     |   |  :__/     :  \\\\   "<<endl;
 cout<<"           |  |   |    \\:   | \\   |   ||  "<<endl;
 cout<<"          /    \\  : :  |:   /  |__|   /|   "<<endl;
 cout<<"      snd |     : : :_/_|  /'._\\  '--|_\\  "<<endl;
 cout<<"          /___.-/_|-'   \\  \\              "<<endl;
 cout<<"                         '-'                "<<endl;
 cout<<"                                            "<<endl;
 cout<<"  'I find your lack of faith disturbing.'   "<<endl;
 cout<<"                             Dart Vader     "<<endl;
 cout<<"                                            "<<endl;
 cout<<"--------------------------------------------"<<endl;
 cout<<"--          Normal termination            --"<<endl;
 cout<<"--------------------------------------------"<<endl;
 return 0;
}

// Function used to read DIRAC output
void read_dirac_out()
{
 bool electronic,first_ao2mos_coefs=true;
 int ishell,iprim,imos,iaos,Nbasis_int;
 int imos1,imos2,imos3,imos4,imos5,imos6,imos7,imos8;
 double occ;
 string line;
 ifstream Dirac_file(dirac_output_file);
 while(getline(Dirac_file,line))
 {
  if(line.length()>8)
  {
   // Nshells info
   if(line.substr(0,8)==" Nshells")
   {
    if(line.substr(8,6)==" Large")
    {
     line=line.substr(14,line.length()-13);
     stringstream ss(line);
     ss>>Nshell_L;
    }
    else if(line.substr(8,6)==" Small")
    {
     line=line.substr(14,line.length()-13);
     stringstream ss(line);
     ss>>Nshell_S;
    }
    else
    {
     line=line.substr(14,line.length()-13);
     stringstream ss(line);
     ss>>Nshell;
     shell2aos=new Shell2AOs[Nshell];
    }
   }
   // Nbasis info
   if(line.substr(0,8)==" Nbasis ")
   {
    if(line.substr(8,6)==" Large")
    {
     line=line.substr(14,line.length()-13);
     stringstream ss(line);
     ss>>Nbasis_L;
    }
    else if(line.substr(8,6)==" Small")
    {
     line=line.substr(14,line.length()-13);
     stringstream ss(line);
     ss>>Nbasis_S;
    }
    else
    {
     line=line.substr(14,line.length()-13);
     stringstream ss(line);
     ss>>Nbasis;
    }
   }
  }
  // Shell types
  if(line.length()==12)
  {
   if(line.substr(0,12)==" Shell type:")
   {
    for(ishell=0;ishell<Nshell;ishell++)
    {
     Dirac_file>>shell2aos[ishell].styp;
     if(shell2aos[ishell].styp== 0)
     {
      shell2aos[ishell].naos=1;
     }
     else if(shell2aos[ishell].styp== 1)
     {
      shell2aos[ishell].naos=3;
     }
     else if(shell2aos[ishell].styp==-2)
     {
      shell2aos[ishell].naos=6;
     }
     else if(shell2aos[ishell].styp==-3)
     {
      shell2aos[ishell].naos=10;
     }
     else if(shell2aos[ishell].styp==-4)
     {
      shell2aos[ishell].naos=15;
     }
     else
     {
      cout<<"Warning! Shell type not supported "<<setw(5)<<shell2aos[ishell].styp<<endl;
     }
    } 
   }
  } 
  // Shell to atom map
  if(line.length()==19)
  {
   if(line.substr(0,19)==" Shell to atom map:")
   {
    for(ishell=0;ishell<Nshell;ishell++)
    {
     Dirac_file>>shell2aos[ishell].atom;
    }
   }
  }
  // Nprim_in per shell
  if(line.length()==32)
  {
   if(line.substr(0,32)==" Number of primitives per shell:")
   {
    for(ishell=0;ishell<Nshell;ishell++)
    {
     Dirac_file>>shell2aos[ishell].nprim;
     shell2aos[ishell].Coef=new double[shell2aos[ishell].nprim];
     shell2aos[ishell].Expon=new double[shell2aos[ishell].nprim];
    }
   }
  }
  // Contraction coefs
  if(line.length()==26)
  {
   if(line.substr(0,26)==" Contraction coefficients:")
   {
    for(ishell=0;ishell<Nshell;ishell++)
    {
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      Dirac_file>>shell2aos[ishell].Coef[iprim];
     }
    }
   }
  }
  // Exponents
  if(line.length()==21)
  {
   if(line.substr(0,21)==" Primitive exponents:")
   {
    for(ishell=0;ishell<Nshell;ishell++)
    {
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      Dirac_file>>shell2aos[ishell].Expon[iprim];
     }
    }
   }
  }
  // Coord shell
  if(line.length()==37)
  {
   if(line.substr(0,37)==" Cartesian coordinates of each shell:")
   {
    for(ishell=0;ishell<Nshell;ishell++)
    {
     Dirac_file>>shell2aos[ishell].Coord[0];
     Dirac_file>>shell2aos[ishell].Coord[1];
     Dirac_file>>shell2aos[ishell].Coord[2];
    }
   }
  }
  if(line.length()>6)
  {
   if(line.substr(0,6)==" Norbs")
   {
    line=line.substr(6,line.length()-6);
    stringstream ss(line);
    ss>>NMOs;
   }
  }
  // Read occ numbers
  if(line.length()>13)
  {
   if(line.substr(0,13)==" Occupancies:")
   {
    line=line.substr(13,line.length()-13);
    stringstream ss(line);
    ss>>NMOs;
    OCCs=new double[8*NMOs];
    for(imos=0;imos<8*NMOs;imos++)
    {
     OCCs[imos]=-TEN;
    }
    electronic=false;
    for(imos=0;imos<NMOs;imos++)
    {
     Dirac_file>>occ;
     if(abs(occ)>pow(TEN,-EIGHT))
     {
      electronic=true;
      imos1=8*imos;
      imos2=imos1+1;imos3=imos2+1;imos4=imos3+1;imos5=imos4+1;imos6=imos5+1;imos7=imos6+1;imos8=imos7+1;
      OCCs[imos1]=occ;OCCs[imos2]=occ;OCCs[imos3]=occ;OCCs[imos4]=occ;OCCs[imos5]=occ;OCCs[imos6]=occ;OCCs[imos7]=occ;OCCs[imos8]=occ;
     }
     if(electronic && abs(occ)<=pow(TEN,-EIGHT))
     {
      occ=ZERO;
      imos1=8*imos;
      imos2=imos1+1;imos3=imos2+1;imos4=imos3+1;imos5=imos4+1;imos6=imos5+1;imos7=imos6+1;imos8=imos7+1;
      OCCs[imos1]=occ;OCCs[imos2]=occ;OCCs[imos3]=occ;OCCs[imos4]=occ;OCCs[imos5]=occ;OCCs[imos6]=occ;OCCs[imos7]=occ;OCCs[imos8]=occ;
     }
    }
   }
  }
  // Read quaternion MO coefs
  if(line.length()>7)
  {
   if(line.substr(0,7)==" NBasis" && first_ao2mos_coefs)
   {
    first_ao2mos_coefs=false;
    line=line.substr(7,line.length()-7);
    stringstream ss(line);
    ss>>Nbasis_int;
    if(Nbasis!=Nbasis_int)
    {
     cout<<"Warning! Nbasis read is not equal to the calculated value."<<endl;
    }
    cout<<"Number of Quat. MOs  : "<<setw(12)<<NMOs<<endl;
    NMOs_LS=NMOs*8;
    cout<<"Number of Scalar MOs : "<<setw(12)<<NMOs_LS<<endl;
    cout<<"Number of AOs (read) : "<<setw(12)<<Nbasis_int<<endl;
    AO2MO_Coef=new complex<double>*[NMOs_LS];
    for(imos=0;imos<NMOs_LS;imos++)
    {
     AO2MO_Coef[imos]=new complex<double>[Nbasis];
     for(iaos=0;iaos<Nbasis;iaos++)
     {
      AO2MO_Coef[imos][iaos]=CZERO;
     }
    }
    getline(Dirac_file,line);
    //cout<<line<<endl; 
    for(imos=0;imos<NMOs;imos++)
    {
     imos1=imos*8;imos2=imos1+1;imos3=imos2+1;imos4=imos3+1;imos5=imos4+1;imos6=imos5+1;imos7=imos6+1;imos8=imos7+1;
     //cout<<"Quaternion to Large component MOs transformation for MO"<<setw(5)<<imos+1<<endl;
     for(iaos=0;iaos<Nbasis_L;iaos++)
     {
      Dirac_file>>Quaternion_coef[0]>>Quaternion_coef[1]>>Quaternion_coef[2]>>Quaternion_coef[3];
      //cout<<Quaternion_coef[0]<<Quaternion_coef[1]<<Quaternion_coef[2]<<Quaternion_coef[3]<<endl;
      complex<double>ztmp4c1( Quaternion_coef[0], Quaternion_coef[1]); 
      complex<double>ztmp4c2(-Quaternion_coef[2], Quaternion_coef[3]); 
      complex<double>ztmp4c3( Quaternion_coef[2], Quaternion_coef[3]); 
      complex<double>ztmp4c4( Quaternion_coef[0],-Quaternion_coef[1]); 
      AO2MO_Coef[imos1][iaos]=ztmp4c1; // unbar alpha L
      AO2MO_Coef[imos2][iaos]=ztmp4c2; // unbar beta L
      AO2MO_Coef[imos5][iaos]=ztmp4c3; // bar alpha L
      AO2MO_Coef[imos6][iaos]=ztmp4c4; // bar beta L
     }
     //cout<<"Quaternion to Small component MOs transformation for MO"<<setw(5)<<imos+1<<endl;
     for(iaos=Nbasis_L;iaos<Nbasis;iaos++)
     {
      Dirac_file>>Quaternion_coef[0]>>Quaternion_coef[1]>>Quaternion_coef[2]>>Quaternion_coef[3];
      //cout<<Quaternion_coef[0]<<Quaternion_coef[1]<<Quaternion_coef[2]<<Quaternion_coef[3]<<endl; 
      complex<double>ztmp4c1( Quaternion_coef[0], Quaternion_coef[1]); 
      complex<double>ztmp4c2(-Quaternion_coef[2], Quaternion_coef[3]); 
      complex<double>ztmp4c3( Quaternion_coef[2], Quaternion_coef[3]); 
      complex<double>ztmp4c4( Quaternion_coef[0],-Quaternion_coef[1]); 
      AO2MO_Coef[imos3][iaos]=ztmp4c1;  // unbar alpha S
      AO2MO_Coef[imos4][iaos]=ztmp4c2;  // unbar beta S
      AO2MO_Coef[imos7][iaos]=ztmp4c3;  // bar alpha S
      AO2MO_Coef[imos8][iaos]=ztmp4c4;  // bar beta S
     } 
    }
    // Complex conjugate the coefficients because now they are stored as rows
    for(imos=0;imos<NMOs_LS;imos++)
    {
     for(iaos=0;iaos<Nbasis;iaos++)
     { 
      AO2MO_Coef[imos][iaos]=conj(AO2MO_Coef[imos][iaos]);
     }
    }
   }
  }
 }
}

// Function used to print the basis file and check how many unique primitives we used in the calc.
void print_basis_file()
{
 bool newcenter;
 int ishell,iprim,icenter;
 vector<double>Coord_x;
 vector<double>Coord_y;
 vector<double>Coord_z;
 ofstream outbasis((dirac_output_name+"basis").c_str());
 outbasis<<"# N              X_A                     Y_A                     Z_A                    Exp                      Quant(nx,ny,nz)"<<endl;
 Nprimitives=0;
 for(ishell=0;ishell<Nshell;ishell++)
 {
  if(shell2aos[ishell].paired2==-1)
  {
   switch(shell2aos[ishell].styp){
    case 0:
     Nprimitives=Nprimitives+shell2aos[ishell].nprim;
     // 0 0 0 
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      outbasis<<setw(3)<<1<<"\t"<<setprecision(15)<<fixed<<scientific<<shell2aos[ishell].Coord[0]<<"\t"<<shell2aos[ishell].Coord[1];
      outbasis<<"\t"<<shell2aos[ishell].Coord[2]<<"\t"<<shell2aos[ishell].Expon[iprim];
      outbasis<<"\t\t"<<0<<"\t"<<0<<"\t"<<0<<endl;
      shell_types.push_back(1);
      prim_exponents.push_back(shell2aos[ishell].Expon[iprim]);
      Coord_x.push_back(shell2aos[ishell].Coord[0]);Coord_y.push_back(shell2aos[ishell].Coord[1]);Coord_z.push_back(shell2aos[ishell].Coord[2]);
     }
     break;
    case 1:
     Nprimitives=Nprimitives+shell2aos[ishell].nprim*3;
     // 1 0 0 
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      outbasis<<setw(3)<<1<<"\t"<<setprecision(15)<<fixed<<scientific<<shell2aos[ishell].Coord[0]<<"\t"<<shell2aos[ishell].Coord[1];
      outbasis<<"\t"<<shell2aos[ishell].Coord[2]<<"\t"<<shell2aos[ishell].Expon[iprim];
      outbasis<<"\t\t"<<1<<"\t"<<0<<"\t"<<0<<endl;
      shell_types.push_back(2);
      prim_exponents.push_back(shell2aos[ishell].Expon[iprim]);
      Coord_x.push_back(shell2aos[ishell].Coord[0]);Coord_y.push_back(shell2aos[ishell].Coord[1]);Coord_z.push_back(shell2aos[ishell].Coord[2]);
     }
     // 0 1 0 
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      outbasis<<setw(3)<<1<<"\t"<<setprecision(15)<<fixed<<scientific<<shell2aos[ishell].Coord[0]<<"\t"<<shell2aos[ishell].Coord[1];
      outbasis<<"\t"<<shell2aos[ishell].Coord[2]<<"\t"<<shell2aos[ishell].Expon[iprim];
      outbasis<<"\t\t"<<0<<"\t"<<1<<"\t"<<0<<endl;
      shell_types.push_back(3);
      prim_exponents.push_back(shell2aos[ishell].Expon[iprim]);
      Coord_x.push_back(shell2aos[ishell].Coord[0]);Coord_y.push_back(shell2aos[ishell].Coord[1]);Coord_z.push_back(shell2aos[ishell].Coord[2]);
     }
     // 0 0 1 
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      outbasis<<setw(3)<<1<<"\t"<<setprecision(15)<<fixed<<scientific<<shell2aos[ishell].Coord[0]<<"\t"<<shell2aos[ishell].Coord[1];
      outbasis<<"\t"<<shell2aos[ishell].Coord[2]<<"\t"<<shell2aos[ishell].Expon[iprim];
      outbasis<<"\t\t"<<0<<"\t"<<0<<"\t"<<1<<endl;
      shell_types.push_back(4);
      prim_exponents.push_back(shell2aos[ishell].Expon[iprim]);
      Coord_x.push_back(shell2aos[ishell].Coord[0]);Coord_y.push_back(shell2aos[ishell].Coord[1]);Coord_z.push_back(shell2aos[ishell].Coord[2]);
     }
     break;
    case -2:
     Nprimitives=Nprimitives+shell2aos[ishell].nprim*6;
     // 2 0 0 
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      outbasis<<setw(3)<<1<<"\t"<<setprecision(15)<<fixed<<scientific<<shell2aos[ishell].Coord[0]<<"\t"<<shell2aos[ishell].Coord[1];
      outbasis<<"\t"<<shell2aos[ishell].Coord[2]<<"\t"<<shell2aos[ishell].Expon[iprim];
      outbasis<<"\t\t"<<2<<"\t"<<0<<"\t"<<0<<endl;
      shell_types.push_back(5);
      prim_exponents.push_back(shell2aos[ishell].Expon[iprim]);
      Coord_x.push_back(shell2aos[ishell].Coord[0]);Coord_y.push_back(shell2aos[ishell].Coord[1]);Coord_z.push_back(shell2aos[ishell].Coord[2]);
     }
     // 1 1 0 
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      outbasis<<setw(3)<<1<<"\t"<<setprecision(15)<<fixed<<scientific<<shell2aos[ishell].Coord[0]<<"\t"<<shell2aos[ishell].Coord[1];
      outbasis<<"\t"<<shell2aos[ishell].Coord[2]<<"\t"<<shell2aos[ishell].Expon[iprim];
      outbasis<<"\t\t"<<1<<"\t"<<1<<"\t"<<0<<endl;
      shell_types.push_back(8);
      prim_exponents.push_back(shell2aos[ishell].Expon[iprim]);
      Coord_x.push_back(shell2aos[ishell].Coord[0]);Coord_y.push_back(shell2aos[ishell].Coord[1]);Coord_z.push_back(shell2aos[ishell].Coord[2]);
     }
     // 1 0 1 
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      outbasis<<setw(3)<<1<<"\t"<<setprecision(15)<<fixed<<scientific<<shell2aos[ishell].Coord[0]<<"\t"<<shell2aos[ishell].Coord[1];
      outbasis<<"\t"<<shell2aos[ishell].Coord[2]<<"\t"<<shell2aos[ishell].Expon[iprim];
      outbasis<<"\t\t"<<1<<"\t"<<0<<"\t"<<1<<endl;
      shell_types.push_back(9);
      prim_exponents.push_back(shell2aos[ishell].Expon[iprim]);
      Coord_x.push_back(shell2aos[ishell].Coord[0]);Coord_y.push_back(shell2aos[ishell].Coord[1]);Coord_z.push_back(shell2aos[ishell].Coord[2]);
     }
     // 0 2 0 
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      outbasis<<setw(3)<<1<<"\t"<<setprecision(15)<<fixed<<scientific<<shell2aos[ishell].Coord[0]<<"\t"<<shell2aos[ishell].Coord[1];
      outbasis<<"\t"<<shell2aos[ishell].Coord[2]<<"\t"<<shell2aos[ishell].Expon[iprim];
      outbasis<<"\t\t"<<0<<"\t"<<2<<"\t"<<0<<endl;
      shell_types.push_back(6);
      prim_exponents.push_back(shell2aos[ishell].Expon[iprim]);
      Coord_x.push_back(shell2aos[ishell].Coord[0]);Coord_y.push_back(shell2aos[ishell].Coord[1]);Coord_z.push_back(shell2aos[ishell].Coord[2]);
     }
     // 0 1 1 
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      outbasis<<setw(3)<<1<<"\t"<<setprecision(15)<<fixed<<scientific<<shell2aos[ishell].Coord[0]<<"\t"<<shell2aos[ishell].Coord[1];
      outbasis<<"\t"<<shell2aos[ishell].Coord[2]<<"\t"<<shell2aos[ishell].Expon[iprim];
      outbasis<<"\t\t"<<0<<"\t"<<1<<"\t"<<1<<endl;
      shell_types.push_back(10);
      prim_exponents.push_back(shell2aos[ishell].Expon[iprim]);
      Coord_x.push_back(shell2aos[ishell].Coord[0]);Coord_y.push_back(shell2aos[ishell].Coord[1]);Coord_z.push_back(shell2aos[ishell].Coord[2]);
     }
     // 0 0 2 
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      outbasis<<setw(3)<<1<<"\t"<<setprecision(15)<<fixed<<scientific<<shell2aos[ishell].Coord[0]<<"\t"<<shell2aos[ishell].Coord[1];
      outbasis<<"\t"<<shell2aos[ishell].Coord[2]<<"\t"<<shell2aos[ishell].Expon[iprim];
      outbasis<<"\t\t"<<0<<"\t"<<0<<"\t"<<2<<endl;
      shell_types.push_back(7);
      prim_exponents.push_back(shell2aos[ishell].Expon[iprim]);
      Coord_x.push_back(shell2aos[ishell].Coord[0]);Coord_y.push_back(shell2aos[ishell].Coord[1]);Coord_z.push_back(shell2aos[ishell].Coord[2]);
     }
     break;
    case -3:
     Nprimitives=Nprimitives+shell2aos[ishell].nprim*10;
     // 3 0 0 
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      outbasis<<setw(3)<<1<<"\t"<<setprecision(15)<<fixed<<scientific<<shell2aos[ishell].Coord[0]<<"\t"<<shell2aos[ishell].Coord[1];
      outbasis<<"\t"<<shell2aos[ishell].Coord[2]<<"\t"<<shell2aos[ishell].Expon[iprim];
      outbasis<<"\t\t"<<3<<"\t"<<0<<"\t"<<0<<endl;
      shell_types.push_back(11);
      prim_exponents.push_back(shell2aos[ishell].Expon[iprim]);
      Coord_x.push_back(shell2aos[ishell].Coord[0]);Coord_y.push_back(shell2aos[ishell].Coord[1]);Coord_z.push_back(shell2aos[ishell].Coord[2]);
     }
     // 2 1 0 
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      outbasis<<setw(3)<<1<<"\t"<<setprecision(15)<<fixed<<scientific<<shell2aos[ishell].Coord[0]<<"\t"<<shell2aos[ishell].Coord[1];
      outbasis<<"\t"<<shell2aos[ishell].Coord[2]<<"\t"<<shell2aos[ishell].Expon[iprim];
      outbasis<<"\t\t"<<2<<"\t"<<1<<"\t"<<0<<endl;
      shell_types.push_back(14);
      prim_exponents.push_back(shell2aos[ishell].Expon[iprim]);
      Coord_x.push_back(shell2aos[ishell].Coord[0]);Coord_y.push_back(shell2aos[ishell].Coord[1]);Coord_z.push_back(shell2aos[ishell].Coord[2]);
     }
     // 2 0 1 
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      outbasis<<setw(3)<<1<<"\t"<<setprecision(15)<<fixed<<scientific<<shell2aos[ishell].Coord[0]<<"\t"<<shell2aos[ishell].Coord[1];
      outbasis<<"\t"<<shell2aos[ishell].Coord[2]<<"\t"<<shell2aos[ishell].Expon[iprim];
      outbasis<<"\t\t"<<2<<"\t"<<0<<"\t"<<1<<endl;
      shell_types.push_back(15);
      prim_exponents.push_back(shell2aos[ishell].Expon[iprim]);
      Coord_x.push_back(shell2aos[ishell].Coord[0]);Coord_y.push_back(shell2aos[ishell].Coord[1]);Coord_z.push_back(shell2aos[ishell].Coord[2]);
     }
     // 1 2 0 
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      outbasis<<setw(3)<<1<<"\t"<<setprecision(15)<<fixed<<scientific<<shell2aos[ishell].Coord[0]<<"\t"<<shell2aos[ishell].Coord[1];
      outbasis<<"\t"<<shell2aos[ishell].Coord[2]<<"\t"<<shell2aos[ishell].Expon[iprim];
      outbasis<<"\t\t"<<1<<"\t"<<2<<"\t"<<0<<endl;
      shell_types.push_back(17);
      prim_exponents.push_back(shell2aos[ishell].Expon[iprim]);
      Coord_x.push_back(shell2aos[ishell].Coord[0]);Coord_y.push_back(shell2aos[ishell].Coord[1]);Coord_z.push_back(shell2aos[ishell].Coord[2]);
     }
     // 1 1 1 
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      outbasis<<setw(3)<<1<<"\t"<<setprecision(15)<<fixed<<scientific<<shell2aos[ishell].Coord[0]<<"\t"<<shell2aos[ishell].Coord[1];
      outbasis<<"\t"<<shell2aos[ishell].Coord[2]<<"\t"<<shell2aos[ishell].Expon[iprim];
      outbasis<<"\t\t"<<1<<"\t"<<1<<"\t"<<1<<endl;
      shell_types.push_back(20);
      prim_exponents.push_back(shell2aos[ishell].Expon[iprim]);
      Coord_x.push_back(shell2aos[ishell].Coord[0]);Coord_y.push_back(shell2aos[ishell].Coord[1]);Coord_z.push_back(shell2aos[ishell].Coord[2]);
     }
     // 1 0 2 
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      outbasis<<setw(3)<<1<<"\t"<<setprecision(15)<<fixed<<scientific<<shell2aos[ishell].Coord[0]<<"\t"<<shell2aos[ishell].Coord[1];
      outbasis<<"\t"<<shell2aos[ishell].Coord[2]<<"\t"<<shell2aos[ishell].Expon[iprim];
      outbasis<<"\t\t"<<1<<"\t"<<0<<"\t"<<2<<endl;
      shell_types.push_back(18);
      prim_exponents.push_back(shell2aos[ishell].Expon[iprim]);
      Coord_x.push_back(shell2aos[ishell].Coord[0]);Coord_y.push_back(shell2aos[ishell].Coord[1]);Coord_z.push_back(shell2aos[ishell].Coord[2]);
     }
     // 0 3 0 
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      outbasis<<setw(3)<<1<<"\t"<<setprecision(15)<<fixed<<scientific<<shell2aos[ishell].Coord[0]<<"\t"<<shell2aos[ishell].Coord[1];
      outbasis<<"\t"<<shell2aos[ishell].Coord[2]<<"\t"<<shell2aos[ishell].Expon[iprim];
      outbasis<<"\t\t"<<0<<"\t"<<3<<"\t"<<0<<endl;
      shell_types.push_back(12);
      prim_exponents.push_back(shell2aos[ishell].Expon[iprim]);
      Coord_x.push_back(shell2aos[ishell].Coord[0]);Coord_y.push_back(shell2aos[ishell].Coord[1]);Coord_z.push_back(shell2aos[ishell].Coord[2]);
     }
     // 0 2 1 
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      outbasis<<setw(3)<<1<<"\t"<<setprecision(15)<<fixed<<scientific<<shell2aos[ishell].Coord[0]<<"\t"<<shell2aos[ishell].Coord[1];
      outbasis<<"\t"<<shell2aos[ishell].Coord[2]<<"\t"<<shell2aos[ishell].Expon[iprim];
      outbasis<<"\t\t"<<0<<"\t"<<2<<"\t"<<1<<endl;
      shell_types.push_back(16);
      prim_exponents.push_back(shell2aos[ishell].Expon[iprim]);
      Coord_x.push_back(shell2aos[ishell].Coord[0]);Coord_y.push_back(shell2aos[ishell].Coord[1]);Coord_z.push_back(shell2aos[ishell].Coord[2]);
     }
     // 0 1 2 
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      outbasis<<setw(3)<<1<<"\t"<<setprecision(15)<<fixed<<scientific<<shell2aos[ishell].Coord[0]<<"\t"<<shell2aos[ishell].Coord[1];
      outbasis<<"\t"<<shell2aos[ishell].Coord[2]<<"\t"<<shell2aos[ishell].Expon[iprim];
      outbasis<<"\t\t"<<0<<"\t"<<1<<"\t"<<2<<endl;
      shell_types.push_back(19);
      prim_exponents.push_back(shell2aos[ishell].Expon[iprim]);
      Coord_x.push_back(shell2aos[ishell].Coord[0]);Coord_y.push_back(shell2aos[ishell].Coord[1]);Coord_z.push_back(shell2aos[ishell].Coord[2]);
     }
     // 0 0 3 
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      outbasis<<setw(3)<<1<<"\t"<<setprecision(15)<<fixed<<scientific<<shell2aos[ishell].Coord[0]<<"\t"<<shell2aos[ishell].Coord[1];
      outbasis<<"\t"<<shell2aos[ishell].Coord[2]<<"\t"<<shell2aos[ishell].Expon[iprim];
      outbasis<<"\t\t"<<0<<"\t"<<0<<"\t"<<3<<endl;
      shell_types.push_back(13);
      prim_exponents.push_back(shell2aos[ishell].Expon[iprim]);
      Coord_x.push_back(shell2aos[ishell].Coord[0]);Coord_y.push_back(shell2aos[ishell].Coord[1]);Coord_z.push_back(shell2aos[ishell].Coord[2]);
     }
     break;
    case -4:
     Nprimitives=Nprimitives+shell2aos[ishell].nprim*15;
     // 4 0 0
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      outbasis<<setw(3)<<1<<"\t"<<setprecision(15)<<fixed<<scientific<<shell2aos[ishell].Coord[0]<<"\t"<<shell2aos[ishell].Coord[1];
      outbasis<<"\t"<<shell2aos[ishell].Coord[2]<<"\t"<<shell2aos[ishell].Expon[iprim];
      outbasis<<"\t\t"<<4<<"\t"<<0<<"\t"<<0<<endl;
      shell_types.push_back(21);
      prim_exponents.push_back(shell2aos[ishell].Expon[iprim]);
      Coord_x.push_back(shell2aos[ishell].Coord[0]);Coord_y.push_back(shell2aos[ishell].Coord[1]);Coord_z.push_back(shell2aos[ishell].Coord[2]);
     }
     // 3 1 0
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      outbasis<<setw(3)<<1<<"\t"<<setprecision(15)<<fixed<<scientific<<shell2aos[ishell].Coord[0]<<"\t"<<shell2aos[ishell].Coord[1];
      outbasis<<"\t"<<shell2aos[ishell].Coord[2]<<"\t"<<shell2aos[ishell].Expon[iprim];
      outbasis<<"\t\t"<<3<<"\t"<<1<<"\t"<<0<<endl;
      shell_types.push_back(24);
      prim_exponents.push_back(shell2aos[ishell].Expon[iprim]);
      Coord_x.push_back(shell2aos[ishell].Coord[0]);Coord_y.push_back(shell2aos[ishell].Coord[1]);Coord_z.push_back(shell2aos[ishell].Coord[2]);
     }
     // 3 0 1
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      outbasis<<setw(3)<<1<<"\t"<<setprecision(15)<<fixed<<scientific<<shell2aos[ishell].Coord[0]<<"\t"<<shell2aos[ishell].Coord[1];
      outbasis<<"\t"<<shell2aos[ishell].Coord[2]<<"\t"<<shell2aos[ishell].Expon[iprim];
      outbasis<<"\t\t"<<3<<"\t"<<0<<"\t"<<1<<endl;
      shell_types.push_back(25);
      prim_exponents.push_back(shell2aos[ishell].Expon[iprim]);
      Coord_x.push_back(shell2aos[ishell].Coord[0]);Coord_y.push_back(shell2aos[ishell].Coord[1]);Coord_z.push_back(shell2aos[ishell].Coord[2]);
     }
     // 2 2 0
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      outbasis<<setw(3)<<1<<"\t"<<setprecision(15)<<fixed<<scientific<<shell2aos[ishell].Coord[0]<<"\t"<<shell2aos[ishell].Coord[1];
      outbasis<<"\t"<<shell2aos[ishell].Coord[2]<<"\t"<<shell2aos[ishell].Expon[iprim];
      outbasis<<"\t\t"<<2<<"\t"<<2<<"\t"<<0<<endl;
      shell_types.push_back(30);
      prim_exponents.push_back(shell2aos[ishell].Expon[iprim]);
      Coord_x.push_back(shell2aos[ishell].Coord[0]);Coord_y.push_back(shell2aos[ishell].Coord[1]);Coord_z.push_back(shell2aos[ishell].Coord[2]);
     }
     // 2 1 1
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      outbasis<<setw(3)<<1<<"\t"<<setprecision(15)<<fixed<<scientific<<shell2aos[ishell].Coord[0]<<"\t"<<shell2aos[ishell].Coord[1];
      outbasis<<"\t"<<shell2aos[ishell].Coord[2]<<"\t"<<shell2aos[ishell].Expon[iprim];
      outbasis<<"\t\t"<<2<<"\t"<<1<<"\t"<<1<<endl;
      shell_types.push_back(33);
      prim_exponents.push_back(shell2aos[ishell].Expon[iprim]);
      Coord_x.push_back(shell2aos[ishell].Coord[0]);Coord_y.push_back(shell2aos[ishell].Coord[1]);Coord_z.push_back(shell2aos[ishell].Coord[2]);
     }
     // 2 0 2
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      outbasis<<setw(3)<<1<<"\t"<<setprecision(15)<<fixed<<scientific<<shell2aos[ishell].Coord[0]<<"\t"<<shell2aos[ishell].Coord[1];
      outbasis<<"\t"<<shell2aos[ishell].Coord[2]<<"\t"<<shell2aos[ishell].Expon[iprim];
      outbasis<<"\t\t"<<2<<"\t"<<0<<"\t"<<2<<endl;
      shell_types.push_back(31);
      prim_exponents.push_back(shell2aos[ishell].Expon[iprim]);
      Coord_x.push_back(shell2aos[ishell].Coord[0]);Coord_y.push_back(shell2aos[ishell].Coord[1]);Coord_z.push_back(shell2aos[ishell].Coord[2]);
     }
     // 1 3 0
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      outbasis<<setw(3)<<1<<"\t"<<setprecision(15)<<fixed<<scientific<<shell2aos[ishell].Coord[0]<<"\t"<<shell2aos[ishell].Coord[1];
      outbasis<<"\t"<<shell2aos[ishell].Coord[2]<<"\t"<<shell2aos[ishell].Expon[iprim];
      outbasis<<"\t\t"<<1<<"\t"<<3<<"\t"<<0<<endl;
      shell_types.push_back(26);
      prim_exponents.push_back(shell2aos[ishell].Expon[iprim]);
      Coord_x.push_back(shell2aos[ishell].Coord[0]);Coord_y.push_back(shell2aos[ishell].Coord[1]);Coord_z.push_back(shell2aos[ishell].Coord[2]);
     }
     // 1 2 1
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      outbasis<<setw(3)<<1<<"\t"<<setprecision(15)<<fixed<<scientific<<shell2aos[ishell].Coord[0]<<"\t"<<shell2aos[ishell].Coord[1];
      outbasis<<"\t"<<shell2aos[ishell].Coord[2]<<"\t"<<shell2aos[ishell].Expon[iprim];
      outbasis<<"\t\t"<<1<<"\t"<<2<<"\t"<<1<<endl;
      shell_types.push_back(34);
      prim_exponents.push_back(shell2aos[ishell].Expon[iprim]);
      Coord_x.push_back(shell2aos[ishell].Coord[0]);Coord_y.push_back(shell2aos[ishell].Coord[1]);Coord_z.push_back(shell2aos[ishell].Coord[2]);
     }
     // 1 1 2
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      outbasis<<setw(3)<<1<<"\t"<<setprecision(15)<<fixed<<scientific<<shell2aos[ishell].Coord[0]<<"\t"<<shell2aos[ishell].Coord[1];
      outbasis<<"\t"<<shell2aos[ishell].Coord[2]<<"\t"<<shell2aos[ishell].Expon[iprim];
      outbasis<<"\t\t"<<1<<"\t"<<1<<"\t"<<2<<endl;
      shell_types.push_back(35);
      prim_exponents.push_back(shell2aos[ishell].Expon[iprim]);
      Coord_x.push_back(shell2aos[ishell].Coord[0]);Coord_y.push_back(shell2aos[ishell].Coord[1]);Coord_z.push_back(shell2aos[ishell].Coord[2]);
     }
     // 1 0 3
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      outbasis<<setw(3)<<1<<"\t"<<setprecision(15)<<fixed<<scientific<<shell2aos[ishell].Coord[0]<<"\t"<<shell2aos[ishell].Coord[1];
      outbasis<<"\t"<<shell2aos[ishell].Coord[2]<<"\t"<<shell2aos[ishell].Expon[iprim];
      outbasis<<"\t\t"<<1<<"\t"<<0<<"\t"<<3<<endl;
      shell_types.push_back(28);
      prim_exponents.push_back(shell2aos[ishell].Expon[iprim]);
      Coord_x.push_back(shell2aos[ishell].Coord[0]);Coord_y.push_back(shell2aos[ishell].Coord[1]);Coord_z.push_back(shell2aos[ishell].Coord[2]);
     }
     // 0 4 0
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      outbasis<<setw(3)<<1<<"\t"<<setprecision(15)<<fixed<<scientific<<shell2aos[ishell].Coord[0]<<"\t"<<shell2aos[ishell].Coord[1];
      outbasis<<"\t"<<shell2aos[ishell].Coord[2]<<"\t"<<shell2aos[ishell].Expon[iprim];
      outbasis<<"\t\t"<<0<<"\t"<<4<<"\t"<<0<<endl;
      shell_types.push_back(22);
      prim_exponents.push_back(shell2aos[ishell].Expon[iprim]);
      Coord_x.push_back(shell2aos[ishell].Coord[0]);Coord_y.push_back(shell2aos[ishell].Coord[1]);Coord_z.push_back(shell2aos[ishell].Coord[2]);
     }
     // 0 3 1
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      outbasis<<setw(3)<<1<<"\t"<<setprecision(15)<<fixed<<scientific<<shell2aos[ishell].Coord[0]<<"\t"<<shell2aos[ishell].Coord[1];
      outbasis<<"\t"<<shell2aos[ishell].Coord[2]<<"\t"<<shell2aos[ishell].Expon[iprim];
      outbasis<<"\t\t"<<0<<"\t"<<3<<"\t"<<1<<endl;
      shell_types.push_back(27);
      prim_exponents.push_back(shell2aos[ishell].Expon[iprim]);
      Coord_x.push_back(shell2aos[ishell].Coord[0]);Coord_y.push_back(shell2aos[ishell].Coord[1]);Coord_z.push_back(shell2aos[ishell].Coord[2]);
     }
     // 0 2 2
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      outbasis<<setw(3)<<1<<"\t"<<setprecision(15)<<fixed<<scientific<<shell2aos[ishell].Coord[0]<<"\t"<<shell2aos[ishell].Coord[1];
      outbasis<<"\t"<<shell2aos[ishell].Coord[2]<<"\t"<<shell2aos[ishell].Expon[iprim];
      outbasis<<"\t\t"<<0<<"\t"<<2<<"\t"<<2<<endl;
      shell_types.push_back(32);
      prim_exponents.push_back(shell2aos[ishell].Expon[iprim]);
      Coord_x.push_back(shell2aos[ishell].Coord[0]);Coord_y.push_back(shell2aos[ishell].Coord[1]);Coord_z.push_back(shell2aos[ishell].Coord[2]);
     }
     // 0 1 3
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      outbasis<<setw(3)<<1<<"\t"<<setprecision(15)<<fixed<<scientific<<shell2aos[ishell].Coord[0]<<"\t"<<shell2aos[ishell].Coord[1];
      outbasis<<"\t"<<shell2aos[ishell].Coord[2]<<"\t"<<shell2aos[ishell].Expon[iprim];
      outbasis<<"\t\t"<<0<<"\t"<<1<<"\t"<<3<<endl;
      shell_types.push_back(29);
      prim_exponents.push_back(shell2aos[ishell].Expon[iprim]);
      Coord_x.push_back(shell2aos[ishell].Coord[0]);Coord_y.push_back(shell2aos[ishell].Coord[1]);Coord_z.push_back(shell2aos[ishell].Coord[2]);
     }
     // 0 0 4
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      outbasis<<setw(3)<<1<<"\t"<<setprecision(15)<<fixed<<scientific<<shell2aos[ishell].Coord[0]<<"\t"<<shell2aos[ishell].Coord[1];
      outbasis<<"\t"<<shell2aos[ishell].Coord[2]<<"\t"<<shell2aos[ishell].Expon[iprim];
      outbasis<<"\t\t"<<0<<"\t"<<0<<"\t"<<4<<endl;
      shell_types.push_back(23);
      prim_exponents.push_back(shell2aos[ishell].Expon[iprim]);
      Coord_x.push_back(shell2aos[ishell].Coord[0]);Coord_y.push_back(shell2aos[ishell].Coord[1]);Coord_z.push_back(shell2aos[ishell].Coord[2]);
     }
     break;
    default:
     cout<<"Warning! Shell-type not supported. "<<setw(5)<<shell2aos[ishell].styp<<endl;
     break;
   }
  }
 }
 outbasis.close();
 // Save coords. of the unique centers (Atoms) 
 Ncenters=1;
 Center_x.push_back(Coord_x[0]);Center_y.push_back(Coord_y[0]);Center_z.push_back(Coord_z[0]);
 for(iprim=0;iprim<Nprimitives;iprim++)
 {
  newcenter=true;
  for(icenter=0;icenter<Ncenters;icenter++)
  {
   if(Coord_x[iprim]==Center_x[icenter] && Coord_y[iprim]==Center_y[icenter] && Coord_z[iprim]==Center_z[icenter])
   {
    newcenter=false;
   }
  }
  if(newcenter)
  {
   Center_x.push_back(Coord_x[iprim]);Center_y.push_back(Coord_y[iprim]);Center_z.push_back(Coord_z[iprim]);
   Ncenters++;
  }
 }
 // Save Primitive to Centers map 
 for(iprim=0;iprim<Nprimitives;iprim++)
 {
  for(icenter=0;icenter<Ncenters;icenter++)
  {
   if(Coord_x[iprim]==Center_x[icenter] && Coord_y[iprim]==Center_y[icenter] && Coord_z[iprim]==Center_z[icenter])
   {
    Prim2Center_map.push_back(icenter+1);
   }
  }
 }
}

void clean_shell2aos()
{
 int ishell,iaos;
 for(ishell=0;ishell<Nshell;ishell++)
 {
  delete[] shell2aos[ishell].Coef; shell2aos[ishell].Coef=NULL;
  delete[] shell2aos[ishell].Expon; shell2aos[ishell].Expon=NULL;
  for(iaos=0;shell2aos[ishell].naos<1;iaos++)
  {
   delete[] shell2aos[ishell].index_min2max[iaos]; shell2aos[ishell].index_min2max[iaos]=NULL;
  }
  delete[] shell2aos[ishell].index_min2max;shell2aos[ishell].index_min2max=NULL;
 }
 delete[] shell2aos;shell2aos=NULL;
}

// Function used to read the 4C 2-RDM, transform it, and reduce it
void read_2rdm4cMO_and_transf()
{
 int index_4cMO[2]={10,10},index_4cMO_prime[2]={10,10};
 long int IMOS; // IMOS used with Scalar MOs (0 to NMOs_LS_4) 
 double Dijkl_4cMO,Trace=ZERO,MEM;
 NMOs_LS_1=NMOs_occ;
 NMOs_LS_2=NMOs_LS_1*NMOs_occ;
 NMOs_LS_3=NMOs_LS_2*NMOs_occ;
 NMOs_LS_4=NMOs_LS_3*NMOs_occ;
 MEM=EIGHT*NMOs_LS_4;
 sumMEM=MEM;
 MEM=MEM/pow(TEN,NINE);
 cout<<setprecision(2)<<fixed;
 cout<<"Memory required ";
 if(MEM>pow(TEN,THREE))
 {
  cout<<setw(10)<<MEM/pow(TEN,THREE)<<" Tb.";
 }
 else
 {
  if(MEM>ONE)
  {
   cout<<setw(10)<<MEM<<" Gb.";
  }
  else
  {
   if(MEM<ONE && MEM>pow(TEN,-THREE))
   {
    cout<<setw(10)<<MEM*pow(TEN,THREE)<<" Mb.";
   }
   else
   {
    cout<<setw(10)<<MEM*pow(TEN,SIX)<<" Kb.";
   }
  }
 }
 cout<<" for storing the Scalar occupied (LS) 2-RDM."<<endl;
 cout<<endl;
 Dijkl_MOsLS=new double[NMOs_LS_4];
 for(IMOS=0;IMOS<NMOs_LS_4;IMOS++){Dijkl_MOsLS[IMOS]=ZERO;} 
 ifstream input_data(dm2_file.c_str(),ios::binary);
 cout<<endl;
 cout<<"Reading the 4c 2-RDM elements"<<endl;
 while(index_4cMO[0]!=0 || index_4cMO[1]!=0 || index_4cMO_prime[0]!=0 || index_4cMO_prime[1]!=0)
 {
  input_data.seekg(RECORD_DELIMITER_LENGTH, ios::cur);
  input_data.read((char*) &index_4cMO[0], sizeof(index_4cMO[0]));
  input_data.read((char*) &index_4cMO[1], sizeof(index_4cMO[1]));
  input_data.read((char*) &index_4cMO_prime[0], sizeof(index_4cMO_prime[0]));
  input_data.read((char*) &index_4cMO_prime[1], sizeof(index_4cMO_prime[1]));
  input_data.read((char*) &Dijkl_4cMO, sizeof(Dijkl_4cMO));
  input_data.seekg(RECORD_DELIMITER_LENGTH, ios::cur);
  if(Dijkl_4cMO!=ZERO && index_4cMO[0]!=0 && index_4cMO[1]!=0 && index_4cMO_prime[0]!=0 && index_4cMO_prime[1]!=0)
  {
   index_4cMO[0]=index_4cMO[0]-1;index_4cMO[1]=index_4cMO[1]-1;index_4cMO_prime[0]=index_4cMO_prime[0]-1;index_4cMO_prime[1]=index_4cMO_prime[1]-1;
   // Note: The function dm2_4c2LS will add contributions to Dijkl_MOsLS if the 2-RDM index_Prim indices are repeated!
   dm2_4c2LS(index_4cMO[0],index_4cMO[1],index_4cMO_prime[0],index_4cMO_prime[1], Dijkl_4cMO); 
   dm2_4c2LS(index_4cMO[1],index_4cMO[0],index_4cMO_prime[0],index_4cMO_prime[1],-Dijkl_4cMO); 
   dm2_4c2LS(index_4cMO[1],index_4cMO[0],index_4cMO_prime[1],index_4cMO_prime[0], Dijkl_4cMO); 
   dm2_4c2LS(index_4cMO[0],index_4cMO[1],index_4cMO_prime[1],index_4cMO_prime[0],-Dijkl_4cMO);
   if(index_4cMO[0]!=index_4cMO_prime[0] || index_4cMO[1]!=index_4cMO_prime[1])
   {
    dm2_4c2LS(index_4cMO_prime[0],index_4cMO_prime[1],index_4cMO[0],index_4cMO[1], Dijkl_4cMO); 
    dm2_4c2LS(index_4cMO_prime[1],index_4cMO_prime[0],index_4cMO[0],index_4cMO[1],-Dijkl_4cMO); 
    dm2_4c2LS(index_4cMO_prime[1],index_4cMO_prime[0],index_4cMO[1],index_4cMO[0], Dijkl_4cMO); 
    dm2_4c2LS(index_4cMO_prime[0],index_4cMO_prime[1],index_4cMO[1],index_4cMO[0],-Dijkl_4cMO);
   }
   index_4cMO[0]=index_4cMO[0]+1;index_4cMO[1]=index_4cMO[1]+1;index_4cMO_prime[0]=index_4cMO_prime[0]+1;index_4cMO_prime[1]=index_4cMO_prime[1]+1;
   if(index_4cMO[0]==index_4cMO_prime[0] && index_4cMO[1]==index_4cMO_prime[1]){Trace=Trace+Dijkl_4cMO;}
  }
 }
 input_data.close();
 cout<<"The 2-RDM elements were read and stored in the Scalar (LS) MO basis"<<endl;
 cout<<"Trace of the 2-RDM stored: "<<setprecision(12)<<fixed<<scientific<<setw(17)<<Trace<<endl;
 cout<<endl;
 if(!transf_cplx)
 {
  // Transform from D_ij,kl (Scalar MO) to D_pq,rs (Primitives) real
  transform_Dijkl2Dpqrs();
  delete[] Dijkl_MOsLS;Dijkl_MOsLS=NULL;
  // Reduce symmetry elements and print it 
  if(!large_mem)
  {
   reduce_print();
  }
 }
 else
 {
  /*
  // Transform from D_ij,kl (Scalar MO) to D_pq,rs (Primitives) complex
  transform_Dijkl2Dpqrs_cplx();
  delete[] Dijkl_MOsLS;Dijkl_MOsLS=NULL;
  // Reduce symmetry elements, transf. to real, and print it 
  reduce_getreal_print();
  */
 }
}

// Function used to transform the 4C 2-RDM to scalar MO 2-RDM
void dm2_4c2LS(int &index_4cMO1,int &index_4cMO2,int &index_4cMO3,int &index_4cMO4, double Dijkl_4cMO)
{
 int imos,imos1;
 long int indices_MO1[4],indices_MO2[4],indices_MO3[4],indices_MO4[4];
 // Transform index from 4cMO to 2(L+S) MO
 indices_MO1[0]=index_4cMO1*4;indices_MO2[0]=index_4cMO2*4;indices_MO3[0]=index_4cMO3*4;indices_MO4[0]=index_4cMO4*4; 
 for(imos=1;imos<4;imos++)
 {
  indices_MO1[imos]=indices_MO1[imos-1]+1; 
  indices_MO2[imos]=indices_MO2[imos-1]+1; 
  indices_MO3[imos]=indices_MO3[imos-1]+1; 
  indices_MO4[imos]=indices_MO4[imos-1]+1; 
 }
 // Produce the 2-RDM in the 2(L+S) scalar MO basis
 for(imos=0;imos<4;imos++)
 {
  for(imos1=0;imos1<4;imos1++)
  {
   Dijkl_MOsLS[indices_MO1[imos]+indices_MO2[imos1]*NMOs_LS_1+indices_MO3[imos]*NMOs_LS_2+indices_MO4[imos1]*NMOs_LS_3]=
   Dijkl_MOsLS[indices_MO1[imos]+indices_MO2[imos1]*NMOs_LS_1+indices_MO3[imos]*NMOs_LS_2+indices_MO4[imos1]*NMOs_LS_3]+Dijkl_4cMO; 
  }
 }
}

// Function used to transform scalar 2-RDM to (real) Primitives
void transform_Dijkl2Dpqrs()
{
 int pivot,index_primitive[2],index_primitive_prime[2];
 long int IMOS,IMOS1,IMOS2,IMOS3;        // IMOS used with Scalar MOs (0 to NMOs_LS_4) 
 long int IPRIM,IPRIM1,IPRIM2,IPRIM3;    // IPRIM (0 to Nprimitives) but can be summed for large Nprimitives number.
 double Dpqrs_re,MEM;
 complex<double>Dpqrs;
 complex<double> *Diqrs_Prims,*Dijrs_Prims,*Dijks_Prims;
 Largest_Prim=-1;
 Nprims1=Nprimitives; 
 Nprims2=Nprims1*Nprimitives; 
 Nprims3=Nprims2*Nprimitives;
 Nprims4=Nprims3*Nprimitives;
 MEM=EIGHT*(TWO*Nprims1+FOUR*Nprims2+EIGHT*Nprims3);
 if(large_mem)
 {
  MEM=(EIGHT*(TWO*Nprims1+FOUR*Nprims2+EIGHT*Nprims3+(Nprims1*(Nprims1+1)*Nprims1*(Nprims1+1))/4));
 }
 MEM=MEM+sumMEM;
 if(MEM>maxmem)
 {
  cout<<endl;
  cout<<"Increase the amount of RAM memory or remove the $large_mem keyword"<<endl;
  MEM=MEM/pow(TEN,NINE);
  cout<<setprecision(2)<<fixed;
  cout<<"Memory required ";
  if(MEM>pow(TEN,THREE))
  {
   cout<<setw(10)<<MEM/pow(TEN,THREE)<<" Tb.";
  }
  else
  {
   if(MEM>ONE)
   {
    cout<<setw(10)<<MEM<<" Gb.";
   }
   else
   {
    if(MEM<ONE && MEM>pow(TEN,-THREE))
    {
     cout<<setw(10)<<MEM*pow(TEN,THREE)<<" Mb.";
    }
    else
    {
     cout<<setw(10)<<MEM*pow(TEN,SIX)<<" Kb.";
    }
   }
  }
  cout<<" for the index transformation and storage of the Scalar occupied (LS) 2-RDM."<<endl; 
  cout<<endl;
 }
 else
 {
  MEM=MEM-sumMEM;
  MEM=MEM/pow(TEN,NINE);
  cout<<setprecision(2)<<fixed;
  cout<<"Memory required ";
  if(MEM>pow(TEN,THREE))
  {
   cout<<setw(10)<<MEM/pow(TEN,THREE)<<" Tb.";
  }
  else
  {
   if(MEM>ONE)
   {
    cout<<setw(10)<<MEM<<" Gb.";
   }
   else
   {
    if(MEM<ONE && MEM>pow(TEN,-THREE))
    {
     cout<<setw(10)<<MEM*pow(TEN,THREE)<<" Mb.";
    }
    else
    {
     cout<<setw(10)<<MEM*pow(TEN,SIX)<<" Kb.";
    }
   }
  }
  cout<<" for the index transformation."<<endl; 
  cout<<endl;
  cout<<"Start writing the transformed-real 2-RDM elements in the primitives basis"<<endl;
  cout<<"[Note: The complex part will be ignored]"<<endl;
  string name="Conv_"+dirac_output_name+"dm2";
  ofstream output_data(name.c_str(),ios::binary);
  Dijks_Prims=new complex<double>[Nprims1];
  Dijrs_Prims=new complex<double>[Nprims2];
  Diqrs_Prims=new complex<double>[Nprims3];
  if(large_mem)
  {
   cout<<"Num. of 2-RDM elem.  : "<<setw(12)<<(Nprims1*(Nprims1+1)*Nprims1*(Nprims1+1))/4<<endl;
   Dpqrs_ALL=new double***[Nprims1];
   for(IPRIM=0;IPRIM<Nprims1;IPRIM++)
   {
    Dpqrs_ALL[IPRIM]=new double**[Nprims1];
    for(IPRIM1=0;IPRIM1<Nprims1;IPRIM1++)
    {
     Dpqrs_ALL[IPRIM][IPRIM1]=new double*[IPRIM+1];
     for(IPRIM2=0;IPRIM2<IPRIM+1;IPRIM2++)
     {
      Dpqrs_ALL[IPRIM][IPRIM1][IPRIM2]=new double[IPRIM1+1];
      for(IPRIM3=0;IPRIM3<IPRIM1+1;IPRIM3++)
      {
       Dpqrs_ALL[IPRIM][IPRIM1][IPRIM2][IPRIM3]=ZERO;
      }
     }
    }
   }
  }
  // Reduce
  cout<<endl;
  cout<<"Automatic reduction of the p <-> r and q <-> s terms"<<endl;
  cout<<endl;
  for(IMOS=0;IMOS<NMOs_LS_1;IMOS++)         // i
  {
   cout<<"Transforming scalar MO "<<setw(5)<<IMOS+1<<"/"<<setw(5)<<NMOs_LS_1<<endl;
   for(IPRIM1=0;IPRIM1<Nprims3;IPRIM1++)    // Initialize qrs
   {
    Diqrs_Prims[IPRIM1]=CZERO;
   }
   // Change j for fixed i
   for(IMOS1=0;IMOS1<NMOs_LS_1;IMOS1++)     // j
   { 
    for(IPRIM2=0;IPRIM2<Nprims2;IPRIM2++)   // Initialize rs
    {
     Dijrs_Prims[IPRIM2]=CZERO;
    }
    // Change k for fixed ij
    for(IMOS2=0;IMOS2<NMOs_LS_1;IMOS2++)     // k
    {
     for(IPRIM3=0;IPRIM3<Nprims1;IPRIM3++)   // Initialize s
     {
      Dijks_Prims[IPRIM3]=CZERO;
     }
     // Change l for fixed ijk
     for(IMOS3=0;IMOS3<NMOs_LS_1;IMOS3++)    // l
     {
      // Change l -> s for fixed ijk
      for(IPRIM3=0;IPRIM3<Nprims1;IPRIM3++)  // s
      {
       Dijks_Prims[IPRIM3]=Dijks_Prims[IPRIM3]+Dijkl_MOsLS[IMOS+IMOS1*NMOs_LS_1+IMOS2*NMOs_LS_2+IMOS3*NMOs_LS_3]*Prim2MO_Coef[IMOS3][IPRIM3]; 
      }
     }
     // Change k -> r for fixed ij
     for(IPRIM3=0;IPRIM3<Nprims1;IPRIM3++)   // s
     {
      for(IPRIM2=0;IPRIM2<Nprims1;IPRIM2++)  // r
      {
       Dijrs_Prims[IPRIM2+IPRIM3*Nprims1]=Dijrs_Prims[IPRIM2+IPRIM3*Nprims1]+Dijks_Prims[IPRIM3]*Prim2MO_Coef[IMOS2][IPRIM2];
      }
     }
    }
    // Change j -> q for fixed i
    for(IPRIM3=0;IPRIM3<Nprims1;IPRIM3++)    // s
    {
     for(IPRIM2=0;IPRIM2<Nprims1;IPRIM2++)   // r
     {
      for(IPRIM1=0;IPRIM1<Nprims1;IPRIM1++)  // q
      {
       Diqrs_Prims[IPRIM1+IPRIM2*Nprims1+IPRIM3*Nprims2]=Diqrs_Prims[IPRIM1+IPRIM2*Nprims1+IPRIM3*Nprims2]
       +Dijrs_Prims[IPRIM2+IPRIM3*Nprims1]*conj(Prim2MO_Coef[IMOS1][IPRIM1]);
      }
     }
    }
   }
   // TODO: Reduce symetry and permutation while we write?
   // Change i -> p
   for(IPRIM3=0;IPRIM3<Nprims1;IPRIM3++)    // s
   {
    for(IPRIM2=0;IPRIM2<Nprims1;IPRIM2++)   // r
    {
     for(IPRIM1=0;IPRIM1<Nprims1;IPRIM1++)  // q
     {
      for(IPRIM=0;IPRIM<Nprims1;IPRIM++)    // p
      {
       Dpqrs=Diqrs_Prims[IPRIM1+IPRIM2*Nprims1+IPRIM3*Nprims2]*conj(Prim2MO_Coef[IMOS][IPRIM]);
       Dpqrs_re=Dpqrs.real();
       if(abs(Dpqrs_re)>threshold)
       {
        index_primitive[0]=IPRIM;index_primitive[1]=IPRIM1;index_primitive_prime[0]=IPRIM2;index_primitive_prime[1]=IPRIM3;
        if(index_primitive[0]>Largest_Prim){Largest_Prim=index_primitive[0];}
        if(index_primitive[1]>Largest_Prim){Largest_Prim=index_primitive[1];}
        if(index_primitive_prime[0]>Largest_Prim){Largest_Prim=index_primitive_prime[0];}
        if(index_primitive_prime[0]>Largest_Prim){Largest_Prim=index_primitive_prime[1];}
        if(index_primitive[0]<index_primitive_prime[0])
        {
         pivot=index_primitive_prime[0];index_primitive_prime[0]=index_primitive[0];index_primitive[0]=pivot;
        }
        if(index_primitive[1]<index_primitive_prime[1])
        {
         pivot=index_primitive_prime[1];index_primitive_prime[1]=index_primitive[1];index_primitive[1]=pivot;
        }
        if(!large_mem)
        {
         Nterms_printed++;
         index_primitive[0]++;index_primitive[1]++;index_primitive_prime[0]++;index_primitive_prime[1]++;
         output_data.seekp(RECORD_DELIMITER_LENGTH, ios::cur);
         output_data.write((char*) &index_primitive[0], sizeof(index_primitive[0]));
         output_data.write((char*) &index_primitive[1], sizeof(index_primitive[1]));
         output_data.write((char*) &index_primitive_prime[0], sizeof(index_primitive_prime[0]));
         output_data.write((char*) &index_primitive_prime[1], sizeof(index_primitive_prime[1]));
         output_data.write((char*) &Dpqrs_re, sizeof(Dpqrs_re));
         output_data.seekp(RECORD_DELIMITER_LENGTH, ios::cur);
        }
        else
        {
         Dpqrs_ALL[index_primitive[0]][index_primitive[1]][index_primitive_prime[0]][index_primitive_prime[1]]=
         Dpqrs_ALL[index_primitive[0]][index_primitive[1]][index_primitive_prime[0]][index_primitive_prime[1]]+Dpqrs_re;
        }
       }
      }
     }
    }
   }
  }
  if(large_mem)
  {
   // Print  it
   for(IPRIM=0;IPRIM<Nprims1;IPRIM++)
   {
    for(IPRIM1=0;IPRIM1<Nprims1;IPRIM1++)
    {
     for(IPRIM2=0;IPRIM2<IPRIM+1;IPRIM2++)
     {
      for(IPRIM3=0;IPRIM3<IPRIM1+1;IPRIM3++)
      {
       Dpqrs_re=Dpqrs_ALL[IPRIM][IPRIM1][IPRIM2][IPRIM3];
       if(abs(Dpqrs_re)>=threshold)
       {
        index_primitive[0]=IPRIM+1;index_primitive[1]=IPRIM1+1;index_primitive_prime[0]=IPRIM2+1;index_primitive_prime[1]=IPRIM3+1;
        output_data.seekp(RECORD_DELIMITER_LENGTH, ios::cur);
        output_data.write((char*) &index_primitive[0], sizeof(index_primitive[0]));
        output_data.write((char*) &index_primitive[1], sizeof(index_primitive[1]));
        output_data.write((char*) &index_primitive_prime[0], sizeof(index_primitive_prime[0]));
        output_data.write((char*) &index_primitive_prime[1], sizeof(index_primitive_prime[1]));
        output_data.write((char*) &Dpqrs_re, sizeof(Dpqrs_re));
        output_data.seekp(RECORD_DELIMITER_LENGTH, ios::cur);
       }
      }
     }
    }
   }
   // Deallocate
   for(IPRIM=0;IPRIM<Nprims1;IPRIM++)
   {
    for(IPRIM1=0;IPRIM1<Nprims1;IPRIM1++)
    {
     for(IPRIM2=0;IPRIM2<IPRIM+1;IPRIM2++)
     {
      delete[] Dpqrs_ALL[IPRIM][IPRIM1][IPRIM2];Dpqrs_ALL[IPRIM][IPRIM1][IPRIM2]=NULL;
     }
     delete[] Dpqrs_ALL[IPRIM][IPRIM1];Dpqrs_ALL[IPRIM][IPRIM1]=NULL;
    }
    delete[] Dpqrs_ALL[IPRIM];Dpqrs_ALL[IPRIM]=NULL;
   }
   delete[] Dpqrs_ALL;Dpqrs_ALL=NULL;
  }
  Dpqrs_re=ZERO;
  index_primitive[0]=0;index_primitive[1]=0;index_primitive_prime[0]=0;index_primitive_prime[1]=0;
  output_data.seekp(RECORD_DELIMITER_LENGTH, ios::cur);
  output_data.write((char*) &index_primitive[0], sizeof(index_primitive[0]));
  output_data.write((char*) &index_primitive[1], sizeof(index_primitive[1]));
  output_data.write((char*) &index_primitive_prime[0], sizeof(index_primitive_prime[0]));
  output_data.write((char*) &index_primitive_prime[1], sizeof(index_primitive_prime[1]));
  output_data.write((char*) &Dpqrs_re, sizeof(Dpqrs_re));
  output_data.seekp(RECORD_DELIMITER_LENGTH, ios::cur);
  output_data.close();
  cout<<endl;
  cout<<"Finished writing the transformed-real 2-RDM elements in the primitives basis"<<endl;
  cout<<"to the file "<<name<<endl;
  if(!large_mem)
  {
   Largest_Prim++;
   cout<<"Num. of printed terms: "<<setw(12)<<Nterms_printed<<endl;
   cout<<"Largest prim. in use : "<<setw(12)<<Largest_Prim<<endl;
  }
  if(Largest_Prim<Nprimitives && !large_mem)
  {
   cout<<"Overwriting Nprimitives variable with largest prim. in use"<<endl;
   Nprimitives=Largest_Prim;
   Nprims1=Nprimitives; 
   Nprims2=Nprims1*Nprimitives; 
   Nprims3=Nprims2*Nprimitives;
   Nprims4=Nprims3*Nprimitives;
  }
  cout<<endl;
  delete[] Dijks_Prims;Dijks_Prims=NULL;
  delete[] Dijrs_Prims;Dijrs_Prims=NULL;
  delete[] Diqrs_Prims;Diqrs_Prims=NULL;
 }
}  

// Read the 2-RDM in primitive basis, reduce it using symmetry, and print it real
void reduce_print()
{
 int index_Prim[2],index_Prim_prime[2];
 long int IPRIM,IPRIM1,IPRIM2,IPRIM3,IPRIM4;    // IPRIM (0 to Nprimitives) but can be summed for large Nprimitives number.
 double Dpqrs,MEM;
 string conv_name="Conv_"+dirac_output_name+"dm2";
 MEM=(EIGHT*((Nprims1*(Nprims1+1)*Nprims1*(Nprims1+1))/4))/pow(TEN,NINE);
 cout<<setprecision(2)<<fixed;
 cout<<"Memory required ";
 if(MEM>pow(TEN,THREE))
 {
  cout<<setw(10)<<MEM/pow(TEN,THREE)<<" Tb.";
 }
 else
 {
  if(MEM>ONE)
  {
   cout<<setw(10)<<MEM<<" Gb.";
  }
  else
  {
   if(MEM<ONE && MEM>pow(TEN,-THREE))
   {
    cout<<setw(10)<<MEM*pow(TEN,THREE)<<" Mb.";
   }
   else
   {
    cout<<setw(10)<<MEM*pow(TEN,SIX)<<" Kb.";
   }
  }
 }
 cout<<" for the reduction."<<endl; 
 cout<<"Num. of 2-RDM elem.  : "<<setw(12)<<(Nprims1*(Nprims1+1)*Nprims1*(Nprims1+1))/4<<endl;
 cout<<endl;
 Dpqrs_ALL=new double***[Nprims1];
 for(IPRIM=0;IPRIM<Nprims1;IPRIM++)
 {
  Dpqrs_ALL[IPRIM]=new double**[Nprims1];
  for(IPRIM1=0;IPRIM1<Nprims1;IPRIM1++)
  {
   Dpqrs_ALL[IPRIM][IPRIM1]=new double*[IPRIM+1];
   for(IPRIM2=0;IPRIM2<IPRIM+1;IPRIM2++)
   {
    Dpqrs_ALL[IPRIM][IPRIM1][IPRIM2]=new double[IPRIM1+1];
    for(IPRIM3=0;IPRIM3<IPRIM1+1;IPRIM3++)
    {
     Dpqrs_ALL[IPRIM][IPRIM1][IPRIM2][IPRIM3]=ZERO;
    }
   }
  }
 }
 ifstream input_data2(conv_name.c_str(),ios::binary);
 cout<<endl;
 index_Prim[0]=10;index_Prim[1]=10;index_Prim_prime[0]=10;index_Prim_prime[1]=10;
 cout<<"Reading the transformed 2-RDM elements from "<<conv_name<<endl;
 while(index_Prim[0]!=0 || index_Prim[1]!=0 || index_Prim_prime[0]!=0 || index_Prim_prime[1]!=0)
 {
  input_data2.seekg(RECORD_DELIMITER_LENGTH, ios::cur);
  input_data2.read((char*) &index_Prim[0], sizeof(index_Prim[0]));
  input_data2.read((char*) &index_Prim[1], sizeof(index_Prim[1]));
  input_data2.read((char*) &index_Prim_prime[0], sizeof(index_Prim_prime[0]));
  input_data2.read((char*) &index_Prim_prime[1], sizeof(index_Prim_prime[1]));
  input_data2.read((char*) &Dpqrs, sizeof(Dpqrs));
  input_data2.seekg(RECORD_DELIMITER_LENGTH, ios::cur);
  if(abs(Dpqrs)!=ZERO && index_Prim[0]!=0 && index_Prim[1]!=0 && index_Prim_prime[0]!=0 && index_Prim_prime[1]!=0)
  {
   index_Prim[0]=index_Prim[0]-1;index_Prim[1]=index_Prim[1]-1;index_Prim_prime[0]=index_Prim_prime[0]-1;index_Prim_prime[1]=index_Prim_prime[1]-1;
   Dpqrs_ALL[index_Prim[0]][index_Prim[1]][index_Prim_prime[0]][index_Prim_prime[1]]=
   Dpqrs_ALL[index_Prim[0]][index_Prim[1]][index_Prim_prime[0]][index_Prim_prime[1]]+Dpqrs;
   index_Prim[0]=index_Prim[0]+1;index_Prim[1]=index_Prim[1]+1;index_Prim_prime[0]=index_Prim_prime[0]+1;index_Prim_prime[1]=index_Prim_prime[1]+1;
  }
 }
 input_data2.close();
 // Print
 cout<<"Printing the reduced transformed 2-RDM"<<endl;
 cout<<"overwriting the file "<<conv_name<<endl;
 cout<<endl;
 ofstream output_data(conv_name.c_str(),ios::binary);
 for(IPRIM=0;IPRIM<Nprims1;IPRIM++)
 {
  for(IPRIM1=0;IPRIM1<Nprims1;IPRIM1++)
  {
   for(IPRIM2=0;IPRIM2<IPRIM+1;IPRIM2++)
   {
    for(IPRIM3=0;IPRIM3<IPRIM1+1;IPRIM3++)
    {
     Dpqrs=Dpqrs_ALL[IPRIM][IPRIM1][IPRIM2][IPRIM3];
     if(abs(Dpqrs)>=threshold)
     {
      index_Prim[0]=IPRIM+1;index_Prim[1]=IPRIM1+1;index_Prim_prime[0]=IPRIM2+1;index_Prim_prime[1]=IPRIM3+1;
      output_data.seekp(RECORD_DELIMITER_LENGTH, ios::cur);
      output_data.write((char*) &index_Prim[0], sizeof(index_Prim[0]));
      output_data.write((char*) &index_Prim[1], sizeof(index_Prim[1]));
      output_data.write((char*) &index_Prim_prime[0], sizeof(index_Prim_prime[0]));
      output_data.write((char*) &index_Prim_prime[1], sizeof(index_Prim_prime[1]));
      output_data.write((char*) &Dpqrs, sizeof(Dpqrs));
      output_data.seekp(RECORD_DELIMITER_LENGTH, ios::cur);
     }
    }
   }
  }
 }
 index_Prim[0]=0;
 index_Prim[1]=0;
 index_Prim_prime[0]=0;
 index_Prim_prime[1]=0;
 Dpqrs=ZERO;
 output_data.seekp(RECORD_DELIMITER_LENGTH, ios::cur);
 output_data.write((char*) &index_Prim[0], sizeof(index_Prim[0]));
 output_data.write((char*) &index_Prim[1], sizeof(index_Prim[1]));
 output_data.write((char*) &index_Prim_prime[0], sizeof(index_Prim_prime[0]));
 output_data.write((char*) &index_Prim_prime[1], sizeof(index_Prim_prime[1]));
 output_data.write((char*) &Dpqrs, sizeof(Dpqrs));
 output_data.seekp(RECORD_DELIMITER_LENGTH, ios::cur);
 output_data.close();
 for(IPRIM=0;IPRIM<Nprims1;IPRIM++)
 {
  for(IPRIM1=0;IPRIM1<Nprims1;IPRIM1++)
  {
   for(IPRIM2=0;IPRIM2<IPRIM+1;IPRIM2++)
   {
    delete[] Dpqrs_ALL[IPRIM][IPRIM1][IPRIM2];Dpqrs_ALL[IPRIM][IPRIM1][IPRIM2]=NULL;
   }
   delete[] Dpqrs_ALL[IPRIM][IPRIM1];Dpqrs_ALL[IPRIM][IPRIM1]=NULL;
  }
  delete[] Dpqrs_ALL[IPRIM];Dpqrs_ALL[IPRIM]=NULL;
 }
 delete[] Dpqrs_ALL;Dpqrs_ALL=NULL;
}

// Currently, only available for ATOMS
void print_wfx()
{
 int iprim,imos,imos1,icenter;
 string line;
 ofstream real_wfx(("dirac_"+dirac_output_name+"RE.wfx").c_str());
 ofstream imag_wfx(("dirac_"+dirac_output_name+"IM.wfx").c_str());
 real_wfx<<setprecision(15)<<fixed<<scientific;
 imag_wfx<<setprecision(15)<<fixed<<scientific;
 line="<Number of Nuclei>";
 real_wfx<<line<<endl; 
 imag_wfx<<line<<endl; 
 real_wfx<<Ncenters<<endl; 
 imag_wfx<<Ncenters<<endl; 
 line="</Number of Nuclei>";
 real_wfx<<line<<endl; 
 imag_wfx<<line<<endl; 
 line="<Number of Occupied Molecular Orbitals>";
 real_wfx<<line<<endl; 
 imag_wfx<<line<<endl; 
 if(OneMO_wfx==-1)
 {
  real_wfx<<NMOs_occ<<endl; 
  imag_wfx<<NMOs_occ<<endl;
 } 
 else
 {
  real_wfx<<4<<endl; 
  imag_wfx<<4<<endl;
 } 
 line="</Number of Occupied Molecular Orbitals>";
 real_wfx<<line<<endl; 
 imag_wfx<<line<<endl; 
 line="<Electronic Spin Multiplicity>"; // Fixed for 'faked' close shell
 real_wfx<<line<<endl;
 imag_wfx<<line<<endl;
 real_wfx<<1<<endl; 
 imag_wfx<<1<<endl; 
 line="</Electronic Spin Multiplicity>";
 real_wfx<<line<<endl;
 imag_wfx<<line<<endl;
 line="<Nuclear Charges>"; // Faked (No need to change them for RHO_OPS!)
 real_wfx<<line<<endl;
 imag_wfx<<line<<endl;
 for(icenter=0;icenter<Ncenters;icenter++)
 {
  real_wfx<<0<<endl; 
  imag_wfx<<0<<endl; 
 }
 line="</Nuclear Charges>"; // Faked
 real_wfx<<line<<endl;
 imag_wfx<<line<<endl;
 line="<Nuclear Cartesian Coordinates>"; // Only C1 symmetry in DIRAC 
 real_wfx<<line<<endl;
 imag_wfx<<line<<endl;
 for(icenter=0;icenter<Ncenters;icenter++)
 {
  real_wfx<<Center_x[icenter]<<endl; // x
  imag_wfx<<Center_x[icenter]<<endl; 
  real_wfx<<Center_y[icenter]<<endl; // y
  imag_wfx<<Center_y[icenter]<<endl; 
  real_wfx<<Center_z[icenter]<<endl; // z
  imag_wfx<<Center_z[icenter]<<endl; 
 }
 line="</Nuclear Cartesian Coordinates>"; // Only C1 symmetry in DIRAC 
 real_wfx<<line<<endl;
 imag_wfx<<line<<endl;
 line="<Number of Primitives>";
 real_wfx<<line<<endl;
 imag_wfx<<line<<endl;
 real_wfx<<Nprimitives<<endl;
 imag_wfx<<Nprimitives<<endl;
 line="</Number of Primitives>";
 real_wfx<<line<<endl;
 imag_wfx<<line<<endl;
 line="<Primitive Centers>"; // Only atomic systems 
 real_wfx<<line<<endl;
 imag_wfx<<line<<endl;
 for(iprim=0;iprim<Nprimitives;iprim++)
 {
  real_wfx<<Prim2Center_map[iprim]<<endl;
  imag_wfx<<Prim2Center_map[iprim]<<endl;
 }
 line="</Primitive Centers>"; 
 real_wfx<<line<<endl;
 imag_wfx<<line<<endl;
 line="<Primitive Types>";
 real_wfx<<line<<endl;
 imag_wfx<<line<<endl;
 for(iprim=0;iprim<Nprimitives;iprim++)
 {
  real_wfx<<shell_types[iprim]<<endl;
  imag_wfx<<shell_types[iprim]<<endl;
 }
 line="</Primitive Types>";
 real_wfx<<line<<endl;
 imag_wfx<<line<<endl;
 line="<Primitive Exponents>";
 real_wfx<<line<<endl;
 imag_wfx<<line<<endl;
 for(iprim=0;iprim<Nprimitives;iprim++)
 {
  real_wfx<<prim_exponents[iprim]<<endl;
  imag_wfx<<prim_exponents[iprim]<<endl;
 }
 line="</Primitive Exponents>";
 real_wfx<<line<<endl;
 imag_wfx<<line<<endl;
 line="<Molecular Orbital Occupation Numbers>";
 real_wfx<<line<<endl;
 imag_wfx<<line<<endl;
 if(OneMO_wfx==-1)
 {
  for(imos=0;imos<NMOs_occ;imos++)
  {
   real_wfx<<MOsLS_occ[imos]<<endl;
   imag_wfx<<MOsLS_occ[imos]<<endl;
  }
 }
 else
 {
  for(imos=0;imos<4;imos++)
  {
   real_wfx<<ONE<<endl;
   imag_wfx<<ONE<<endl;
  }
 }
 line="</Molecular Orbital Occupation Numbers>";
 real_wfx<<line<<endl;
 imag_wfx<<line<<endl;
 line="<Molecular Orbital Spin Types>"; // Faked spin
 real_wfx<<line<<endl;
 imag_wfx<<line<<endl;
 if(OneMO_wfx==-1)
 {
  for(imos=0;imos<NMOs_occ;imos++)
  {
   line=" Alpha";
   real_wfx<<line<<endl;
   imag_wfx<<line<<endl;
  }
 }
 else
 {
  for(imos=0;imos<4;imos++)
  {
   line=" Alpha";
   real_wfx<<line<<endl;
   imag_wfx<<line<<endl;
  }
 }
 line="</Molecular Orbital Spin Types>";
 real_wfx<<line<<endl;
 imag_wfx<<line<<endl;
 line="<Molecular Orbital Primitive Coefficients>"; 
 real_wfx<<line<<endl;
 imag_wfx<<line<<endl;
 if(OneMO_wfx==-1)
 {
  for(imos=0;imos<NMOs_occ;imos++)
  {
   line="</MO Number>";
   real_wfx<<line<<endl;
   imag_wfx<<line<<endl;
   for(iprim=0;iprim<Nprimitives;iprim++)
   {
    real_wfx<<Prim2MO_Coef[imos][iprim].real()<<endl;
    imag_wfx<<Prim2MO_Coef[imos][iprim].imag()<<endl;
   }
  }
 }
 else
 {
  imos1=0;
  for(imos=0;imos<4;imos++)
  {
   line="</MO Number>";
   real_wfx<<line<<endl;
   imag_wfx<<line<<endl;
   for(iprim=0;iprim<Nprimitives;iprim++)
   {
    real_wfx<<One_Prim2MO_Coef_RE[imos1]<<endl;
    imag_wfx<<One_Prim2MO_Coef_IM[imos1]<<endl;
    imos1++;
   }
  }
 }
 line="</Molecular Orbital Primitive Coefficients>"; 
 real_wfx<<line<<endl;
 imag_wfx<<line<<endl;
 imag_wfx.close(); 
 real_wfx.close(); 
}

//
/* MRM: Switched off because we do not need the imaginary part of the 2-RDM in primitives
// 
// Function used to transform scalar 2-RDM to (complex) Primitives
void transform_Dijkl2Dpqrs_cplx()
{
 int index_primitive[2],index_primitive_prime[2];
 long int IMOS,IMOS1,IMOS2,IMOS3;        // IMOS used with Scalar MOs (0 to NMOs_LS_4) 
 long int IPRIM,IPRIM1,IPRIM2,IPRIM3;    // IPRIM (0 to Nprimitives) but can be summed for large Nprimitives number.
 double Dpqrs_re,Dpqrs_im,MEM;
 complex<double>Dpqrs;
 complex<double> *Diqrs_Prims,*Dijrs_Prims,*Dijks_Prims;
 Largest_Prim=-1;
 Nprims1=Nprimitives; 
 Nprims2=Nprims1*Nprimitives; 
 Nprims3=Nprims2*Nprimitives;
 Nprims4=Nprims3*Nprimitives;
 MEM=EIGHT*(TWO*Nprims1+FOUR*Nprims2+EIGHT*Nprims3)/pow(TEN,NINE);
 cout<<setprecision(2)<<fixed;
 cout<<"Memory required ";
 if(MEM>pow(TEN,THREE))
 {
  cout<<setw(10)<<MEM/pow(TEN,THREE)<<" Tb.";
 }
 else
 {
  if(MEM>ONE)
  {
   cout<<setw(10)<<MEM<<" Gb.";
  }
  else
  {
   if(MEM<ONE && MEM>pow(TEN,-THREE))
   {
    cout<<setw(10)<<MEM*pow(TEN,THREE)<<" Mb.";
   }
   else
   {
    cout<<setw(10)<<MEM*pow(TEN,SIX)<<" Kb.";
   }
  }
 }
 cout<<" for the index transformation."<<endl; 
 cout<<endl;
 cout<<"Start writing the transformed-complex 2-RDM elements in the primitives basis"<<endl;
 string name="Conv_cplx_"+dirac_output_name+"dm2";
 ofstream output_data(name.c_str(),ios::binary);
 Dijks_Prims=new complex<double>[Nprims1];
 Dijrs_Prims=new complex<double>[Nprims2];
 Diqrs_Prims=new complex<double>[Nprims3];
 cout<<endl;
 for(IMOS=0;IMOS<NMOs_LS_1;IMOS++)         // i
 {
  cout<<"Transforming scalar MO "<<setw(5)<<IMOS+1<<"/"<<setw(5)<<NMOs_LS_1<<endl;
  for(IPRIM1=0;IPRIM1<Nprims3;IPRIM1++)    // Initialize qrs
  {
   Diqrs_Prims[IPRIM1]=CZERO;
  }
  // Change j for fixed i
  for(IMOS1=0;IMOS1<NMOs_LS_1;IMOS1++)     // j
  { 
   for(IPRIM2=0;IPRIM2<Nprims2;IPRIM2++)   // Initialize rs
   {
    Dijrs_Prims[IPRIM2]=CZERO;
   }
   // Change k for fixed ij
   for(IMOS2=0;IMOS2<NMOs_LS_1;IMOS2++)     // k
   {
    for(IPRIM3=0;IPRIM3<Nprims1;IPRIM3++)   // Initialize s
    {
     Dijks_Prims[IPRIM3]=CZERO;
    }
    // Change l for fixed ijk
    for(IMOS3=0;IMOS3<NMOs_LS_1;IMOS3++)    // l
    {
     // Change l -> s for fixed ijk
     for(IPRIM3=0;IPRIM3<Nprims1;IPRIM3++)  // s
     {
      Dijks_Prims[IPRIM3]=Dijks_Prims[IPRIM3]+Dijkl_MOsLS[IMOS+IMOS1*NMOs_LS_1+IMOS2*NMOs_LS_2+IMOS3*NMOs_LS_3]*Prim2MO_Coef[IMOS3][IPRIM3]; 
     }
    }
    // Change k -> r for fixed ij
    for(IPRIM3=0;IPRIM3<Nprims1;IPRIM3++)   // s
    {
     for(IPRIM2=0;IPRIM2<Nprims1;IPRIM2++)  // r
     {
      Dijrs_Prims[IPRIM2+IPRIM3*Nprims1]=Dijrs_Prims[IPRIM2+IPRIM3*Nprims1]+Dijks_Prims[IPRIM3]*Prim2MO_Coef[IMOS2][IPRIM2];
     }
    }
   }
   // Change j -> q for fixed i
   for(IPRIM3=0;IPRIM3<Nprims1;IPRIM3++)    // s
   {
    for(IPRIM2=0;IPRIM2<Nprims1;IPRIM2++)   // r
    {
     for(IPRIM1=0;IPRIM1<Nprims1;IPRIM1++)  // q
     {
      Diqrs_Prims[IPRIM1+IPRIM2*Nprims1+IPRIM3*Nprims2]=Diqrs_Prims[IPRIM1+IPRIM2*Nprims1+IPRIM3*Nprims2]
      +Dijrs_Prims[IPRIM2+IPRIM3*Nprims1]*conj(Prim2MO_Coef[IMOS1][IPRIM1]);
     }
    }
   }
  }
  // Change i -> p
  for(IPRIM3=0;IPRIM3<Nprims1;IPRIM3++)    // s
  {
   for(IPRIM2=0;IPRIM2<Nprims1;IPRIM2++)   // r
   {
    for(IPRIM1=0;IPRIM1<Nprims1;IPRIM1++)  // q
    {
     for(IPRIM=0;IPRIM<Nprims1;IPRIM++)    // p
     {
      Dpqrs=Diqrs_Prims[IPRIM1+IPRIM2*Nprims1+IPRIM3*Nprims2]*conj(Prim2MO_Coef[IMOS][IPRIM]);
      if(abs(Dpqrs)>threshold)
      {
       Nterms_printed++;
       Dpqrs_re=Dpqrs.real();
       Dpqrs_im=Dpqrs.imag();
       index_primitive[0]=IPRIM+1;index_primitive[1]=IPRIM1+1;index_primitive_prime[0]=IPRIM2+1;index_primitive_prime[1]=IPRIM3+1;
       if(index_primitive[0]>Largest_Prim){Largest_Prim=index_primitive[0];}
       if(index_primitive[1]>Largest_Prim){Largest_Prim=index_primitive[1];}
       if(index_primitive_prime[0]>Largest_Prim){Largest_Prim=index_primitive_prime[0];}
       if(index_primitive_prime[0]>Largest_Prim){Largest_Prim=index_primitive_prime[1];}
       output_data.seekp(RECORD_DELIMITER_LENGTH, ios::cur);
       output_data.write((char*) &index_primitive[0], sizeof(index_primitive[0]));
       output_data.write((char*) &index_primitive[1], sizeof(index_primitive[1]));
       output_data.write((char*) &index_primitive_prime[0], sizeof(index_primitive_prime[0]));
       output_data.write((char*) &index_primitive_prime[1], sizeof(index_primitive_prime[1]));
       output_data.write((char*) &Dpqrs_re, sizeof(Dpqrs_re));
       output_data.write((char*) &Dpqrs_im, sizeof(Dpqrs_im));
       output_data.seekp(RECORD_DELIMITER_LENGTH, ios::cur);
      }
     }
    }
   }
  }
 } 
 Dpqrs_re=ZERO;
 Dpqrs_im=ZERO;
 index_primitive[0]=0;index_primitive[1]=0;index_primitive_prime[0]=0;index_primitive_prime[1]=0;
 output_data.seekp(RECORD_DELIMITER_LENGTH, ios::cur);
 output_data.write((char*) &index_primitive[0], sizeof(index_primitive[0]));
 output_data.write((char*) &index_primitive[1], sizeof(index_primitive[1]));
 output_data.write((char*) &index_primitive_prime[0], sizeof(index_primitive_prime[0]));
 output_data.write((char*) &index_primitive_prime[1], sizeof(index_primitive_prime[1]));
 output_data.write((char*) &Dpqrs_re, sizeof(Dpqrs_re));
 output_data.write((char*) &Dpqrs_im, sizeof(Dpqrs_im));
 output_data.seekp(RECORD_DELIMITER_LENGTH, ios::cur);
 output_data.close();
 cout<<endl;
 cout<<"Finished writing the transformed-complex 2-RDM elements in the primitives basis"<<endl;
 cout<<"Num. of printed terms: "<<setw(12)<<Nterms_printed<<endl;
 cout<<"Largest prim. in use : "<<setw(12)<<Largest_Prim<<endl;
 if(Largest_Prim<Nprimitives)
 {
  cout<<"Overwriting Nprimitives variable with largest prim. in use"<<endl;
  Nprimitives=Largest_Prim;
  Nprims1=Nprimitives; 
  Nprims2=Nprims1*Nprimitives; 
  Nprims3=Nprims2*Nprimitives;
  Nprims4=Nprims3*Nprimitives;
 }
 cout<<endl;
 delete[] Dijks_Prims;Dijks_Prims=NULL;
 delete[] Dijrs_Prims;Dijrs_Prims=NULL;
 delete[] Diqrs_Prims;Diqrs_Prims=NULL;
}  

// Read the 2-RDM in primitive basis, reduce it using symmetry, and print it real
void reduce_getreal_print()
{
 int index_Prim[2],index_Prim_prime[2];
 long int IPRIM,IPRIM1,IPRIM2,IPRIM3,IPRIM4;    // IPRIM (0 to Nprimitives) but can be summed for large Nprimitives number.
 double Dpqrs_Prim_RE,Dpqrs_Prim_IM,Dpqrs,MEM;
 string conv_name="Conv_cplx_"+dirac_output_name+"dm2";
 MEM=EIGHT*((TWO*EIGHT)*Nprims4)/pow(TEN,NINE);
 cout<<setprecision(2)<<fixed;
 cout<<"Memory required ";
 if(MEM>pow(TEN,THREE))
 {
  cout<<setw(10)<<MEM/pow(TEN,THREE)<<" Tb.";
 }
 else
 {
  if(MEM>ONE)
  {
   cout<<setw(10)<<MEM<<" Gb.";
  }
  else
  {
   if(MEM<ONE && MEM>pow(TEN,-THREE))
   {
    cout<<setw(10)<<MEM*pow(TEN,THREE)<<" Mb.";
   }
   else
   {
    cout<<setw(10)<<MEM*pow(TEN,SIX)<<" Kb.";
   }
  }
 }
 cout<<" for the reduction."<<endl; 
 cout<<"Num. of 2-RDM elem.  : "<<setw(12)<<Nprims4<<endl;
 cout<<endl;
 Dpqrs_ALL_cplx=new complex<double>[Nprims4];
 ifstream input_data2(conv_name.c_str(),ios::binary);
 cout<<endl;
 index_Prim[0]=10;index_Prim[1]=10;index_Prim_prime[0]=10;index_Prim_prime[1]=10;
 cout<<"Reading the transformed 2-RDM elements from "<<conv_name<<endl;
 while(index_Prim[0]!=0 || index_Prim[1]!=0 || index_Prim_prime[0]!=0 || index_Prim_prime[1]!=0)
 {
  input_data2.seekg(RECORD_DELIMITER_LENGTH, ios::cur);
  input_data2.read((char*) &index_Prim[0], sizeof(index_Prim[0]));
  input_data2.read((char*) &index_Prim[1], sizeof(index_Prim[1]));
  input_data2.read((char*) &index_Prim_prime[0], sizeof(index_Prim_prime[0]));
  input_data2.read((char*) &index_Prim_prime[1], sizeof(index_Prim_prime[1]));
  input_data2.read((char*) &Dpqrs_Prim_RE, sizeof(Dpqrs_Prim_RE));
  input_data2.read((char*) &Dpqrs_Prim_IM, sizeof(Dpqrs_Prim_IM));
  input_data2.seekg(RECORD_DELIMITER_LENGTH, ios::cur);
  complex<double>ztmpPRIM(Dpqrs_Prim_RE,Dpqrs_Prim_IM);
  if(abs(ztmpPRIM)!=ZERO && index_Prim[0]!=0 && index_Prim[1]!=0 && index_Prim_prime[0]!=0 && index_Prim_prime[1]!=0)
  {
   index_Prim[0]=index_Prim[0]-1;index_Prim[1]=index_Prim[1]-1;index_Prim_prime[0]=index_Prim_prime[0]-1;index_Prim_prime[1]=index_Prim_prime[1]-1;
   Dpqrs_ALL_cplx[index_Prim[0]+index_Prim[1]*Nprims1+index_Prim_prime[0]*Nprims2+index_Prim_prime[1]*Nprims3]= 
   Dpqrs_ALL_cplx[index_Prim[0]+index_Prim[1]*Nprims1+index_Prim_prime[0]*Nprims2+index_Prim_prime[1]*Nprims3]+ztmpPRIM;
   index_Prim[0]=index_Prim[0]+1;index_Prim[1]=index_Prim[1]+1;index_Prim_prime[0]=index_Prim_prime[0]+1;index_Prim_prime[1]=index_Prim_prime[1]+1;
  }
 }
 input_data2.close();
 // Reduce
 cout<<"Reducing the p <-> r and q <-> s terms"<<endl;
 for(IPRIM=0;IPRIM<Nprims4;IPRIM++)
 {
  IPRIM1=(IPRIM/Nprims3);
  IPRIM2=((IPRIM-IPRIM1*Nprims3)/Nprims2);
  IPRIM3=((IPRIM-IPRIM2*Nprims2-IPRIM1*Nprims3)/Nprims1);
  IPRIM4=(IPRIM-IPRIM3*Nprims1-IPRIM2*Nprims2-IPRIM1*Nprims3);
  if(IPRIM4!=IPRIM2 && IPRIM3!=IPRIM1)
  {
   Dpqrs_ALL_cplx[IPRIM]=Dpqrs_ALL_cplx[IPRIM]
                 +Dpqrs_ALL_cplx[IPRIM2+IPRIM3*Nprims1+IPRIM4*Nprims2+IPRIM1*Nprims3]
                 +Dpqrs_ALL_cplx[IPRIM4+IPRIM1*Nprims1+IPRIM2*Nprims2+IPRIM3*Nprims3]
                 +Dpqrs_ALL_cplx[IPRIM2+IPRIM1*Nprims1+IPRIM4*Nprims2+IPRIM3*Nprims3];
   Dpqrs_ALL_cplx[IPRIM2+IPRIM3*Nprims1+IPRIM4*Nprims2+IPRIM1*Nprims3]=CZERO;
   Dpqrs_ALL_cplx[IPRIM4+IPRIM1*Nprims1+IPRIM2*Nprims2+IPRIM3*Nprims3]=CZERO;
   Dpqrs_ALL_cplx[IPRIM2+IPRIM1*Nprims1+IPRIM4*Nprims2+IPRIM3*Nprims3]=CZERO;
  }
  if(IPRIM4!=IPRIM2 && IPRIM3==IPRIM1)
  {
   Dpqrs_ALL_cplx[IPRIM]=Dpqrs_ALL_cplx[IPRIM]
                 +Dpqrs_ALL_cplx[IPRIM2+IPRIM3*Nprims1+IPRIM4*Nprims2+IPRIM1*Nprims3];
   Dpqrs_ALL_cplx[IPRIM2+IPRIM3*Nprims1+IPRIM4*Nprims2+IPRIM1*Nprims3]=CZERO;
  }
  if(IPRIM4==IPRIM2 && IPRIM3!=IPRIM1)
  {
   Dpqrs_ALL_cplx[IPRIM]=Dpqrs_ALL_cplx[IPRIM]
                 +Dpqrs_ALL_cplx[IPRIM4+IPRIM1*Nprims1+IPRIM2*Nprims2+IPRIM3*Nprims3];
   Dpqrs_ALL_cplx[IPRIM4+IPRIM1*Nprims1+IPRIM2*Nprims2+IPRIM3*Nprims3]=CZERO;
  }
 }
 // Print
 cout<<"Printing the reduced (real) transformed 2-RDM"<<endl;
 // TODO: add this symmetry
 //
 //if(symmrr_prime){cout<<"keep only one term among prim_p(1) prim_q(2) prim_r(1) prim_s(2) = prim_q(2) prim_p(1) prim_s(2) prim_r(1)"<<endl;}
 //
 conv_name="Conv_"+dirac_output_name+"dm2";
 ofstream output_data(conv_name.c_str(),ios::binary);
 for(IPRIM=0;IPRIM<Nprims4;IPRIM++)
 {
  Dpqrs=Dpqrs_ALL_cplx[IPRIM].real();
  if(abs(Dpqrs)>=threshold)
  {
   index_Prim_prime[1]=(int)(IPRIM/Nprims3);
   index_Prim_prime[0]=(int)((IPRIM-index_Prim_prime[1]*Nprims3)/Nprims2);
   index_Prim[1]=(int)((IPRIM-index_Prim_prime[0]*Nprims2-index_Prim_prime[1]*Nprims3)/Nprims1);
   index_Prim[0]=(int)(IPRIM-index_Prim[1]*Nprims1-index_Prim_prime[0]*Nprims2-index_Prim_prime[1]*Nprims3);
   // Symmetry reduction?
   //
   //if(!(index_Prim[0]==index_Prim[1] && index_Prim[1]==index_Prim_prime[0] && index_Prim_prime[0]==index_Prim_prime[1]) && symmrr_prime)
   //{
   // Dpqrs_ALL_cplx[index_Prim[0]+index_Prim[1]*Nprims1+index_Prim_prime[0]*Nprims2+index_Prim_prime[1]*Nprims3]=CZERO;
   // Dpqrs_ALL_cplx[index_Prim[1]+index_Prim[0]*Nprims1+index_Prim_prime[1]*Nprims2+index_Prim_prime[0]*Nprims3]=CZERO;
   //}
   //
   index_Prim[0]++;index_Prim[1]++;index_Prim_prime[0]++;index_Prim_prime[1]++;
   output_data.seekp(RECORD_DELIMITER_LENGTH, ios::cur);
   output_data.write((char*) &index_Prim[0], sizeof(index_Prim[0]));
   output_data.write((char*) &index_Prim[1], sizeof(index_Prim[1]));
   output_data.write((char*) &index_Prim_prime[0], sizeof(index_Prim_prime[0]));
   output_data.write((char*) &index_Prim_prime[1], sizeof(index_Prim_prime[1]));
   output_data.write((char*) &Dpqrs, sizeof(Dpqrs));
   output_data.seekp(RECORD_DELIMITER_LENGTH, ios::cur);
  }
 }
 index_Prim[0]=0;
 index_Prim[1]=0;
 index_Prim_prime[0]=0;
 index_Prim_prime[1]=0;
 Dpqrs=ZERO;
 output_data.seekp(RECORD_DELIMITER_LENGTH, ios::cur);
 output_data.write((char*) &index_Prim[0], sizeof(index_Prim[0]));
 output_data.write((char*) &index_Prim[1], sizeof(index_Prim[1]));
 output_data.write((char*) &index_Prim_prime[0], sizeof(index_Prim_prime[0]));
 output_data.write((char*) &index_Prim_prime[1], sizeof(index_Prim_prime[1]));
 output_data.write((char*) &Dpqrs, sizeof(Dpqrs));
 output_data.seekp(RECORD_DELIMITER_LENGTH, ios::cur);
 output_data.close();
 delete[] Dpqrs_ALL_cplx;Dpqrs_ALL_cplx=NULL;
}
*/

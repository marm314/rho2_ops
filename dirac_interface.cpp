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

using namespace std;

void print_basis_file();
void read_dirac_out();
void clean_shell2aos();
void print_wfx();

int Nprimitives,Nbasis,Nbasis_L,Nbasis_S,Nshell,Nshell_L,Nshell_S,NMOs,NMOs_LS,NMOs_occ,OneMO_wfx=-1;
struct Shell2AOs
{
 int styp,atom,nprim,naos,paired2=-1;
 int **index_min2max; 
 double *Coef,*Expon,Coord[3];
};
Shell2AOs *shell2aos;
double Quaternion_coef[4];
double *OCCs,**Prim2AO_Coef;
complex<double> **AO2MO_Coef,**Prim2MO_Coef;
string dirac_output_name,dirac_output_file;
vector<int>shell_types;
vector<double>prim_exponents;
vector<double>MOsLS_occ;
vector<double>One_Prim2MO_Coef_RE;
vector<double>One_Prim2MO_Coef_IM;

int main(int argc, char *argv[])
{
 cout<<"----------------------------------------"<<endl;
 cout<<"----------------------------------------"<<endl;
 cout<<"---  DIRAC INTERFACE FOR RHO2_OPS    ---"<<endl;
 cout<<"----------------------------------------"<<endl;
 cout<<"----------------------------------------"<<endl;
 cout<<"-- Developed by: M. Rodriguez-Mayorga --"<<endl;
 cout<<"--   email: marm3.14@gmail.com        --"<<endl;
 cout<<"----------------------------------------"<<endl;
 cout<<"----------------------------------------"<<endl;
 cout<<endl;
 if(argc!=2 && argc!=3)
 {
  cout<<endl;
  cout<<"Please, Include the name of the DIRAC program as an argument."<<endl;
  cout<<endl;
  cout<<"name_name.out"<<endl;
  cout<<"or"<<endl;
  cout<<"name_name.out one_mo_wfx(integer optional)"<<endl;
  cout<<endl;
  cout<<"----------------------------------------"<<endl;
  cout<<"--        Normal termination          --"<<endl;
  cout<<"----------------------------------------"<<endl;
  cout<<"----------------------------------------"<<endl;
  return -1;
 }
 bool repeated_prims;
 int ishell,ishell1,iprim,iprim1,iaos,iaos1,imos,imos1,imos2;
 int naos;
 dirac_output_file=argv[1];
 dirac_output_name=dirac_output_file.substr(0,dirac_output_file.length()-3);
 if(argc==3)
 {
  OneMO_wfx=atoi(argv[2]);
  OneMO_wfx=(OneMO_wfx-1)*4; // even->unbar, odd->bar
  cout<<"Orbital selection is swittched on. Scalar (LS) orbitals to print in the WFX file: "<<setw(5)<<OneMO_wfx+1<<" to "<<setw(5)<<OneMO_wfx+4<<endl;
 }
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
 cout<<"Building the matrix Prim2MO_Coefs[MO][Primitives] ((2x4MO) x 2(L+S)) = AO2MO_Coef[MO][AO]*Prim2AO_Coef[AO][Primitives]"<<endl;
 Prim2MO_Coef=new complex<double>*[NMOs_LS]; 
 for(imos=0;imos<NMOs_LS;imos++)
 {
  Prim2MO_Coef[imos]=new complex<double>[Nprimitives];
  for(iprim=0;iprim<Nprimitives;iprim++)
  {
   complex<double>ztmp0=(ZERO,ZERO);
   Prim2MO_Coef[imos][iprim]=ztmp0;
   for(iaos=0;iaos<Nbasis;iaos++)
   {
    Prim2MO_Coef[imos][iprim]=Prim2MO_Coef[imos][iprim]+AO2MO_Coef[imos][iaos]*Prim2AO_Coef[iaos][iprim];
   }
  }
 }
 cout<<"Matrix Prim2MO_Coefs[MO][Primitives] for bar and unbar orbs. built."<<endl;
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
 ofstream coefs_file((dirac_output_name+"coef").c_str());
 ofstream coefs_file_pos((dirac_output_name+"coef_pos").c_str());
 imos1=0;
 coefs_file<<setprecision(12)<<fixed<<scientific;
 coefs_file_pos<<setprecision(12)<<fixed<<scientific;
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
  if(OCCs[imos]==-TEN) // Write positronic MO coefs to file
  {
   for(iprim=0;iprim<Nprimitives;iprim++)
   {
    coefs_file_pos<<setw(20)<<Prim2MO_Coef[imos][iprim].real()<<setw(20)<<Prim2MO_Coef[imos][iprim].imag();
   }
  }
  if(OCCs[imos]>pow(TEN,-EIGHT)) // Overwrite positronic MO coefs with electronic ones and print the electronic ones
  {
   for(iprim=0;iprim<Nprimitives;iprim++)
   {
    Prim2MO_Coef[imos1][iprim]=Prim2MO_Coef[imos][iprim];
    coefs_file<<setw(20)<<Prim2MO_Coef[imos1][iprim].real()<<setw(20)<<Prim2MO_Coef[imos1][iprim].imag();
   }
   MOsLS_occ.push_back(OCCs[imos]);
   imos1++;
  }
 }
 NMOs_occ=imos1;
 coefs_file.close();
 coefs_file_pos.close();
 cout<<"Num. of occ MO 2(L+S): "<<setw(12)<<NMOs_occ<<endl;
 // Print a WFX file for RHO_OPS (currently only available for ATOMS)
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
 // Read 2-RDM (binary)
 // TODO: Read the 2-RDM and (anti)symmetrize it. Then, transform it to 4component (expand indices to scalars).
 

 // Deallocate arrays Primitive to MO coefs (rows).
 for(imos=0;imos<NMOs_LS;imos++)
 {
  delete [] Prim2MO_Coef[imos];Prim2MO_Coef[imos]=NULL;
 }
 delete[] Prim2MO_Coef;Prim2MO_Coef=NULL;
 delete[] OCCs;OCCs=NULL;
 cout<<endl;
 cout<<"----------------------------------------"<<endl;
 cout<<"--        Normal termination          --"<<endl;
 cout<<"----------------------------------------"<<endl;
 cout<<"----------------------------------------"<<endl;
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
    cout<<"Number of 2(L+S) MOs : "<<setw(12)<<NMOs_LS<<endl;
    cout<<"Number of AOs (read) : "<<setw(12)<<Nbasis_int<<endl;
    AO2MO_Coef=new complex<double>*[NMOs_LS];
    for(imos=0;imos<NMOs_LS;imos++)
    {
     AO2MO_Coef[imos]=new complex<double>[Nbasis];
     for(iaos=0;iaos<Nbasis;iaos++)
     {
      complex<double>ztmp0=(ZERO,ZERO);
      AO2MO_Coef[imos][iaos]=ztmp0;
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
      complex<double>ztmp0( Quaternion_coef[0], Quaternion_coef[1]); 
      complex<double>ztmp1(-Quaternion_coef[2], Quaternion_coef[3]); 
      complex<double>ztmp2( Quaternion_coef[2], Quaternion_coef[3]); 
      complex<double>ztmp3( Quaternion_coef[0],-Quaternion_coef[1]); 
      AO2MO_Coef[imos1][iaos]=ztmp0; // unbar alpha L
      AO2MO_Coef[imos2][iaos]=ztmp1; // unbar beta L
      AO2MO_Coef[imos5][iaos]=ztmp2; // bar alpha L
      AO2MO_Coef[imos6][iaos]=ztmp3; // bar beta L
     }
     //cout<<"Quaternion to Small component MOs transformation for MO"<<setw(5)<<imos+1<<endl;
     for(iaos=Nbasis_L;iaos<Nbasis;iaos++)
     {
      Dirac_file>>Quaternion_coef[0]>>Quaternion_coef[1]>>Quaternion_coef[2]>>Quaternion_coef[3];
      //cout<<Quaternion_coef[0]<<Quaternion_coef[1]<<Quaternion_coef[2]<<Quaternion_coef[3]<<endl; 
      complex<double>ztmp0( Quaternion_coef[0], Quaternion_coef[1]); 
      complex<double>ztmp1(-Quaternion_coef[2], Quaternion_coef[3]); 
      complex<double>ztmp2( Quaternion_coef[2], Quaternion_coef[3]); 
      complex<double>ztmp3( Quaternion_coef[0],-Quaternion_coef[1]); 
      AO2MO_Coef[imos3][iaos]=ztmp0;  // unbar alpha S
      AO2MO_Coef[imos4][iaos]=ztmp1;  // unbar beta S
      AO2MO_Coef[imos7][iaos]=ztmp2;  // bar alpha S
      AO2MO_Coef[imos8][iaos]=ztmp3;  // bar beta S
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
 int ishell,iprim;
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
     }
     // 0 1 0 
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      outbasis<<setw(3)<<1<<"\t"<<setprecision(15)<<fixed<<scientific<<shell2aos[ishell].Coord[0]<<"\t"<<shell2aos[ishell].Coord[1];
      outbasis<<"\t"<<shell2aos[ishell].Coord[2]<<"\t"<<shell2aos[ishell].Expon[iprim];
      outbasis<<"\t\t"<<0<<"\t"<<1<<"\t"<<0<<endl;
      shell_types.push_back(3);
      prim_exponents.push_back(shell2aos[ishell].Expon[iprim]);
     }
     // 0 0 1 
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      outbasis<<setw(3)<<1<<"\t"<<setprecision(15)<<fixed<<scientific<<shell2aos[ishell].Coord[0]<<"\t"<<shell2aos[ishell].Coord[1];
      outbasis<<"\t"<<shell2aos[ishell].Coord[2]<<"\t"<<shell2aos[ishell].Expon[iprim];
      outbasis<<"\t\t"<<0<<"\t"<<0<<"\t"<<1<<endl;
      shell_types.push_back(4);
      prim_exponents.push_back(shell2aos[ishell].Expon[iprim]);
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
     }
     // 1 1 0 
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      outbasis<<setw(3)<<1<<"\t"<<setprecision(15)<<fixed<<scientific<<shell2aos[ishell].Coord[0]<<"\t"<<shell2aos[ishell].Coord[1];
      outbasis<<"\t"<<shell2aos[ishell].Coord[2]<<"\t"<<shell2aos[ishell].Expon[iprim];
      outbasis<<"\t\t"<<1<<"\t"<<1<<"\t"<<0<<endl;
      shell_types.push_back(8);
      prim_exponents.push_back(shell2aos[ishell].Expon[iprim]);
     }
     // 1 0 1 
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      outbasis<<setw(3)<<1<<"\t"<<setprecision(15)<<fixed<<scientific<<shell2aos[ishell].Coord[0]<<"\t"<<shell2aos[ishell].Coord[1];
      outbasis<<"\t"<<shell2aos[ishell].Coord[2]<<"\t"<<shell2aos[ishell].Expon[iprim];
      outbasis<<"\t\t"<<1<<"\t"<<0<<"\t"<<1<<endl;
      shell_types.push_back(9);
      prim_exponents.push_back(shell2aos[ishell].Expon[iprim]);
     }
     // 0 2 0 
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      outbasis<<setw(3)<<1<<"\t"<<setprecision(15)<<fixed<<scientific<<shell2aos[ishell].Coord[0]<<"\t"<<shell2aos[ishell].Coord[1];
      outbasis<<"\t"<<shell2aos[ishell].Coord[2]<<"\t"<<shell2aos[ishell].Expon[iprim];
      outbasis<<"\t\t"<<0<<"\t"<<2<<"\t"<<0<<endl;
      shell_types.push_back(6);
      prim_exponents.push_back(shell2aos[ishell].Expon[iprim]);
     }
     // 0 1 1 
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      outbasis<<setw(3)<<1<<"\t"<<setprecision(15)<<fixed<<scientific<<shell2aos[ishell].Coord[0]<<"\t"<<shell2aos[ishell].Coord[1];
      outbasis<<"\t"<<shell2aos[ishell].Coord[2]<<"\t"<<shell2aos[ishell].Expon[iprim];
      outbasis<<"\t\t"<<0<<"\t"<<1<<"\t"<<1<<endl;
      shell_types.push_back(10);
      prim_exponents.push_back(shell2aos[ishell].Expon[iprim]);
     }
     // 0 0 2 
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      outbasis<<setw(3)<<1<<"\t"<<setprecision(15)<<fixed<<scientific<<shell2aos[ishell].Coord[0]<<"\t"<<shell2aos[ishell].Coord[1];
      outbasis<<"\t"<<shell2aos[ishell].Coord[2]<<"\t"<<shell2aos[ishell].Expon[iprim];
      outbasis<<"\t\t"<<0<<"\t"<<0<<"\t"<<2<<endl;
      shell_types.push_back(7);
      prim_exponents.push_back(shell2aos[ishell].Expon[iprim]);
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
     }
     // 2 1 0 
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      outbasis<<setw(3)<<1<<"\t"<<setprecision(15)<<fixed<<scientific<<shell2aos[ishell].Coord[0]<<"\t"<<shell2aos[ishell].Coord[1];
      outbasis<<"\t"<<shell2aos[ishell].Coord[2]<<"\t"<<shell2aos[ishell].Expon[iprim];
      outbasis<<"\t\t"<<2<<"\t"<<1<<"\t"<<0<<endl;
      shell_types.push_back(14);
      prim_exponents.push_back(shell2aos[ishell].Expon[iprim]);
     }
     // 2 0 1 
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      outbasis<<setw(3)<<1<<"\t"<<setprecision(15)<<fixed<<scientific<<shell2aos[ishell].Coord[0]<<"\t"<<shell2aos[ishell].Coord[1];
      outbasis<<"\t"<<shell2aos[ishell].Coord[2]<<"\t"<<shell2aos[ishell].Expon[iprim];
      outbasis<<"\t\t"<<2<<"\t"<<0<<"\t"<<1<<endl;
      shell_types.push_back(15);
      prim_exponents.push_back(shell2aos[ishell].Expon[iprim]);
     }
     // 1 2 0 
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      outbasis<<setw(3)<<1<<"\t"<<setprecision(15)<<fixed<<scientific<<shell2aos[ishell].Coord[0]<<"\t"<<shell2aos[ishell].Coord[1];
      outbasis<<"\t"<<shell2aos[ishell].Coord[2]<<"\t"<<shell2aos[ishell].Expon[iprim];
      outbasis<<"\t\t"<<1<<"\t"<<2<<"\t"<<0<<endl;
      shell_types.push_back(17);
      prim_exponents.push_back(shell2aos[ishell].Expon[iprim]);
     }
     // 1 1 1 
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      outbasis<<setw(3)<<1<<"\t"<<setprecision(15)<<fixed<<scientific<<shell2aos[ishell].Coord[0]<<"\t"<<shell2aos[ishell].Coord[1];
      outbasis<<"\t"<<shell2aos[ishell].Coord[2]<<"\t"<<shell2aos[ishell].Expon[iprim];
      outbasis<<"\t\t"<<1<<"\t"<<1<<"\t"<<1<<endl;
      shell_types.push_back(20);
      prim_exponents.push_back(shell2aos[ishell].Expon[iprim]);
     }
     // 1 0 2 
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      outbasis<<setw(3)<<1<<"\t"<<setprecision(15)<<fixed<<scientific<<shell2aos[ishell].Coord[0]<<"\t"<<shell2aos[ishell].Coord[1];
      outbasis<<"\t"<<shell2aos[ishell].Coord[2]<<"\t"<<shell2aos[ishell].Expon[iprim];
      outbasis<<"\t\t"<<1<<"\t"<<0<<"\t"<<2<<endl;
      shell_types.push_back(18);
      prim_exponents.push_back(shell2aos[ishell].Expon[iprim]);
     }
     // 0 3 0 
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      outbasis<<setw(3)<<1<<"\t"<<setprecision(15)<<fixed<<scientific<<shell2aos[ishell].Coord[0]<<"\t"<<shell2aos[ishell].Coord[1];
      outbasis<<"\t"<<shell2aos[ishell].Coord[2]<<"\t"<<shell2aos[ishell].Expon[iprim];
      outbasis<<"\t\t"<<0<<"\t"<<3<<"\t"<<0<<endl;
      shell_types.push_back(12);
      prim_exponents.push_back(shell2aos[ishell].Expon[iprim]);
     }
     // 0 2 1 
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      outbasis<<setw(3)<<1<<"\t"<<setprecision(15)<<fixed<<scientific<<shell2aos[ishell].Coord[0]<<"\t"<<shell2aos[ishell].Coord[1];
      outbasis<<"\t"<<shell2aos[ishell].Coord[2]<<"\t"<<shell2aos[ishell].Expon[iprim];
      outbasis<<"\t\t"<<0<<"\t"<<2<<"\t"<<1<<endl;
      shell_types.push_back(16);
      prim_exponents.push_back(shell2aos[ishell].Expon[iprim]);
     }
     // 0 1 2 
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      outbasis<<setw(3)<<1<<"\t"<<setprecision(15)<<fixed<<scientific<<shell2aos[ishell].Coord[0]<<"\t"<<shell2aos[ishell].Coord[1];
      outbasis<<"\t"<<shell2aos[ishell].Coord[2]<<"\t"<<shell2aos[ishell].Expon[iprim];
      outbasis<<"\t\t"<<0<<"\t"<<1<<"\t"<<2<<endl;
      shell_types.push_back(19);
      prim_exponents.push_back(shell2aos[ishell].Expon[iprim]);
     }
     // 0 0 3 
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      outbasis<<setw(3)<<1<<"\t"<<setprecision(15)<<fixed<<scientific<<shell2aos[ishell].Coord[0]<<"\t"<<shell2aos[ishell].Coord[1];
      outbasis<<"\t"<<shell2aos[ishell].Coord[2]<<"\t"<<shell2aos[ishell].Expon[iprim];
      outbasis<<"\t\t"<<0<<"\t"<<0<<"\t"<<3<<endl;
      shell_types.push_back(13);
      prim_exponents.push_back(shell2aos[ishell].Expon[iprim]);
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
     }
     // 3 1 0
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      outbasis<<setw(3)<<1<<"\t"<<setprecision(15)<<fixed<<scientific<<shell2aos[ishell].Coord[0]<<"\t"<<shell2aos[ishell].Coord[1];
      outbasis<<"\t"<<shell2aos[ishell].Coord[2]<<"\t"<<shell2aos[ishell].Expon[iprim];
      outbasis<<"\t\t"<<3<<"\t"<<1<<"\t"<<0<<endl;
      shell_types.push_back(24);
      prim_exponents.push_back(shell2aos[ishell].Expon[iprim]);
     }
     // 3 0 1
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      outbasis<<setw(3)<<1<<"\t"<<setprecision(15)<<fixed<<scientific<<shell2aos[ishell].Coord[0]<<"\t"<<shell2aos[ishell].Coord[1];
      outbasis<<"\t"<<shell2aos[ishell].Coord[2]<<"\t"<<shell2aos[ishell].Expon[iprim];
      outbasis<<"\t\t"<<3<<"\t"<<0<<"\t"<<1<<endl;
      shell_types.push_back(25);
      prim_exponents.push_back(shell2aos[ishell].Expon[iprim]);
     }
     // 2 2 0
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      outbasis<<setw(3)<<1<<"\t"<<setprecision(15)<<fixed<<scientific<<shell2aos[ishell].Coord[0]<<"\t"<<shell2aos[ishell].Coord[1];
      outbasis<<"\t"<<shell2aos[ishell].Coord[2]<<"\t"<<shell2aos[ishell].Expon[iprim];
      outbasis<<"\t\t"<<2<<"\t"<<2<<"\t"<<0<<endl;
      shell_types.push_back(30);
      prim_exponents.push_back(shell2aos[ishell].Expon[iprim]);
     }
     // 2 1 1
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      outbasis<<setw(3)<<1<<"\t"<<setprecision(15)<<fixed<<scientific<<shell2aos[ishell].Coord[0]<<"\t"<<shell2aos[ishell].Coord[1];
      outbasis<<"\t"<<shell2aos[ishell].Coord[2]<<"\t"<<shell2aos[ishell].Expon[iprim];
      outbasis<<"\t\t"<<2<<"\t"<<1<<"\t"<<1<<endl;
      shell_types.push_back(33);
      prim_exponents.push_back(shell2aos[ishell].Expon[iprim]);
     }
     // 2 0 2
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      outbasis<<setw(3)<<1<<"\t"<<setprecision(15)<<fixed<<scientific<<shell2aos[ishell].Coord[0]<<"\t"<<shell2aos[ishell].Coord[1];
      outbasis<<"\t"<<shell2aos[ishell].Coord[2]<<"\t"<<shell2aos[ishell].Expon[iprim];
      outbasis<<"\t\t"<<2<<"\t"<<0<<"\t"<<2<<endl;
      shell_types.push_back(31);
      prim_exponents.push_back(shell2aos[ishell].Expon[iprim]);
     }
     // 1 3 0
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      outbasis<<setw(3)<<1<<"\t"<<setprecision(15)<<fixed<<scientific<<shell2aos[ishell].Coord[0]<<"\t"<<shell2aos[ishell].Coord[1];
      outbasis<<"\t"<<shell2aos[ishell].Coord[2]<<"\t"<<shell2aos[ishell].Expon[iprim];
      outbasis<<"\t\t"<<1<<"\t"<<3<<"\t"<<0<<endl;
      shell_types.push_back(26);
      prim_exponents.push_back(shell2aos[ishell].Expon[iprim]);
     }
     // 1 2 1
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      outbasis<<setw(3)<<1<<"\t"<<setprecision(15)<<fixed<<scientific<<shell2aos[ishell].Coord[0]<<"\t"<<shell2aos[ishell].Coord[1];
      outbasis<<"\t"<<shell2aos[ishell].Coord[2]<<"\t"<<shell2aos[ishell].Expon[iprim];
      outbasis<<"\t\t"<<1<<"\t"<<2<<"\t"<<1<<endl;
      shell_types.push_back(34);
      prim_exponents.push_back(shell2aos[ishell].Expon[iprim]);
     }
     // 1 1 2
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      outbasis<<setw(3)<<1<<"\t"<<setprecision(15)<<fixed<<scientific<<shell2aos[ishell].Coord[0]<<"\t"<<shell2aos[ishell].Coord[1];
      outbasis<<"\t"<<shell2aos[ishell].Coord[2]<<"\t"<<shell2aos[ishell].Expon[iprim];
      outbasis<<"\t\t"<<1<<"\t"<<1<<"\t"<<2<<endl;
      shell_types.push_back(35);
      prim_exponents.push_back(shell2aos[ishell].Expon[iprim]);
     }
     // 1 0 3
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      outbasis<<setw(3)<<1<<"\t"<<setprecision(15)<<fixed<<scientific<<shell2aos[ishell].Coord[0]<<"\t"<<shell2aos[ishell].Coord[1];
      outbasis<<"\t"<<shell2aos[ishell].Coord[2]<<"\t"<<shell2aos[ishell].Expon[iprim];
      outbasis<<"\t\t"<<1<<"\t"<<0<<"\t"<<3<<endl;
      shell_types.push_back(28);
      prim_exponents.push_back(shell2aos[ishell].Expon[iprim]);
     }
     // 0 4 0
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      outbasis<<setw(3)<<1<<"\t"<<setprecision(15)<<fixed<<scientific<<shell2aos[ishell].Coord[0]<<"\t"<<shell2aos[ishell].Coord[1];
      outbasis<<"\t"<<shell2aos[ishell].Coord[2]<<"\t"<<shell2aos[ishell].Expon[iprim];
      outbasis<<"\t\t"<<0<<"\t"<<4<<"\t"<<0<<endl;
      shell_types.push_back(22);
      prim_exponents.push_back(shell2aos[ishell].Expon[iprim]);
     }
     // 0 3 1
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      outbasis<<setw(3)<<1<<"\t"<<setprecision(15)<<fixed<<scientific<<shell2aos[ishell].Coord[0]<<"\t"<<shell2aos[ishell].Coord[1];
      outbasis<<"\t"<<shell2aos[ishell].Coord[2]<<"\t"<<shell2aos[ishell].Expon[iprim];
      outbasis<<"\t\t"<<0<<"\t"<<3<<"\t"<<1<<endl;
      shell_types.push_back(27);
      prim_exponents.push_back(shell2aos[ishell].Expon[iprim]);
     }
     // 0 2 2
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      outbasis<<setw(3)<<1<<"\t"<<setprecision(15)<<fixed<<scientific<<shell2aos[ishell].Coord[0]<<"\t"<<shell2aos[ishell].Coord[1];
      outbasis<<"\t"<<shell2aos[ishell].Coord[2]<<"\t"<<shell2aos[ishell].Expon[iprim];
      outbasis<<"\t\t"<<0<<"\t"<<2<<"\t"<<2<<endl;
      shell_types.push_back(32);
      prim_exponents.push_back(shell2aos[ishell].Expon[iprim]);
     }
     // 0 1 3
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      outbasis<<setw(3)<<1<<"\t"<<setprecision(15)<<fixed<<scientific<<shell2aos[ishell].Coord[0]<<"\t"<<shell2aos[ishell].Coord[1];
      outbasis<<"\t"<<shell2aos[ishell].Coord[2]<<"\t"<<shell2aos[ishell].Expon[iprim];
      outbasis<<"\t\t"<<0<<"\t"<<1<<"\t"<<3<<endl;
      shell_types.push_back(29);
      prim_exponents.push_back(shell2aos[ishell].Expon[iprim]);
     }
     // 0 0 4
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      outbasis<<setw(3)<<1<<"\t"<<setprecision(15)<<fixed<<scientific<<shell2aos[ishell].Coord[0]<<"\t"<<shell2aos[ishell].Coord[1];
      outbasis<<"\t"<<shell2aos[ishell].Coord[2]<<"\t"<<shell2aos[ishell].Expon[iprim];
      outbasis<<"\t\t"<<0<<"\t"<<0<<"\t"<<4<<endl;
      shell_types.push_back(23);
      prim_exponents.push_back(shell2aos[ishell].Expon[iprim]);
     }
     break;
    default:
     cout<<"Warning! Shell-type not supported. "<<setw(5)<<shell2aos[ishell].styp<<endl;
     break;
   }
  }
 }
 outbasis.close();
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

// Currently, only available for ATOMS
void print_wfx()
{
 int iprim,imos,imos1;
 string line;
 ofstream real_wfx(("dirac_"+dirac_output_name+"RE.wfx").c_str());
 ofstream imag_wfx(("dirac_"+dirac_output_name+"IM.wfx").c_str());
 real_wfx<<setprecision(15)<<fixed<<scientific;
 imag_wfx<<setprecision(15)<<fixed<<scientific;
 line="<Number of Nuclei>";
 real_wfx<<line<<endl; 
 imag_wfx<<line<<endl; 
 real_wfx<<1<<endl; 
 imag_wfx<<1<<endl; 
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
 line="<Electronic Spin Multiplicity>";
 real_wfx<<line<<endl;
 imag_wfx<<line<<endl;
 real_wfx<<1<<endl; 
 imag_wfx<<1<<endl; 
 line="</Electronic Spin Multiplicity>";
 real_wfx<<line<<endl;
 imag_wfx<<line<<endl;
 line="<Nuclear Charges>"; // Faked
 real_wfx<<line<<endl;
 imag_wfx<<line<<endl;
 real_wfx<<0<<endl; 
 imag_wfx<<0<<endl; 
 line="</Nuclear Charges>"; // Faked
 real_wfx<<line<<endl;
 imag_wfx<<line<<endl;
 line="<Nuclear Cartesian Coordinates>"; // Only atomic systems 
 real_wfx<<line<<endl;
 imag_wfx<<line<<endl;
 real_wfx<<ZERO<<endl; // x
 imag_wfx<<ZERO<<endl; 
 real_wfx<<ZERO<<endl; // y
 imag_wfx<<ZERO<<endl; 
 real_wfx<<ZERO<<endl; // z
 imag_wfx<<ZERO<<endl; 
 line="</Nuclear Cartesian Coordinates>"; // Only atomic systems 
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
  real_wfx<<1<<endl;
  imag_wfx<<1<<endl;
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

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

int Nprimitives,Nbasis,Nbasis_L,Nbasis_S,Nshell,Nshell_L,Nshell_S,NMOs,NMOs_LS,NMOs_occ;
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
 if(argc!=2)
 {
  cout<<endl;
  cout<<"Please, Include the name of the DIRAC program as an argument."<<endl;
  cout<<endl;
  cout<<"----------------------------------------"<<endl;
  cout<<"--        Normal termination          --"<<endl;
  cout<<"----------------------------------------"<<endl;
  cout<<"----------------------------------------"<<endl;
  return -1;
 }
 bool repeated_prims;
 int ishell,ishell1,iprim,iprim1,iaos,iaos1,imos,imos1;
 int naos;
 dirac_output_file=argv[1];
 dirac_output_name=dirac_output_file.substr(0,dirac_output_file.length()-3);
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
 cout<<"Number of primitives : "<<setw(12)<<Nprimitives<<endl;
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
 cout<<"AO to primitives map (per shell):"<<endl;
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
 cout<<"Building Matrix Prim2MO_Coefs[imos][iprim] (4MO x 2(L+S)) = AO2MO_Coef[imos][iaos]*Prim2AO_Coef[iaos][iprim]"<<endl;
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
 cout<<"Matrix Prim2MO_Coefs[imos][iprim] built."<<endl;
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
 // WARNING! WARNING! Overwrite the positronic states with the electronic states ("sort electronic states")!
 ofstream coefs_file((dirac_output_name+"coef").c_str());
 ofstream coefs_file_pos((dirac_output_name+"coef_pos").c_str());
 imos1=0;
 coefs_file<<setprecision(12)<<fixed<<scientific;
 coefs_file_pos<<setprecision(12)<<fixed<<scientific;
 for(imos=0;imos<NMOs_LS;imos++)
 {
  if(OCCs[imos]==-TEN) // Write positronic states
  {
   for(iprim=0;iprim<Nprimitives;iprim++)
   {
    coefs_file_pos<<setw(20)<<Prim2MO_Coef[imos][iprim].real()<<setw(20)<<Prim2MO_Coef[imos][iprim].imag();
   }
  }
  if(OCCs[imos]>pow(TEN,-EIGHT)) // Overwrite positronic with electronic states
  {
   for(iprim=0;iprim<Nprimitives;iprim++)
   {
    Prim2MO_Coef[imos1][iprim]=Prim2MO_Coef[imos][iprim];
    coefs_file<<setw(20)<<Prim2MO_Coef[imos1][iprim].real()<<setw(20)<<Prim2MO_Coef[imos1][iprim].imag();
   }
   imos1++;
  }
 }
 NMOs_occ=imos1;
 coefs_file.close();
 coefs_file_pos.close();
 cout<<"Num. of occ MO 2(L+S): "<<setw(12)<<NMOs_occ<<endl;
 // Read 2-RDM (binary)
 // TODO: Check that 1st MOs states are electronic! read the 2-RDM and symmetrize it AND transform it to 4component (expand indices!) and print it
 print_wfx();
 

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
 bool electronic;
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
   if(line.substr(0,7)==" NBasis")
   {
    line=line.substr(7,line.length()-7);
    stringstream ss(line);
    ss>>Nbasis_int;
    if(Nbasis!=Nbasis_int)
    {
     cout<<"Warning! Nbasis read is not equal to the calculated value."<<endl;
    }
    cout<<"Number of quat. MOs  : "<<setw(12)<<NMOs<<endl;
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
     cout<<"Quaternion to Large component MOs transformation for MO"<<setw(5)<<imos+1<<endl;
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
     cout<<"Quaternion to Small component MOs transformation for MO"<<setw(5)<<imos+1<<endl;
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
     }
     // 0 1 0 
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      outbasis<<setw(3)<<1<<"\t"<<setprecision(15)<<fixed<<scientific<<shell2aos[ishell].Coord[0]<<"\t"<<shell2aos[ishell].Coord[1];
      outbasis<<"\t"<<shell2aos[ishell].Coord[2]<<"\t"<<shell2aos[ishell].Expon[iprim];
      outbasis<<"\t\t"<<0<<"\t"<<1<<"\t"<<0<<endl;
     }
     // 0 0 1 
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      outbasis<<setw(3)<<1<<"\t"<<setprecision(15)<<fixed<<scientific<<shell2aos[ishell].Coord[0]<<"\t"<<shell2aos[ishell].Coord[1];
      outbasis<<"\t"<<shell2aos[ishell].Coord[2]<<"\t"<<shell2aos[ishell].Expon[iprim];
      outbasis<<"\t\t"<<0<<"\t"<<0<<"\t"<<1<<endl;
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
     }
     // 1 1 0 
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      outbasis<<setw(3)<<1<<"\t"<<setprecision(15)<<fixed<<scientific<<shell2aos[ishell].Coord[0]<<"\t"<<shell2aos[ishell].Coord[1];
      outbasis<<"\t"<<shell2aos[ishell].Coord[2]<<"\t"<<shell2aos[ishell].Expon[iprim];
      outbasis<<"\t\t"<<1<<"\t"<<1<<"\t"<<0<<endl;
     }
     // 1 0 1 
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      outbasis<<setw(3)<<1<<"\t"<<setprecision(15)<<fixed<<scientific<<shell2aos[ishell].Coord[0]<<"\t"<<shell2aos[ishell].Coord[1];
      outbasis<<"\t"<<shell2aos[ishell].Coord[2]<<"\t"<<shell2aos[ishell].Expon[iprim];
      outbasis<<"\t\t"<<1<<"\t"<<0<<"\t"<<1<<endl;
     }
     // 0 2 0 
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      outbasis<<setw(3)<<1<<"\t"<<setprecision(15)<<fixed<<scientific<<shell2aos[ishell].Coord[0]<<"\t"<<shell2aos[ishell].Coord[1];
      outbasis<<"\t"<<shell2aos[ishell].Coord[2]<<"\t"<<shell2aos[ishell].Expon[iprim];
      outbasis<<"\t\t"<<0<<"\t"<<2<<"\t"<<0<<endl;
     }
     // 0 1 1 
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      outbasis<<setw(3)<<1<<"\t"<<setprecision(15)<<fixed<<scientific<<shell2aos[ishell].Coord[0]<<"\t"<<shell2aos[ishell].Coord[1];
      outbasis<<"\t"<<shell2aos[ishell].Coord[2]<<"\t"<<shell2aos[ishell].Expon[iprim];
      outbasis<<"\t\t"<<0<<"\t"<<1<<"\t"<<1<<endl;
     }
     // 0 0 2 
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      outbasis<<setw(3)<<1<<"\t"<<setprecision(15)<<fixed<<scientific<<shell2aos[ishell].Coord[0]<<"\t"<<shell2aos[ishell].Coord[1];
      outbasis<<"\t"<<shell2aos[ishell].Coord[2]<<"\t"<<shell2aos[ishell].Expon[iprim];
      outbasis<<"\t\t"<<0<<"\t"<<0<<"\t"<<2<<endl;
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
     }
     // 2 1 0 
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      outbasis<<setw(3)<<1<<"\t"<<setprecision(15)<<fixed<<scientific<<shell2aos[ishell].Coord[0]<<"\t"<<shell2aos[ishell].Coord[1];
      outbasis<<"\t"<<shell2aos[ishell].Coord[2]<<"\t"<<shell2aos[ishell].Expon[iprim];
      outbasis<<"\t\t"<<2<<"\t"<<1<<"\t"<<0<<endl;
     }
     // 2 0 1 
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      outbasis<<setw(3)<<1<<"\t"<<setprecision(15)<<fixed<<scientific<<shell2aos[ishell].Coord[0]<<"\t"<<shell2aos[ishell].Coord[1];
      outbasis<<"\t"<<shell2aos[ishell].Coord[2]<<"\t"<<shell2aos[ishell].Expon[iprim];
      outbasis<<"\t\t"<<2<<"\t"<<0<<"\t"<<1<<endl;
     }
     // 1 2 0 
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      outbasis<<setw(3)<<1<<"\t"<<setprecision(15)<<fixed<<scientific<<shell2aos[ishell].Coord[0]<<"\t"<<shell2aos[ishell].Coord[1];
      outbasis<<"\t"<<shell2aos[ishell].Coord[2]<<"\t"<<shell2aos[ishell].Expon[iprim];
      outbasis<<"\t\t"<<1<<"\t"<<2<<"\t"<<0<<endl;
     }
     // 1 1 1 
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      outbasis<<setw(3)<<1<<"\t"<<setprecision(15)<<fixed<<scientific<<shell2aos[ishell].Coord[0]<<"\t"<<shell2aos[ishell].Coord[1];
      outbasis<<"\t"<<shell2aos[ishell].Coord[2]<<"\t"<<shell2aos[ishell].Expon[iprim];
      outbasis<<"\t\t"<<1<<"\t"<<1<<"\t"<<1<<endl;
     }
     // 1 0 2 
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      outbasis<<setw(3)<<1<<"\t"<<setprecision(15)<<fixed<<scientific<<shell2aos[ishell].Coord[0]<<"\t"<<shell2aos[ishell].Coord[1];
      outbasis<<"\t"<<shell2aos[ishell].Coord[2]<<"\t"<<shell2aos[ishell].Expon[iprim];
      outbasis<<"\t\t"<<1<<"\t"<<0<<"\t"<<2<<endl;
     }
     // 0 3 0 
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      outbasis<<setw(3)<<1<<"\t"<<setprecision(15)<<fixed<<scientific<<shell2aos[ishell].Coord[0]<<"\t"<<shell2aos[ishell].Coord[1];
      outbasis<<"\t"<<shell2aos[ishell].Coord[2]<<"\t"<<shell2aos[ishell].Expon[iprim];
      outbasis<<"\t\t"<<0<<"\t"<<3<<"\t"<<0<<endl;
     }
     // 0 2 1 
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      outbasis<<setw(3)<<1<<"\t"<<setprecision(15)<<fixed<<scientific<<shell2aos[ishell].Coord[0]<<"\t"<<shell2aos[ishell].Coord[1];
      outbasis<<"\t"<<shell2aos[ishell].Coord[2]<<"\t"<<shell2aos[ishell].Expon[iprim];
      outbasis<<"\t\t"<<0<<"\t"<<2<<"\t"<<1<<endl;
     }
     // 0 1 2 
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      outbasis<<setw(3)<<1<<"\t"<<setprecision(15)<<fixed<<scientific<<shell2aos[ishell].Coord[0]<<"\t"<<shell2aos[ishell].Coord[1];
      outbasis<<"\t"<<shell2aos[ishell].Coord[2]<<"\t"<<shell2aos[ishell].Expon[iprim];
      outbasis<<"\t\t"<<0<<"\t"<<1<<"\t"<<2<<endl;
     }
     // 0 0 3 
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      outbasis<<setw(3)<<1<<"\t"<<setprecision(15)<<fixed<<scientific<<shell2aos[ishell].Coord[0]<<"\t"<<shell2aos[ishell].Coord[1];
      outbasis<<"\t"<<shell2aos[ishell].Coord[2]<<"\t"<<shell2aos[ishell].Expon[iprim];
      outbasis<<"\t\t"<<0<<"\t"<<0<<"\t"<<3<<endl;
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
     }
     // 3 1 0
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      outbasis<<setw(3)<<1<<"\t"<<setprecision(15)<<fixed<<scientific<<shell2aos[ishell].Coord[0]<<"\t"<<shell2aos[ishell].Coord[1];
      outbasis<<"\t"<<shell2aos[ishell].Coord[2]<<"\t"<<shell2aos[ishell].Expon[iprim];
      outbasis<<"\t\t"<<3<<"\t"<<1<<"\t"<<0<<endl;
     }
     // 3 0 1
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      outbasis<<setw(3)<<1<<"\t"<<setprecision(15)<<fixed<<scientific<<shell2aos[ishell].Coord[0]<<"\t"<<shell2aos[ishell].Coord[1];
      outbasis<<"\t"<<shell2aos[ishell].Coord[2]<<"\t"<<shell2aos[ishell].Expon[iprim];
      outbasis<<"\t\t"<<3<<"\t"<<0<<"\t"<<1<<endl;
     }
     // 2 2 0
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      outbasis<<setw(3)<<1<<"\t"<<setprecision(15)<<fixed<<scientific<<shell2aos[ishell].Coord[0]<<"\t"<<shell2aos[ishell].Coord[1];
      outbasis<<"\t"<<shell2aos[ishell].Coord[2]<<"\t"<<shell2aos[ishell].Expon[iprim];
      outbasis<<"\t\t"<<2<<"\t"<<2<<"\t"<<0<<endl;
     }
     // 2 1 1
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      outbasis<<setw(3)<<1<<"\t"<<setprecision(15)<<fixed<<scientific<<shell2aos[ishell].Coord[0]<<"\t"<<shell2aos[ishell].Coord[1];
      outbasis<<"\t"<<shell2aos[ishell].Coord[2]<<"\t"<<shell2aos[ishell].Expon[iprim];
      outbasis<<"\t\t"<<2<<"\t"<<1<<"\t"<<1<<endl;
     }
     // 2 0 2
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      outbasis<<setw(3)<<1<<"\t"<<setprecision(15)<<fixed<<scientific<<shell2aos[ishell].Coord[0]<<"\t"<<shell2aos[ishell].Coord[1];
      outbasis<<"\t"<<shell2aos[ishell].Coord[2]<<"\t"<<shell2aos[ishell].Expon[iprim];
      outbasis<<"\t\t"<<2<<"\t"<<0<<"\t"<<2<<endl;
     }
     // 1 3 0
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      outbasis<<setw(3)<<1<<"\t"<<setprecision(15)<<fixed<<scientific<<shell2aos[ishell].Coord[0]<<"\t"<<shell2aos[ishell].Coord[1];
      outbasis<<"\t"<<shell2aos[ishell].Coord[2]<<"\t"<<shell2aos[ishell].Expon[iprim];
      outbasis<<"\t\t"<<1<<"\t"<<3<<"\t"<<0<<endl;
     }
     // 1 2 1
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      outbasis<<setw(3)<<1<<"\t"<<setprecision(15)<<fixed<<scientific<<shell2aos[ishell].Coord[0]<<"\t"<<shell2aos[ishell].Coord[1];
      outbasis<<"\t"<<shell2aos[ishell].Coord[2]<<"\t"<<shell2aos[ishell].Expon[iprim];
      outbasis<<"\t\t"<<1<<"\t"<<2<<"\t"<<1<<endl;
     }
     // 1 1 2
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      outbasis<<setw(3)<<1<<"\t"<<setprecision(15)<<fixed<<scientific<<shell2aos[ishell].Coord[0]<<"\t"<<shell2aos[ishell].Coord[1];
      outbasis<<"\t"<<shell2aos[ishell].Coord[2]<<"\t"<<shell2aos[ishell].Expon[iprim];
      outbasis<<"\t\t"<<1<<"\t"<<1<<"\t"<<2<<endl;
     }
     // 1 0 3
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      outbasis<<setw(3)<<1<<"\t"<<setprecision(15)<<fixed<<scientific<<shell2aos[ishell].Coord[0]<<"\t"<<shell2aos[ishell].Coord[1];
      outbasis<<"\t"<<shell2aos[ishell].Coord[2]<<"\t"<<shell2aos[ishell].Expon[iprim];
      outbasis<<"\t\t"<<1<<"\t"<<0<<"\t"<<3<<endl;
     }
     // 0 4 0
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      outbasis<<setw(3)<<1<<"\t"<<setprecision(15)<<fixed<<scientific<<shell2aos[ishell].Coord[0]<<"\t"<<shell2aos[ishell].Coord[1];
      outbasis<<"\t"<<shell2aos[ishell].Coord[2]<<"\t"<<shell2aos[ishell].Expon[iprim];
      outbasis<<"\t\t"<<0<<"\t"<<4<<"\t"<<0<<endl;
     }
     // 0 3 1
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      outbasis<<setw(3)<<1<<"\t"<<setprecision(15)<<fixed<<scientific<<shell2aos[ishell].Coord[0]<<"\t"<<shell2aos[ishell].Coord[1];
      outbasis<<"\t"<<shell2aos[ishell].Coord[2]<<"\t"<<shell2aos[ishell].Expon[iprim];
      outbasis<<"\t\t"<<0<<"\t"<<3<<"\t"<<1<<endl;
     }
     // 0 2 2
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      outbasis<<setw(3)<<1<<"\t"<<setprecision(15)<<fixed<<scientific<<shell2aos[ishell].Coord[0]<<"\t"<<shell2aos[ishell].Coord[1];
      outbasis<<"\t"<<shell2aos[ishell].Coord[2]<<"\t"<<shell2aos[ishell].Expon[iprim];
      outbasis<<"\t\t"<<0<<"\t"<<2<<"\t"<<2<<endl;
     }
     // 0 1 3
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      outbasis<<setw(3)<<1<<"\t"<<setprecision(15)<<fixed<<scientific<<shell2aos[ishell].Coord[0]<<"\t"<<shell2aos[ishell].Coord[1];
      outbasis<<"\t"<<shell2aos[ishell].Coord[2]<<"\t"<<shell2aos[ishell].Expon[iprim];
      outbasis<<"\t\t"<<0<<"\t"<<1<<"\t"<<3<<endl;
     }
     // 0 0 4
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      outbasis<<setw(3)<<1<<"\t"<<setprecision(15)<<fixed<<scientific<<shell2aos[ishell].Coord[0]<<"\t"<<shell2aos[ishell].Coord[1];
      outbasis<<"\t"<<shell2aos[ishell].Coord[2]<<"\t"<<shell2aos[ishell].Expon[iprim];
      outbasis<<"\t\t"<<0<<"\t"<<0<<"\t"<<4<<endl;
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

void print_wfx()
{
 string line;
 ofstream real_wfx((dirac_output_name+"RE.wfx").c_str());
 ofstream imag_wfx((dirac_output_name+"IM.wfx").c_str());
 line="<Number of Nuclei>";
 real_wfx<<line<<endl; 
 imag_wfx<<line<<endl; 

 imag_wfx.close(); 
 real_wfx.close(); 
}

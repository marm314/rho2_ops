#include<iostream>
#include<iomanip>
#include<stdlib.h>
#include<cstring>
#include<fstream>
#include<vector>
#include<sstream>
#include"Numbers.h"
#include"String_ops.h"

using namespace std;


void line_fill(int &number, string line_read);
void Quant_fill(int **Quant,int styp);
void write_basis(FILE *pFile,string name_dm2,int &nucleous,double Atom_coord[3],double Exp,int nx,int ny,int nz);
//Global variables:
int RECORD_DELIMITER_LENGTH=4,sys;
bool donofdm2=false; 
bool int8alldm2=false; 

int main(int argc, char *argv[])
{
 cout<<"##########################################################################"<<endl;
 cout<<"##########################################################################"<<endl;
 cout<<"# Conversion of DM2 from MO to Primitives (Program in C++)               #"<<endl;
 cout<<"# written by Dr. Mauricio Rodriguez M email: marm3.14@gmail.com          #"<<endl;
 cout<<"##########################################################################"<<endl;
 cout<<"##########################################################################"<<endl;
 cout<<endl;
 cout<<endl;
 cout<<"#*****************************************************************************#";
 cout<<endl;
 cout<<"# Copyright (C) 2016 Dr. Mauricio A. Rodriguez Mayorga                        #";
 cout<<endl;
 cout<<"# Ph.D. student at Donostia International Physics Center (DIPC)               #";
 cout<<endl;
 cout<<"# for support and comments send an email to: marm314@gmail.com                #";
 cout<<endl;
 cout<<"#*****************************************************************************#";
 cout<<endl;
 cout<<"#  This program is free software: you can redistribute it and/or modify       #";
 cout<<endl;
 cout<<"#  it under the terms of the GNU General Public License as published by       #";
 cout<<endl;
 cout<<"#  the Free Software Foundation, either version 3 of the License, or          #";
 cout<<endl;
 cout<<"#  (at your option) any later version.                                        #";
 cout<<endl;
 cout<<"#  This program is distributed in the hope that it will be useful,            #";
 cout<<endl;
 cout<<"#  but WITHOUT ANY WARRANTY; without even the implied warranty of             #";
 cout<<endl;
 cout<<"#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              #";
 cout<<endl;
 cout<<"#  GNU General Public License for more details.                               #";
 cout<<endl;
 cout<<"#  You should have received a copy of the GNU General Public License          #";
 cout<<endl;
 cout<<"#  along with this program.  If not, see <http://www.gnu.org/licenses/>.      #";
 cout<<endl;
 cout<<"#*****************************************************************************#";
 cout<<endl;
 if(argc==5 || argc==6 || argc==7)
 {
  bool contractSP=false,use_dm2_iiii_aa=false,two_dm2_mat=false,reduce=false,all_dm2_are_given=false,red_sym=false;
  int i,j,k,l,AO,nbasis,nbasis2,nprim_shell,stype,contr_coef,spcontr_coef,iprim,iprim_tot,smap,prim_exp,natoms;
  int *shelltype,*nprim_per_shell,*shell_to_atom,*Atomic_Z,nucleous,**Quant;
  int pivot,element[2]={10},element_prime[2]={10};
  long int ii,jj,kk,ll,new_index,N2,N3,N4,N5,elementL[2]={10},element_primeL[2]={10},Ired,Jred,Kred,Lred;
  double threshold,Dijkl,Dijkl_change,trace=ZERO,*contrac_coeff,*SPcontrac_coeff,*Exponents,**MO_COEF,**L,**B,**Cartes_coord;
  double *C_coef_send, *SPC_coef_send,*Exponents_send,Atom_coord[3];
  string name_dm2(argv[1]);
  string name_fchk(argv[2]);
  threshold=atof(argv[3]);
  string y_or_n(argv[4]);
  if(y_or_n=="y" || y_or_n=="Y")
  {
   two_dm2_mat=true;
  }
  if(argc==6)
  {
   string y_or_n2(argv[5]);
   lowercase(y_or_n2);
   if(y_or_n2=="i4")
   {
    use_dm2_iiii_aa=true;
   }
   if(y_or_n2=="r")
   {
    reduce=true;
   }
   if(y_or_n2=="rs")
   {
    reduce=true;
    red_sym=true;
   }
   if(y_or_n2=="i4r")
   {
    use_dm2_iiii_aa=true;
    reduce=true;
   }
   if(y_or_n2=="i4rs")
   {
    use_dm2_iiii_aa=true;
    reduce=true;
    red_sym=true;
   }
   if(y_or_n2=="a")
   {
    all_dm2_are_given=true;
   }
   if(y_or_n2=="ar")
   {
    all_dm2_are_given=true;
    reduce=true;
   }
   if(y_or_n2=="ars")
   {
    all_dm2_are_given=true;
    reduce=true;
    red_sym=true;
   }
   if(y_or_n2=="donof")
   {
    donofdm2=true;
   }
   if(y_or_n2=="donofr")
   {
    donofdm2=true;
    reduce=true;
   }
   if(y_or_n2=="donofrs")
   {
    donofdm2=true;
    reduce=true;
    red_sym=true;
   }
   if(y_or_n2=="i8all")
   {
    int8alldm2=true; 
   }
   if(y_or_n2=="i8allr")
   {
    int8alldm2=true; 
    reduce=true;
   }
   if(y_or_n2=="i8allrs")
   {
    int8alldm2=true; 
    reduce=true;
    red_sym=true;
   }
  }
  string line;
  ifstream open_fchk;
  ifstream check_dm2;
  open_fchk.open((name_fchk).c_str());
  if(open_fchk.good())
  {
   cout<<endl;
   cout<<"Computing the product of coefficients matrix."<<endl;
   cout<<"[Coefficients for passing from MO -> Primitives (B matrix in the notes)]"<<endl;
   cout<<endl;
   while(getline(open_fchk,line))
   {
    if(line.substr(0,15)=="Number of atoms")
    {
     line_fill(natoms,line);
    }
    else if(line.substr(0,25)=="Number of basis functions")
    {
     line_fill(nbasis,line);
    }
    else if(line.substr(0,14)=="Atomic numbers")
    {
     Atomic_Z=new int[natoms];
     for(i=0;i<natoms;i++)
     {
      open_fchk>>Atomic_Z[i];
     }
    }
    else if(line.substr(0,29)=="Current cartesian coordinates")
    {
     Cartes_coord=new double*[natoms];
     for(i=0;i<natoms;i++)
     {
      Cartes_coord[i]=new double[3];
     }
     for(i=0;i<natoms;i++)
     {
      for(j=0;j<3;j++)
      {
       open_fchk>>Cartes_coord[i][j];
      }
     }
    }
    else if(line.substr(0,11)=="Shell types")
    {
     line_fill(stype,line);
     shelltype=new int[stype];
     for(i=0;i<stype;i++)
     {
      open_fchk>>shelltype[i];
     }
    }
    else if(line.substr(0,30)=="Number of primitives per shell")
    {
     line_fill(nprim_shell,line);
     nprim_per_shell=new int[nprim_shell];
     for(i=0;i<nprim_shell;i++)
     {
      open_fchk>>nprim_per_shell[i];
     }
    }
    else if(line.substr(0,17)=="Shell to atom map")
    {
     line_fill(smap,line);
     shell_to_atom=new int[smap];
     for(i=0;i<smap;i++)
     {
      open_fchk>>shell_to_atom[i];
     }
    }
    else if(line.substr(0,19)=="Primitive exponents")
    {
     line_fill(prim_exp,line);
     Exponents=new double[prim_exp];
     for(i=0;i<prim_exp;i++)
     {
      open_fchk>>Exponents[i];
     }
    }
    else if(line.substr(0,24)=="Contraction coefficients")
    {
     line_fill(contr_coef,line);
     contrac_coeff=new double[contr_coef];
     for(i=0;i<contr_coef;i++)
     {
      open_fchk>>contrac_coeff[i];
     }
    }
    else if(line.substr(0,31)=="P(S=P) Contraction coefficients")
    {
     line_fill(spcontr_coef,line);
     contractSP=true;
     SPcontrac_coeff=new double[spcontr_coef];
     for(i=0;i<spcontr_coef;i++)
     {
      open_fchk>>SPcontrac_coeff[i];
     }
    }
    else if(line.substr(0,21)=="Alpha MO coefficients")
    {
     MO_COEF=new double*[nbasis];
     for(i=0;i<nbasis;i++)
     {
      MO_COEF[i]=new double[nbasis];
     }
     for(i=0;i<nbasis;i++)
     {
      for(j=0;j<nbasis;j++)
      {
       open_fchk>>MO_COEF[i][j];
      }
     }
    }
    else if(line.substr(0,20)=="Beta MO coefficients" || line.substr(0,20)=="B3TA MO coefficients")
    {cout<<"Warning! 'Beta MO coefficient' found. Check only for HF calculations that UHF is not being used!"<<endl;}
    else{}
   }
   open_fchk.close();
   //Compute the size of the basis in primitives (nbasis2):
   nbasis2=0;
   for(i=0;i<stype;i++)
   {
    if((shelltype[i]==0 || shelltype[i]==1) || (shelltype[i]==2 ||shelltype[i]==3))
    {
     nbasis2=nbasis2+(shelltype[i]+1)*(shelltype[i]+2)*nprim_per_shell[i]/2;
    }
    if(shelltype[i]==-1)
    {
     nbasis2=nbasis2+4*nprim_per_shell[i]; // One s plus three p
    }
    if(shelltype[i]==-2)
    {
     nbasis2=nbasis2+6*nprim_per_shell[i]; // Build using the 6D
    }
    if(shelltype[i]==-3)
    {
     nbasis2=nbasis2+10*nprim_per_shell[i]; // Build using the 10F
    }
    if(shelltype[i]==-4)
    {
     nbasis2=nbasis2+15*nprim_per_shell[i]; // Build using the 15G
    }
   }
   cout<<"The size of the basis in cartesian primitives is: "<<nbasis2<<endl;
   L=new double*[nbasis];
   for(i=0;i<nbasis;i++)
   {
    L[i]=new double[nbasis2];
   }
   for(i=0;i<nbasis;i++)
   {
    for(j=0;j<nbasis2;j++)
    {
     L[i][j]=ZERO;
    }
   }
   //Create basis file
   FILE *pFile;
   pFile=fopen((name_dm2.substr(0,name_dm2.length()-4)+".basis").c_str(),"w");
   fputs("# Z  \t\t X_A \t\t Y_A \t\t Z_A \t\t Exp \t\t Quant(nx,ny,nz)\n",pFile);
   for(i=0;i<stype;i++)
   {
    nucleous=Atomic_Z[shell_to_atom[i]-1];
    Atom_coord[0]=Cartes_coord[shell_to_atom[i]-1][0];
    Atom_coord[1]=Cartes_coord[shell_to_atom[i]-1][1];
    Atom_coord[2]=Cartes_coord[shell_to_atom[i]-1][2];
    iprim=nprim_per_shell[i];
    iprim_tot=0;
    for(j=0;j<i;j++)
    {iprim_tot=iprim_tot+nprim_per_shell[j];}
    Exponents_send=new double[iprim];
    k=0;
    for(j=iprim_tot;j<iprim+iprim_tot;j++)
    {
     Exponents_send[k]=Exponents[j];
     k++;
    }
    Quant=new int*[(abs(shelltype[i])+2)*(abs(shelltype[i])+1)/2];
    for(j=0;j<(abs(shelltype[i])+2)*(abs(shelltype[i])+1)/2;j++)
    {
     Quant[j]=new int[3];
    }
    Quant_fill(Quant,abs(shelltype[i]));
    if((shelltype[i]==0 || shelltype[i]==1) || (shelltype[i]==2 || shelltype[i]==3) || shelltype[i]==4)
    {
     for(j=0;j<(shelltype[i]+2)*(shelltype[i]+1)/2;j++)
     {
      for(k=0;k<iprim;k++)
      {
       write_basis(pFile,name_dm2,nucleous,Atom_coord,Exponents_send[k],Quant[j][0],Quant[j][1],Quant[j][2]);
      }
     }
    }
    else if(shelltype[i]==-1)
    {
     //s Gaussian primitives
     for(k=0;k<iprim;k++)
     {
      write_basis(pFile,name_dm2,nucleous,Atom_coord,Exponents_send[k],0,0,0);
     }
     //p Gaussian primitives
     for(j=0;j<(abs(shelltype[i])+2)*(abs(shelltype[i])+1)/2;j++)
     {
      for(k=0;k<iprim;k++)
      {
       write_basis(pFile,name_dm2,nucleous,Atom_coord,Exponents_send[k],Quant[j][0],Quant[j][1],Quant[j][2]);
      }
     }
    }
    else if(shelltype[i]==-2 || shelltype[i]==-3 || shelltype[i]==-4)
    {
     for(k=0;k<iprim;k++)
     {
      for(j=0;j<(abs(shelltype[i])+2)*(abs(shelltype[i])+1)/2;j++)
      {
       write_basis(pFile,name_dm2,nucleous,Atom_coord,Exponents_send[k],Quant[j][0],Quant[j][1],Quant[j][2]);
      }
     }
    }
    else{}
    for(j=0;j<(abs(shelltype[i])+2)*(abs(shelltype[i])+1)/2;j++)
    {
     delete[] Quant[j];Quant[j]=NULL;
    }
    delete[] Quant;Quant=NULL;
    delete[] Exponents_send;Exponents_send=NULL;
   }
   fclose(pFile);
   //Create L matrix
   AO=0;
   l=0;
   for(i=0;i<stype;i++)
   {
    iprim=nprim_per_shell[i];
    iprim_tot=0;
    for(j=0;j<i;j++)
    {iprim_tot=iprim_tot+nprim_per_shell[j];}
    C_coef_send=new double[iprim];
    if(contractSP)
    {
     SPC_coef_send=new double[iprim];
    }
    k=0;
    for(j=iprim_tot;j<iprim+iprim_tot;j++)
    {
     C_coef_send[k]=contrac_coeff[j];
     if(contractSP)
     {SPC_coef_send[k]=SPcontrac_coeff[j];}
     k++;
    }
    if((shelltype[i]==0 || shelltype[i]==1) || (shelltype[i]==2 || shelltype[i]==3) || shelltype[i]==4)
    {
     for(j=0;j<(shelltype[i]+2)*(shelltype[i]+1)/2;j++)
     {
      for(k=0;k<iprim;k++)
      {
       L[AO][l]=C_coef_send[k];
       l++;
      }
      AO++;
     }
    }
    if(shelltype[i]==-1)
    {
     //The s orbital
     for(k=0;k<iprim;k++)
     {
      L[AO][l]=C_coef_send[k];
      l++;
     }
     AO++;
     //The p orbitals
     for(j=0;j<3;j++)
     {
      for(k=0;k<iprim;k++)
      {
       L[AO][l]=SPC_coef_send[k];
       l++;
      }
      AO++;
     }
    }
   //See H. Bernhard Schlegel and Michael J. Frisch
   //"Transformation Between Cartesian and Pure Spherical Harmonic Gaussians"
   //Int. J. Quant. Chem., 54, 83-87 (1995).
    if(shelltype[i]==-2)
    {
     for(k=0;k<iprim;k++)
     {
      //(2,0)
      L[AO][l]=-HALF*C_coef_send[k];L[AO][l+1]=-HALF*C_coef_send[k];L[AO][l+2]=C_coef_send[k];
      //2^(-1/2){(2,1)+(2,-1)}
      L[AO+1][l+4]=C_coef_send[k];
      //(-2)^(-1/2){(2,1)-(2,-1)}
      L[AO+2][l+5]=C_coef_send[k];
      //2^(-1/2){(2,2)+(2,-2)}
      L[AO+3][l]=pow(THREE/FOUR,HALF)*C_coef_send[k];L[AO+3][l+1]=-pow(THREE/FOUR,HALF)*C_coef_send[k];
      //(-2)^(-1/2){(2,2)-(2,-2)}
      L[AO+4][l+3]=C_coef_send[k];
      l=l+6;
     }
     AO=AO+5;
    }
    if(shelltype[i]==-3)
    {
     for(k=0;k<iprim;k++)
     {
      //(3,0)
      L[AO][l+2]=C_coef_send[k];L[AO][l+4]=-(THREE/(TWO*pow(FIVE,HALF)))*C_coef_send[k];L[AO][l+5]=-(THREE/(TWO*pow(FIVE,HALF)))*C_coef_send[k];
      //2^(-1/2){(3,1)+(3,-1)}
      L[AO+1][l+7]=pow(SIX/FIVE,HALF)*C_coef_send[k];L[AO+1][l]=-(pow(SIX,HALF)/FOUR)*C_coef_send[k];
      L[AO+1][l+6]=-(pow(SIX/FIVE,HALF)/FOUR)*C_coef_send[k];
      //(-2)^(-1/2){(3,1)-(3,-1)}
      L[AO+2][l+8]=pow(SIX/FIVE,HALF)*C_coef_send[k];L[AO+2][l+1]=-(pow(SIX,HALF)/FOUR)*C_coef_send[k];
      L[AO+2][l+3]=-(pow(SIX/FIVE,HALF)/FOUR)*C_coef_send[k];
      //2^(-1/2){(3,2)+(3,-2)}
      L[AO+3][l+4]=pow(SIX/EIGHT,HALF)*C_coef_send[k];L[AO+3][l+5]=-pow(SIX/EIGHT,HALF)*C_coef_send[k];
      //(-2)^(-1/2){(3,2)-(3,-2)}
      L[AO+4][l+9]=C_coef_send[k];
      //2^(-1/2){(3,3)+(3,-3)}
      L[AO+5][l]=(pow(TEN,HALF)/FOUR)*C_coef_send[k];L[AO+5][l+6]=-(THREE/(TWO*pow(TWO,HALF)))*C_coef_send[k];
      //(-2)^(-1/2){(3,3)-(3,-3)}
      L[AO+6][l+1]=-(pow(TEN,HALF)/FOUR)*C_coef_send[k];L[AO+6][l+3]=(THREE/(TWO*pow(TWO,HALF)))*C_coef_send[k];
      l=l+10;
     }
     AO=AO+7;
    }
    if(shelltype[i]==-4)
    {
     for(k=0;k<iprim;k++)
     {
      //(4,0)
      L[AO][l+2]=C_coef_send[k];
      L[AO][l]=(THREE/EIGHT)*C_coef_send[k];
      L[AO][l+1]=(THREE/EIGHT)*C_coef_send[k];
      L[AO][l+10]=-(THREE*pow(THREE,HALF)/pow(SEVEN*FIVE,HALF))*C_coef_send[k];
      L[AO][l+11]=-(THREE*pow(THREE,HALF)/pow(SEVEN*FIVE,HALF))*C_coef_send[k];
      L[AO][l+9]=(THREE*pow(THREE,HALF)/pow(SEVEN*FIVE,HALF))*HALF*HALF*C_coef_send[k];
      //2^(-1/2){(4,1)+(4,-1)}
      L[AO+1][l+7]=pow(TEN/SEVEN,HALF)*C_coef_send[k];
      L[AO+1][l+4]=-(THREE*pow(TEN/SEVEN,HALF)/FOUR)*C_coef_send[k];
      L[AO+1][l+13]=-(THREE*pow(TWO/SEVEN,HALF)/FOUR)*C_coef_send[k];
      //(-2)^(-1/2){(4,1)-(4,-1)}
      L[AO+2][l+8]=pow(TEN/SEVEN,HALF)*C_coef_send[k];
      L[AO+2][l+6]=-(THREE*pow(TEN/SEVEN,HALF)/FOUR)*C_coef_send[k];
      L[AO+2][l+12]=-(THREE*pow(TWO/SEVEN,HALF)/FOUR)*C_coef_send[k];
      //2^(-1/2){(4,2)+(4,-2)}
      L[AO+3][l+10]=(THREE*pow(THREE/SEVEN,HALF)/TWO)*C_coef_send[k];
      L[AO+3][l+11]=-(THREE*pow(THREE/SEVEN,HALF)/TWO)*C_coef_send[k];
      L[AO+3][l]=-HALF*HALF*pow(FIVE,HALF)*C_coef_send[k];
      L[AO+3][l+1]=HALF*HALF*pow(FIVE,HALF)*C_coef_send[k];
      //(-2)^(-1/2){(4,2)-(4,-2)}
      L[AO+4][l+14]=(THREE/pow(SEVEN,HALF))*C_coef_send[k];
      L[AO+4][l+3]=-HALF*pow(FIVE/SEVEN,HALF)*C_coef_send[k];
      L[AO+4][l+5]=-HALF*pow(FIVE/SEVEN,HALF)*C_coef_send[k];
      //2^(-1/2){(4,3)+(4,-3)}
      L[AO+5][l+4]=(pow(TEN,HALF)/FOUR)*C_coef_send[k];
      L[AO+5][l+13]=-THREE*(pow(TWO,HALF)/FOUR)*C_coef_send[k];
      //(-2)^(-1/2){(4,3)-(4,-3)}
      L[AO+6][l+6]=(-pow(TEN,HALF)/FOUR)*C_coef_send[k];
      L[AO+6][l+12]=THREE*(pow(TWO,HALF)/FOUR)*C_coef_send[k];
      //2^(-1/2){(4,4)+(4,-4)}
      L[AO+7][l]=(pow(SEVEN*FIVE,HALF)/EIGHT)*C_coef_send[k];
      L[AO+7][l+1]=(pow(SEVEN*FIVE,HALF)/EIGHT)*C_coef_send[k];
      L[AO+7][l+9]=-(THREE*pow(THREE,HALF)/FOUR)*C_coef_send[k];
      //(-2)^(-1/2){(4,4)-(4,-4)}
      L[AO+8][l+3]=pow(FIVE/FOUR,HALF)*C_coef_send[k];
      L[AO+8][l+5]=-pow(FIVE/FOUR,HALF)*C_coef_send[k];
      l=l+15;
     }
     AO=AO+9;
    }
    delete[] C_coef_send;
    C_coef_send=NULL;
    if(contractSP)
    {delete[] SPC_coef_send;SPC_coef_send=NULL;}
   }
   B=new double*[nbasis];
   for(i=0;i<nbasis;i++)
   {
    B[i]=new double[nbasis2];
   }
   for(i=0;i<nbasis;i++)
   {
    for(j=0;j<nbasis2;j++)
    {
     B[i][j]=ZERO;
    }
   }
   //B = MO_COEF x  L
   for(i=0;i<nbasis;i++)
   {
    for(j=0;j<nbasis2;j++)
    {
     for(k=0;k<nbasis;k++)
     {
      B[i][j]=B[i][j]+MO_COEF[i][k]*L[k][j];
     }
    }
   }
   for(i=0;i<nbasis;i++)
   {
    delete[] L[i];L[i]=NULL;
   }
   delete[] L; L=NULL;
   for(i=0;i<natoms;i++)
   {
    delete[] Cartes_coord[i];Cartes_coord[i]=NULL;
   }
   delete[] Cartes_coord;Cartes_coord=NULL;
   delete[] Atomic_Z;Atomic_Z=NULL;
   delete[] shelltype;shelltype=NULL;
   delete[] nprim_per_shell;nprim_per_shell=NULL;
   delete[] shell_to_atom;shell_to_atom=NULL;
   delete[] contrac_coeff; contrac_coeff=NULL;
   delete[] Exponents;Exponents=NULL;
   if(contractSP)
   {
    delete[] SPcontrac_coeff;SPcontrac_coeff=NULL;
   }
   for(i=0;i<nbasis;i++)
   {
    delete[] MO_COEF[i]; MO_COEF[i]=NULL;
   }
   delete[] MO_COEF; MO_COEF=NULL;
   check_dm2.open((name_dm2).c_str());
   if(check_dm2.good())
   {
    check_dm2.close();
    cout<<endl;
    N2=nbasis2;
    N3=N2*N2;
    N4=N3*N2;
    N5=N4*N2;
    vector <double> Dijkl_term(N5,ZERO);
    cout<<"Storing of the spinless 2RDM..."<<endl;
    {
     ifstream input_data(name_dm2.c_str(),ios::binary);
     if(!donofdm2 && !int8alldm2)
     {
      while(element[0]!=0 || element[1]!=0 || element_prime[0]!=0 || element_prime[1]!=0)
      {
       input_data.seekg(RECORD_DELIMITER_LENGTH, ios::cur);
       input_data.read((char*) &element[0], sizeof(element[0]));
       input_data.read((char*) &element[1], sizeof(element[1]));
       input_data.read((char*) &element_prime[0], sizeof(element_prime[0]));
       input_data.read((char*) &element_prime[1], sizeof(element_prime[1]));
       input_data.read((char*) &Dijkl, sizeof(Dijkl));
       input_data.seekg(RECORD_DELIMITER_LENGTH, ios::cur);
       if(Dijkl!=ZERO && element[0]!=0 && element[1]!=0 && element_prime[0]!=0 && element_prime[1]!=0)
       {
        if(!all_dm2_are_given)
        {
         if(!use_dm2_iiii_aa)
         {
          if(element[0]%2!=element_prime[0]%2 && element[0]%2==element_prime[1]%2)
          {
           pivot=element_prime[0];
           element_prime[0]=element_prime[1];
           element_prime[1]=pivot;
           Dijkl=-Dijkl;
          }
          if(element[0]%2==element_prime[0]%2 && element[1]%2==element_prime[1]%2)
          {
           if(element[0]%2==element[1]%2)
           {
            if(element[0]%2==0)
            {
             ii=(element[0]/2)-1;
             jj=(element[1]/2)-1;
             kk=(element_prime[0]/2)-1;
             ll=(element_prime[1]/2)-1;
            }
            else
            {
             ii=(element[0]-1)/2;
             jj=(element[1]-1)/2;
             kk=(element_prime[0]-1)/2;
             ll=(element_prime[1]-1)/2;
            }
            Dijkl_term[ii+jj*N2+kk*N3+ll*N4]=Dijkl_term[ii+jj*N2+kk*N3+ll*N4]+Dijkl;
            Dijkl_term[jj+ii*N2+kk*N3+ll*N4]=Dijkl_term[jj+ii*N2+kk*N3+ll*N4]-Dijkl;
            Dijkl_term[jj+ii*N2+ll*N3+kk*N4]=Dijkl_term[jj+ii*N2+ll*N3+kk*N4]+Dijkl;
            Dijkl_term[ii+jj*N2+ll*N3+kk*N4]=Dijkl_term[ii+jj*N2+ll*N3+kk*N4]-Dijkl;
            if(ii!=kk || jj!=ll)
            {
             Dijkl_term[kk+ll*N2+ii*N3+jj*N4]=Dijkl_term[kk+ll*N2+ii*N3+jj*N4]+Dijkl;
             Dijkl_term[ll+kk*N2+ii*N3+jj*N4]=Dijkl_term[ll+kk*N2+ii*N3+jj*N4]-Dijkl;
             Dijkl_term[ll+kk*N2+jj*N3+ii*N4]=Dijkl_term[ll+kk*N2+jj*N3+ii*N4]+Dijkl;
             Dijkl_term[kk+ll*N2+jj*N3+ii*N4]=Dijkl_term[kk+ll*N2+jj*N3+ii*N4]-Dijkl;
            }
           }
           else
           {
            if(element[0]%2==0)
            {
             ii=(element[0]/2)-1;
             jj=(element[1]-1)/2;
             kk=(element_prime[0]/2)-1;
             ll=(element_prime[1]-1)/2;
            }
            else
            {
             ii=(element[0]-1)/2;
             jj=(element[1]/2)-1;
             kk=(element_prime[0]-1)/2;
             ll=(element_prime[1]/2)-1;
            }
            Dijkl_term[ii+jj*N2+kk*N3+ll*N4]=Dijkl_term[ii+jj*N2+kk*N3+ll*N4]+Dijkl;
            Dijkl_term[jj+ii*N2+ll*N3+kk*N4]=Dijkl_term[jj+ii*N2+ll*N3+kk*N4]+Dijkl;
            if(ii!=kk || jj!=ll)
            {
             Dijkl_term[kk+ll*N2+ii*N3+jj*N4]=Dijkl_term[kk+ll*N2+ii*N3+jj*N4]+Dijkl;
             Dijkl_term[ll+kk*N2+jj*N3+ii*N4]=Dijkl_term[ll+kk*N2+jj*N3+ii*N4]+Dijkl;
            }
           }
          }
         }
         else
         {
          if(element[0]%2!=element_prime[0]%2 && element[0]%2==element_prime[1]%2)
          {
           pivot=element_prime[0];
           element_prime[0]=element_prime[1];
           element_prime[1]=pivot;
           Dijkl=-Dijkl;
          }
          if(element[0]%2==element_prime[0]%2 && element[1]%2==element_prime[1]%2)
          {
           if(element[0]%2==element[1]%2)
           {
            if(element[0]%2==0)
            {
             ii=(element[0]/2)-1;
             jj=(element[1]/2)-1;
             kk=(element_prime[0]/2)-1;
             ll=(element_prime[1]/2)-1;
            }
            else
            {
             ii=(element[0]-1)/2;
             jj=(element[1]-1)/2;
             kk=(element_prime[0]-1)/2;
             ll=(element_prime[1]-1)/2;
            }
            Dijkl_term[ii+jj*N2+kk*N3+ll*N4]=Dijkl_term[ii+jj*N2+kk*N3+ll*N4]+Dijkl;
           }
           else
           {
            if(element[0]%2==0)
            {
             ii=(element[0]/2)-1;
             jj=(element[1]-1)/2;
             kk=(element_prime[0]/2)-1;
             ll=(element_prime[1]-1)/2;
            }
            else
            {
             ii=(element[0]-1)/2;
             jj=(element[1]/2)-1;
             kk=(element_prime[0]-1)/2;
             ll=(element_prime[1]/2)-1;
            }
            Dijkl_term[ii+jj*N2+kk*N3+ll*N4]=Dijkl_term[ii+jj*N2+kk*N3+ll*N4]+Dijkl;
           }
          }
         }
        }
        else
        {
         if(element[0]%2==0)
         {
          ii=(element[0]/2)-1;
         }
         else
         {
          ii=(element[0]-1)/2;
         }
         if(element[1]%2==0)
         {
          jj=(element[1]/2)-1;
         } 
         else
         {
          jj=(element[1]-1)/2;
         }
         if(element_prime[0]%2==0)
         {
          kk=(element_prime[0]/2)-1;
         }
         else
         {
          kk=(element_prime[0]-1)/2;
         }
         if(element_prime[1]%2==0)
         {
          ll=(element_prime[1])/2-1;
         }
         else
         {
          ll=(element_prime[1]-1)/2;
         }
         Dijkl_term[ii+jj*N2+kk*N3+ll*N4]=Dijkl_term[ii+jj*N2+kk*N3+ll*N4]+Dijkl;
        }
       }
      }
     }
     else
     {
      element[0]=10;
      while(element[0]!=0 || element[1]!=0 || element_prime[0]!=0 || element_prime[1]!=0)
      {
       input_data.seekg(RECORD_DELIMITER_LENGTH, ios::cur);
       input_data.read((char*) &element[0], sizeof(element[0]));
       input_data.read((char*) &i, sizeof(i));
       input_data.read((char*) &element[1], sizeof(element[1]));
       input_data.read((char*) &i, sizeof(i));
       input_data.read((char*) &element_prime[0], sizeof(element_prime[0]));
       input_data.read((char*) &i, sizeof(i));
       input_data.read((char*) &element_prime[1], sizeof(element_prime[1]));
       input_data.read((char*) &i, sizeof(i));
       input_data.read((char*) &Dijkl, sizeof(Dijkl));
       input_data.seekg(RECORD_DELIMITER_LENGTH, ios::cur);
       if(Dijkl!=ZERO && element[0]!=0 && element[1]!=0 && element_prime[0]!=0 && element_prime[1]!=0)
       {
        ii=element[0]-1;
        jj=element[1]-1;
        kk=element_prime[0]-1;
        ll=element_prime[1]-1;
        Dijkl_term[ii+jj*N2+kk*N3+ll*N4]=Dijkl_term[ii+jj*N2+kk*N3+ll*N4]+TWO*Dijkl;
       }
      }
     }
     input_data.close();
    }
///////////////////////////////////////////////////////////////////////////////
/// Print spinless DM2                                                       //
///////////////////////////////////////////////////////////////////////////////
    if(argc==7)
    {
     ofstream spilss_2rdm("sl-2RDM");
     spilss_2rdm<<setprecision(10)<<fixed<<scientific;
     for(ii=0;ii<N2;ii++)
     {
      for(jj=0;jj<N2;jj++)
      {   
       for(kk=0;kk<N2;kk++)
       {
        for(ll=0;ll<N2;ll++)
        {  
         if(abs(Dijkl_term[ii+jj*N2+kk*N3+ll*N4])>=threshold)
         {  
          spilss_2rdm<<setw(4)<<ii+1<<setw(4)<<jj+1<<setw(4)<<kk+1<<setw(4)<<ll+1<<setw(21)<<Dijkl_term[ii+jj*N2+kk*N3+ll*N4]<<endl; 
         }
        }
       }
      }
     }
     spilss_2rdm.close();
    }
//////////////////////////////////////////////////////////////////////////////// 
    for(ii=0;ii<N2;ii++)
    {
     for(jj=0;jj<N2;jj++)
     {
      if(abs(Dijkl_term[ii*(1+N3)+jj*(N2+N4)])>=threshold)
      {
       trace=trace+Dijkl_term[ii*(1+N3)+jj*(N2+N4)];
      }
     }
    }
    cout<<"Storing done!"<<endl;
    cout<<endl;
    cout<<"Trace of the 2RDM stored:  "<<setprecision(10)<<fixed<<scientific<<trace<<endl;
    cout<<endl;
    //Change basis
    ofstream temp_dm2;
    ifstream temp_dm2_2;
    cout<<"Change indexes: ^2 Dij,kl -> ^2 Dpq,rs"<<endl;
    cout<<endl;
    //Change l:
    if(!two_dm2_mat)
    {
     temp_dm2.open("Temp.dm2",ios::binary);
     for(ii=0;ii<N5;ii++)
     {
      Dijkl=Dijkl_term[ii];
      if(abs(Dijkl)>=threshold)
      {
       element_primeL[1]=ii/N4;
       for(jj=0;jj<N2;jj++)
       {
        Dijkl_change=Dijkl*B[element_primeL[1]][jj];
        if(abs(Dijkl_change)>=threshold)
        {
         new_index=ii-element_primeL[1]*N4+jj*N4;
         temp_dm2.seekp(RECORD_DELIMITER_LENGTH, ios::cur);
         temp_dm2.write((char*) &new_index, sizeof(new_index));
         temp_dm2.write((char*) &Dijkl_change, sizeof(Dijkl_change));
         temp_dm2.seekp(RECORD_DELIMITER_LENGTH, ios::cur);
        }
       }
      }
      Dijkl_term[ii]=ZERO;
     }
     new_index=-1;
     Dijkl_change=ZERO;
     temp_dm2.seekp(RECORD_DELIMITER_LENGTH, ios::cur);
     temp_dm2.write((char*) &new_index, sizeof(new_index));
     temp_dm2.write((char*) &Dijkl_change, sizeof(Dijkl_change));
     temp_dm2.seekp(RECORD_DELIMITER_LENGTH, ios::cur);
     temp_dm2.close();
     temp_dm2_2.open("Temp.dm2",ios::binary);
     temp_dm2_2.seekg(RECORD_DELIMITER_LENGTH, ios::cur);
     temp_dm2_2.read((char*) &new_index, sizeof(new_index));
     temp_dm2_2.read((char*) &Dijkl_change, sizeof(Dijkl_change));
     temp_dm2_2.seekg(RECORD_DELIMITER_LENGTH, ios::cur);
     while(new_index!=-1)
     {
      Dijkl_term[new_index]=Dijkl_term[new_index]+Dijkl_change;
      temp_dm2_2.seekg(RECORD_DELIMITER_LENGTH, ios::cur);
      temp_dm2_2.read((char*) &new_index, sizeof(new_index));
      temp_dm2_2.read((char*) &Dijkl_change, sizeof(Dijkl_change));
      temp_dm2_2.seekg(RECORD_DELIMITER_LENGTH, ios::cur);
     }
     temp_dm2_2.close();
     sys=system("/bin/rm Temp.dm2");
     cout<<"l coefficient changed"<<endl;
    }
    //Change k:
    if(!two_dm2_mat)
    {
     temp_dm2.open("Temp.dm2",ios::binary);
     for(ii=0;ii<N5;ii++)
     {
      Dijkl=Dijkl_term[ii];
      if(abs(Dijkl)>=threshold)
      {
       element_primeL[1]=ii/N4;
       element_primeL[0]=(ii-element_primeL[1]*N4)/N3;
       for(jj=0;jj<N2;jj++)
       {
        Dijkl_change=Dijkl*B[element_primeL[0]][jj];
        if(abs(Dijkl_change)>=threshold)
        {
         new_index=ii-element_primeL[0]*N3+jj*N3;
         temp_dm2.seekp(RECORD_DELIMITER_LENGTH, ios::cur);
         temp_dm2.write((char*) &new_index, sizeof(new_index));
         temp_dm2.write((char*) &Dijkl_change, sizeof(Dijkl_change));
         temp_dm2.seekp(RECORD_DELIMITER_LENGTH, ios::cur);
        }
       }
      }
      Dijkl_term[ii]=ZERO;
     }
     new_index=-1;
     Dijkl_change=ZERO;
     temp_dm2.seekp(RECORD_DELIMITER_LENGTH, ios::cur);
     temp_dm2.write((char*) &new_index, sizeof(new_index));
     temp_dm2.write((char*) &Dijkl_change, sizeof(Dijkl_change));
     temp_dm2.seekp(RECORD_DELIMITER_LENGTH, ios::cur);
     temp_dm2.close();
     temp_dm2_2.open("Temp.dm2",ios::binary);
     temp_dm2_2.seekg(RECORD_DELIMITER_LENGTH, ios::cur);
     temp_dm2_2.read((char*) &new_index, sizeof(new_index));
     temp_dm2_2.read((char*) &Dijkl_change, sizeof(Dijkl_change));
     temp_dm2_2.seekg(RECORD_DELIMITER_LENGTH, ios::cur);
     while(new_index!=-1)
     {
      Dijkl_term[new_index]=Dijkl_term[new_index]+Dijkl_change;
      temp_dm2_2.seekg(RECORD_DELIMITER_LENGTH, ios::cur);
      temp_dm2_2.read((char*) &new_index, sizeof(new_index));
      temp_dm2_2.read((char*) &Dijkl_change, sizeof(Dijkl_change));
      temp_dm2_2.seekg(RECORD_DELIMITER_LENGTH, ios::cur);
     }
     temp_dm2_2.close();
     sys=system("/bin/rm Temp.dm2");
     cout<<"k coefficient changed"<<endl;
    }
    //Change j:
    if(!two_dm2_mat)
    {
     temp_dm2.open("Temp.dm2",ios::binary);
     for(ii=0;ii<N5;ii++)
     {
      Dijkl=Dijkl_term[ii];
      if(abs(Dijkl)>=threshold)
      {
       element_primeL[1]=ii/N4;
       element_primeL[0]=(ii-element_primeL[1]*N4)/N3;
       elementL[1]=(ii-element_primeL[0]*N3-element_primeL[1]*N4)/N2;
       for(jj=0;jj<N2;jj++)
       {
        Dijkl_change=Dijkl*B[elementL[1]][jj];
        if(abs(Dijkl_change)>=threshold)
        {
         new_index=ii-elementL[1]*N2+jj*N2;
         temp_dm2.seekp(RECORD_DELIMITER_LENGTH, ios::cur);
         temp_dm2.write((char*) &new_index, sizeof(new_index));
         temp_dm2.write((char*) &Dijkl_change, sizeof(Dijkl_change));
         temp_dm2.seekp(RECORD_DELIMITER_LENGTH, ios::cur);
        }
       }
      }
      Dijkl_term[ii]=ZERO;
     }
     new_index=-1;
     Dijkl_change=ZERO;
     temp_dm2.seekp(RECORD_DELIMITER_LENGTH, ios::cur);
     temp_dm2.write((char*) &new_index, sizeof(new_index));
     temp_dm2.write((char*) &Dijkl_change, sizeof(Dijkl_change));
     temp_dm2.seekp(RECORD_DELIMITER_LENGTH, ios::cur);
     temp_dm2.close();
     temp_dm2_2.open("Temp.dm2",ios::binary);
     temp_dm2_2.seekg(RECORD_DELIMITER_LENGTH, ios::cur);
     temp_dm2_2.read((char*) &new_index, sizeof(new_index));
     temp_dm2_2.read((char*) &Dijkl_change, sizeof(Dijkl_change));
     temp_dm2_2.seekg(RECORD_DELIMITER_LENGTH, ios::cur);
     while(new_index!=-1)
     {
      Dijkl_term[new_index]=Dijkl_term[new_index]+Dijkl_change;
      temp_dm2_2.seekg(RECORD_DELIMITER_LENGTH, ios::cur);
      temp_dm2_2.read((char*) &new_index, sizeof(new_index));
      temp_dm2_2.read((char*) &Dijkl_change, sizeof(Dijkl_change));
      temp_dm2_2.seekg(RECORD_DELIMITER_LENGTH, ios::cur);
     }
     temp_dm2_2.close();
     sys=system("/bin/rm Temp.dm2");
     cout<<"j coefficient changed"<<endl;
    }
    //Change i:
    if(!two_dm2_mat)
    {
     temp_dm2.open("Temp.dm2",ios::binary);
     for(ii=0;ii<N5;ii++)
     {
      Dijkl=Dijkl_term[ii];
      if(abs(Dijkl)>=threshold)
      {
       element_primeL[1]=ii/N4;
       element_primeL[0]=(ii-element_primeL[1]*N4)/N3;
       elementL[1]=(ii-element_primeL[0]*N3-element_primeL[1]*N4)/N2;
       elementL[0]=(ii-elementL[1]*N2-element_primeL[0]*N3-element_primeL[1]*N4);
       for(jj=0;jj<N2;jj++)
       {
        Dijkl_change=Dijkl*B[elementL[0]][jj];
        if(abs(Dijkl_change)>=threshold)
        {
         new_index=ii-elementL[0]+jj;
         temp_dm2.seekp(RECORD_DELIMITER_LENGTH, ios::cur);
         temp_dm2.write((char*) &new_index, sizeof(new_index));
         temp_dm2.write((char*) &Dijkl_change, sizeof(Dijkl_change));
         temp_dm2.seekp(RECORD_DELIMITER_LENGTH, ios::cur);
        }
       }
      }
      Dijkl_term[ii]=ZERO;
     }
     new_index=-1;
     Dijkl_change=ZERO;
     temp_dm2.seekp(RECORD_DELIMITER_LENGTH, ios::cur);
     temp_dm2.write((char*) &new_index, sizeof(new_index));
     temp_dm2.write((char*) &Dijkl_change, sizeof(Dijkl_change));
     temp_dm2.seekp(RECORD_DELIMITER_LENGTH, ios::cur);
     temp_dm2.close();
     temp_dm2_2.open("Temp.dm2",ios::binary);
     temp_dm2_2.seekg(RECORD_DELIMITER_LENGTH, ios::cur);
     temp_dm2_2.read((char*) &new_index, sizeof(new_index));
     temp_dm2_2.read((char*) &Dijkl_change, sizeof(Dijkl_change));
     temp_dm2_2.seekg(RECORD_DELIMITER_LENGTH, ios::cur);
     while(new_index!=-1)
     {
      Dijkl_term[new_index]=Dijkl_term[new_index]+Dijkl_change;
      temp_dm2_2.seekg(RECORD_DELIMITER_LENGTH, ios::cur);
      temp_dm2_2.read((char*) &new_index, sizeof(new_index));
      temp_dm2_2.read((char*) &Dijkl_change, sizeof(Dijkl_change));
      temp_dm2_2.seekg(RECORD_DELIMITER_LENGTH, ios::cur);
     }
     temp_dm2_2.close();
     sys=system("/bin/rm Temp.dm2");
     cout<<"i coefficient changed"<<endl;
    }
    if(two_dm2_mat)
    {
     vector <double> Dijkl_term2(N5,ZERO);
     for(ii=0;ii<N5;ii++)
     {
      Dijkl=Dijkl_term[ii];
      if(abs(Dijkl)>=threshold)
      {
       element_primeL[1]=ii/N4;
       for(jj=0;jj<N2;jj++)
       {
        Dijkl_change=Dijkl*B[element_primeL[1]][jj];
        if(abs(Dijkl_change)>=threshold)
        {
         new_index=ii-element_primeL[1]*N4+jj*N4;
         Dijkl_term2[new_index]=Dijkl_term2[new_index]+Dijkl_change;
        }
       }
      }
      Dijkl_term[ii]=ZERO;
     }
     cout<<"l coefficient changed"<<endl;
     for(ii=0;ii<N5;ii++)
     {
      Dijkl=Dijkl_term2[ii];
      if(abs(Dijkl)>=threshold)
      {
       element_primeL[1]=ii/N4;
       element_primeL[0]=(ii-element_primeL[1]*N4)/N3;
       for(jj=0;jj<N2;jj++)
       {
        Dijkl_change=Dijkl*B[element_primeL[0]][jj];
        if(abs(Dijkl_change)>=threshold)
        {
         new_index=ii-element_primeL[0]*N3+jj*N3;
         Dijkl_term[new_index]=Dijkl_term[new_index]+Dijkl_change;
        }
       }
      }
      Dijkl_term2[ii]=ZERO;
     }
     cout<<"k coefficient changed"<<endl;
     for(ii=0;ii<N5;ii++)
     {
      Dijkl=Dijkl_term[ii];
      if(abs(Dijkl)>=threshold)
      {
       element_primeL[1]=ii/N4;
       element_primeL[0]=(ii-element_primeL[1]*N4)/N3;
       elementL[1]=(ii-element_primeL[0]*N3-element_primeL[1]*N4)/N2;
       for(jj=0;jj<N2;jj++)
       {
        Dijkl_change=Dijkl*B[elementL[1]][jj];
        if(abs(Dijkl_change)>=threshold)
        {
         new_index=ii-elementL[1]*N2+jj*N2;
         Dijkl_term2[new_index]=Dijkl_term2[new_index]+Dijkl_change;
        }
       }
      }
      Dijkl_term[ii]=ZERO;
     }
     cout<<"j coefficient changed"<<endl;
     for(ii=0;ii<N5;ii++)
     {
      Dijkl=Dijkl_term2[ii];
      if(abs(Dijkl)>=threshold)
      {
       element_primeL[1]=ii/N4;
       element_primeL[0]=(ii-element_primeL[1]*N4)/N3;
       elementL[1]=(ii-element_primeL[0]*N3-element_primeL[1]*N4)/N2;
       elementL[0]=(ii-elementL[1]*N2-element_primeL[0]*N3-element_primeL[1]*N4);
       for(jj=0;jj<N2;jj++)
       {
        Dijkl_change=Dijkl*B[elementL[0]][jj];
        if(abs(Dijkl_change)>=threshold)
        {
         new_index=ii-elementL[0]+jj;
         Dijkl_term[new_index]=Dijkl_term[new_index]+Dijkl_change;
        }
       }
      }
      Dijkl_term2[ii]=ZERO;
     }
     cout<<"i coefficient changed"<<endl;
    }
    //Change done, proceed to print:
    cout<<"Change done!"<<endl;
    cout<<endl;
    trace=ZERO;
    for(ii=0;ii<N2;ii++)
    {
     for(jj=0;jj<N2;jj++)
     {
      if(abs(Dijkl_term[ii*(1+N3)+jj*(N2+N4)])>=threshold)
      {
       trace=trace+Dijkl_term[ii*(1+N3)+jj*(N2+N4)];
      }
     }
    }
    cout<<endl;
    cout<<"Trace of the transformed 2RDM:  "<<setprecision(10)<<fixed<<scientific<<trace<<endl;
    cout<<endl;
    if(reduce)
    {
     cout<<endl;
     cout<<"Reduce the p <-> r and q <-> s terms"<<endl;
     if(red_sym)
     {
      cout<<"also keep only one term among orb_p(1) orb_q(2) orb_r(1) orb_s(2) = orb_q(2) orb_p(1) orb_s(2) orb_r(1)"<<endl;
     }
     cout<<endl;
     for(ii=0;ii<N5;ii++)
     {
      if(abs(Dijkl_term[ii])>=threshold)
      {
       Lred=(ii/N4);
       Kred=((ii-Lred*N4)/N3);
       Jred=((ii-Kred*N3-Lred*N4)/N2);
       Ired=(ii-Jred*N2-Kred*N3-Lred*N4);
       if(Ired!=Kred && Jred!=Lred)
       {
        Dijkl_term[ii]=Dijkl_term[ii]
                      +Dijkl_term[Kred+Jred*N2+Ired*N3+Lred*N4]
                      +Dijkl_term[Ired+Lred*N2+Kred*N3+Jred*N4]
                      +Dijkl_term[Kred+Lred*N2+Ired*N3+Jred*N4];
        Dijkl_term[Kred+Jred*N2+Ired*N3+Lred*N4]=ZERO;
        Dijkl_term[Ired+Lred*N2+Kred*N3+Jred*N4]=ZERO;
        Dijkl_term[Kred+Lred*N2+Ired*N3+Jred*N4]=ZERO;
       }
       if(Ired!=Kred && Jred==Lred)
       {
        Dijkl_term[ii]=Dijkl_term[ii]
                      +Dijkl_term[Kred+Jred*N2+Ired*N3+Lred*N4];
        Dijkl_term[Kred+Jred*N2+Ired*N3+Lred*N4]=ZERO;
       }
       if(Ired==Kred && Jred!=Lred)
       {
        Dijkl_term[ii]=Dijkl_term[ii]
                      +Dijkl_term[Ired+Lred*N2+Kred*N3+Jred*N4];
        Dijkl_term[Ired+Lred*N2+Kred*N3+Jred*N4]=ZERO;
       }
      }
     }
    }
    cout<<"Printing the matrix:"<<endl;
    cout<<endl;
    ofstream dm2_file;
    dm2_file.open(("Conv_"+name_dm2).c_str(), ios::binary);
    for(ii=0;ii<N5;ii++)
    {
     Dijkl=Dijkl_term[ii];
     if(abs(Dijkl)>=threshold)
     {
      element_prime[1]=(int)(ii/N4);
      element_prime[0]=(int)((ii-element_prime[1]*N4)/N3);
      element[1]=(int)((ii-element_prime[0]*N3-element_prime[1]*N4)/N2);
      element[0]=(int)(ii-element[1]*N2-element_prime[0]*N3-element_prime[1]*N4);
      if(!(element[0]==element[1] && element[1]==element_prime[0] && element_prime[0]==element_prime[1]) && (reduce && red_sym))
      {
       Dijkl_term[element[0]+element[1]*N2+element_prime[0]*N3+element_prime[1]*N4]=ZERO;
       Dijkl_term[element[1]+element[0]*N2+element_prime[1]*N3+element_prime[0]*N4]=ZERO;
      } 
      element[0]++;element[1]++;element_prime[0]++;element_prime[1]++;
      dm2_file.seekp(RECORD_DELIMITER_LENGTH, ios::cur);
      dm2_file.write((char*) &element[0], sizeof(element[0]));
      dm2_file.write((char*) &element[1], sizeof(element[1]));
      dm2_file.write((char*) &element_prime[0], sizeof(element_prime[0]));
      dm2_file.write((char*) &element_prime[1], sizeof(element_prime[1]));
      dm2_file.write((char*) &Dijkl, sizeof(Dijkl));
      dm2_file.seekp(RECORD_DELIMITER_LENGTH, ios::cur);
     }
    }
    element[0]=0;
    element[1]=0;
    element_prime[0]=0;
    element_prime[1]=0;
    Dijkl=ZERO;
    dm2_file.seekp(RECORD_DELIMITER_LENGTH, ios::cur);
    dm2_file.write((char*) &element[0], sizeof(element[0]));
    dm2_file.write((char*) &element[1], sizeof(element[1]));
    dm2_file.write((char*) &element_prime[0], sizeof(element_prime[0]));
    dm2_file.write((char*) &element_prime[1], sizeof(element_prime[1]));
    dm2_file.write((char*) &Dijkl, sizeof(Dijkl));
    dm2_file.seekp(RECORD_DELIMITER_LENGTH, ios::cur);
    dm2_file.close();
   }
   else
   {
    cout<<endl;
    check_dm2.close();
    cout<<"Error! Could not open file: "<<name_dm2<<endl;
    cout<<endl;
   }
   for(i=0;i<nbasis;i++)
   {
    delete[] B[i];B[i]=NULL;
   }
   delete[] B; B=NULL;
  }
  else
  {
   cout<<endl;
   open_fchk.close();
   cout<<"Error! Could not open file: "<<name_fchk<<endl;
   cout<<endl;
  }
 }
 else
 {
  cout<<endl;
  cout<<"Please, include the parameters needed"<<endl;
  cout<<"name_file.dm2 name_file.fchk threshold(e.g. 1e-10) y.or.n(store 2 DM2)"<<endl;
  cout<<"or"<<endl;
  cout<<"name_file.dm2 name_file.fchk threshold(e.g. 1e-10) y.or.n(store 2 DM2) i4(use aa_2D^ii _ii terms)"<<endl;
  cout<<"or"<<endl;
  cout<<"name_file.dm2 name_file.fchk threshold(e.g. 1e-10) y.or.n(store 2 DM2) r(reduce p <-> r and q <-> s terms)"<<endl;
  cout<<"or"<<endl;
  cout<<"name_file.dm2 name_file.fchk threshold(e.g. 1e-10) y.or.n(store 2 DM2) rs(reduce p <-> r, q <-> s terms and symmetry)"<<endl;
  cout<<"or"<<endl;
  cout<<"name_file.dm2 name_file.fchk threshold(e.g. 1e-10) y.or.n(store 2 DM2) a(all 2D elements are given)"<<endl;
  cout<<"or"<<endl;
  cout<<"name_file.dm2 name_file.fchk threshold(e.g. 1e-10) y.or.n(store 2 DM2) donof(2D elements are from DoNOF code)"<<endl;
  cout<<"or"<<endl;
  cout<<"name_file.dm2 name_file.fchk threshold(e.g. 1e-10) y.or.n(store 2 DM2) i8all(2D elements are from PSI4 code)"<<endl;
  cout<<"or"<<endl;
  cout<<"name_file.dm2 name_file.fchk threshold(e.g. 1e-10) y.or.n(store 2 DM2) i4r(aa_2D^ii _ii terms and reduce p <-> r and q <-> s terms)"<<endl;
  cout<<"or"<<endl;
  cout<<"name_file.dm2 name_file.fchk threshold(e.g. 1e-10) y.or.n(store 2 DM2) i4rs(aa_2D^ii _ii terms and reduce p <-> r, q <-> s terms and symm)"<<endl;
  cout<<"or"<<endl;
  cout<<"name_file.dm2 name_file.fchk threshold(e.g. 1e-10) y.or.n(store 2 DM2) ar(all 2D elements are given and reduce p <-> r and q <-> s terms)"<<endl;
  cout<<"or"<<endl;
  cout<<"name_file.dm2 name_file.fchk threshold(e.g. 1e-10) y.or.n(store 2 DM2) ars(all 2D elements are given and reduce p <-> r, q <-> s terms and symm)"<<endl;
  cout<<"or"<<endl;
  cout<<"name_file.dm2 name_file.fchk threshold(e.g. 1e-10) y.or.n(store 2 DM2) donofr(2D elements are from DoNOF code and reduce p <-> r and q <-> s terms)"<<endl;
  cout<<endl;
  cout<<"name_file.dm2 name_file.fchk threshold(e.g. 1e-10) y.or.n(store 2 DM2) donofrs(2D elements are from DoNOF code and reduce p <-> r, q <-> s terms and symm)"<<endl;
  cout<<endl;
  cout<<"name_file.dm2 name_file.fchk threshold(e.g. 1e-10) y.or.n(store 2 DM2) i8allr(2D elements are from PSI4 code and reduce p <-> r and q <-> s terms)"<<endl;
  cout<<endl;
  cout<<"name_file.dm2 name_file.fchk threshold(e.g. 1e-10) y.or.n(store 2 DM2) i8allrs(2D elements are from PSI4 code and reduce p <-> r, q <-> s terms and symm)"<<endl;
  cout<<endl;
  cout<<endl;
 }
 cout<<"#################################################"<<endl;
 cout<<endl;
 cout<<" Normal termination of chimpanC code          "<<endl;
 cout<<endl;
 cout<<"                   __,_ __                    "<<endl;
 cout<<"             .--.  .-"<<"     "<<"-.  .--.    "<<endl;
 cout<<"            / .. \\/  .-. .-.  \\/ .. \\         "<<endl;
 cout<<"           | |  '|  /   Y   \\  |'  | |        "<<endl;
 cout<<"           | \\   \\  \\ 0 | 0 /  /   / |        "<<endl;
 cout<<"            \\ '- ,\\.-"<<"`` ``"<<"-./, -' /         "<<endl;
 cout<<"             `'-' /_   ^ ^   _\\ '-'`          "<<endl;
 cout<<"             .--'|  \\._ _ _./  |'--.          "<<endl;
 cout<<"           /`    \\   \\.-.  /   /    `\\        "<<endl;
 cout<<"          /       '._/  |-' _.'       \\       "<<endl;
 cout<<"         /          ;  /--~'   |       \\      "<<endl;
 cout<<"        /        .'\\|.-\\--.     \\       \\     "<<endl;
 cout<<"       /   .'-. /.-.;\\  |\\|'~'-.|\\       \\    "<<endl;
 cout<<"       \\       `-./`|_\\_/ `     `\\'.      \\   "<<endl;
 cout<<"        '.      ;     ___)        '.`;    /   "<<endl;
 cout<<"          '-.,_ ;     ___)          \\/   /    "<<endl;
 cout<<"           \\   ``'------'\\       \\   `  /     "<<endl;
 cout<<"            '.    \\       '.      |   ;/_     "<<endl;
 cout<<"          ___>     '.       \\_ _ _/   ,  '--. "<<endl;
 cout<<"        .'   '.   .-~~~~~-. /     |--'`~~-.  \\ "<<endl;
 cout<<"       // / .---'/  .-~~-._/ / / /---..__.'  /"<<endl;
 cout<<"      ((_(_/    /  /      (_(_(_(---.__    .' "<<endl;
 cout<<"                | |     _              `~~`   "<<endl;
 cout<<"                | |     \\'.                   "<<endl;
 cout<<"                 \\ '....' |                   "<<endl;
 cout<<"                  '.,___.'                    "<<endl;
 cout<<"                                              "<<endl;
 cout<<endl;
 cout<<"#################################################"<<endl;
 return 0;
}

//Functions for FCHK
void line_fill(int &number, string line_read)
{
  line_read=line_read.substr(56,line_read.length());
  stringstream ss(line_read); //string to int in ss
  ss>>number;
}

//Write basis file:
void write_basis(FILE *pFile,string name_dm2,int &nucleous,double Atom_coord[3],double Exp,int nx,int ny,int nz)
{
 string aux,temp;
 temp=name_dm2.substr(0,name_dm2.length()-4)+"_basis.temp";
 ofstream temporal((temp).c_str());
 temporal<<setw(3)<<nucleous<<"\t"<<setprecision(8)<<fixed<<scientific<<Atom_coord[0]<<"\t"<<Atom_coord[1]<<"\t"<<Atom_coord[2]<<"\t"<<Exp<<"\t\t"<<nx<<"\t"<<ny<<"\t"<<nz<<endl;
 temporal.close();
 ifstream temporal2;
 temporal2.open((temp).c_str());
 sys=system(("rm "+temp).c_str());
 getline(temporal2,aux);
 fputs((aux+"\n").c_str(),pFile);
}

//Generate n, l and m
void Quant_fill(int **Quant,int styp)
{
 int i,j,k,qfil,permute[3];
 qfil=0;
 for(i=styp;i>-1;i=i-1)
 {for(j=i;j>-1;j=j-1)
   {for(k=j;k>-1;k=k-1)
    {
     if((i+j+k)==styp)
     {
       Quant[qfil][0]=i;
       Quant[qfil][1]=j;
       Quant[qfil][2]=k;
       qfil++;
       if(j!=k)
       {
        Quant[qfil][0]=i;
        Quant[qfil][1]=k;
        Quant[qfil][2]=j;
        qfil++;
       }
       if(i!=j)
       {
        Quant[qfil][0]=j;
        Quant[qfil][1]=i;
        Quant[qfil][2]=k;
        qfil++;
       }
       if((i!=j) && (j!=k))
       {
        Quant[qfil][0]=j;
        Quant[qfil][1]=k;
        Quant[qfil][2]=i;
        qfil++;
        Quant[qfil][0]=k;
        Quant[qfil][1]=i;
        Quant[qfil][2]=j;
        qfil++;
       }
       if(k!=i)
       {
        Quant[qfil][0]=k;
        Quant[qfil][1]=j;
        Quant[qfil][2]=i;
        qfil++;
       }
     }
   }
  }
 }
 if(styp==3)
 {
  if(!donofdm2)
  {
   permute[0]=Quant[5][0];permute[1]=Quant[5][1];permute[2]=Quant[5][2];
   Quant[5][0]=Quant[4][0];Quant[5][1]=Quant[4][1];Quant[5][2]=Quant[4][2];
   Quant[4][0]=Quant[3][0];Quant[4][1]=Quant[3][1];Quant[4][2]=Quant[3][2];
   Quant[3][0]=permute[0];Quant[3][1]=permute[1];Quant[3][2]=permute[2];
   permute[0]=Quant[7][0];permute[1]=Quant[7][1];permute[2]=Quant[7][2];
   Quant[7][0]=Quant[8][0];Quant[7][1]=Quant[8][1];Quant[7][2]=Quant[8][2];
   Quant[8][0]=permute[0];Quant[8][1]=permute[1];Quant[8][2]=permute[2];
  }
  else
  {
   permute[0]=Quant[3][0];permute[1]=Quant[3][1];permute[2]=Quant[3][2];
   Quant[3][0]=Quant[4][0];Quant[3][1]=Quant[4][1];Quant[3][2]=Quant[4][2];
   Quant[4][0]=Quant[5][0];Quant[4][1]=Quant[5][1];Quant[4][2]=Quant[5][2];
   Quant[5][0]=permute[0];Quant[5][1]=permute[1];Quant[5][2]=permute[2];
   permute[0]=Quant[8][0];permute[1]=Quant[8][1];permute[2]=Quant[8][2];
   Quant[8][0]=Quant[7][0];Quant[8][1]=Quant[7][1];Quant[8][2]=Quant[7][2];
   Quant[7][0]=Quant[6][0];Quant[7][1]=Quant[6][1];Quant[7][2]=Quant[6][2];
   Quant[6][0]=permute[0];Quant[6][1]=permute[1];Quant[6][2]=permute[2];
  }
 }
 if(styp==4)
 {
  permute[0]=Quant[6][0];permute[1]=Quant[6][1];permute[2]=Quant[6][2];
  Quant[6][0]=Quant[7][0];Quant[6][1]=Quant[7][1];Quant[6][2]=Quant[7][2];
  Quant[7][0]=permute[0];Quant[7][1]=permute[1];Quant[7][2]=permute[2];
 }
}

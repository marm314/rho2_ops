////////////////////////////////////////////
//----------------------------------------//
//----------------------------------------//
//---      DM2 Hartree-Fock            ---//
//---      DM2 Hartree-Fock (like)     ---//
//----------------------------------------//
//----------------------------------------//
//-- Developed by: M. Rodriguez-Mayorga --//
//--   email: marm3.14@gmail.com        --//
//----------------------------------------//
//----------------------------------------//
////////////////////////////////////////////
#include<iostream>
#include<fstream>
#include<cstring>
#include<cmath>
#include<iomanip> 
#include<stdlib.h>
#include<stdio.h>
#define ZERO 0.0e0
#define TWO 2.0e0

using namespace std;
////////////////
// Functions  //
////////////////
void trace_dmn(string name_file);
void Create_HF_dm2_from_dm1(int *dm1,int nbasis, string name_dmn);
void Create_HF_dm1_file(int *dm1,int nbasis, string name_dmn);
void Create_HFl_dm2_from_dm1(double **dm1,int nbasis, string name_dmn);
double dm1_to_dm2(int &i, int &j, int &k, int &l, int *dm1);
double dm1_to_dm2_hfl(int &i, int &j, int &k, int &l, double **dm1);
void fill_in_dm1(double **dm1, string name_file);
//////////////////////
// Global variables //
//////////////////////
const int RECORD_DELIMITER_LENGTH=4;
double threshold;
///////////////////
// Main Function //
///////////////////
int main(int argc, char *argv[])
{
 if(argc==5)
 {
  //Argument 0 name_program, 1 first.dm2, 2 nelects, 3 nbasis, 4 multiplicity.  
  int i,nalpha,nbeta,nexcess;
  string name_file(argv[1]);
  int  nelectrons=(int)atof(argv[2]);
  int  nbasis=atof(argv[3]);
  int  multiplicity=atof(argv[4]);
  int *aux,*dm1;
  threshold=ZERO;
  nexcess=(multiplicity-1);
  nbeta=(nelectrons-nexcess)/2;
  nalpha=nelectrons-nbeta;
  aux=new int[nbasis];
  dm1=new int[2*nbasis];
  for(i=0;i<nbasis;i++)
  {
   aux[i]=0;
   dm1[2*i]=ZERO;
   dm1[2*i+1]=ZERO;
  }
  for(i=0;i<nalpha;i++)
  {
   aux[i]=aux[i]+1;
  }
  for(i=0;i<nbeta;i++)
  {
   aux[i]=aux[i]+1;
  }
  for(i=0;i<nbasis;i++)
  {
   if(aux[i]!=0)
   {
    if(aux[i]%2==0)
    {
     dm1[2*i]=1;
     dm1[2*i+1]=1;
    }
    else
    {
     dm1[2*i]=1;
    }
   }
   else
   {
    break;
   }
  } 
  delete[] aux; aux=NULL; 
  nbasis=2*nbasis;  
  Create_HF_dm1_file(dm1,nbasis,name_file);
  Create_HF_dm2_from_dm1(dm1,nbasis,name_file);
  delete[] dm1; dm1=NULL;
  name_file="HF_"+name_file+".dm2"; 
  trace_dmn(name_file);
 }
 else if(argc==4)
 {
  int i,j;
  string name_file_dm1(argv[1]);
  int  nbasis=atof(argv[2]);
  threshold=atof(argv[3]);
  double **dm1;
  nbasis=2*nbasis;
  dm1=new double*[nbasis];
  for(i=0;i<nbasis;i++)
  {
   dm1[i]=new double[nbasis]; 
  }
  for(i=0;i<nbasis;i++)
  {
   for(j=0;j<nbasis;j++)
   {
    dm1[i][j]=ZERO; 
   }
  }
  fill_in_dm1(dm1,name_file_dm1); 
  Create_HFl_dm2_from_dm1(dm1,nbasis,name_file_dm1); 
  name_file_dm1="HFl_"+name_file_dm1.substr(0,name_file_dm1.length()-4)+".dm2";
  trace_dmn(name_file_dm1);
  for(i=0;i<nbasis;i++)
  {
   delete[] dm1[i];dm1[i]=NULL;
  } 
  delete[] dm1; dm1=NULL;   
 }
 else
 {
  cout<<"Write the parameters:"<<endl;
  cout<<"a) name.dm1(in) nbasis(num AOs) threshold > To build a  HFlike-DM2"<<endl;
  cout<<"b) name.dm2(out) nelectrons nbasis(num AOs) multiplicity(1 for singlet)  -> To build a HF-DM2"<<endl;
 }
 return 0;
}
////////////////////////////
// Functions declarations //
////////////////////////////
void fill_in_dm1(double **dm1,string name_file)
{
 int element[1],element_prime[1];
 double Dij,Trace=ZERO;
 element[0]=10;element_prime[0]=10;
 ifstream input_data(name_file.c_str(), ios::binary);
 if(input_data.good())
 {
  while((element[0]!=0 || element_prime[0]!=0))
  {
   input_data.seekg(RECORD_DELIMITER_LENGTH, ios::cur);
   input_data.read((char*) &element[0], sizeof(element[0]));
   input_data.read((char*) &element_prime[0], sizeof(element_prime[0]));
   input_data.read((char*) &Dij, sizeof(Dij));
   input_data.seekg(RECORD_DELIMITER_LENGTH, ios::cur);
   if(abs(Dij)>=threshold)
   {
    if(element[0]==element_prime[0])
    {Trace=Trace+Dij;}
    if(dm1[element[0]-1][element_prime[0]-1]==ZERO)
    {
     dm1[element[0]-1][element_prime[0]-1]=Dij;
    }
    else
    {
     dm1[element[0]-1][element_prime[0]-1]=dm1[element[0]-1][element_prime[0]-1]+Dij;
    }
    dm1[element_prime[0]-1][element[0]-1]=dm1[element[0]-1][element_prime[0]-1];
   }
  }
  cout<<"Trace(DM1 "<<name_file<<" in): "<<setprecision(12)<<fixed<<Trace<<endl;
 }
 input_data.close();
}

void trace_dmn(string name_file)
{
 int element[2],element_prime[2];
 double Dijkl,Trace=ZERO;
 element[0]=10;element[1]=10;element_prime[0]=10;element_prime[1]=10;
 ifstream input_data(name_file.c_str(), ios::binary);
 if(input_data.good())
 {
  while((element[0]!=0 || element_prime[0]!=0) || (element[1]!=0 || element_prime[1]!=0))
  {
   input_data.seekg(RECORD_DELIMITER_LENGTH, ios::cur);
   input_data.read((char*) &element[0], sizeof(element[0]));
   input_data.read((char*) &element[1], sizeof(element[1]));
   input_data.read((char*) &element_prime[0], sizeof(element_prime[0]));
   input_data.read((char*) &element_prime[1], sizeof(element_prime[1]));
   input_data.read((char*) &Dijkl, sizeof(Dijkl));
   input_data.seekg(RECORD_DELIMITER_LENGTH, ios::cur);
   if(abs(Dijkl)>=ZERO)
   {
    if(element[0]==element_prime[0] && element[1]==element_prime[1])
    {Trace=Trace+Dijkl;}
   }
  }
  cout<<"Trace(DM2 "<<name_file<<" out): "<<setprecision(12)<<fixed<<Trace*TWO<<endl;
 }
 input_data.close();
}

void Create_HFl_dm2_from_dm1(double **dm1,int nbasis, string name_dmn)
{
 int i,j,k,l,elements[4];
 double DijklHFl;
 ofstream output_data(("HFl_"+name_dmn.substr(0,name_dmn.length()-4)+".dm2").c_str(),ios::out | ios::binary);
 for(i=0;i<nbasis-1;i++)
 {
  for(j=i+1;j<nbasis;j++)
  {  
   for(k=i;k<nbasis-1;k++)
   {
    for(l=k+1;l<nbasis;l++)
    {
     if(i==k)
     { 
      if(j<=l)
      {
       DijklHFl=dm1_to_dm2_hfl(i,j,k,l,dm1);
       if(abs(DijklHFl)>threshold)
       {
        elements[0]=i+1;
        elements[1]=j+1;
        elements[2]=k+1;
        elements[3]=l+1;
        output_data.seekp(RECORD_DELIMITER_LENGTH, ios::cur);
        output_data.write((char*) &elements[0], sizeof(elements[0]));
        output_data.write((char*) &elements[1], sizeof(elements[1]));
        output_data.write((char*) &elements[2], sizeof(elements[2]));
        output_data.write((char*) &elements[3], sizeof(elements[3]));
        output_data.write((char*) &DijklHFl, sizeof(DijklHFl));
        output_data.seekp(RECORD_DELIMITER_LENGTH, ios::cur);
       }
      } 
     } 
     else
     {
      DijklHFl=dm1_to_dm2_hfl(i,j,k,l,dm1);
      if(abs(DijklHFl)>threshold)
      {
       elements[0]=i+1;
       elements[1]=j+1;
       elements[2]=k+1;
       elements[3]=l+1;
       output_data.seekp(RECORD_DELIMITER_LENGTH, ios::cur);
       output_data.write((char*) &elements[0], sizeof(elements[0]));
       output_data.write((char*) &elements[1], sizeof(elements[1]));
       output_data.write((char*) &elements[2], sizeof(elements[2]));
       output_data.write((char*) &elements[3], sizeof(elements[3]));
       output_data.write((char*) &DijklHFl, sizeof(DijklHFl));
       output_data.seekp(RECORD_DELIMITER_LENGTH, ios::cur);
      }
     }
    } 
   }
  }
 } 
 elements[0]=0;
 elements[1]=0;
 elements[2]=0;
 elements[3]=0;
 DijklHFl=ZERO;
 output_data.seekp(RECORD_DELIMITER_LENGTH, ios::cur);
 output_data.write((char*) &elements[0], sizeof(elements[0]));
 output_data.write((char*) &elements[1], sizeof(elements[1]));
 output_data.write((char*) &elements[2], sizeof(elements[2]));
 output_data.write((char*) &elements[3], sizeof(elements[3]));
 output_data.write((char*) &DijklHFl, sizeof(DijklHFl));
 output_data.seekp(RECORD_DELIMITER_LENGTH, ios::cur);
 output_data.close();
}

void Create_HF_dm1_file(int *dm1,int nbasis, string name_dmn)
{
 int i,elements[1];
 double niHF;
 ofstream output_data(("HF_"+name_dmn+".dm1").c_str(),ios::out | ios::binary);
 for(i=0;i<nbasis;i++)
 {
  elements[0]=i+1;
  output_data.seekp(RECORD_DELIMITER_LENGTH, ios::cur);
  output_data.write((char*) &elements[0], sizeof(elements[0]));
  output_data.write((char*) &elements[0], sizeof(elements[0]));
  output_data.write((char*) &niHF, sizeof(niHF));
  output_data.seekp(RECORD_DELIMITER_LENGTH, ios::cur);
 }
 elements[0]=0;niHF=ZERO;
 output_data.seekp(RECORD_DELIMITER_LENGTH, ios::cur);
 output_data.write((char*) &elements[0], sizeof(elements[0]));
 output_data.write((char*) &elements[0], sizeof(elements[0]));
 output_data.write((char*) &niHF, sizeof(niHF));
 output_data.seekp(RECORD_DELIMITER_LENGTH, ios::cur);
 output_data.close();
}

void Create_HF_dm2_from_dm1(int *dm1,int nbasis, string name_dmn)
{
 int i,j,k,l,elements[4];
 double DijklHF;
 ofstream output_data(("HF_"+name_dmn+".dm2").c_str(),ios::out | ios::binary);
 for(i=0;i<nbasis-1;i++)
 {
  for(j=i+1;j<nbasis;j++)
  {  
   for(k=i;k<nbasis-1;k++)
   {
    for(l=k+1;l<nbasis;l++)
    {
     if(i==k)
     { 
      if(j<=l)
      {
       DijklHF=dm1_to_dm2(i,j,k,l,dm1);
       if(abs(DijklHF)>threshold)
       {
        elements[0]=i+1;
        elements[1]=j+1;
        elements[2]=k+1;
        elements[3]=l+1;
        output_data.seekp(RECORD_DELIMITER_LENGTH, ios::cur);
        output_data.write((char*) &elements[0], sizeof(elements[0]));
        output_data.write((char*) &elements[1], sizeof(elements[1]));
        output_data.write((char*) &elements[2], sizeof(elements[2]));
        output_data.write((char*) &elements[3], sizeof(elements[3]));
        output_data.write((char*) &DijklHF, sizeof(DijklHF));
        output_data.seekp(RECORD_DELIMITER_LENGTH, ios::cur);
       }
      } 
     } 
     else
     {
      DijklHF=dm1_to_dm2(i,j,k,l,dm1);
      if(abs(DijklHF)>threshold)
      {
       elements[0]=i+1;
       elements[1]=j+1;
       elements[2]=k+1;
       elements[3]=l+1;
       output_data.seekp(RECORD_DELIMITER_LENGTH, ios::cur);
       output_data.write((char*) &elements[0], sizeof(elements[0]));
       output_data.write((char*) &elements[1], sizeof(elements[1]));
       output_data.write((char*) &elements[2], sizeof(elements[2]));
       output_data.write((char*) &elements[3], sizeof(elements[3]));
       output_data.write((char*) &DijklHF, sizeof(DijklHF));
       output_data.seekp(RECORD_DELIMITER_LENGTH, ios::cur);
      }
     }
    } 
   }
  }
 } 
 elements[0]=0;
 elements[1]=0;
 elements[2]=0;
 elements[3]=0;
 DijklHF=ZERO;
 output_data.seekp(RECORD_DELIMITER_LENGTH, ios::cur);
 output_data.write((char*) &elements[0], sizeof(elements[0]));
 output_data.write((char*) &elements[1], sizeof(elements[1]));
 output_data.write((char*) &elements[2], sizeof(elements[2]));
 output_data.write((char*) &elements[3], sizeof(elements[3]));
 output_data.write((char*) &DijklHF, sizeof(DijklHF));
 output_data.seekp(RECORD_DELIMITER_LENGTH, ios::cur);
 output_data.close();
}

double dm1_to_dm2(int &i, int &j, int &k, int &l, int *dm1)
{
 double Dijkl=ZERO;
 if(i==k && j==l)
 {
  Dijkl=(double)dm1[i]*(double)dm1[j];
 }
 if(i==l && j==k)
 {
  Dijkl=Dijkl-(double)dm1[i]*(double)dm1[j];
 }
 return Dijkl; 
}

double dm1_to_dm2_hfl(int &i, int &j, int &k, int &l, double **dm1)
{
 double Dijkl=ZERO;
 Dijkl=dm1[i][k]*dm1[j][l]-dm1[i][l]*dm1[j][k];
 return Dijkl; 
}

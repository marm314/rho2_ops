#include<iostream>
#include<iomanip>
#include<stdlib.h>
#include<fstream>
#include<cstring>
#include<vector>
#include<cmath>
#include<sstream>
#include<algorithm>
#define ZERO 0.0e0

using namespace std;

void find_indices(int i,int j,int k,int l,int &element0,int &element1,int &element_prime0,int &element_prime1,int *order,int nbasis);
int RECORD_DELIMITER_LENGTH=4;

vector <string> MOsyms;
vector <string> Calcsyms;

int main(int argc,char *argv[])
{
 cout<<"----------------------------------------"<<endl;
 cout<<"----------------------------------------"<<endl;
 cout<<"---   PSI4 INTERFACE FOR RHO2_OPS    ---"<<endl;
 cout<<"----------------------------------------"<<endl;
 cout<<"----------------------------------------"<<endl;
 cout<<"-- Developed by: M. Rodriguez-Mayorga --"<<endl;
 cout<<"--   email: marm3.14@gmail.com        --"<<endl;
 cout<<"----------------------------------------"<<endl;
 cout<<"----------------------------------------"<<endl;
 cout<<endl;
 if(argc!=3)
 {
  cout<<"Include the parameters:"<<endl;
  cout<<"name_psi4.out threshold(e.g. 1e-12)"<<endl;
  cout<<endl;
  cout<<"----------------------------------------"<<endl;
  cout<<"--        Normal termination          --"<<endl;
  cout<<"----------------------------------------"<<endl;
  cout<<"----------------------------------------"<<endl;
  return 0;
 }
 bool read=true;
 int i,j,k,l,nth,nbasis,nsym=0;
 int N1,N2,N5col,Ndocc=0,Nuocc=0;
 int *order_orbs,element[2],element_prime[2];
 long int ii,jj;
 double **DM2,Dijkl;
 nth=0;
 string line,subline;
 string name_out=argv[1];
 double threshold=atof(argv[2]);
 string name_dm2;
 name_dm2=name_out.substr(0,(name_out.length()-3))+"dm2";
 ofstream tmp_file("tmp_psi4.txt");
 ifstream read_out(name_out);
 while(getline(read_out,line))
 {
  if(line.length()>30)
  {
   if(line.substr(16,4)=="DOCC" || line.substr(16,4)=="UOCC")
   {
    if(line.substr(16,4)=="DOCC")
    {
     stringstream ss1(line.substr(20,9)); 
     ss1>>Ndocc;
     cout<<"Number of frozen orbitals   : "<<setw(12)<<Ndocc<<endl;
    }
    else
    {
     stringstream ss2(line.substr(20,9)); 
     ss2>>Nuocc;
     cout<<"Number of secondary orbitals: "<<setw(12)<<Nuocc<<endl;
    }
   }
  }
  if(line.length()>29)
  {
   if(line.substr(0,28)=="    Number of basis function")
   {
    stringstream ss(line.substr(29,line.length())); 
    ss>>nbasis;
    N1=nbasis;
    N2=N1*N1; 
    cout<<"Size of the basis read: "<<setw(12)<<nbasis<<endl;
   }
  }
  if(line.length()>23)
  {
   if(line.substr(6,14)=="Pre-Iterations")
   {
    for(i=0;i<4;i++)
    {
     getline(read_out,line);
    }
    do
    {
     getline(read_out,line);
     if(line.substr(3,7)!="-------")
     {
      subline=line.substr(4,5);
      Calcsyms.push_back(subline);
     }
     nsym++;
    }while(line.substr(3,7)!="-------");
    nsym--;
    cout<<"Size of the nsym  read: "<<setw(12)<<nsym<<endl;
   }
  }
  if(line.length()>18)
  {
   if(line.substr(12,3)=="Sym")
   {
    i=14;
    do
    {
     i++;
     if(line[i]!=' ')
     {
      j=0;
      subline="   ";
      do
      {
       subline[j]=line[i];
       j++;
       i++;
      }while(line[i]!=' ' && i<line.length());
      MOsyms.push_back(subline);
     }
    }while(i<line.length());
   }
   if(line.substr(5,13)=="MO-basis TPDM")
   {
    cout<<"TPDM found, proceed to read it"<<endl;
    getline(read_out,line);
    do
    {
     for(i=0;i<4;i++)
     {
      getline(read_out,line);
     }
     if(line.length()>0)
     {
      tmp_file<<line.substr(5,line.length()-5)<<endl;
      for(i=0;i<N2-1;i++)
      {
       getline(read_out,line);
       tmp_file<<line.substr(5,line.length()-5)<<endl;
      }
     }
     else
     {
      read=false;
     }
    }while(read);
   }
  }
 } 
 read_out.close();
 tmp_file.close();
 // Orb sym information
 cout<<endl;
 cout<<"Orbital symmetries used in the PSI4 calculation:"<<endl;
 cout<<endl;
 j=0;
 for(i=0;i<Calcsyms.size();i++)
 {
  Calcsyms[i].erase(std::remove_if(Calcsyms[i].begin(),Calcsyms[i].end(),::isspace),Calcsyms[i].end());
  cout<<setw(5)<<Calcsyms[i];
  j++;
  if(j==9){cout<<endl;j=0;}
 }
 if(j!=0){cout<<endl;}
 cout<<endl;
 cout<<"Orbital symmetries found in PSI4 output:"<<endl;
 cout<<endl;
 j=0;
 for(i=0;i<MOsyms.size();i++)
 {
  MOsyms[i].erase(std::remove_if(MOsyms[i].begin(),MOsyms[i].end(),::isspace),MOsyms[i].end());
  cout<<setw(5)<<MOsyms[i];
  j++;
  if(j==9){cout<<endl;j=0;}
 }
 if(j!=0){cout<<endl;}
 cout<<endl;
 // Allocate order of orbs
 N1=nbasis; 
 order_orbs=new int[N1];
 k=0;
 for(i=0;i<Calcsyms.size();i++)
 {
  for(j=0;j<MOsyms.size();j++)
  {
   if(MOsyms[j]==Calcsyms[i])
   {
    order_orbs[k]=j;
    k++;
   }
  }
 }
 cout<<endl;
 cout<<"Order of the orbitals attending to the symmetry (as used in DETCI):"<<endl;
 cout<<endl;
 j=0;
 for(i=0;i<N1;i++)
 {
  j++;
  cout<<setw(5)<<order_orbs[i]+1;
  if(j==9){cout<<endl;j=0;}
 }
 if(j!=0){cout<<endl;}
 cout<<endl;
 k=0;
 for(i=0;i<Calcsyms.size();i++)
 {
  for(j=0;j<N1;j++)
  {
   if(MOsyms[j]==Calcsyms[i])
   {
    order_orbs[j]=k;
    k++;
   }
  }
 }
 cout<<endl;
 cout<<"Numbering of the orbitals attending to the symmetry (as used in DETCI):"<<endl;
 cout<<endl;
 j=0;
 for(i=0;i<N1;i++)
 {
  j++;
  cout<<setw(5)<<order_orbs[i]+1;
  if(j==9){cout<<endl;j=0;}
 }
 if(j!=0){cout<<endl;}
 cout<<endl;
 // Allocate DM2 matrix
 DM2=new double*[N2];
 for(i=0;i<N2;i++)
 {
  DM2[i]=new double[N2];
  for(j=0;j<N2;j++)
  {
   DM2[i][j]=ZERO;
  }
 }
 k=0;N5col=N2/5;
 ifstream read_dm2("tmp_psi4.txt");
 for(l=0;l<N5col;l++)
 {
  for(j=0;j<N2;j++)
  {
   for(i=k;i<k+5;i++)
   {
    read_dm2>>DM2[j][i];
   }
  }
  k=k+5;
 }
 for(j=0;j<N2;j++)
 {
  for(i=k;i<N2;i++)
  {
   read_dm2>>DM2[j][i];
  }
 }
 read_dm2.close();
 system(("rm tmp_psi4.txt"));
 // Print the unformatted DM2 matrix
 ofstream dm2_out(name_dm2,ios::out | ios::binary);  
 for(i=0;i<N1;i++)
 {
  for(k=0;k<N1;k++)
  {
   ii=k+i*N1;
   for(j=0;j<N1;j++)
   {
    for(l=0;l<N1;l++)
    {
     jj=l+j*N1;
     if(abs(DM2[ii][jj])>threshold)
     {
      Dijkl=DM2[ii][jj];
      find_indices(i,j,k,l,element[0],element[1],element_prime[0],element_prime[1],order_orbs,N1);
      //cout<<"MRM: "<<element[0]+1<<" "<<element[1]+1<<" "<<element_prime[0]+1<<" "<<element_prime[1]+1<<" "<<Dijkl<<endl;
      //cout<<"MRM: "<<i<<" "<<j<<" "<<k<<" "<<l<<" "<<Dijkl<<endl;
      dm2_out.seekp(RECORD_DELIMITER_LENGTH, ios::cur);
      dm2_out.write((char*) &element[0], sizeof(element[0]));
      dm2_out.write((char*) &nth, sizeof(nth));
      dm2_out.write((char*) &element[1], sizeof(element[1]));
      dm2_out.write((char*) &nth, sizeof(nth));
      dm2_out.write((char*) &element_prime[0], sizeof(element_prime[0]));
      dm2_out.write((char*) &nth, sizeof(nth));
      dm2_out.write((char*) &element_prime[1], sizeof(element_prime[1]));
      dm2_out.write((char*) &nth, sizeof(nth));
      dm2_out.write((char*) &Dijkl, sizeof(Dijkl));
      dm2_out.seekp(RECORD_DELIMITER_LENGTH, ios::cur);
     }
    }
   }
  }
 }
 element[0]=0;
 element[1]=0;
 element_prime[0]=0;
 element_prime[1]=0;
 Dijkl=ZERO;
 dm2_out.seekp(RECORD_DELIMITER_LENGTH, ios::cur);
 dm2_out.write((char*) &element[0], sizeof(element[0]));
 dm2_out.write((char*) &nth, sizeof(nth));
 dm2_out.write((char*) &element[1], sizeof(element[1]));
 dm2_out.write((char*) &nth, sizeof(nth));
 dm2_out.write((char*) &element_prime[0], sizeof(element_prime[0]));
 dm2_out.write((char*) &nth, sizeof(nth));
 dm2_out.write((char*) &element_prime[1], sizeof(element_prime[1]));
 dm2_out.write((char*) &nth, sizeof(nth));
 dm2_out.write((char*) &Dijkl, sizeof(Dijkl));
 dm2_out.seekp(RECORD_DELIMITER_LENGTH, ios::cur);
 dm2_out.close();
 cout<<endl;
 cout<<"File "<<name_dm2<<" prepared for chimpanC."<<endl;
 cout<<"Note: use the donof option in chimpanC."<<endl;
 cout<<endl;
 // Deallocate all dynamic arrays
 for(i=0;i<N2;i++)
 {
  delete[] DM2[i];DM2[i]=NULL;
 } 
 delete[] DM2;DM2=NULL;
 delete[] order_orbs;order_orbs=NULL;
 // Print end of file
 cout<<endl;
 cout<<"----------------------------------------"<<endl;
 cout<<"--        Normal termination          --"<<endl;
 cout<<"----------------------------------------"<<endl;
 cout<<"----------------------------------------"<<endl;
 return 0;
}

void find_indices(int i,int j,int k,int l,int &element0,int &element1,int &element_prime0,int &element_prime1,int *order,int nbasis)
{
 int index;
 for(index=0;index<nbasis;index++)
 {
  if(i==order[index]){element0=index+1;} 
  if(j==order[index]){element1=index+1;} 
  if(k==order[index]){element_prime0=index+1;} 
  if(l==order[index]){element_prime1=index+1;} 
 }
}

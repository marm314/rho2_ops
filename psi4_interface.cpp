#include<iostream>
#include<iomanip>
#include<stdlib.h>
#include<fstream>
#include<cstring>
#include<vector>
#include<cmath>
#include<sstream>
#include<algorithm>
#include"Numbers.h"

using namespace std;

void find_indices(int i,int j,int k,int l,int &element0,int &element1,int &element_prime0,int &element_prime1,int *order,int nbasis);
double dm1_to_dm2_hfl(int &i, int &j, int &k, int &l, double **dm1);
int RECORD_DELIMITER_LENGTH=4;

vector<string>MOsyms;
vector<string>Calcsyms;

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
 bool read=true,*active_map;
 int i,j,k,l,m,nth,nbasis,nsym=0;
 int N1,N2,N5col,Ndocc=0,Nuocc=0,Norb_act=0,Norb_act2;
 int *order_orbs,element[2],element_prime[2],*docc_orbs,*uocc_orbs,*active_orbs;
 long int ii,jj;
 double **DM2,**DM1,**DM1spin,Dijkl,Trace,Trace1dm,Nelect;
 string line,subline;
 string name_out=argv[1];
 double threshold=atof(argv[2]);
 string name_dm2,name_dm1;
 name_dm2=name_out.substr(0,(name_out.length()-3))+"dm2";
 name_dm1=name_out.substr(0,(name_out.length()-3))+"dm1"; // TODO
 ofstream tmp_file("tmp_psi4.txt");
 ifstream read_out(name_out.c_str());
 ofstream prnt_sym("tmp1.txt");
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
     docc_orbs=new int[Calcsyms.size()];
     prnt_sym<<line.substr(29,line.length()-29)<<endl;
    }
    else
    {
     stringstream ss2(line.substr(20,9)); 
     ss2>>Nuocc;
     cout<<"Number of secondary orbitals: "<<setw(12)<<Nuocc<<endl;
     uocc_orbs=new int[Calcsyms.size()];
     prnt_sym<<line.substr(29,line.length()-29)<<endl;
    }
   }
   if(line.substr(6,14)=="Active (total)")
   {
    stringstream ss3(line.substr(20,9)); 
    ss3>>Norb_act;
    cout<<"Number of active orbitals   : "<<setw(12)<<Norb_act<<endl;
    active_orbs=new int[Calcsyms.size()];
    prnt_sym<<line.substr(29,line.length()-29)<<endl;
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
    Norb_act2=Norb_act*Norb_act;
    cout<<"Size of the DM2 matrix "<<setw(10)<<Norb_act2<<" x "<<setw(10)<<Norb_act2<<endl;
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
      for(i=0;i<Norb_act2-1;i++)
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
 prnt_sym.close();
 read_out.close();
 tmp_file.close();
 // Frozen, active, and secondary orbitals
 active_map=new bool[N2];
 for(i=0;i<N2;i++)
 {
  active_map[i]=true;
 }
 ifstream read_sym("tmp1.txt");
 for(i=0;i<Calcsyms.size();i++)
 {
  read_sym>>docc_orbs[i];
 }
 for(i=0;i<Calcsyms.size();i++)
 {
  read_sym>>active_orbs[i];
 }
 for(i=0;i<Calcsyms.size();i++)
 {
  read_sym>>uocc_orbs[i];
 }
 read_sym.close();
 cout<<endl;
 cout<<"Frozen orbitals per-symmetry groups:"<<endl;
 cout<<endl;
 j=0;
 for(i=0;i<Calcsyms.size();i++)
 {
  cout<<setw(5)<<docc_orbs[i];
  j++;
  if(j==9){cout<<endl;j=0;}
 } 
 if(j!=0){cout<<endl;}
 cout<<endl;
 cout<<"Active orbitals per-symmetry groups:"<<endl;
 cout<<endl;
 j=0;
 for(i=0;i<Calcsyms.size();i++)
 {
  cout<<setw(5)<<active_orbs[i];
  j++;
  if(j==9){cout<<endl;j=0;}
 } 
 if(j!=0){cout<<endl;}
 cout<<endl;
 cout<<"Secondary orbitals per-symmetry groups:"<<endl;
 cout<<endl;
 j=0;
 for(i=0;i<Calcsyms.size();i++)
 {
  cout<<setw(5)<<uocc_orbs[i];
  j++;
  if(j==9){cout<<endl;j=0;}
 } 
 if(j!=0){cout<<endl;}
 cout<<endl;
 nth=system("rm tmp1.txt");
 l=0;
 for(i=0;i<Calcsyms.size();i++)
 {
  for(j=0;j<docc_orbs[i];j++)
  {
   for(k=0;k<N1;k++)
   {
    active_map[l]=false;
    l++;
   }
  }
  for(j=0;j<active_orbs[i];j++)
  {
   for(k=0;k<Calcsyms.size();k++)
   {
    for(m=0;m<docc_orbs[k];m++)
    {
     active_map[l]=false;
     l++;
    }
    for(m=0;m<active_orbs[k];m++)
    {
     l++;
    }
    for(m=0;m<uocc_orbs[k];m++)
    {
     active_map[l]=false;
     l++;
    }
   }
  }
  for(j=0;j<uocc_orbs[i];j++)
  {
   for(k=0;k<N1;k++)
   {
    active_map[l]=false;
    l++;
   }
  }
 }
 /*
 // Print for debug the map with info about frozen, active, and secondary orbitals 
 l=0;
 for(i=0;i<N1;i++)
 {
  for(j=0;j<N1;j++)
  {
   if(active_map[l])
   {
    cout<<setw(5)<<" T ";
   }
   else
   {
    cout<<setw(5)<<" F ";
   }
   l++;
  }
  cout<<endl;
 }
 */  
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
 k=0;N5col=Norb_act2/5;
 ifstream read_dm2("tmp_psi4.txt");
 for(l=0;l<N5col;l++)
 {
  for(j=0;j<Norb_act2;j++)
  {
   for(i=k;i<k+5;i++)
   {
    read_dm2>>DM2[j][i];
   }
  }
  k=k+5;
 }
 for(j=0;j<Norb_act2;j++)
 {
  for(i=k;i<Norb_act2;i++)
  {
   read_dm2>>DM2[j][i];
  }
 }
 read_dm2.close(); 
 nth=system("rm tmp_psi4.txt");
 // Sort the 2-RDM
 if(N2!=Norb_act2)
 {
  //Columns
  k=Norb_act2;
  for(i=N2-1;i>-1;i--)
  {
   if(active_map[i])
   {
    k--;
    for(j=0;j<N2;j++)
    {
     DM2[j][i]=DM2[j][k];
     DM2[j][k]=ZERO;
    }
   }
  }
  //Rows
  k=Norb_act2;
  for(i=N2-1;i>-1;i--)
  {
   if(active_map[i])
   {
    k--;
    for(j=0;j<N2;j++)
    {
     DM2[i][j]=DM2[k][j];
     DM2[k][j]=ZERO;
    }
   }
  }
 }
 // Print the unformatted DM2 matrix and produce the spinless 1-RDM
 DM1=new double*[N1];
 for(i=0;i<N1;i++)
 {
  DM1[i]=new double[N1];
  for(j=0;j<N1;j++)
  {
   DM1[i][j]=ZERO;
  }
 }
 Trace=ZERO;
  //Print the stored 2-RDM
 nth=0;
 ofstream dm2_out(name_dm2.c_str(),ios::out | ios::binary);  
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
      if(element[1]==element_prime[1])
      {
       DM1[element[0]-1][element_prime[0]-1]=DM1[element[0]-1][element_prime[0]-1]+Dijkl;
       if(element[0]==element_prime[0])
       {
        Trace=Trace+Dijkl;
       }
      }
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
 for(i=0;i<N2;i++)
 {
  delete[] DM2[i];DM2[i]=NULL;
 } 
 delete[] DM2;DM2=NULL;
  //Build 1-RDM and for closed-shell systems the spin-with 1-RDM 
 if(N1!=Norb_act)
 {
  DM1spin=new double*[2*N1];
  for(i=0;i<2*N1;i++)
  {
   DM1spin[i]=new double[2*N1];
   for(j=0;j<2*N1;j++)
   {
    DM1spin[i][j]=ZERO;
   }
  }
 }
 Nelect=HALF*(ONE+sqrt(ONE+EIGHT*Trace));
 Trace1dm=ZERO;
 for(i=0;i<N1;i++)
 {
  for(j=0;j<N1;j++)
  {
   DM1[i][j]=TWO*DM1[i][j]/(Nelect-ONE); 
   if(i==j)
   {
    if(N1!=Norb_act)
    {
     if(i<Ndocc)
     {
      DM1[i][i]=TWO; // Prepared for closed-shell
     }
    }
    Trace1dm=Trace1dm+DM1[i][i]; 
   }
  }
  if(N1!=Norb_act)
  {
   for(j=0;j<N1;j++)
   {
    k=2*i;
    l=2*j;
    DM1spin[k][l]=HALF*DM1[i][j];
    k=2*i+1;
    l=2*j+1;
    DM1spin[k][l]=HALF*DM1[i][j];
   }
  }
 }
 if(N1!=Norb_act)
 {
  for(i=0;i<2*N1;i++)
  {
   for(j=0;j<2*N1;j++)
   {
    for(k=0;k<2*N1;k++)
    {
     for(l=0;l<2*N1;l++)
     {
      if(i%2==k%2 && j%2==l%2)
      {
       if(i<2*Ndocc || j<2*Ndocc || k<2*Ndocc || l<2*Ndocc) // Do HFL if at least one of them is fully occ (frozen with occ=1).
       {
        Dijkl=dm1_to_dm2_hfl(i,j,k,l,DM1spin);
        if(abs(Dijkl)>threshold)
        {
         if(i==k && j==l)
         {
          Trace=Trace+Dijkl;
         }
         if(i%2==0)
         {
          element[0]=i/2; 
          element_prime[0]=k/2; 
         }
         else
         {
          element[0]=(i-1)/2; 
          element_prime[0]=(k-1)/2; 
         }
         if(j%2==0)
         {
          element[1]=j/2; 
          element_prime[1]=l/2; 
         }
         else
         {
          element[1]=(j-1)/2; 
          element_prime[1]=(l-1)/2; 
         }
         element[0]=element[0]+1;
         element[1]=element[1]+1;
         element_prime[0]=element_prime[0]+1;
         element_prime[1]=element_prime[1]+1;
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
   }
  }
 }
 cout<<endl;
 cout<<endl;
 cout<<"Trace of the printed 2-RDM (Lowdin): "<<setprecision(8)<<fixed<<setw(17)<<Trace<<endl;
 cout<<"Trace of the contracted 1-RDM      : "<<setprecision(8)<<fixed<<setw(17)<<Trace1dm<<endl;
 cout<<endl;
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
 cout<<"File "<<name_dm2<<" prepared for chimpanC."<<endl;
 cout<<endl;
 cout<<"Note: Use the i8all option in chimpanC."<<endl;
 cout<<endl;
 // Deallocate all dynamic arrays
 if(N1!=Norb_act)
 {
  for(i=0;i<2*N1;i++)
  {
   delete[] DM1spin[i];DM1spin[i]=NULL;
  }
  delete[] DM1spin;DM1spin=NULL;
 }
 for(i=0;i<N1;i++)
 {
  delete[] DM1[i];DM1[i]=NULL;
 }
 delete[] DM1;DM1=NULL;
 delete[] order_orbs;order_orbs=NULL;
 delete[] docc_orbs;docc_orbs=NULL;
 delete[] uocc_orbs;uocc_orbs=NULL;
 delete[] active_map;active_map=NULL;
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

double dm1_to_dm2_hfl(int &i, int &j, int &k, int &l, double **dm1)
{
 double Dijkl=ZERO;
 Dijkl=dm1[i][k]*dm1[j][l]-dm1[i][l]*dm1[j][k];
 return HALF*Dijkl; // Lowdin normalization
}

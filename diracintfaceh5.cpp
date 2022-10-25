#include<iostream>
#include<iomanip>
#include<vector>
#include<cstring>
#include<cmath>
#include<fstream>
#include"hdf5.h"
#define zero 0.0e0
#define two  2.0e0
#define tol6 1.0e-6
#define tol8 1.0e-8

using namespace std;

void read_int(hid_t file_id,string path,string var, vector<int> *to_store);
void read_dou(hid_t file_id,string path,string var, vector<double> *to_store);
void read_str(hid_t file_id, string path, string var, vector<string> *to_store);
void asign(vector<int> *s2atom_map,vector<double> *Shell_coord_L,vector<double> *Shell_coord_S);
bool coor_atoms_shells(vector<double> *Shell_coord_L,vector<double> *Atom_pos);
void h5_to_temp(string h5_file);

int main(int argc, char *argv[]) 
{
 string h5_file;
 if(argc!=2)
 {
  cout<<"Include the name of the H5 file"<<endl;
  return -1;
 }
 else
 {
  h5_file=argv[1];
 }
 h5_to_temp(h5_file);  
 return 0;
}

void h5_to_temp(string h5_file)
{
 bool atom_shell_coord;
 int i,j,k,nbasis_tot,coefs_quat;
 double ang2au=1.8897259886;
 double **MO2AO_coefs_quat;
 hid_t file_id, dset_id, dspace_id; /* identifiers */
 herr_t status;
 vector<int> Nbasis_L,Shell_type_L,Nprim_shell_L,Nshells_L;
 vector<double> Shell_coord_L,Prim_exponents_L,Prim_coefs_L;
 vector<int> Nbasis_S,Shell_type_S,Nprim_shell_S,Nshells_S;
 vector<double> Shell_coord_S,Prim_exponents_S,Prim_coefs_S;
 vector<int> Nbasis_MO,Nz,s2atom_map;
 vector<double> MO_occs,MO_coefs,Nuc_charge,Atom_pos;
 string path;
 ofstream tmp_h5((h5_file.substr(0,h5_file.length()-3)+".tmp").c_str());
 tmp_h5<<endl;
 tmp_h5<<"******************"<<endl;
 tmp_h5<<"Open DIRAC H5 file"<<endl;
 tmp_h5<<"******************"<<endl;
 tmp_h5<<endl;

 file_id = H5Fopen(h5_file.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
 // Reading charges
  path="/input/molecule/";
  // Retrieve nuclear charges
  read_dou(file_id,path,"nuc_charge",&Nuc_charge);
  tmp_h5<<" Nuc charges:"<<setw(9)<<Nuc_charge.size()<<endl;
  j=0;tmp_h5<<" ";
  for(i=0;i<Nuc_charge.size();i++)
  {
   tmp_h5<<setprecision(8)<<fixed<<setw(24)<<Nuc_charge[i];
   j++;
   if(j==6){tmp_h5<<endl;j=0;tmp_h5<<" ";}
  }
  tmp_h5<<endl;
  // Retrieve nuclear charges
  read_dou(file_id,path,"geometry",&Atom_pos);
  tmp_h5<<" Atomic pos :"<<setw(9)<<Atom_pos.size()<<endl;
  j=0;tmp_h5<<" ";
  for(i=0;i<Atom_pos.size();i++)
  {
   tmp_h5<<setprecision(8)<<fixed<<setw(24)<<Atom_pos[i]*ang2au;
   j++;
   if(j==6){tmp_h5<<endl;j=0;tmp_h5<<" ";}
  }
  tmp_h5<<endl;
 // Reading Large component basis
  path="/input/aobasis/1/";
  // Retrieve Num. AO Large
  read_int(file_id,path,"n_ao",&Nbasis_L);
  tmp_h5<<" Nbasis  Large "<<setw(12)<<Nbasis_L[0]<<endl;
  // Retrieve Num. large shells
  read_int(file_id,path,"n_shells",&Nshells_L);
  tmp_h5<<" Nshells Large "<<setw(12)<<Nshells_L[0]<<endl;
  // Retrieve Shell_type_L
  read_int(file_id,path,"orbmom",&Shell_type_L);
  // Retrieve Nprim per shell L
  read_int(file_id,path,"n_prim",&Nprim_shell_L);
  // double 1D vector
  // Retrieve Shell_coord_L
  read_dou(file_id,path,"center",&Shell_coord_L);
  // Retrieve Prim_exponents_L
  read_dou(file_id,path,"exponents",&Prim_exponents_L);
  // Retrieve Prim_coefs_L
  read_dou(file_id,path,"contractions",&Prim_coefs_L);
 // Reading Small component basis
  path="/input/aobasis/2/";
  // Retrieve Num. AO Small
  read_int(file_id,path,"n_ao",&Nbasis_S);
  tmp_h5<<" Nbasis  Small "<<setw(12)<<Nbasis_S[0]<<endl;
  tmp_h5<<" Nbasis        "<<setw(12)<<Nbasis_L[0]+Nbasis_S[0]<<endl;
  // Retrieve Num. small shells
  read_int(file_id,path,"n_shells",&Nshells_S);
  tmp_h5<<" Nshells Small "<<setw(12)<<Nshells_S[0]<<endl;
  tmp_h5<<" Nshells       "<<setw(12)<<Nshells_L[0]+Nshells_S[0]<<endl;
  // Retrieve Shell_type_S
  read_int(file_id,path,"orbmom",&Shell_type_S);
  tmp_h5<<endl;
  tmp_h5<<" Shell type:"<<endl;
  j=0;tmp_h5<<" ";
  for(i=0;i<Shell_type_L.size();i++)
  {
   Shell_type_L[i]=Shell_type_L[i]-1;
   tmp_h5<<setw(10)<<Shell_type_L[i];
   j++;
   if(j==10){tmp_h5<<endl;j=0;tmp_h5<<" ";}
  }
  for(i=0;i<Shell_type_S.size();i++)
  {
   Shell_type_S[i]=Shell_type_S[i]-1;
   tmp_h5<<setw(10)<<Shell_type_S[i];
   j++;
   if(j==10){tmp_h5<<endl;j=0;tmp_h5<<" ";}
  }
  tmp_h5<<endl;
  // Retrieve Nprim per shell S
  read_dou(file_id,path,"center",&Shell_coord_S);
  tmp_h5<<endl;
  tmp_h5<<" Shell to atom map:"<<endl;
  asign(&s2atom_map,&Shell_coord_L,&Shell_coord_S);
  j=0;tmp_h5<<" ";
  for(i=0;i<s2atom_map.size();i++)
  {
   tmp_h5<<setw(10)<<s2atom_map[i];
   j++;
   if(j==10){tmp_h5<<endl;j=0;tmp_h5<<" ";}
  }
  tmp_h5<<endl;
  // Retrieve Shell_coord_S
  tmp_h5<<endl;
  tmp_h5<<" Number of primitives per shell:"<<endl;
  read_int(file_id,path,"n_prim",&Nprim_shell_S);
  j=0;tmp_h5<<" ";
  for(i=0;i<Nprim_shell_L.size();i++)
  {
   tmp_h5<<setw(10)<<Nprim_shell_L[i];
   j++;
   if(j==10){tmp_h5<<endl;j=0;tmp_h5<<" ";}
  }
  for(i=0;i<Nprim_shell_S.size();i++)
  {
   tmp_h5<<setw(10)<<Nprim_shell_S[i];
   j++;
   if(j==10){tmp_h5<<endl;j=0;tmp_h5<<" ";}
  }
  tmp_h5<<endl;
  // Retrieve Prim_coefs_S
  read_dou(file_id,path,"contractions",&Prim_coefs_S);
  tmp_h5<<endl;
  tmp_h5<<" Contraction coefficients:"<<endl;
  j=0;tmp_h5<<" ";
  for(i=0;i<Prim_coefs_L.size();i++)
  {
   tmp_h5<<setprecision(8)<<fixed<<setw(24)<<Prim_coefs_L[i];
   j++;
   if(j==6){tmp_h5<<endl;j=0;tmp_h5<<" ";}
  }
  for(i=0;i<Prim_coefs_S.size();i++)
  {
   tmp_h5<<setprecision(8)<<fixed<<setw(24)<<Prim_coefs_S[i];
   j++;
   if(j==6){tmp_h5<<endl;j=0;tmp_h5<<" ";}
  }
  tmp_h5<<endl;
  // Retrieve Prim_exponents_S
  read_dou(file_id,path,"exponents",&Prim_exponents_S);
  tmp_h5<<endl;
  tmp_h5<<" Primitive exponents:"<<endl;
  j=0;tmp_h5<<" ";
  for(i=0;i<Prim_exponents_L.size();i++)
  {
   tmp_h5<<setprecision(8)<<fixed<<setw(24)<<Prim_exponents_L[i];
   j++;
   if(j==6){tmp_h5<<endl;j=0;tmp_h5<<" ";}
  }
  for(i=0;i<Prim_exponents_S.size();i++)
  {
   tmp_h5<<setprecision(8)<<fixed<<setw(24)<<Prim_exponents_S[i];
   j++;
   if(j==6){tmp_h5<<endl;j=0;tmp_h5<<" ";}
  }
  tmp_h5<<endl;
  // Print shell cartes coordinates
  tmp_h5<<endl;
  tmp_h5<<" Cartesian coordinates of each shell:"<<endl;
  j=0;tmp_h5<<" ";
  for(i=0;i<Shell_coord_L.size();i++)
  {
   tmp_h5<<setprecision(8)<<fixed<<setw(24)<<Shell_coord_L[i];
   j++;
   if(j==6){tmp_h5<<endl;j=0;tmp_h5<<" ";}
  }
  for(i=0;i<Shell_coord_S.size();i++)
  {
   tmp_h5<<setprecision(8)<<fixed<<setw(24)<<Shell_coord_S[i];
   j++;
   if(j==6){tmp_h5<<endl;j=0;tmp_h5<<" ";}
  }
  tmp_h5<<endl;
 // Read results information
  path="/result/wavefunctions/scf/mobasis/";
  // Retrieve MO_occs
  read_dou(file_id,path,"occupations",&MO_occs);
  tmp_h5<<" Occupancies:"<<setw(9)<<MO_occs.size()<<endl;
  j=0;tmp_h5<<" ";
  for(i=0;i<MO_occs.size();i++)
  {
   tmp_h5<<setprecision(8)<<fixed<<setw(24)<<MO_occs[i];
   j++;
   if(j==6){tmp_h5<<endl;j=0;tmp_h5<<" ";}
  }
  tmp_h5<<endl;
  // Retrieve Num. MO
  read_int(file_id,path,"n_mo",&Nbasis_MO);
  tmp_h5<<endl;
  tmp_h5<<" Norbs         "<<setw(12)<<Nbasis_MO[0]<<endl;
  // Retrieve MO_coefs
  read_dou(file_id,path,"orbitals",&MO_coefs);
  coefs_quat=MO_coefs.size()/4;
  nbasis_tot=coefs_quat/Nbasis_MO[0];
  tmp_h5<<" NBasis        "<<setw(12)<<nbasis_tot<<endl;
  // Retrieve Nz
  read_int(file_id,path,"nz",&Nz);
  tmp_h5<<" Nz            "<<setw(12)<<Nz[0]<<endl;
  tmp_h5<<endl;
  MO2AO_coefs_quat=new double*[coefs_quat];
  for(i=0;i<coefs_quat;i++)
  {
   MO2AO_coefs_quat[i]=new double[4];
  } 
  tmp_h5<<endl;
  k=0;
  for(i=0;i<4;i++)
  {
   for(j=0;j<coefs_quat;j++)
   {
    MO2AO_coefs_quat[j][i]=MO_coefs[k];
    k++;
   }
  }
  for(i=0;i<Nbasis_MO[0];i++)
  {
   for(j=0;j<nbasis_tot;j++)
   {
    for(k=0;k<4;k++)
    {
     tmp_h5<<setprecision(8)<<fixed<<setw(20)<<MO2AO_coefs_quat[j+i*nbasis_tot][k];
    }
    tmp_h5<<endl;
   }
   tmp_h5<<endl;
  }
  tmp_h5<<endl;
  for(i=0;i<coefs_quat;i++)
  {
   delete[] MO2AO_coefs_quat[i];MO2AO_coefs_quat[i]=NULL;
  } 
  delete[] MO2AO_coefs_quat;MO2AO_coefs_quat=NULL;

 // Close the file
 status = H5Fclose(file_id);
 tmp_h5.close();
 atom_shell_coord=coor_atoms_shells(&Shell_coord_L,&Atom_pos);
 if(!atom_shell_coord)
 {
  cout<<"Warning! Atomic coordinates and shell coordinates do not coincide"<<endl;
 }
}

void read_int(hid_t file_id, string path, string var, vector<int> *to_store)
{
 int idims,ndims,res_sz;
 hid_t dset_id, dspace_id;
 herr_t status;
 string string_dset=path+var; 
 dset_id = H5Dopen(file_id,string_dset.c_str(),H5P_DEFAULT);
 dspace_id = H5Dget_space(dset_id);
 ndims = H5Sget_simple_extent_ndims(dspace_id);
 hsize_t res_dims[ndims];
 status = H5Sget_simple_extent_dims(dspace_id, res_dims, NULL);
 res_sz = 1;
 for(idims=0;idims<ndims;idims++)
 {
  res_sz *= res_dims[idims];
 }
 to_store[0].resize(res_sz);
 status = H5Dread(dset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, to_store[0].data());
 status = H5Dclose(dset_id);
}

void read_dou(hid_t file_id, string path, string var, vector<double> *to_store)
{
 int idims,ndims,res_sz;
 hid_t dset_id, dspace_id;
 herr_t status;
 string string_dset=path+var; 
 dset_id = H5Dopen(file_id,string_dset.c_str(),H5P_DEFAULT);
 dspace_id = H5Dget_space(dset_id);
 ndims = H5Sget_simple_extent_ndims(dspace_id);
 hsize_t res_dims[ndims];
 status = H5Sget_simple_extent_dims(dspace_id, res_dims, NULL);
 res_sz = 1;
 for(idims=0;idims<ndims;idims++)
 {
  res_sz *= res_dims[idims];
 }
 to_store[0].resize(res_sz);
 status = H5Dread(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, to_store[0].data());
 status = H5Dclose(dset_id);
}

void read_str(hid_t file_id, string path, string var, vector<string> *to_store)
{
 int idims,ndims,res_sz;
 vector<int>first_read;
 hid_t dset_id, dspace_id, h5t_string;
 herr_t status;
 string string_dset=path+var; 
 dset_id = H5Dopen(file_id,string_dset.c_str(),H5P_DEFAULT);
 dspace_id = H5Dget_space(dset_id);
 ndims = H5Sget_simple_extent_ndims(dspace_id);
 hsize_t res_dims[ndims];
 status = H5Sget_simple_extent_dims(dspace_id, res_dims, NULL);
 res_sz = 1;
 for(idims=0;idims<ndims;idims++)
 {
  res_sz *= res_dims[idims];
 }
 to_store[0].resize(res_sz);
 first_read.resize(res_sz);
 h5t_string=H5Dget_type(dset_id);
 status = H5Dread(dset_id, h5t_string, H5S_ALL, H5S_ALL, H5P_DEFAULT, first_read.data());
 for(idims=0;idims<first_read.size();idims++)
 {
  cout<<first_read[idims]<<endl;
 }

 status = H5Dclose(dset_id);
}

void asign(vector<int> *s2atom_map,vector<double> *Shell_coord_L,vector<double> *Shell_coord_S)
{
 int i,j,k,ind,nshells=(Shell_coord_L[0].size()+Shell_coord_S[0].size())/3;
 bool newcent;
 double diff;
 struct NewCenter
 {
  int index;
  double Coord[3];
 };
 vector<NewCenter>newcenter;
 NewCenter centertemp;
 ind=1;
 centertemp.index=ind;
 centertemp.Coord[0]=Shell_coord_L[0].at(0);
 centertemp.Coord[1]=Shell_coord_L[0].at(1);
 centertemp.Coord[2]=Shell_coord_L[0].at(2);
 newcenter.push_back(centertemp);
 ind++;
 for(i=1;i<Shell_coord_L[0].size()/3;i++)
 {
  newcent=true;
  for(j=0;j<newcenter.size();j++)
  {
   diff=zero;
   for(k=0;k<3;k++)
   {
    diff=diff+pow(newcenter[j].Coord[k]-Shell_coord_L[0].at(3*i+k),two);
   }
   if(sqrt(diff)<tol8)
   {
    newcent=false;j=newcenter.size();
   }
  }
  if(newcent)
  {
   centertemp.index=ind;
   for(k=0;k<3;k++)
   {
    centertemp.Coord[k]=Shell_coord_L[0].at(3*i+k);
   }
   newcenter.push_back(centertemp);
   ind++;
  } 
 }
 for(i=0;i<Shell_coord_S[0].size()/3;i++)
 {
  newcent=true;
  for(j=0;j<newcenter.size();j++)
  {
   diff=zero;
   for(k=0;k<3;k++)
   {
    diff=diff+pow(newcenter[j].Coord[k]-Shell_coord_S[0].at(3*i+k),two);
   }
   if(sqrt(diff)<tol8)
   {
    newcent=false;j=newcenter.size();
   }
  }
  if(newcent)
  {
   cout<<"Warning! New center found for Small component w.r.t. Large component"<<endl;
   centertemp.index=ind;
   for(k=0;k<3;k++)
   {
    centertemp.Coord[k]=Shell_coord_S[0].at(3*i+k);
   }
   newcenter.push_back(centertemp);
   ind++;
  } 
 }
 for(i=0;i<Shell_coord_L[0].size()/3;i++)
 {
  for(j=0;j<newcenter.size();j++)
  {
   diff=zero;
   for(k=0;k<3;k++)
   {
    diff=diff+pow(newcenter[j].Coord[k]-Shell_coord_L[0].at(3*i+k),two);
   }
   if(sqrt(diff)<tol8)
   {
    s2atom_map[0].push_back(newcenter[j].index);
    j=newcenter.size();
   }
  }
 }
 for(i=0;i<Shell_coord_S[0].size()/3;i++)
 {
  for(j=0;j<newcenter.size();j++)
  {
   diff=zero;
   for(k=0;k<3;k++)
   {
    diff=diff+pow(newcenter[j].Coord[k]-Shell_coord_S[0].at(3*i+k),two);
   }
   if(sqrt(diff)<tol8)
   {
    s2atom_map[0].push_back(newcenter[j].index);
    j=newcenter.size();
   }
  }
 }
}

bool coor_atoms_shells(vector<double> *Shell_coord_L,vector<double> *Atom_pos)
{
 bool atshell_coord=true,found=true;
 int i,j,k;
 double diff,ang2au=1.8897259886;
 for(i=0;i<Atom_pos[0].size()/3;i++)
 {
  found=false;
  for(j=0;j<Shell_coord_L[0].size()/3;j++)
  {
   diff=zero;
   for(k=0;k<3;k++)
   {
    diff=diff+pow(ang2au*Atom_pos[0].at(3*i+k)-Shell_coord_L[0].at(3*j+k),two);
   } 
   if(sqrt(diff)<tol6)
   {
    found=true;
    j=Shell_coord_L[0].size();
   }
  }
  if(!found) 
  {
   atshell_coord=false;
   i=Atom_pos[0].size();
  }
 }  
 return atshell_coord;
}

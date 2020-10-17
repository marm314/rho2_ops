#include<iostream>
#include<fstream>
#include<omp.h>
#include"sphere_lebedev_rule.h"
#include"gauss_quad.h"
#include"legendre_quadrature.h"
#include"Input_commands.h"
#include"Numbers.h"
#include"Mathematical_Functions.h"
#ifdef HAVE_MPI
#include<mpi.h>
#endif

using namespace std;

void initialize(double Rmax[5][5]);
void grid_avail(int &Order);
void user_defined_angular(int &nang,int &nang2,double *x,double *y,double *z,double *w_theta_phi,string name_basis);
double normGauss(double &expon, int &n, int &l,int &m);
void terms_dm2(string name_file);
void fill_in_dm2(string name_file);
void calc_time(double DATE[2][4]);
//double Vint(double &S, int &order, double *s, double *w,double &Sijkl, double &alpha_ijkl,double &zeta_ijkl,
//double &Si,double &Sj,double &Sk, double &Sl, int &nxi,int &nxj,int &nxk,int &nxl);
double Sab(int &order, double *s, double *w,double &S_ab, double &Sa,double &Sb,double &exp_a,double &exp_b, int &nxa,int &nxb);
double Mab(int &order, double *s, double *w,double &S_ab, double &Sa,double &Sb,double &exp_a,double &exp_b, int &nxa,int &nxb);
double M2ab(int &order, double *s, double *w,double &S_ab, double &Sa,double &Sb,double &exp_a,double &exp_b, int &nxa,int &nxb);
double Vint3D(double S_intra[3], int &order, double *s, double *w,double Sijkl[3], double &alpha_ijkl,double &zeta_ijkl,
double Primitive_coords[4][3], int nsx_ijkl[4],int nsy_ijkl[4],int nsz_ijkl[4]);
double Vint3D_2(double S_extra[3], int &order, double *s, double *w,double Sijkl[3], double &alpha_ijkl,double &zeta_ijkl,
double Primitive_coords[4][3], int nsx_ijkl[4],int nsy_ijkl[4],int nsz_ijkl[4]);
double Jab(double &expa, double &expb,int &lxa,int &lya,int &lza,int &lxb,
int &lyb,int &lzb,double &Xa,double &Ya,double &Za,double &Xb,double &Yb,double &Zb);

//Global variables
long int nterms;
int ID=0,nproc=1;
double threshold,Rmax[5][5]={ZERO},DATE[2][4],localIr,globalIr;
const int RECORD_DELIMITER_LENGTH=4;
ifstream date_file;
struct basis_info
{
 int Z;
 double Cart_cord[3];
 double Expon;
 bool newterm;
 int nx,ny,nz;
};
struct DM2
{
 int indexes[4];
 double Dijkl;
};
DM2 *dm2;

int main(int argc, char *argv[])
{
#ifdef HAVE_MPI
MPI::Init();
int id=MPI::COMM_WORLD.Get_rank();
ID=id;
int npr[1];
MPI_Comm_size(MPI_COMM_WORLD,npr);
nproc=npr[0];
#endif
if(ID==0)
{
 cout<<"##########################################################################"<<endl;
 cout<<"##########################################################################"<<endl;
 cout<<"# Evaluation of the Intracule, Extracule and Second Moments using the    #"<<endl;
 cout<<"# DM2 in Primitives (in C++) written by M.Sc. Mauricio Rodriguez Mayorga #"<<endl;
 cout<<"# email: marm3.14@gmail.com                                              #"<<endl;
 cout<<"##########################################################################"<<endl;
 cout<<"##########################################################################"<<endl;
 cout<<"#*****************************************************************************#";
 cout<<endl;
 cout<<"# Copyright (C) 2016 M.Sc. Mauricio A. Rodriguez Mayorga                      #";
 cout<<endl;
 cout<<"# Ph.D. student at University of Girona (Girona) and DIPC                     #";
 cout<<endl;
 cout<<"# for support and comments send an email to: marm3.14@gmail.com               #";
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
 cout<<"# See. J. Cioslowski and G. Liu, J. Chem. Phys., 105, 4151 (1996)             #";
 cout<<endl;
 cout<<"#*****************************************************************************#";
 cout<<endl;
 /////////////////////
 //Begin clock      //
 /////////////////////
 system("date +%j' '%H' '%M' '%S > date_intrac.date");
 date_file.open("date_intrac.date");
 date_file>>DATE[0][0]>>DATE[0][1]>>DATE[0][2]>>DATE[0][3];
 date_file.close();
 system("rm date_intrac.date");
}
 if(argc==2)
 {
  int i,j,k,nprims,Lmax,Nroot_2Lmax_plus_1,Nroot_Lmax_plus_2,nrad,order_ijkl,max_exp_ijkl,nx_sum,ny_sum,nz_sum,nang,nang2,counter=0;
  int nx_exp[4],ny_exp[4],nz_exp[4];
  double Integral,Integral2,Init,Step,Last,rscan,DMNfact,alpha,a,b,zeta_ik,zeta_jl,zeta_ijkl,e_ijkl,alpha_ijkl,symmetry_omitted_terms=ZERO;
  double Nelectrons,threshold2,Jik,Jjl,Cnorm_4gauss,Cnorm_ang,Aijkl,Xprime,Yprime,Zprime,lambda_rs,lambda_scr;
  double TPS[3][3],Coord_atom[4][3],res[16],Rik[3],Rjl[3],Rijkl[3],Sik[3],Sjl[3],Mik[3],Mjl[3],M2ik[3],Point_intra[3]={ZERO},Point_extra[3]={ZERO};
  double **r_intrac,**w_intrac,**r_extrac,**w_extrac,*r_mom,*w_mom,**Tot_rad,*x,*y,*z,*w_theta_phi,*r_legendre,*w_legendre,**Jacobian_legendre,*Shannon_rad;
  bool legendre=false,parallel=false;
  Step=ZERO;
  string aux(argv[1]);
  Input Input_commands(aux);
  /////////////////////
  //Compute Intracule//
  /////////////////////
  if(Input_commands.intracule)
  {
if(ID==0)
{
   cout<<"#******************************************#"<<endl;
   cout<<"#          Compute the Intracule           #"<<endl;
   cout<<"#******************************************#"<<endl;
}
   Init=Input_commands.radial_Init;
   Step=Input_commands.radial_Step;
   Last=Input_commands.radial_Last;
   lambda_rs=Input_commands.lambda_rs;
   lambda_scr=Input_commands.lambda_scr;
   string name_dm2=Input_commands.name_dm2;
   string name_basis=Input_commands.name_basis;
   threshold=Input_commands.threshold_in;
   parallel=Input_commands.parallel;
   ////////////////
   //Check times //
   ////////////////
if(ID==0)
{
   if(Input_commands.time_intra){cout<<endl;cout<<"Time 1 (D H M S): "<<endl;system("date +%j' '%H' '%M' '%S ");cout<<endl;}
}
   ////////////////
   //End time    //
   ////////////////
if(ID==0)
{
   if(parallel)
   {
    cout<<"OpenMP parallelized evaluations of the intracule will be done."<<endl;
   }
   else
   {
    cout<<"Serial evaluations of the intracule will be done."<<endl;
   }
   cout<<"Threshold for DM2 terms:"<<setprecision(10)<<fixed<<scientific<<setw(18)<<Input_commands.threshold_in<<endl;
}
   //We define the DMNfactor to correct from McWeeney Normalization [N(N-1)/2]2! to Lowdin N(N-1)/2 for the DMN matrix.
   DMNfact=HALF;
   //Create angular grid
   nang=Input_commands.order_ang;
   nang2=Input_commands.order_ang2;
   if(nang>0 && nang2==0)
   {
    grid_avail(nang);
    Cnorm_ang=TWO;
   }
   else if(nang>0 && nang2>0)
   {
    nang=nang*nang2;    
    Cnorm_ang=ONE;
   }
   else
   {
if(ID==0)
{
    cout<<"None angular grid chosen for spherical average"<<setw(17)<<endl;
}
    nang=1;
    Cnorm_ang=ONE;
   }
   x=new double[nang];
   y=new double[nang];
   z=new double[nang];
   w_theta_phi=new double[nang];
   if(nang>1)
   {
    for(i=0;i<nang;i++)
    {
     x[i]=ZERO;
     y[i]=ZERO;
     z[i]=ZERO;
     w_theta_phi[i]=ZERO;
    }
    if(nang2==0)
    {
     ld_by_order(nang,x,y,z,w_theta_phi);
    }
    else
    {
     user_defined_angular(Input_commands.order_ang,Input_commands.order_ang2,x,y,z,w_theta_phi,Input_commands.name_basis);
    }
   }
   else
   {
    // We allow to scan only in z12 axis!
    x[0]=ZERO;
    y[0]=ZERO;
    z[0]=ONE;
    w_theta_phi[0]=ONE;
   }
   //Create radial grid
   legendre=Input_commands.legendre;
   if(!legendre)
   {
    nrad=0;
    for(rscan=Init;rscan<Last+Step;rscan=rscan+Step){nrad++;}
    Tot_rad=new double*[nrad];
    Shannon_rad=new double[nrad];
    for(i=0;i<nrad;i++)
    {
     Shannon_rad[i]=ZERO;
     Tot_rad[i]=new double[nang];
     // Initialize I(R) for scheme 3 in the paper
     for(j=0;j<nang;j++)
     {
      Tot_rad[i][j]=ZERO;
     }
    }
   }
   else
   {
    nrad=Input_commands.radial_grid;
    Tot_rad=new double*[nrad];
    Shannon_rad=new double[nrad];
    for(i=0;i<nrad;i++)
    {
     Shannon_rad[i]=ZERO;
     Tot_rad[i]=new double[nang];
     // Initialize I(R) for scheme 3 in the paper
     for(j=0;j<nang;j++)
     {
      Tot_rad[i][j]=ZERO;
     }
    }
    r_legendre=new double[nrad];
    w_legendre=new double[nrad];
    Jacobian_legendre=new double*[nrad];
    for(i=0;i<nrad;i++)
    {
     Jacobian_legendre[i]=new double[5];
    }
    Init=ZERO;
    if(Input_commands.radial_Last!=1e99)
    {
     Last=Input_commands.radial_Last;
    }
    else
    {
     Last=ONE;
    }
#ifdef HAVE_MPI
MPI_Barrier(MPI_COMM_WORLD);
#endif
if(ID==0)
{
    legendre_quadrature(name_basis.substr(0,name_basis.length()-6),nrad,Init,Last);
    //Read quadrature info
    ifstream read_rad_quad;
    // Read weights
    read_rad_quad.open((name_basis.substr(0,name_basis.length()-6)+"_w.txt").c_str());
    for(i=0;i<nrad;i++)
    {
     read_rad_quad>>w_legendre[i];
    }
    read_rad_quad.close();
    // Read roots
    read_rad_quad.open((name_basis.substr(0,name_basis.length()-6)+"_x.txt").c_str());
    for(i=0;i<nrad;i++)
    {
     read_rad_quad>>r_legendre[i];
    }
    read_rad_quad.close();
    // Remove quadrature files
    system(("rm "+name_basis.substr(0,name_basis.length()-6)+"_r.txt").c_str());
    system(("rm "+name_basis.substr(0,name_basis.length()-6)+"_w.txt").c_str());
    system(("rm "+name_basis.substr(0,name_basis.length()-6)+"_x.txt").c_str());
}
#ifdef HAVE_MPI
MPI_Bcast(w_legendre,nrad,MPI_DOUBLE,0,MPI_COMM_WORLD);
MPI_Bcast(r_legendre,nrad,MPI_DOUBLE,0,MPI_COMM_WORLD);
#endif
    if(Input_commands.radial_Last==1e99)
    {
if(ID==0)
{
     cout<<"Superior lim. Legendre :         Infinity"<<endl;
}
     for(i=0;i<nrad;i++)
     {
      for(j=0;j<5;j++)
      {
       Jacobian_legendre[i][j]=pow(r_legendre[i]/(ONE-r_legendre[i]),(double)j)/pow(ONE-r_legendre[i],TWO);
      }
      //We convert the variable r_legendre[i] from 0->1 to 0->Infinity (to be used later, where real [not from 0 to 1] evals of intracule must be performed!)
      r_legendre[i]=r_legendre[i]/(ONE-r_legendre[i]);
     }
    }
    else
    {
if(ID==0)
{
     cout<<"Superior lim. Legendre :"<<setprecision(10)<<fixed<<scientific<<setw(18)<<Last<<endl;
}
     for(i=0;i<nrad;i++)
     {
      for(j=0;j<5;j++)
      {
       Jacobian_legendre[i][j]=pow(r_legendre[i],(double)j);
      }
     }
    }
   }
   // Read basis info
   ifstream open_basis;
   open_basis.open((name_basis).c_str());
   if(open_basis.good())
   {
    nprims=0;
    while(getline(open_basis,aux))
    {
     nprims++;
    }
    open_basis.close();
    nprims=nprims-1;
    basis_info *BASIS_INFO=new basis_info[nprims];
    open_basis.open((name_basis).c_str());
    getline(open_basis,aux);
    for(i=0;i<nprims;i++)
    {
     open_basis>>BASIS_INFO[i].Z>>BASIS_INFO[i].Cart_cord[0]>>BASIS_INFO[i].Cart_cord[1]>>BASIS_INFO[i].Cart_cord[2];
     open_basis>>BASIS_INFO[i].Expon>>BASIS_INFO[i].nx>>BASIS_INFO[i].ny>>BASIS_INFO[i].nz;
    }
    open_basis.close();
    //We use the number of primitives to compute threshold 2 and we also initialize Rmax matrix here
if(ID==0)
{
    cout<<"tau                    :"<<setprecision(10)<<fixed<<scientific<<setw(18)<<Input_commands.tau<<endl;
}
    threshold2=TWO*Input_commands.tau/((double)nprims*((double)nprims+ONE));
if(ID==0)
{
    cout<<"2*tau/[Np(Np+1)]       :"<<setprecision(10)<<fixed<<scientific<<setw(18)<<threshold2<<endl;
    cout<<"[Np = Number of primitives]"<<endl;
}
    initialize(Rmax);
    //Compute Nroot_2Lmax_plus_1
    Lmax=0;
    for(i=0;i<nprims;i++)
    {
     if(BASIS_INFO[i].nx>Lmax)
     {
      Lmax=BASIS_INFO[i].nx;
     }
    }
if(ID==0)
{
    cout<<"Lmax                   :"<<setw(18)<<Lmax<<endl;
}
    //Create quadrature rule order
    Nroot_2Lmax_plus_1=2*Lmax+1;
if(ID==0)
{
    cout<<"Quadrature rule order  :"<<setw(18)<<Nroot_2Lmax_plus_1<<endl;
    cout<<"[Quadrature for primitive integrals. Order = 2 Lmax + 1 ]"<<endl;
}
    alpha=ZERO;
    a=ZERO;
    b=ONE;
    r_intrac=new double*[Nroot_2Lmax_plus_1];
    w_intrac=new double*[Nroot_2Lmax_plus_1];
    for(i=0;i<Nroot_2Lmax_plus_1;i++)
    {
     r_intrac[i]=new double[i+1];
     w_intrac[i]=new double[i+1];
     for(j=0;j<i+1;j++)
     {
      r_intrac[i][j]=ZERO;
      w_intrac[i][j]=ZERO;
     }
    }
    //Read quadrature info
    ifstream read_quad_r,read_quad_w;
    for(i=0;i<Nroot_2Lmax_plus_1;i++)
    {
     j=i+1;
#ifdef HAVE_MPI
double *wtmp,*rtmp;
wtmp=new double[j];
rtmp=new double[j];
MPI_Barrier(MPI_COMM_WORLD);
#endif
if(ID==0)
{
     gauss_hermite_rule(name_basis.substr(0,name_basis.length()-6),alpha,a,b,j);
     // Read weights
     read_quad_w.open((name_basis.substr(0,name_basis.length()-6)+"_w.txt").c_str());
     // Read roots
     read_quad_r.open((name_basis.substr(0,name_basis.length()-6)+"_x.txt").c_str());
     for(j=0;j<i+1;j++)
     {
      read_quad_w>>w_intrac[i][j];
      read_quad_r>>r_intrac[i][j];
#ifdef HAVE_MPI
wtmp[j]=w_intrac[i][j];
rtmp[j]=r_intrac[i][j];
#endif
     }
     read_quad_r.close();
     read_quad_w.close();
     // Remove quadrature files
     system(("rm "+name_basis.substr(0,name_basis.length()-6)+"_r.txt").c_str());
     system(("rm "+name_basis.substr(0,name_basis.length()-6)+"_w.txt").c_str());
     system(("rm "+name_basis.substr(0,name_basis.length()-6)+"_x.txt").c_str());
}
#ifdef HAVE_MPI
MPI_Bcast(wtmp,j,MPI_DOUBLE,0,MPI_COMM_WORLD);
MPI_Bcast(rtmp,j,MPI_DOUBLE,0,MPI_COMM_WORLD);
for(j=0;j<i+1;j++)
{
w_intrac[i][j]=wtmp[j];
r_intrac[i][j]=rtmp[j];
}
delete[] wtmp;wtmp=NULL;
delete[] rtmp;rtmp=NULL;
#endif
    }
    //Read DM2  file
    ifstream open_dm2;
    open_dm2.open((name_dm2).c_str());
    if(open_dm2.good())
    {
     open_dm2.close();
     nterms=0;
     //Check number of terms
     terms_dm2(name_dm2);
#ifdef HAVE_MPI
if(ID==0)
{
cout<<"Reducing the number of DM2 terms to "<<setw(8)<<nterms/nproc+1<<" per proc."<<endl;
}
nterms=nterms/nproc+1;
#endif
     dm2=new DM2[nterms];
     //Store the DM2
if(ID==0)
{
     cout<<"Storage of the DM2"<<endl;
}
     fill_in_dm2(name_dm2);
if(ID==0)
{
     cout<<"Storage done!"<<endl;
     //Number of terms to compute
     cout<<"Total intracule points : "<<setw(17)<<nang*nrad<<endl;
     cout<<"[Total intracule points = Order angular grid x Num radial points |u|.]"<<endl;
     cout<<"Intracule points x DM2 : "<<setw(17)<<(double)nterms*(double)nang*(double)nrad<<endl;
     cout<<"[Intracule points x DM2 terms = Nterms DM2 x Order angular grid x Num radial points |u|.]"<<endl;
     cout<<"Note: If no angular grid was chosen, Order angular grid is equal to ONE."<<endl;
     //Evaluate intracule
     ////////////////
     //Check times //
     ////////////////
     if(Input_commands.time_intra){cout<<endl;cout<<"Time 2 (D H M S): "<<endl;system("date +%j' '%H' '%M' '%S ");cout<<endl;}
     ////////////////
     //End time    //
     ////////////////
     cout<<"#*****************************************************************************#";
     cout<<endl;
}
     if(parallel)
     {
      if(Input_commands.nthreads>omp_get_max_threads())
      {Input_commands.nthreads=omp_get_max_threads();}
if(ID==0)
{
      cout<<"Running"<<setw(4)<<Input_commands.nthreads<<" threads."<<endl;
}
      #pragma omp parallel num_threads(Input_commands.nthreads) \
       private(i,j,k,Point_intra,Xprime,Yprime,Zprime,Jik,Jjl,Cnorm_4gauss,Coord_atom, \
       nx_sum,ny_sum,nz_sum,order_ijkl,max_exp_ijkl,nx_exp,ny_exp,nz_exp,zeta_ik,zeta_jl, \
       zeta_ijkl,Aijkl,e_ijkl,Rik,Rjl,Rijkl,alpha_ijkl,rscan,Integral,Integral2) \
       shared(BASIS_INFO,dm2,r_intrac,w_intrac,r_legendre,x,y,z,w_theta_phi,Tot_rad,Shannon_rad,nterms,nrad, \
       nang,nang2) \
       reduction(+:symmetry_omitted_terms,counter)
      {
       //Variables that belong to each thread and are only needed here
       int nth=omp_get_num_threads();
       int ith=omp_get_thread_num();
       double **Tot_rad_th;
       Tot_rad_th=new double*[nrad];
       for(i=0;i<nrad;i++)
       {
        Tot_rad_th[i]=new double[nang];
        for(j=0;j<nang;j++)
        {
         Tot_rad_th[i][j]=ZERO;
        }
       }
       for(i=ith;i<nterms;i=i+nth)
       {
        Cnorm_4gauss=ONE;
        for(j=0;j<4;j++)
        {
         //Compute the normalization factor of the gaussians (ijkl)
         Cnorm_4gauss=Cnorm_4gauss*
                     normGauss(BASIS_INFO[dm2[i].indexes[j]].Expon,
                               BASIS_INFO[dm2[i].indexes[j]].nx,
                               BASIS_INFO[dm2[i].indexes[j]].ny,
                               BASIS_INFO[dm2[i].indexes[j]].nz);
         //Take the coord. of this 4 gaussians
         Coord_atom[j][0]=BASIS_INFO[dm2[i].indexes[j]].Cart_cord[0];
         Coord_atom[j][1]=BASIS_INFO[dm2[i].indexes[j]].Cart_cord[1];
         Coord_atom[j][2]=BASIS_INFO[dm2[i].indexes[j]].Cart_cord[2];
         // Take the quantum numbers of the primitives
         nx_exp[j]=BASIS_INFO[dm2[i].indexes[j]].nx;
         ny_exp[j]=BASIS_INFO[dm2[i].indexes[j]].ny;
         nz_exp[j]=BASIS_INFO[dm2[i].indexes[j]].nz;
        }
        //Jab is the upper bound fucntion of Eq. 36
        Jik=Jab(BASIS_INFO[dm2[i].indexes[0]].Expon,
                BASIS_INFO[dm2[i].indexes[2]].Expon,
                BASIS_INFO[dm2[i].indexes[0]].nx,
                BASIS_INFO[dm2[i].indexes[0]].ny,
                BASIS_INFO[dm2[i].indexes[0]].nz,
                BASIS_INFO[dm2[i].indexes[2]].nx,
                BASIS_INFO[dm2[i].indexes[2]].ny,
                BASIS_INFO[dm2[i].indexes[2]].nz,
                BASIS_INFO[dm2[i].indexes[0]].Cart_cord[0],
                BASIS_INFO[dm2[i].indexes[0]].Cart_cord[1],
                BASIS_INFO[dm2[i].indexes[0]].Cart_cord[2],
                BASIS_INFO[dm2[i].indexes[2]].Cart_cord[0],
                BASIS_INFO[dm2[i].indexes[2]].Cart_cord[1],
                BASIS_INFO[dm2[i].indexes[2]].Cart_cord[2]);
        Jjl=Jab(BASIS_INFO[dm2[i].indexes[1]].Expon,
                BASIS_INFO[dm2[i].indexes[3]].Expon,
                BASIS_INFO[dm2[i].indexes[1]].nx,
                BASIS_INFO[dm2[i].indexes[1]].ny,
                BASIS_INFO[dm2[i].indexes[1]].nz,
                BASIS_INFO[dm2[i].indexes[3]].nx,
                BASIS_INFO[dm2[i].indexes[3]].ny,
                BASIS_INFO[dm2[i].indexes[3]].nz,
                BASIS_INFO[dm2[i].indexes[1]].Cart_cord[0],
                BASIS_INFO[dm2[i].indexes[1]].Cart_cord[1],
                BASIS_INFO[dm2[i].indexes[1]].Cart_cord[2],
                BASIS_INFO[dm2[i].indexes[3]].Cart_cord[0],
                BASIS_INFO[dm2[i].indexes[3]].Cart_cord[1],
                BASIS_INFO[dm2[i].indexes[3]].Cart_cord[2]);
        if((abs(dm2[i].Dijkl)*Cnorm_4gauss*pow(Jik*Jjl,HALF))>=threshold2)
        {
         counter++;
         //Non intracule coord. dependent part
         zeta_ik=BASIS_INFO[dm2[i].indexes[0]].Expon+BASIS_INFO[dm2[i].indexes[2]].Expon;
         zeta_jl=BASIS_INFO[dm2[i].indexes[1]].Expon+BASIS_INFO[dm2[i].indexes[3]].Expon;
         zeta_ijkl=zeta_ik+zeta_jl;
         Aijkl=dm2[i].Dijkl*pow(zeta_ijkl,-HALF*THREE)*Cnorm_4gauss;
         for(j=0;j<3;j++)
         {
          Rik[j]=Coord_atom[0][j]-Coord_atom[2][j];
          Rjl[j]=Coord_atom[1][j]-Coord_atom[3][j];
         }
         Aijkl=Aijkl
              *exp(-(BASIS_INFO[dm2[i].indexes[0]].Expon*BASIS_INFO[dm2[i].indexes[2]].Expon*scalarP(Rik)/zeta_ik)
              -(BASIS_INFO[dm2[i].indexes[1]].Expon*BASIS_INFO[dm2[i].indexes[3]].Expon*scalarP(Rjl)/zeta_jl));
         e_ijkl=zeta_ik*zeta_jl/zeta_ijkl;
         //Redifine and recycle Rik and Rjl since we alredy used them:
         for(j=0;j<3;j++)
         {
          Rik[j]=(BASIS_INFO[dm2[i].indexes[0]].Expon*Coord_atom[0][j]+BASIS_INFO[dm2[i].indexes[2]].Expon*Coord_atom[2][j])/zeta_ik;
          Rjl[j]=(BASIS_INFO[dm2[i].indexes[1]].Expon*Coord_atom[1][j]+BASIS_INFO[dm2[i].indexes[3]].Expon*Coord_atom[3][j])/zeta_jl;
          Rijkl[j]=(Rik[j]*zeta_ik+Rjl[j]*zeta_jl)/zeta_ijkl;
         }
         alpha_ijkl=(zeta_ik-zeta_jl)/(zeta_ijkl*TWO);
         nx_sum=nx_exp[0]+nx_exp[1]+nx_exp[2]+nx_exp[3];
         ny_sum=ny_exp[0]+ny_exp[1]+ny_exp[2]+ny_exp[3];
         nz_sum=nz_exp[0]+nz_exp[1]+nz_exp[2]+nz_exp[3];
         max_exp_ijkl=nx_sum;
         if(max_exp_ijkl<ny_sum)
         {max_exp_ijkl=ny_sum;}
         if(max_exp_ijkl<nz_sum)
         {max_exp_ijkl=nz_sum;}
         if(max_exp_ijkl%2==0)
         {
          order_ijkl=(max_exp_ijkl/2)+1;
         }
         else
         {
          order_ijkl=((max_exp_ijkl-1)/2)+1;
         }
         //Intracule coord. dependent part
         rscan=Init;
         for(j=0;j<nrad;j++)
         {
          for(k=0;k<nang;k++)
          {
           if(z[k]>=ZERO)
           {
            if(legendre)
            {
             Point_intra[0]=r_legendre[j]*x[k];
             Point_intra[1]=r_legendre[j]*y[k];
             Point_intra[2]=r_legendre[j]*z[k];
            }
            else
            {
             Point_intra[0]=rscan*x[k];
             Point_intra[1]=rscan*y[k];
             Point_intra[2]=rscan*z[k];
            }
            Xprime=pow(e_ijkl,HALF)*(Point_intra[0]+Rik[0]-Rjl[0]);
            Yprime=pow(e_ijkl,HALF)*(Point_intra[1]+Rik[1]-Rjl[1]);
            Zprime=pow(e_ijkl,HALF)*(Point_intra[2]+Rik[2]-Rjl[2]);
            if(z[k]>ZERO)
            {
             //If there is no angular grid Cnorm_ang=1 otherwise is 2 since z12<0 are not evaluated and the points are
             //symmetrically distributed on the unit sphere
             Tot_rad_th[j][k]=Tot_rad_th[j][k]+Cnorm_ang*Aijkl*exp(-(pow(Xprime,TWO)+pow(Yprime,TWO)+pow(Zprime,TWO)))
                     *Vint3D(Point_intra,order_ijkl,r_intrac[order_ijkl-1],w_intrac[order_ijkl-1],
                      Rijkl,alpha_ijkl,zeta_ijkl,Coord_atom,nx_exp,ny_exp,nz_exp);
            }
            else
            {
             Tot_rad_th[j][k]=Tot_rad_th[j][k]+Aijkl*exp(-(pow(Xprime,TWO)+pow(Yprime,TWO)+pow(Zprime,TWO)))
                     *Vint3D(Point_intra,order_ijkl,r_intrac[order_ijkl-1],w_intrac[order_ijkl-1],
                      Rijkl,alpha_ijkl,zeta_ijkl,Coord_atom,nx_exp,ny_exp,nz_exp);
            }
           }
           else
           {
            if(nang2==0)
            {
             symmetry_omitted_terms++;
            }
            else
            {
             if(legendre)
             {
              Point_intra[0]=r_legendre[j]*x[k];
              Point_intra[1]=r_legendre[j]*y[k];
              Point_intra[2]=r_legendre[j]*z[k];
             }
             else
             {
              Point_intra[0]=rscan*x[k];
              Point_intra[1]=rscan*y[k];
              Point_intra[2]=rscan*z[k];
             }
             Xprime=pow(e_ijkl,HALF)*(Point_intra[0]+Rik[0]-Rjl[0]);
             Yprime=pow(e_ijkl,HALF)*(Point_intra[1]+Rik[1]-Rjl[1]);
             Zprime=pow(e_ijkl,HALF)*(Point_intra[2]+Rik[2]-Rjl[2]);
             Tot_rad_th[j][k]=Tot_rad_th[j][k]+Aijkl*exp(-(pow(Xprime,TWO)+pow(Yprime,TWO)+pow(Zprime,TWO)))
                     *Vint3D(Point_intra,order_ijkl,r_intrac[order_ijkl-1],w_intrac[order_ijkl-1],
                      Rijkl,alpha_ijkl,zeta_ijkl,Coord_atom,nx_exp,ny_exp,nz_exp);
            }
           }
          }
          rscan=rscan+Step;
         }
        }
       }
       #pragma omp critical
       {
        for(i=0;i<nrad;i++)
        {
         for(j=0;j<nang;j++)
         {
          Tot_rad[i][j]=Tot_rad[i][j]+Tot_rad_th[i][j];
         }
        }
       }
       for(i=0;i<nrad;i++)
       {
        delete[] Tot_rad_th[i];Tot_rad_th[i]=NULL;
       }
       delete[] Tot_rad_th; Tot_rad_th=NULL;
       #pragma omp barrier
       #pragma omp master
       {
        for(i=0;i<nrad;i++)
        {
         Integral=ZERO;
         Integral2=ZERO;
         for(j=0;j<nang;j++)
         {
          Integral=Integral+w_theta_phi[j]*Tot_rad[i][j];
          if(abs(Tot_rad[i][j])>=pow(TEN,-TWO*SIX))
          {
           if(nang2==0)
           {
            if(z[j]>ZERO)
            {
             Integral2=Integral2-w_theta_phi[j]*TWO*abs(HALF*Tot_rad[i][j])*log(abs(HALF*Tot_rad[i][j]));
            }
            else
            {
             Integral2=Integral2-w_theta_phi[j]*abs(Tot_rad[i][j])*log(abs(Tot_rad[i][j]));
            }
           } 
           else
           {
            Integral2=Integral2-w_theta_phi[j]*abs(Tot_rad[i][j])*log(abs(Tot_rad[i][j]));
           }
          }
         }
         Tot_rad[i][0]=DMNfact*FOUR*PI*Integral;
         Shannon_rad[i]=DMNfact*FOUR*PI*Integral2;
#ifdef HAVE_MPI
globalIr=ZERO;
localIr=Tot_rad[i][0];
MPI_Reduce(&localIr,&globalIr,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
Tot_rad[i][0]=globalIr;
Shannon_rad[i]=ZERO;
#endif
        }
       }
      }
     }
     else
     {
if(ID==0)
{
      cout<<"Running"<<setw(4)<<1<<" thread."<<endl;
}
      for(i=0;i<nterms;i++)
      {
       Cnorm_4gauss=ONE;
       for(j=0;j<4;j++)
       {
        //Compute the normalization factor of the gaussians (ijkl)
        Cnorm_4gauss=Cnorm_4gauss*
                    normGauss(BASIS_INFO[dm2[i].indexes[j]].Expon,
                              BASIS_INFO[dm2[i].indexes[j]].nx,
                              BASIS_INFO[dm2[i].indexes[j]].ny,
                              BASIS_INFO[dm2[i].indexes[j]].nz);
        //Take the coord. of this 4 gaussians
        Coord_atom[j][0]=BASIS_INFO[dm2[i].indexes[j]].Cart_cord[0];
        Coord_atom[j][1]=BASIS_INFO[dm2[i].indexes[j]].Cart_cord[1];
        Coord_atom[j][2]=BASIS_INFO[dm2[i].indexes[j]].Cart_cord[2];
        // Take the quantum numbers of the primitives
        nx_exp[j]=BASIS_INFO[dm2[i].indexes[j]].nx;
        ny_exp[j]=BASIS_INFO[dm2[i].indexes[j]].ny;
        nz_exp[j]=BASIS_INFO[dm2[i].indexes[j]].nz;
       }
       //Jab is the upper bound fucntion of Eq. 36
       Jik=Jab(BASIS_INFO[dm2[i].indexes[0]].Expon,
               BASIS_INFO[dm2[i].indexes[2]].Expon,
               BASIS_INFO[dm2[i].indexes[0]].nx,
               BASIS_INFO[dm2[i].indexes[0]].ny,
               BASIS_INFO[dm2[i].indexes[0]].nz,
               BASIS_INFO[dm2[i].indexes[2]].nx,
               BASIS_INFO[dm2[i].indexes[2]].ny,
               BASIS_INFO[dm2[i].indexes[2]].nz,
               BASIS_INFO[dm2[i].indexes[0]].Cart_cord[0],
               BASIS_INFO[dm2[i].indexes[0]].Cart_cord[1],
               BASIS_INFO[dm2[i].indexes[0]].Cart_cord[2],
               BASIS_INFO[dm2[i].indexes[2]].Cart_cord[0],
               BASIS_INFO[dm2[i].indexes[2]].Cart_cord[1],
               BASIS_INFO[dm2[i].indexes[2]].Cart_cord[2]);
       Jjl=Jab(BASIS_INFO[dm2[i].indexes[1]].Expon,
               BASIS_INFO[dm2[i].indexes[3]].Expon,
               BASIS_INFO[dm2[i].indexes[1]].nx,
               BASIS_INFO[dm2[i].indexes[1]].ny,
               BASIS_INFO[dm2[i].indexes[1]].nz,
               BASIS_INFO[dm2[i].indexes[3]].nx,
               BASIS_INFO[dm2[i].indexes[3]].ny,
               BASIS_INFO[dm2[i].indexes[3]].nz,
               BASIS_INFO[dm2[i].indexes[1]].Cart_cord[0],
               BASIS_INFO[dm2[i].indexes[1]].Cart_cord[1],
               BASIS_INFO[dm2[i].indexes[1]].Cart_cord[2],
               BASIS_INFO[dm2[i].indexes[3]].Cart_cord[0],
               BASIS_INFO[dm2[i].indexes[3]].Cart_cord[1],
               BASIS_INFO[dm2[i].indexes[3]].Cart_cord[2]);
       if((abs(dm2[i].Dijkl)*Cnorm_4gauss*pow(Jik*Jjl,HALF))>=threshold2)
       {
        counter++;
        //Non intracule coord. dependent part
        zeta_ik=BASIS_INFO[dm2[i].indexes[0]].Expon+BASIS_INFO[dm2[i].indexes[2]].Expon;
        zeta_jl=BASIS_INFO[dm2[i].indexes[1]].Expon+BASIS_INFO[dm2[i].indexes[3]].Expon;
        zeta_ijkl=zeta_ik+zeta_jl;
        Aijkl=dm2[i].Dijkl*pow(zeta_ijkl,-HALF*THREE)*Cnorm_4gauss;
        for(j=0;j<3;j++)
        {
         Rik[j]=Coord_atom[0][j]-Coord_atom[2][j];
         Rjl[j]=Coord_atom[1][j]-Coord_atom[3][j];
        }
        Aijkl=Aijkl
             *exp(-(BASIS_INFO[dm2[i].indexes[0]].Expon*BASIS_INFO[dm2[i].indexes[2]].Expon*scalarP(Rik)/zeta_ik)
             -(BASIS_INFO[dm2[i].indexes[1]].Expon*BASIS_INFO[dm2[i].indexes[3]].Expon*scalarP(Rjl)/zeta_jl));
        e_ijkl=zeta_ik*zeta_jl/zeta_ijkl;
        //Redifine and recycle Rik and Rjl since we alredy used them:
        for(j=0;j<3;j++)
        {
         Rik[j]=(BASIS_INFO[dm2[i].indexes[0]].Expon*Coord_atom[0][j]+BASIS_INFO[dm2[i].indexes[2]].Expon*Coord_atom[2][j])/zeta_ik;
         Rjl[j]=(BASIS_INFO[dm2[i].indexes[1]].Expon*Coord_atom[1][j]+BASIS_INFO[dm2[i].indexes[3]].Expon*Coord_atom[3][j])/zeta_jl;
         Rijkl[j]=(Rik[j]*zeta_ik+Rjl[j]*zeta_jl)/zeta_ijkl;
        }
        alpha_ijkl=(zeta_ik-zeta_jl)/(zeta_ijkl*TWO);
        nx_sum=nx_exp[0]+nx_exp[1]+nx_exp[2]+nx_exp[3];
        ny_sum=ny_exp[0]+ny_exp[1]+ny_exp[2]+ny_exp[3];
        nz_sum=nz_exp[0]+nz_exp[1]+nz_exp[2]+nz_exp[3];
        max_exp_ijkl=nx_sum;
        if(max_exp_ijkl<ny_sum)
        {max_exp_ijkl=ny_sum;}
        if(max_exp_ijkl<nz_sum)
        {max_exp_ijkl=nz_sum;}
        if(max_exp_ijkl%2==0)
        {
         order_ijkl=(max_exp_ijkl/2)+1;
        }
        else
        {
         order_ijkl=((max_exp_ijkl-1)/2)+1;
        }
        //Intracule coord. dependent part
        rscan=Init;
        for(j=0;j<nrad;j++)
        {
         for(k=0;k<nang;k++)
         {
          if(z[k]>=ZERO)
          {
           if(legendre)
           {
            Point_intra[0]=r_legendre[j]*x[k];
            Point_intra[1]=r_legendre[j]*y[k];
            Point_intra[2]=r_legendre[j]*z[k];
           }
           else
           {
            Point_intra[0]=rscan*x[k];
            Point_intra[1]=rscan*y[k];
            Point_intra[2]=rscan*z[k];
           }
           Xprime=pow(e_ijkl,HALF)*(Point_intra[0]+Rik[0]-Rjl[0]);
           Yprime=pow(e_ijkl,HALF)*(Point_intra[1]+Rik[1]-Rjl[1]);
           Zprime=pow(e_ijkl,HALF)*(Point_intra[2]+Rik[2]-Rjl[2]);
           if(z[k]>ZERO)
           {
            //If there is no angular grid Cnorm_ang=1 otherwise is 2 since z12<0 are not evaluated and the points are
            //symmetrically distributed in the unit sphere
            Tot_rad[j][k]=Tot_rad[j][k]+Cnorm_ang*Aijkl*exp(-(pow(Xprime,TWO)+pow(Yprime,TWO)+pow(Zprime,TWO)))*
            Vint3D(Point_intra,order_ijkl,r_intrac[order_ijkl-1],w_intrac[order_ijkl-1],Rijkl,alpha_ijkl,zeta_ijkl,
            Coord_atom,nx_exp,ny_exp,nz_exp);
           }
           else
           {
            Tot_rad[j][k]=Tot_rad[j][k]+Aijkl*exp(-(pow(Xprime,TWO)+pow(Yprime,TWO)+pow(Zprime,TWO)))*
            Vint3D(Point_intra,order_ijkl,r_intrac[order_ijkl-1],w_intrac[order_ijkl-1],Rijkl,alpha_ijkl,zeta_ijkl,
            Coord_atom,nx_exp,ny_exp,nz_exp);
           }
          }
          else
          {
           if(nang2==0)
           {
            symmetry_omitted_terms++;
           }
           else
           {
            if(legendre)
            {
             Point_intra[0]=r_legendre[j]*x[k];
             Point_intra[1]=r_legendre[j]*y[k];
             Point_intra[2]=r_legendre[j]*z[k];
            }
            else
            {
             Point_intra[0]=rscan*x[k];
             Point_intra[1]=rscan*y[k];
             Point_intra[2]=rscan*z[k];
            }
            Xprime=pow(e_ijkl,HALF)*(Point_intra[0]+Rik[0]-Rjl[0]);
            Yprime=pow(e_ijkl,HALF)*(Point_intra[1]+Rik[1]-Rjl[1]);
            Zprime=pow(e_ijkl,HALF)*(Point_intra[2]+Rik[2]-Rjl[2]);
            Tot_rad[j][k]=Tot_rad[j][k]+Aijkl*exp(-(pow(Xprime,TWO)+pow(Yprime,TWO)+pow(Zprime,TWO)))*
            Vint3D(Point_intra,order_ijkl,r_intrac[order_ijkl-1],w_intrac[order_ijkl-1],Rijkl,alpha_ijkl,zeta_ijkl,
            Coord_atom,nx_exp,ny_exp,nz_exp);
           }
          }
         }
         rscan=rscan+Step;
        }
       }
      }
     }
     ////////////////
     //Check times //
     ////////////////
if(ID==0)
{
     if(Input_commands.time_intra){cout<<"Time 3 (D H M S): "<<endl;system("date +%j' '%H' '%M' '%S ");cout<<endl;}
     ////////////////
     //End time    //
     ////////////////
     cout<<"#*****************************************************************************#";
     cout<<endl;
     cout<<"Computed terms from DM2:"<<setw(17)<<counter<<endl;
     cout<<"[Terms that passed |^2Dijkl Iijkl(R)| >= 2*tau/[Np(Np+1)]]"<<endl;
     if(nang!=0)
     {
      cout<<"Symmetry omit. terms   :"<<setw(17)<<symmetry_omitted_terms<<endl;
      cout<<"[Omitted terms are intracule points with z12 < 0. They are omitted because I(R)=I(-R)]"<<endl;
     }
     cout<<"#*****************************************************************************#";
     cout<<endl;
     cout<<"#*****************************************************************************#";
     cout<<endl;
     cout<<"Let u=|r2-r1| and I(u) = int I_3D(u,Omega_12) d Omega_12 "<<endl;
     cout<<"#*****************************************************************************#";
     cout<<endl;
     cout<<"#*****************************************************************************#";
     cout<<endl;
#ifdef HAVE_MPI
     cout<<"#    u              I(u)                I(u)u^2           I(u)*erf[l*u]*u^2   I(u)*exp[-l*u]*u^2"<<endl;
#else
     cout<<"#    u              I(u)                I(u)u^2           I(u)*erf[l*u]*u^2   I(u)*exp[-l*u]*u^2  S_M[I_3D(u,Omega_12)]"<<endl;
#endif
}
     for(i=0;i<10;i++){res[i]=ZERO;}
     rscan=Init;
     for(i=0;i<nrad;i++)
     {
      if(!parallel)
      {
       Integral=ZERO;
       Integral2=ZERO;
       for(j=0;j<nang;j++)
       {
        Integral=Integral+w_theta_phi[j]*Tot_rad[i][j];
        if(abs(Tot_rad[i][j])>=pow(TEN,-TWO*SIX))
        {
         if(nang2==0)
         {
          if(z[j]>ZERO)
          {
           Integral2=Integral2-w_theta_phi[j]*TWO*abs(HALF*Tot_rad[i][j])*log(abs(HALF*Tot_rad[i][j]));
          }
          else
          {
           Integral2=Integral2-w_theta_phi[j]*abs(Tot_rad[i][j])*log(abs(Tot_rad[i][j]));
          }
         } 
         else
         {
          Integral2=Integral2-w_theta_phi[j]*abs(Tot_rad[i][j])*log(abs(Tot_rad[i][j]));
         }
        }
       }
       Tot_rad[i][0]=DMNfact*FOUR*PI*Integral;
       Shannon_rad[i]=DMNfact*FOUR*PI*Integral2;
#ifdef HAVE_MPI
globalIr=ZERO;
localIr=Tot_rad[i][0];
MPI_Reduce(&localIr,&globalIr,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
Tot_rad[i][0]=globalIr;
Shannon_rad[i]=ZERO;
#endif
      }
      if(abs(Tot_rad[i][0])<pow(TEN,-THREE*FIVE)){Tot_rad[i][0]=ZERO;}
      if(abs(Tot_rad[i][0])>ZERO)
      {
       if(legendre)
       {
if(ID==0)
{
        cout<<setprecision(8)<<fixed<<scientific<<r_legendre[i]<<setw(20)<<Tot_rad[i][0]<<setw(20);
        cout<<Tot_rad[i][0]*pow(r_legendre[i],TWO)<<setw(20)<<Tot_rad[i][0]*pow(r_legendre[i],TWO)*erf(r_legendre[i]*lambda_rs)<<setw(20);
#ifdef HAVE_MPI
        cout<<Tot_rad[i][0]*pow(r_legendre[i],TWO)*exp(-r_legendre[i]*lambda_scr);
#else
        cout<<Tot_rad[i][0]*pow(r_legendre[i],TWO)*exp(-r_legendre[i]*lambda_scr)<<setw(20)<<Shannon_rad[i];
#endif
        cout<<endl;
}
       }
       else
       {
if(ID==0)
{
        cout<<setprecision(8)<<fixed<<scientific<<rscan<<setw(20)<<Tot_rad[i][0]<<setw(20);
        cout<<Tot_rad[i][0]*pow(rscan,TWO)<<setw(20)<<Tot_rad[i][0]*pow(rscan,TWO)*erf(rscan*lambda_rs)<<setw(20);
#ifdef HAVE_MPI
        cout<<Tot_rad[i][0]*pow(rscan,TWO)*exp(-rscan*lambda_scr);
#else
        cout<<Tot_rad[i][0]*pow(rscan,TWO)*exp(-rscan*lambda_scr)<<setw(20)<<Shannon_rad[i];
#endif
        cout<<endl;
}
       }
      }
      if(legendre)
      {
       res[0]=res[0]+Tot_rad[i][0]*Jacobian_legendre[i][0]*w_legendre[i];
       res[1]=res[1]+Tot_rad[i][0]*Jacobian_legendre[i][1]*w_legendre[i];
       res[2]=res[2]+Tot_rad[i][0]*Jacobian_legendre[i][2]*w_legendre[i];
       res[3]=res[3]+Tot_rad[i][0]*Jacobian_legendre[i][3]*w_legendre[i];
       res[4]=res[4]+Tot_rad[i][0]*Jacobian_legendre[i][4]*w_legendre[i];
       res[5]=res[5]+Shannon_rad[i]*Jacobian_legendre[i][2]*w_legendre[i];
       res[6]=res[6]+Tot_rad[i][0]*erf(r_legendre[i]*lambda_rs)*Jacobian_legendre[i][1]*w_legendre[i];
       res[7]=res[7]+Tot_rad[i][0]*erf(r_legendre[i]*lambda_rs)*Jacobian_legendre[i][2]*w_legendre[i];
       res[8]=res[8]+Tot_rad[i][0]*exp(-r_legendre[i]*lambda_scr)*Jacobian_legendre[i][1]*w_legendre[i];
       res[9]=res[9]+Tot_rad[i][0]*exp(-r_legendre[i]*lambda_scr)*Jacobian_legendre[i][2]*w_legendre[i];
      }
      else
      {
       if(i>0)
       {
        res[0]=res[0]+HALF*(Tot_rad[i-1][0]+Tot_rad[i][0])*Step;
        res[1]=res[1]+HALF*(Tot_rad[i-1][0]*(rscan-Step)+Tot_rad[i][0]*rscan)*Step;
        res[2]=res[2]+HALF*(Tot_rad[i-1][0]*pow(rscan-Step,TWO)+Tot_rad[i][0]*pow(rscan,TWO))*Step;
        res[3]=res[3]+HALF*(Tot_rad[i-1][0]*pow(rscan-Step,THREE)+Tot_rad[i][0]*pow(rscan,THREE))*Step;
        res[4]=res[4]+HALF*(Tot_rad[i-1][0]*pow(rscan-Step,FOUR)+Tot_rad[i][0]*pow(rscan,FOUR))*Step;
        res[5]=res[5]+HALF*(Shannon_rad[i-1]*pow(rscan-Step,TWO)+Shannon_rad[i]*pow(rscan,TWO))*Step;
        res[6]=res[6]+HALF*(Tot_rad[i-1][0]*erf((rscan-Step)*lambda_rs)*(rscan-Step)+Tot_rad[i][0]*erf(rscan*lambda_rs)*(rscan))*Step;
        res[7]=res[7]+HALF*(Tot_rad[i-1][0]*erf((rscan-Step)*lambda_rs)*pow(rscan-Step,TWO)
              +Tot_rad[i][0]*erf(rscan*lambda_rs)*pow(rscan,TWO))*Step;
        res[8]=res[8]+HALF*(Tot_rad[i-1][0]*exp(-(rscan-Step)*lambda_scr)*(rscan-Step)+Tot_rad[i][0]*exp(-rscan*lambda_scr)*(rscan))*Step;
        res[7]=res[7]+HALF*(Tot_rad[i-1][0]*exp(-(rscan-Step)*lambda_scr)*pow(rscan-Step,TWO)
              +Tot_rad[i][0]*exp(-rscan*lambda_scr)*pow(rscan,TWO))*Step;
       }
      }
      rscan=rscan+Step;
     }
if(ID==0)
{    
     cout<<endl;
     cout<<"#################################################"<<endl;
     cout<<endl;
     if(legendre)
     {
      cout<<"Integration by Legendre Quadrature rule:"<<endl;
     }
     else
     {
      cout<<"Integration by Trapezoidal rule:"<<endl;
     }
     cout<<"int I(u)           du  :"<<setw(18)<<res[0]<<endl;
     cout<<"int I(u) u         du  :"<<setw(18)<<res[1]<<endl;
     cout<<"int I(u) u^2       du  :"<<setw(18)<<res[2]<<endl;
     cout<<"int I(u) u^3       du  :"<<setw(18)<<res[3]<<endl;
     cout<<"int I(u) u^4       du  :"<<setw(18)<<res[4]<<endl;
     res[2]=round(res[2]);
     cout<<"Tr[^2Dijkl]            :"<<setw(18)<<res[2]<<endl;
     cout<<"Vee                    :"<<setw(18)<<res[1]<<endl;
     cout<<"[Trace is N(N-1)/2]"<<endl;
     cout<<"<  u  > per pair of e- :"<<setw(18)<<res[3]/res[2]<<endl;
     cout<<"< u^2 > per pair of e- :"<<setw(18)<<res[4]/res[2]<<endl;
#ifndef HAVE_MPI
     cout<<"S_Norm[I(u,Omega_12)]  :"<<setw(18)<<(res[2]*log(res[2])+res[5])/res[2]<<endl;
#endif
     cout<<"Range sep lambda (l)   :"<<setw(18)<<lambda_rs<<endl;
     cout<<"Range sep Vee          :"<<setw(18)<<res[6]<<endl;
     cout<<"Range sep Tr[^2Dijkl]  :"<<setw(18)<<res[7]<<endl;
     cout<<"Screening lambda (l)   :"<<setw(18)<<lambda_scr<<endl;
     cout<<"Screened Vee           :"<<setw(18)<<res[8]<<endl;
     cout<<"Screened Tr[^2Dijkl]   :"<<setw(18)<<res[9]<<endl;
     cout<<"<u> and <u^2> are per pair of electrons to ensure that a PROBABILITY is used!"<<endl;
#ifndef HAVE_MPI
     cout<<"S_Norm[I_3D(u,Omega_12)] = (M*Log[M]+S_M[I_3D(u,Omega_12)])/M where M = Tr[^2Dijkl] and S_M is normalized to M"<<endl;
#endif
     cout<<"Range sep = Integrals of the intracule including erf[lambda_rs*u] "<<endl;
     cout<<"Screening = Integrals of the intracule including exp[-lambda_scr*u] "<<endl;
     cout<<endl;
}
     delete[] dm2;
if(ID==0)
{
     if(Input_commands.rweight)
     {
      if(legendre)
      {
       ofstream r_weights;
       r_weights.open("r_weights.txt");
       r_weights<<"#Radial_coordinate (u) \t   I(u)  \t\t Jacobian (for u^2) \t Weight"<<endl;
       for(i=0;i<nrad;i++)
       {
        if(abs(Tot_rad[i][0])>ZERO)
        {
         r_weights<<setprecision(10)<<fixed<<scientific<<r_legendre[i]<<"\t"<<Tot_rad[i][0]<<"\t"<<Jacobian_legendre[i][2]<<"\t"<<w_legendre[i]<<endl;
        }
       }
       r_weights.close();
       cout<<"r_weights.txt file was generated for the present job"<<endl;
       cout<<endl;
      }
      else
      {
       cout<<"No weights were generated for the radial scan"<<endl;
      }
     }
}
    }
    else
    {
if(ID==0)
{
     cout<<endl;
}
     open_dm2.close();
if(ID==0)
{
     cout<<"Error! Could not open file: "<<name_dm2<<endl;
     cout<<endl;
}
    }
    delete[] BASIS_INFO;
    for(i=0;i<Nroot_2Lmax_plus_1;i++)
    {
     delete[] r_intrac[i];r_intrac[i]=NULL;
     delete[] w_intrac[i];w_intrac[i]=NULL;
    }
    delete[] r_intrac;r_intrac=NULL;
    delete[] w_intrac;w_intrac=NULL;
   }
   else
   {
if(ID==0)
{
    cout<<endl;
}
    open_basis.close();
if(ID==0)
{
    cout<<"Error! Could not open file: "<<name_basis<<endl;
    cout<<endl;
}
   }
   delete[] x;
   delete[] y;
   delete[] z;
   delete[] w_theta_phi;
   for(i=0;i<nrad;i++)
   {
    delete[] Tot_rad[i];Tot_rad[i]=NULL;
   }
   delete[] Tot_rad; Tot_rad=NULL;
   delete[] Shannon_rad; Shannon_rad=NULL;
   if(legendre)
   {
    delete[] r_legendre; r_legendre=NULL;
    delete[] w_legendre; w_legendre=NULL;
    for(i=0;i<nrad;i++)
    {
     delete[] Jacobian_legendre[i]; //Jacobian_legendre[i]=NULL;
    }
    delete[] Jacobian_legendre; //Jacobian_legendre=NULL;
   }
  }
  /////////////////////
  //Compute Extracule//
  /////////////////////
  if(Input_commands.extracule)
  {
if(ID==0)
{
   cout<<"#******************************************#"<<endl;
   cout<<"#          Compute the Extracule           #"<<endl;
   cout<<"#******************************************#"<<endl;
}
   Init=Input_commands.radial_Init;
   Step=Input_commands.radial_Step;
   Last=Input_commands.radial_Last;
   string name_dm2=Input_commands.name_dm2;
   string name_basis=Input_commands.name_basis;
   threshold=Input_commands.threshold_in;
   parallel=Input_commands.parallel;
   ////////////////
   //Check times //
   ////////////////
if(ID==0)
{
   if(Input_commands.time_extra){cout<<endl;cout<<"Time 1 (D H M S): "<<endl;system("date +%j' '%H' '%M' '%S ");cout<<endl;}
}
   ////////////////
   //End time    //
   ////////////////
if(ID==0)
{
   if(parallel)
   {
    cout<<"OpenMP parallelized evaluations of the extracule will be done."<<endl;
   }
   else
   {
    cout<<"Serial evaluations of the extracule will be done."<<endl;
   }
   cout<<"Threshold for DM2 terms:"<<setprecision(10)<<fixed<<scientific<<setw(18)<<Input_commands.threshold_in<<endl;
}
   //We define the DMNfactor to correct from McWeeney Normalization [N(N-1)/2]2! to Lowdin N(N-1)/2 for the DMN matrix.
   DMNfact=HALF;
   //Create angular grid
   nang=Input_commands.order_ang;
   nang2=Input_commands.order_ang2;
   if(nang>0 && nang2==0)
   {
    grid_avail(nang);
   }
   else if(nang>0 && nang2>0)
   {
    nang=nang*nang2;
   }
   else
   {
if(ID==0)
{
    cout<<"None angular grid chosen for spherical average"<<setw(17)<<endl;
}
    nang=1;
   }
   x=new double[nang];
   y=new double[nang];
   z=new double[nang];
   w_theta_phi=new double[nang];
   if(nang>1)
   {
    for(i=0;i<nang;i++)
    {
     x[i]=ZERO;
     y[i]=ZERO;
     z[i]=ZERO;
     w_theta_phi[i]=ZERO;
    }
    if(nang2==0)
    {
     ld_by_order(nang,x,y,z,w_theta_phi);
    }
    else
    {
     user_defined_angular(Input_commands.order_ang,Input_commands.order_ang2,x,y,z,w_theta_phi,Input_commands.name_basis);
    }
   }
   else
   {
    // We allow to scan only in z12 axis!
    x[0]=ZERO;
    y[0]=ZERO;
    z[0]=ONE;
    w_theta_phi[0]=ONE;
   }
   //Create radial grid
   legendre=Input_commands.legendre;
   if(!legendre)
   {
    nrad=0;
    for(rscan=Init;rscan<Last+Step;rscan=rscan+Step){nrad++;}
    Tot_rad=new double*[nrad];
    Shannon_rad=new double[nrad];
    for(i=0;i<nrad;i++)
    {
     Shannon_rad[i]=ZERO;
     Tot_rad[i]=new double[nang];
     // Initialize I(R) for scheme 3 in the paper
     for(j=0;j<nang;j++)
     {
      Tot_rad[i][j]=ZERO;
     }
    }
   }
   else
   {
    nrad=Input_commands.radial_grid;
    Tot_rad=new double*[nrad];
    Shannon_rad=new double[nrad];
    for(i=0;i<nrad;i++)
    {
     Shannon_rad[i]=ZERO;
     Tot_rad[i]=new double[nang];
     // Initialize I(R) for scheme 3 in the paper
     for(j=0;j<nang;j++)
     {
      Tot_rad[i][j]=ZERO;
     }
    }
    r_legendre=new double[nrad];
    w_legendre=new double[nrad];
    Jacobian_legendre=new double*[nrad];
    for(i=0;i<nrad;i++)
    {
     Jacobian_legendre[i]=new double[5];
    }
    Init=ZERO;
    if(Input_commands.radial_Last!=1e99)
    {
     Last=Input_commands.radial_Last;
    }
    else
    {
     Last=ONE;
    }
#ifdef HAVE_MPI
MPI_Barrier(MPI_COMM_WORLD);
#endif
if(ID==0)
{
    legendre_quadrature(name_basis.substr(0,name_basis.length()-6),nrad,Init,Last);
    //Read quadrature info
    ifstream read_rad_quad;
    // Read weights
    read_rad_quad.open((name_basis.substr(0,name_basis.length()-6)+"_w.txt").c_str());
    for(i=0;i<nrad;i++)
    {
     read_rad_quad>>w_legendre[i];
    }
    read_rad_quad.close();
    // Read roots
    read_rad_quad.open((name_basis.substr(0,name_basis.length()-6)+"_x.txt").c_str());
    for(i=0;i<nrad;i++)
    {
     read_rad_quad>>r_legendre[i];
    }
    read_rad_quad.close();
    // Remove quadrature files
    system(("rm "+name_basis.substr(0,name_basis.length()-6)+"_r.txt").c_str());
    system(("rm "+name_basis.substr(0,name_basis.length()-6)+"_w.txt").c_str());
    system(("rm "+name_basis.substr(0,name_basis.length()-6)+"_x.txt").c_str());
}
#ifdef HAVE_MPI
MPI_Bcast(w_legendre,nrad,MPI_DOUBLE,0,MPI_COMM_WORLD);
MPI_Bcast(r_legendre,nrad,MPI_DOUBLE,0,MPI_COMM_WORLD);
#endif
    if(Input_commands.radial_Last==1e99)
    {
if(ID==0)
{
     cout<<"Superior lim. Legendre :         Infinity"<<endl;
}
     for(i=0;i<nrad;i++)
     {
      for(j=0;j<5;j++)
      {
       Jacobian_legendre[i][j]=pow(r_legendre[i]/(ONE-r_legendre[i]),(double)j)/pow(ONE-r_legendre[i],TWO);
      }
      //We convert the variable r_legendre[i] from 0->1 to 0->Infinity (to be used later, where real [not from 0 to 1] evals of intracule must be performed!)
      r_legendre[i]=r_legendre[i]/(ONE-r_legendre[i]);
     }
    }
    else
    {
if(ID==0)
{
     cout<<"Superior lim. Legendre :"<<setprecision(10)<<fixed<<scientific<<setw(18)<<Last<<endl;
}
     for(i=0;i<nrad;i++)
     {
      for(j=0;j<5;j++)
      {
       Jacobian_legendre[i][j]=pow(r_legendre[i],(double)j);
      }
     }
    }
   }
   // Read basis info
   ifstream open_basis;
   open_basis.open((name_basis).c_str());
   if(open_basis.good())
   {
    nprims=0;
    while(getline(open_basis,aux))
    {
     nprims++;
    }
    open_basis.close();
    nprims=nprims-1;
    basis_info *BASIS_INFO=new basis_info[nprims];
    open_basis.open((name_basis).c_str());
    getline(open_basis,aux);
    for(i=0;i<nprims;i++)
    {
     open_basis>>BASIS_INFO[i].Z>>BASIS_INFO[i].Cart_cord[0]>>BASIS_INFO[i].Cart_cord[1]>>BASIS_INFO[i].Cart_cord[2];
     open_basis>>BASIS_INFO[i].Expon>>BASIS_INFO[i].nx>>BASIS_INFO[i].ny>>BASIS_INFO[i].nz;
    }
    open_basis.close();
    //We use the number of primitives to compute threshold 2 and we also initialize Rmax matrix here
if(ID==0)
{
    cout<<"tau                    :"<<setprecision(10)<<fixed<<scientific<<setw(18)<<Input_commands.tau<<endl;
}
    threshold2=TWO*Input_commands.tau/((double)nprims*((double)nprims+ONE));
if(ID==0)
{
    cout<<"2*tau/[Np(Np+1)]       :"<<setprecision(10)<<fixed<<scientific<<setw(18)<<threshold2<<endl;
    cout<<"[Np = Number of primitives]"<<endl;
}
    initialize(Rmax);
    //Compute Nroot_2Lmax_plus_1
    Lmax=0;
    for(i=0;i<nprims;i++)
    {
     if(BASIS_INFO[i].nx>Lmax)
     {
      Lmax=BASIS_INFO[i].nx;
     }
    }
if(ID==0)
{
    cout<<"Lmax                   :"<<setw(18)<<Lmax<<endl;
}
    //Create quadrature rule order
    Nroot_2Lmax_plus_1=2*Lmax+1;
if(ID==0)
{
    cout<<"Quadrature rule order  :"<<setw(18)<<Nroot_2Lmax_plus_1<<endl;
    cout<<"[Quadrature for primitive integrals. Order = 2 Lmax + 1 ]"<<endl;
}
    alpha=ZERO;
    a=ZERO;
    b=ONE;
    r_extrac=new double*[Nroot_2Lmax_plus_1];
    w_extrac=new double*[Nroot_2Lmax_plus_1];
    for(i=0;i<Nroot_2Lmax_plus_1;i++)
    {
     r_extrac[i]=new double[i+1];
     w_extrac[i]=new double[i+1];
     for(j=0;j<i+1;j++)
     {
      r_extrac[i][j]=ZERO;
      w_extrac[i][j]=ZERO;
     }
    }
    //Read quadrature info
    ifstream read_quad_r,read_quad_w;
    for(i=0;i<Nroot_2Lmax_plus_1;i++)
    {
     j=i+1;
#ifdef HAVE_MPI
double *wtmp,*rtmp;
wtmp=new double[j];
rtmp=new double[j];
MPI_Barrier(MPI_COMM_WORLD);
#endif
if(ID==0)
{
     gauss_hermite_rule(name_basis.substr(0,name_basis.length()-6),alpha,a,b,j);
     // Read weights
     read_quad_w.open((name_basis.substr(0,name_basis.length()-6)+"_w.txt").c_str());
     // Read roots
     read_quad_r.open((name_basis.substr(0,name_basis.length()-6)+"_x.txt").c_str());
     for(j=0;j<i+1;j++)
     {
      read_quad_w>>w_extrac[i][j];
      read_quad_r>>r_extrac[i][j];
#ifdef HAVE_MPI
wtmp[j]=w_extrac[i][j];
rtmp[j]=r_extrac[i][j];
#endif
     }
     read_quad_r.close();
     read_quad_w.close();
     // Remove quadrature files
     system(("rm "+name_basis.substr(0,name_basis.length()-6)+"_r.txt").c_str());
     system(("rm "+name_basis.substr(0,name_basis.length()-6)+"_w.txt").c_str());
     system(("rm "+name_basis.substr(0,name_basis.length()-6)+"_x.txt").c_str());
}
#ifdef HAVE_MPI
MPI_Bcast(wtmp,j,MPI_DOUBLE,0,MPI_COMM_WORLD);
MPI_Bcast(rtmp,j,MPI_DOUBLE,0,MPI_COMM_WORLD);
for(j=0;j<i+1;j++)
{
w_extrac[i][j]=wtmp[j];
r_extrac[i][j]=rtmp[j];
}
delete[] wtmp;wtmp=NULL;
delete[] rtmp;rtmp=NULL;
#endif
    }
    //Read DM2  file
    ifstream open_dm2;
    open_dm2.open((name_dm2).c_str());
    if(open_dm2.good())
    {
     open_dm2.close();
     nterms=0;
     //Check number of terms
     terms_dm2(name_dm2);
#ifdef HAVE_MPI
if(ID==0)
{
cout<<"Reducing the number of DM2 terms to "<<setw(8)<<nterms/nproc+1<<" per proc."<<endl;
}
nterms=nterms/nproc+1;
#endif
     dm2=new DM2[nterms];
     //Store the DM2
if(ID==0)
{
     cout<<"Storage of the DM2"<<endl;
}
     fill_in_dm2(name_dm2);
if(ID==0)
{
     cout<<"Storage done!"<<endl;
     //Number of terms to compute
     cout<<"Total extracule points : "<<setw(17)<<nang*nrad<<endl;
     cout<<"[Total extracule points = Order angular grid x Num radial points |U|.]"<<endl;
     cout<<"Extracule points x DM2 : "<<setw(17)<<(double)nterms*(double)nang*(double)nrad<<endl;
     cout<<"[Extracule points x DM2 terms = Nterms DM2 x Order angular grid x Num radial points |U|.]"<<endl;
     cout<<"Note: If no angular grid was chosen, Order angular grid is equal to ONE."<<endl;
     //Evaluate intracule
     ////////////////
     //Check times //
     ////////////////
     if(Input_commands.time_extra){cout<<endl;cout<<"Time 2 (D H M S): "<<endl;system("date +%j' '%H' '%M' '%S ");cout<<endl;}
     ////////////////
     //End time    //
     ////////////////
     cout<<"#*****************************************************************************#";
     cout<<endl;
}
     if(parallel)
     {
      if(Input_commands.nthreads>omp_get_max_threads())
      {Input_commands.nthreads=omp_get_max_threads();}
if(ID==0)
{
      cout<<"Running"<<setw(4)<<Input_commands.nthreads<<" threads."<<endl;
}
      #pragma omp parallel num_threads(Input_commands.nthreads) \
       private(i,j,k,Point_extra,Xprime,Yprime,Zprime,Jik,Jjl,Cnorm_4gauss,Coord_atom, \
       nx_sum,ny_sum,nz_sum,order_ijkl,max_exp_ijkl,nx_exp,ny_exp,nz_exp,zeta_ik,zeta_jl, \
       zeta_ijkl,Aijkl,e_ijkl,Rik,Rjl,Rijkl,alpha_ijkl,rscan,Integral,Integral2) \
       shared(BASIS_INFO,dm2,r_extrac,w_extrac,r_legendre,x,y,z,w_theta_phi,Tot_rad,Shannon_rad, \
       nterms,nrad,nang) \
       reduction(+:symmetry_omitted_terms,counter)
      {
       //Variables that belong to each thread and are only needed here
       int nth=omp_get_num_threads();
       int ith=omp_get_thread_num();
       double **Tot_rad_th;
       Tot_rad_th=new double*[nrad];
       for(i=0;i<nrad;i++)
       {
        Tot_rad_th[i]=new double[nang];
        for(j=0;j<nang;j++)
        {
         Tot_rad_th[i][j]=ZERO;
        }
       }
       for(i=ith;i<nterms;i=i+nth)
       {
        Cnorm_4gauss=ONE;
        for(j=0;j<4;j++)
        {
         //Compute the normalization factor of the gaussians (ijkl)
         Cnorm_4gauss=Cnorm_4gauss*
                     normGauss(BASIS_INFO[dm2[i].indexes[j]].Expon,
                               BASIS_INFO[dm2[i].indexes[j]].nx,
                               BASIS_INFO[dm2[i].indexes[j]].ny,
                               BASIS_INFO[dm2[i].indexes[j]].nz);
         //Take the coord. of this 4 gaussians
         Coord_atom[j][0]=BASIS_INFO[dm2[i].indexes[j]].Cart_cord[0];
         Coord_atom[j][1]=BASIS_INFO[dm2[i].indexes[j]].Cart_cord[1];
         Coord_atom[j][2]=BASIS_INFO[dm2[i].indexes[j]].Cart_cord[2];
         // Take the quantum numbers of the primitives
         nx_exp[j]=BASIS_INFO[dm2[i].indexes[j]].nx;
         ny_exp[j]=BASIS_INFO[dm2[i].indexes[j]].ny;
         nz_exp[j]=BASIS_INFO[dm2[i].indexes[j]].nz;
        }
        //Jab is the upper bound fucntion of Eq. 39
        Jik=Jab(BASIS_INFO[dm2[i].indexes[0]].Expon,
                BASIS_INFO[dm2[i].indexes[2]].Expon,
                BASIS_INFO[dm2[i].indexes[0]].nx,
                BASIS_INFO[dm2[i].indexes[0]].ny,
                BASIS_INFO[dm2[i].indexes[0]].nz,
                BASIS_INFO[dm2[i].indexes[2]].nx,
                BASIS_INFO[dm2[i].indexes[2]].ny,
                BASIS_INFO[dm2[i].indexes[2]].nz,
                BASIS_INFO[dm2[i].indexes[0]].Cart_cord[0],
                BASIS_INFO[dm2[i].indexes[0]].Cart_cord[1],
                BASIS_INFO[dm2[i].indexes[0]].Cart_cord[2],
                BASIS_INFO[dm2[i].indexes[2]].Cart_cord[0],
                BASIS_INFO[dm2[i].indexes[2]].Cart_cord[1],
                BASIS_INFO[dm2[i].indexes[2]].Cart_cord[2]);
        Jjl=Jab(BASIS_INFO[dm2[i].indexes[1]].Expon,
                BASIS_INFO[dm2[i].indexes[3]].Expon,
                BASIS_INFO[dm2[i].indexes[1]].nx,
                BASIS_INFO[dm2[i].indexes[1]].ny,
                BASIS_INFO[dm2[i].indexes[1]].nz,
                BASIS_INFO[dm2[i].indexes[3]].nx,
                BASIS_INFO[dm2[i].indexes[3]].ny,
                BASIS_INFO[dm2[i].indexes[3]].nz,
                BASIS_INFO[dm2[i].indexes[1]].Cart_cord[0],
                BASIS_INFO[dm2[i].indexes[1]].Cart_cord[1],
                BASIS_INFO[dm2[i].indexes[1]].Cart_cord[2],
                BASIS_INFO[dm2[i].indexes[3]].Cart_cord[0],
                BASIS_INFO[dm2[i].indexes[3]].Cart_cord[1],
                BASIS_INFO[dm2[i].indexes[3]].Cart_cord[2]);
        if(EIGHT*(abs(dm2[i].Dijkl)*Cnorm_4gauss*pow(Jik*Jjl,HALF))>=threshold2)
        {
         counter++;
         //Non extracule coord. dependent part
         zeta_ik=BASIS_INFO[dm2[i].indexes[0]].Expon+BASIS_INFO[dm2[i].indexes[2]].Expon;
         zeta_jl=BASIS_INFO[dm2[i].indexes[1]].Expon+BASIS_INFO[dm2[i].indexes[3]].Expon;
         zeta_ijkl=zeta_ik+zeta_jl;
         Aijkl=EIGHT*dm2[i].Dijkl*pow(zeta_ijkl,-HALF*THREE)*Cnorm_4gauss;
         for(j=0;j<3;j++)
         {
          Rik[j]=Coord_atom[0][j]-Coord_atom[2][j];
          Rjl[j]=Coord_atom[1][j]-Coord_atom[3][j];
         }
         Aijkl=Aijkl
              *exp(-(BASIS_INFO[dm2[i].indexes[0]].Expon*BASIS_INFO[dm2[i].indexes[2]].Expon*scalarP(Rik)/zeta_ik)
              -(BASIS_INFO[dm2[i].indexes[1]].Expon*BASIS_INFO[dm2[i].indexes[3]].Expon*scalarP(Rjl)/zeta_jl));
         e_ijkl=zeta_ik*zeta_jl/zeta_ijkl;
         //Redifine and recycle Rik and Rjl since we alredy used them:
         for(j=0;j<3;j++)
         {
          Rik[j]=(BASIS_INFO[dm2[i].indexes[0]].Expon*Coord_atom[0][j]+BASIS_INFO[dm2[i].indexes[2]].Expon*Coord_atom[2][j])/zeta_ik;
          Rjl[j]=(BASIS_INFO[dm2[i].indexes[1]].Expon*Coord_atom[1][j]+BASIS_INFO[dm2[i].indexes[3]].Expon*Coord_atom[3][j])/zeta_jl;
          Rijkl[j]=(-Rik[j]*zeta_ik+Rjl[j]*zeta_jl)/zeta_ijkl;
         }
         alpha_ijkl=(zeta_ik-zeta_jl)/(zeta_ijkl*TWO);
         nx_sum=nx_exp[0]+nx_exp[1]+nx_exp[2]+nx_exp[3];
         ny_sum=ny_exp[0]+ny_exp[1]+ny_exp[2]+ny_exp[3];
         nz_sum=nz_exp[0]+nz_exp[1]+nz_exp[2]+nz_exp[3];
         max_exp_ijkl=nx_sum;
         if(max_exp_ijkl<ny_sum)
         {max_exp_ijkl=ny_sum;}
         if(max_exp_ijkl<nz_sum)
         {max_exp_ijkl=nz_sum;}
         if(max_exp_ijkl%2==0)
         {
          order_ijkl=(max_exp_ijkl/2)+1;
         }
         else
         {
          order_ijkl=((max_exp_ijkl-1)/2)+1;
         }
         //Intracule coord. dependent part
         rscan=Init;
         for(j=0;j<nrad;j++)
         {
          for(k=0;k<nang;k++)
          {
           if(legendre)
           {
            Point_extra[0]=r_legendre[j]*x[k];
            Point_extra[1]=r_legendre[j]*y[k];
            Point_extra[2]=r_legendre[j]*z[k];
           }
           else
           {
            Point_extra[0]=rscan*x[k];
            Point_extra[1]=rscan*y[k];
            Point_extra[2]=rscan*z[k];
           }
           Xprime=TWO*pow(e_ijkl,HALF)*(Point_extra[0]-HALF*(Rik[0]+Rjl[0]));
           Yprime=TWO*pow(e_ijkl,HALF)*(Point_extra[1]-HALF*(Rik[1]+Rjl[1]));
           Zprime=TWO*pow(e_ijkl,HALF)*(Point_extra[2]-HALF*(Rik[2]+Rjl[2]));
           Tot_rad_th[j][k]=Tot_rad_th[j][k]+Aijkl*exp(-(pow(Xprime,TWO)+pow(Yprime,TWO)+pow(Zprime,TWO)))
                     *Vint3D_2(Point_extra,order_ijkl,r_extrac[order_ijkl-1],w_extrac[order_ijkl-1],
                      Rijkl,alpha_ijkl,zeta_ijkl,Coord_atom,nx_exp,ny_exp,nz_exp);
          }
          rscan=rscan+Step;
         }
        }
       }
       #pragma omp critical
       {
        for(i=0;i<nrad;i++)
        {
         for(j=0;j<nang;j++)
         {
          Tot_rad[i][j]=Tot_rad[i][j]+Tot_rad_th[i][j];
         }
        }
       }
       for(i=0;i<nrad;i++)
       {
        delete[] Tot_rad_th[i];Tot_rad_th[i]=NULL;
       }
       delete[] Tot_rad_th; Tot_rad_th=NULL;
       #pragma omp barrier
       #pragma omp master
       {
        for(i=0;i<nrad;i++)
        {
         Integral=ZERO;
         Integral2=ZERO;
         for(j=0;j<nang;j++)
         {
          Integral=Integral+w_theta_phi[j]*Tot_rad[i][j];
          if(abs(Tot_rad[i][j])>=pow(TEN,-TWO*SIX))
          {
           Integral2=Integral2-w_theta_phi[j]*abs(Tot_rad[i][j])*log(abs(Tot_rad[i][j]));
          }
         }
         Tot_rad[i][0]=DMNfact*FOUR*PI*Integral;
         Shannon_rad[i]=DMNfact*FOUR*PI*Integral2;
#ifdef HAVE_MPI
globalIr=ZERO;
localIr=Tot_rad[i][0];
MPI_Reduce(&localIr,&globalIr,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
Tot_rad[i][0]=globalIr;
Shannon_rad[i]=ZERO;
#endif
        }
       }
      }
     }
     else
     {
if(ID==0)
{
      cout<<"Running"<<setw(4)<<1<<" thread."<<endl;
}
      for(i=0;i<nterms;i++)
      {
       Cnorm_4gauss=ONE;
       for(j=0;j<4;j++)
       {
        //Compute the normalization factor of the gaussians (ijkl)
        Cnorm_4gauss=Cnorm_4gauss*
                    normGauss(BASIS_INFO[dm2[i].indexes[j]].Expon,
                              BASIS_INFO[dm2[i].indexes[j]].nx,
                              BASIS_INFO[dm2[i].indexes[j]].ny,
                              BASIS_INFO[dm2[i].indexes[j]].nz);
        //Take the coord. of this 4 gaussians
        Coord_atom[j][0]=BASIS_INFO[dm2[i].indexes[j]].Cart_cord[0];
        Coord_atom[j][1]=BASIS_INFO[dm2[i].indexes[j]].Cart_cord[1];
        Coord_atom[j][2]=BASIS_INFO[dm2[i].indexes[j]].Cart_cord[2];
        // Take the quantum numbers of the primitives
        nx_exp[j]=BASIS_INFO[dm2[i].indexes[j]].nx;
        ny_exp[j]=BASIS_INFO[dm2[i].indexes[j]].ny;
        nz_exp[j]=BASIS_INFO[dm2[i].indexes[j]].nz;
       }
       //Jab is the upper bound fucntion of Eq. 39
       Jik=Jab(BASIS_INFO[dm2[i].indexes[0]].Expon,
               BASIS_INFO[dm2[i].indexes[2]].Expon,
               BASIS_INFO[dm2[i].indexes[0]].nx,
               BASIS_INFO[dm2[i].indexes[0]].ny,
               BASIS_INFO[dm2[i].indexes[0]].nz,
               BASIS_INFO[dm2[i].indexes[2]].nx,
               BASIS_INFO[dm2[i].indexes[2]].ny,
               BASIS_INFO[dm2[i].indexes[2]].nz,
               BASIS_INFO[dm2[i].indexes[0]].Cart_cord[0],
               BASIS_INFO[dm2[i].indexes[0]].Cart_cord[1],
               BASIS_INFO[dm2[i].indexes[0]].Cart_cord[2],
               BASIS_INFO[dm2[i].indexes[2]].Cart_cord[0],
               BASIS_INFO[dm2[i].indexes[2]].Cart_cord[1],
               BASIS_INFO[dm2[i].indexes[2]].Cart_cord[2]);
       Jjl=Jab(BASIS_INFO[dm2[i].indexes[1]].Expon,
               BASIS_INFO[dm2[i].indexes[3]].Expon,
               BASIS_INFO[dm2[i].indexes[1]].nx,
               BASIS_INFO[dm2[i].indexes[1]].ny,
               BASIS_INFO[dm2[i].indexes[1]].nz,
               BASIS_INFO[dm2[i].indexes[3]].nx,
               BASIS_INFO[dm2[i].indexes[3]].ny,
               BASIS_INFO[dm2[i].indexes[3]].nz,
               BASIS_INFO[dm2[i].indexes[1]].Cart_cord[0],
               BASIS_INFO[dm2[i].indexes[1]].Cart_cord[1],
               BASIS_INFO[dm2[i].indexes[1]].Cart_cord[2],
               BASIS_INFO[dm2[i].indexes[3]].Cart_cord[0],
               BASIS_INFO[dm2[i].indexes[3]].Cart_cord[1],
               BASIS_INFO[dm2[i].indexes[3]].Cart_cord[2]);
       if(EIGHT*(abs(dm2[i].Dijkl)*Cnorm_4gauss*pow(Jik*Jjl,HALF))>=threshold2)
       {
        counter++;
        //Non extracule coord. dependent part
        zeta_ik=BASIS_INFO[dm2[i].indexes[0]].Expon+BASIS_INFO[dm2[i].indexes[2]].Expon;
        zeta_jl=BASIS_INFO[dm2[i].indexes[1]].Expon+BASIS_INFO[dm2[i].indexes[3]].Expon;
        zeta_ijkl=zeta_ik+zeta_jl;
        Aijkl=EIGHT*dm2[i].Dijkl*pow(zeta_ijkl,-HALF*THREE)*Cnorm_4gauss;
        for(j=0;j<3;j++)
        {
         Rik[j]=Coord_atom[0][j]-Coord_atom[2][j];
         Rjl[j]=Coord_atom[1][j]-Coord_atom[3][j];
        }
        Aijkl=Aijkl
             *exp(-(BASIS_INFO[dm2[i].indexes[0]].Expon*BASIS_INFO[dm2[i].indexes[2]].Expon*scalarP(Rik)/zeta_ik)
             -(BASIS_INFO[dm2[i].indexes[1]].Expon*BASIS_INFO[dm2[i].indexes[3]].Expon*scalarP(Rjl)/zeta_jl));
        e_ijkl=zeta_ik*zeta_jl/zeta_ijkl;
        //Redifine and recycle Rik and Rjl since we alredy used them:
        for(j=0;j<3;j++)
        {
         Rik[j]=(BASIS_INFO[dm2[i].indexes[0]].Expon*Coord_atom[0][j]+BASIS_INFO[dm2[i].indexes[2]].Expon*Coord_atom[2][j])/zeta_ik;
         Rjl[j]=(BASIS_INFO[dm2[i].indexes[1]].Expon*Coord_atom[1][j]+BASIS_INFO[dm2[i].indexes[3]].Expon*Coord_atom[3][j])/zeta_jl;
         Rijkl[j]=(-Rik[j]*zeta_ik+Rjl[j]*zeta_jl)/zeta_ijkl;
        }
        alpha_ijkl=(zeta_ik-zeta_jl)/(zeta_ijkl*TWO);
        nx_sum=nx_exp[0]+nx_exp[1]+nx_exp[2]+nx_exp[3];
        ny_sum=ny_exp[0]+ny_exp[1]+ny_exp[2]+ny_exp[3];
        nz_sum=nz_exp[0]+nz_exp[1]+nz_exp[2]+nz_exp[3];
        max_exp_ijkl=nx_sum;
        if(max_exp_ijkl<ny_sum)
        {max_exp_ijkl=ny_sum;}
        if(max_exp_ijkl<nz_sum)
        {max_exp_ijkl=nz_sum;}
        if(max_exp_ijkl%2==0)
        {
         order_ijkl=(max_exp_ijkl/2)+1;
        }
        else
        {
         order_ijkl=((max_exp_ijkl-1)/2)+1;
        }
        //Intracule coord. dependent part
        rscan=Init;
        for(j=0;j<nrad;j++)
        {
         for(k=0;k<nang;k++)
         {
          if(legendre)
          {
           Point_extra[0]=r_legendre[j]*x[k];
           Point_extra[1]=r_legendre[j]*y[k];
           Point_extra[2]=r_legendre[j]*z[k];
          }
          else
          {
           Point_extra[0]=rscan*x[k];
           Point_extra[1]=rscan*y[k];
           Point_extra[2]=rscan*z[k];
          }
          Xprime=TWO*pow(e_ijkl,HALF)*(Point_extra[0]-HALF*(Rik[0]+Rjl[0]));
          Yprime=TWO*pow(e_ijkl,HALF)*(Point_extra[1]-HALF*(Rik[1]+Rjl[1]));
          Zprime=TWO*pow(e_ijkl,HALF)*(Point_extra[2]-HALF*(Rik[2]+Rjl[2]));
          Tot_rad[j][k]=Tot_rad[j][k]+Aijkl*exp(-(pow(Xprime,TWO)+pow(Yprime,TWO)+pow(Zprime,TWO)))
                     *Vint3D_2(Point_extra,order_ijkl,r_extrac[order_ijkl-1],w_extrac[order_ijkl-1],
                      Rijkl,alpha_ijkl,zeta_ijkl,Coord_atom,nx_exp,ny_exp,nz_exp);
         }
         rscan=rscan+Step;
        }
       }
      }
     }
     ////////////////
     //Check times //
     ////////////////
if(ID==0)
{
     if(Input_commands.time_extra){cout<<"Time 3 (D H M S): "<<endl;system("date +%j' '%H' '%M' '%S ");cout<<endl;}
     ////////////////
     //End time    //
     ////////////////
     cout<<"#*****************************************************************************#";
     cout<<endl;
     cout<<"Computed terms from DM2:"<<setw(17)<<counter<<endl;
     cout<<"[Terms that passed |^2Dijkl Eijkl(R)| >= 2*tau/[Np(Np+1)]]"<<endl;
     cout<<"#*****************************************************************************#";
     cout<<endl;
     cout<<"#*****************************************************************************#";
     cout<<endl;
     cout<<"Let U=|r2+r1|/2 and E(U) = int E_3D(U,Omega_12) d Omega_12 "<<endl;
     cout<<"#*****************************************************************************#";
     cout<<endl;
     cout<<"#*****************************************************************************#";
     cout<<endl;
#ifdef HAVE_MPI
     cout<<"#    U              E(U)                E(U)U^2"<<endl;
#else
     cout<<"#    U              E(U)                E(U)U^2           S_M[E_3D(U,Omega_12)]"<<endl;
#endif
}
     for(i=0;i<6;i++){res[i]=ZERO;}
     rscan=Init;
     for(i=0;i<nrad;i++)
     {
      if(!parallel)
      {
       Integral=ZERO;
       Integral2=ZERO;
       for(j=0;j<nang;j++)
       {
        Integral=Integral+w_theta_phi[j]*Tot_rad[i][j];
        if(abs(Tot_rad[i][j])>=pow(TEN,-TWO*SIX))
        {
         Integral2=Integral2-w_theta_phi[j]*abs(Tot_rad[i][j])*log(abs(Tot_rad[i][j]));
        }
       }
       Tot_rad[i][0]=DMNfact*FOUR*PI*Integral;
       Shannon_rad[i]=DMNfact*FOUR*PI*Integral2;
#ifdef HAVE_MPI
globalIr=ZERO;
localIr=Tot_rad[i][0];
MPI_Reduce(&localIr,&globalIr,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
Tot_rad[i][0]=globalIr;
Shannon_rad[i]=ZERO;
#endif
      }
      if(abs(Tot_rad[i][0])<pow(TEN,-THREE*FIVE)){Tot_rad[i][0]=ZERO;}
      if(abs(Tot_rad[i][0])>ZERO)
      {
       if(legendre)
       {
if(ID==0)
{
        cout<<setprecision(8)<<fixed<<scientific<<r_legendre[i]<<setw(20)<<Tot_rad[i][0]<<setw(20);
#ifdef HAVE_MPI
        cout<<Tot_rad[i][0]*pow(r_legendre[i],TWO);
#else
        cout<<Tot_rad[i][0]*pow(r_legendre[i],TWO)<<setw(20)<<Shannon_rad[i];
#endif
        cout<<endl;
}
       }
       else
       {
if(ID==0)
{
        cout<<setprecision(8)<<fixed<<scientific<<rscan<<setw(20)<<Tot_rad[i][0]<<setw(20);
#ifdef HAVE_MPI
        cout<<Tot_rad[i][0]*pow(rscan,TWO);
#else
        cout<<Tot_rad[i][0]*pow(rscan,TWO)<<setw(20)<<Shannon_rad[i];
#endif
        cout<<endl;
}
       }
      }
      if(legendre)
      {
       res[0]=res[0]+Tot_rad[i][0]*Jacobian_legendre[i][0]*w_legendre[i];
       res[1]=res[1]+Tot_rad[i][0]*Jacobian_legendre[i][1]*w_legendre[i];
       res[2]=res[2]+Tot_rad[i][0]*Jacobian_legendre[i][2]*w_legendre[i];
       res[3]=res[3]+Tot_rad[i][0]*Jacobian_legendre[i][3]*w_legendre[i];
       res[4]=res[4]+Tot_rad[i][0]*Jacobian_legendre[i][4]*w_legendre[i];
       res[5]=res[5]+Shannon_rad[i]*Jacobian_legendre[i][2]*w_legendre[i];
      }
      else
      {
       if(i>0)
       {
        res[0]=res[0]+HALF*(Tot_rad[i-1][0]+Tot_rad[i][0])*Step;
        res[1]=res[1]+HALF*(Tot_rad[i-1][0]*(rscan-Step)+Tot_rad[i][0]*rscan)*Step;
        res[2]=res[2]+HALF*(Tot_rad[i-1][0]*pow(rscan-Step,TWO)+Tot_rad[i][0]*pow(rscan,TWO))*Step;
        res[3]=res[3]+HALF*(Tot_rad[i-1][0]*pow(rscan-Step,THREE)+Tot_rad[i][0]*pow(rscan,THREE))*Step;
        res[4]=res[4]+HALF*(Tot_rad[i-1][0]*pow(rscan-Step,FOUR)+Tot_rad[i][0]*pow(rscan,FOUR))*Step;
        res[5]=res[5]+HALF*(Shannon_rad[i-1]*pow(rscan-Step,TWO)+Shannon_rad[i]*pow(rscan,TWO))*Step;
       }
      }
      rscan=rscan+Step;
     }
if(ID==0)
{
     cout<<endl;
     cout<<"#################################################"<<endl;
     cout<<endl;
     if(legendre)
     {
      cout<<"Integration by Legendre Quadrature rule:"<<endl;
     }
     else
     {
      cout<<"Integration by Trapezoidal rule:"<<endl;
     }
     cout<<"int E(U)      dU         :"<<setw(18)<<res[0]<<endl;
     cout<<"int E(U) U    dU         :"<<setw(18)<<res[1]<<endl;
     cout<<"int E(U) U^2  dU         :"<<setw(18)<<res[2]<<endl;
     cout<<"int E(U) U^3  dU         :"<<setw(18)<<res[3]<<endl;
     cout<<"int E(U) U^4  dU         :"<<setw(18)<<res[4]<<endl;
     res[2]=round(res[2]);
     cout<<"Tr[^2Dijkl]              :"<<setw(18)<<res[2]<<endl;
     cout<<"[Trace is N(N-1)/2]"<<endl;
     cout<<"<  U  > per pair of e-   :"<<setw(18)<<res[3]/res[2]<<endl;
     cout<<"< U^2 > per pair of e-   :"<<setw(18)<<res[4]/res[2]<<endl;
#ifndef HAVE_MPI
     cout<<"S_Norm[E_3D(U,Omega_12)] :"<<setw(18)<<(res[2]*log(res[2])+res[5])/res[2]<<endl;
#endif
     cout<<"<U> and <U^2> are per pair of electrons to ensure that a PROBABILITY is used!"<<endl;
#ifndef HAVE_MPI
     cout<<"S_Norm[E_3D(U,Omega_12)] = (M*Log[M]+S_M[E_3D(U,Omega_12)])/M where M = Tr[^2Dijkl] and S_M is normalized to M"<<endl;
#endif
     cout<<endl;
}
     delete[] dm2;
if(ID==0)
{
     if(Input_commands.rweight)
     {
      if(legendre)
      {
       ofstream r_weights;
       r_weights.open("r_weights.txt");
       r_weights<<"#Radial_coordinate (U) \t   E(U)  \t\t Jacobian (for U^2) \t Weight"<<endl;
       for(i=0;i<nrad;i++)
       {
        if(abs(Tot_rad[i][0])>ZERO)
        {
         r_weights<<setprecision(10)<<fixed<<scientific<<r_legendre[i]<<"\t"<<Tot_rad[i][0]<<"\t"<<Jacobian_legendre[i][2]<<"\t"<<w_legendre[i]<<endl;
        }
       }
       r_weights.close();
       cout<<"r_weights.txt file was generated for the present job"<<endl;
       cout<<endl;
      }
      else
      {
       cout<<"No weights were generated for the radial scan"<<endl;
      }
     }
}
    }
    else
    {
if(ID==0)
{
     cout<<endl;
}
     open_dm2.close();
if(ID==0)
{
     cout<<"Error! Could not open file: "<<name_dm2<<endl;
     cout<<endl;
}
    }
    delete[] BASIS_INFO;
    for(i=0;i<Nroot_2Lmax_plus_1;i++)
    {
     delete[] r_extrac[i];r_extrac[i]=NULL;
     delete[] w_extrac[i];w_extrac[i]=NULL;
    }
    delete[] r_extrac;r_extrac=NULL;
    delete[] w_extrac;w_extrac=NULL;
   }
   else
   {
if(ID==0)
{
    cout<<endl;
}
    open_basis.close();
if(ID==0)
{
    cout<<"Error! Could not open file: "<<name_basis<<endl;
    cout<<endl;
}
   }
   delete[] x;
   delete[] y;
   delete[] z;
   delete[] w_theta_phi;
   for(i=0;i<nrad;i++)
   {
    delete[] Tot_rad[i];Tot_rad[i]=NULL;
   }
   delete[] Tot_rad; Tot_rad=NULL;
   delete[] Shannon_rad; Shannon_rad=NULL;
   if(legendre)
   {
    delete[] r_legendre; r_legendre=NULL;
    delete[] w_legendre; w_legendre=NULL;
    for(i=0;i<nrad;i++)
    {
     delete[] Jacobian_legendre[i]; //Jacobian_legendre[i]=NULL;
    }
    delete[] Jacobian_legendre; //Jacobian_legendre=NULL;
   }
  }
  //////////////////
  //Second moments//
  //////////////////
  if(Input_commands.second_moments)
  {
   cout<<"#******************************************#"<<endl;
   cout<<"#          Compute Moments                 #"<<endl;
   cout<<"#******************************************#"<<endl;
   string name_dm2=Input_commands.name_dm2;
   string name_basis=Input_commands.name_basis;
   threshold=Input_commands.threshold_in;
   parallel=Input_commands.parallel;
   // Read basis info
   ifstream open_basis;
   open_basis.open((name_basis).c_str());
   if(open_basis.good())
   {
    nprims=0;
    while(getline(open_basis,aux))
    {
     nprims++;
    }
    open_basis.close();
    nprims=nprims-1;
    basis_info *BASIS_INFO=new basis_info[nprims];
    open_basis.open((name_basis).c_str());
    getline(open_basis,aux);
    for(i=0;i<nprims;i++)
    {
     open_basis>>BASIS_INFO[i].Z>>BASIS_INFO[i].Cart_cord[0]>>BASIS_INFO[i].Cart_cord[1]>>BASIS_INFO[i].Cart_cord[2];
     open_basis>>BASIS_INFO[i].Expon>>BASIS_INFO[i].nx>>BASIS_INFO[i].ny>>BASIS_INFO[i].nz;
    }
    open_basis.close();
    //Compute Nroot_2Lmax_plus_1
    Lmax=0;
    for(i=0;i<nprims;i++)
    {
     if(BASIS_INFO[i].nx>Lmax)
     {
      Lmax=BASIS_INFO[i].nx;
     }
    }
    cout<<"Lmax                   :"<<setw(17)<<Lmax<<endl;
    //Create quadrature rule order
    Nroot_Lmax_plus_2=Lmax+2;
    cout<<"Quadrature rule order  :"<<setw(17)<<Nroot_Lmax_plus_2<<endl;
    cout<<"[Quadrature for primitive integrals. Order = Lmax + 2 ]"<<endl;
    alpha=ZERO;
    a=ZERO;
    b=ONE;
    gauss_hermite_rule(name_basis.substr(0,name_basis.length()-6),alpha,a,b,Nroot_Lmax_plus_2);
    r_mom=new double[Nroot_Lmax_plus_2];
    w_mom=new double[Nroot_Lmax_plus_2];
    //Read quadrature info
    ifstream read_quad;
    // Read weights
    read_quad.open((name_basis.substr(0,name_basis.length()-6)+"_w.txt").c_str());
    for(i=0;i<Nroot_Lmax_plus_2;i++)
    {
     read_quad>>w_mom[i];
    }
    read_quad.close();
    // Read roots
    read_quad.open((name_basis.substr(0,name_basis.length()-6)+"_x.txt").c_str());
    for(i=0;i<Nroot_Lmax_plus_2;i++)
    {
     read_quad>>r_mom[i];
    }
    read_quad.close();
    // Remove quadrature files
    system(("rm "+name_basis.substr(0,name_basis.length()-6)+"_r.txt").c_str());
    system(("rm "+name_basis.substr(0,name_basis.length()-6)+"_w.txt").c_str());
    system(("rm "+name_basis.substr(0,name_basis.length()-6)+"_x.txt").c_str());
    ifstream open_dm2;
    open_dm2.open((name_dm2).c_str());
    if(open_dm2.good())
    {
     open_dm2.close();
     nterms=0;
     //Check number of terms
     terms_dm2(name_dm2);
     dm2=new DM2[nterms];
     //Store the DM2
     cout<<"Storage of the DM2"<<endl;
     fill_in_dm2(name_dm2);
     cout<<"Storage done!"<<endl;
     //Compute moments
     cout<<"#*****************************************************************************#";
     cout<<endl;
     for(i=0;i<16;i++){res[i]=ZERO;}
     cout<<"Running"<<setw(4)<<1<<" thread."<<endl;
     for(i=0;i<nterms;i++)
     {
      Cnorm_4gauss=ONE;
      for(j=0;j<4;j++)
      {
       //Compute the normalization factor of the gaussians (ijkl)
       Cnorm_4gauss=Cnorm_4gauss*
                    normGauss(BASIS_INFO[dm2[i].indexes[j]].Expon,
                              BASIS_INFO[dm2[i].indexes[j]].nx,
                              BASIS_INFO[dm2[i].indexes[j]].ny,
                              BASIS_INFO[dm2[i].indexes[j]].nz);
        //Take the coord. of this 4 gaussians
        Coord_atom[j][0]=BASIS_INFO[dm2[i].indexes[j]].Cart_cord[0];
        Coord_atom[j][1]=BASIS_INFO[dm2[i].indexes[j]].Cart_cord[1];
        Coord_atom[j][2]=BASIS_INFO[dm2[i].indexes[j]].Cart_cord[2];
        // Take the quantum numbers of the primitives
        nx_exp[j]=BASIS_INFO[dm2[i].indexes[j]].nx;
        ny_exp[j]=BASIS_INFO[dm2[i].indexes[j]].ny;
        nz_exp[j]=BASIS_INFO[dm2[i].indexes[j]].nz;
       }
      // Here zeta_ijkl = (z_ik * z_jl)^-3/2
      zeta_ijkl=pow((BASIS_INFO[dm2[i].indexes[0]].Expon+BASIS_INFO[dm2[i].indexes[2]].Expon)*(BASIS_INFO[dm2[i].indexes[1]].Expon+BASIS_INFO[dm2[i].indexes[3]].Expon),-HALF*THREE);
      // Ri-Rk and Rj-Rl
      for(j=0;j<3;j++)
      {
       Rik[j]=Coord_atom[0][j]-Coord_atom[2][j];
       Rjl[j]=Coord_atom[1][j]-Coord_atom[3][j];
      }
      // Compute Aijkl -> Aijkl x Integrals is the primitive contribution for the moments.
      Aijkl=dm2[i].Dijkl*Cnorm_4gauss*zeta_ijkl
            *exp(-((BASIS_INFO[dm2[i].indexes[0]].Expon*BASIS_INFO[dm2[i].indexes[2]].Expon)/(BASIS_INFO[dm2[i].indexes[0]].Expon+BASIS_INFO[dm2[i].indexes[2]].Expon))*scalarP(Rik))
            *exp(-((BASIS_INFO[dm2[i].indexes[1]].Expon*BASIS_INFO[dm2[i].indexes[3]].Expon)/(BASIS_INFO[dm2[i].indexes[1]].Expon+BASIS_INFO[dm2[i].indexes[3]].Expon))*scalarP(Rjl));
      //Redifine and recycle Rik and Rjl since we alredy used them:
      for(j=0;j<3;j++)
      {
       Rik[j]=BASIS_INFO[dm2[i].indexes[0]].Expon*Coord_atom[0][j]+BASIS_INFO[dm2[i].indexes[2]].Expon*Coord_atom[2][j];
       Rjl[j]=BASIS_INFO[dm2[i].indexes[1]].Expon*Coord_atom[1][j]+BASIS_INFO[dm2[i].indexes[3]].Expon*Coord_atom[3][j];
       Rik[j]=Rik[j]/(BASIS_INFO[dm2[i].indexes[0]].Expon+BASIS_INFO[dm2[i].indexes[2]].Expon);
       Rjl[j]=Rjl[j]/(BASIS_INFO[dm2[i].indexes[1]].Expon+BASIS_INFO[dm2[i].indexes[3]].Expon);
      }
      //r1
      Sik[0]=Sab(Nroot_Lmax_plus_2,r_mom,w_mom,Rik[0],Coord_atom[0][0],Coord_atom[2][0],BASIS_INFO[dm2[i].indexes[0]].Expon,BASIS_INFO[dm2[i].indexes[2]].Expon,nx_exp[0],nx_exp[2]);
      Sik[1]=Sab(Nroot_Lmax_plus_2,r_mom,w_mom,Rik[1],Coord_atom[0][1],Coord_atom[2][1],BASIS_INFO[dm2[i].indexes[0]].Expon,BASIS_INFO[dm2[i].indexes[2]].Expon,ny_exp[0],ny_exp[2]);
      Sik[2]=Sab(Nroot_Lmax_plus_2,r_mom,w_mom,Rik[2],Coord_atom[0][2],Coord_atom[2][2],BASIS_INFO[dm2[i].indexes[0]].Expon,BASIS_INFO[dm2[i].indexes[2]].Expon,nz_exp[0],nz_exp[2]);
      Mik[0]=Mab(Nroot_Lmax_plus_2,r_mom,w_mom,Rik[0],Coord_atom[0][0],Coord_atom[2][0],BASIS_INFO[dm2[i].indexes[0]].Expon,BASIS_INFO[dm2[i].indexes[2]].Expon,nx_exp[0],nx_exp[2]);
      Mik[1]=Mab(Nroot_Lmax_plus_2,r_mom,w_mom,Rik[1],Coord_atom[0][1],Coord_atom[2][1],BASIS_INFO[dm2[i].indexes[0]].Expon,BASIS_INFO[dm2[i].indexes[2]].Expon,ny_exp[0],ny_exp[2]);
      Mik[2]=Mab(Nroot_Lmax_plus_2,r_mom,w_mom,Rik[2],Coord_atom[0][2],Coord_atom[2][2],BASIS_INFO[dm2[i].indexes[0]].Expon,BASIS_INFO[dm2[i].indexes[2]].Expon,nz_exp[0],nz_exp[2]);
      M2ik[0]=M2ab(Nroot_Lmax_plus_2,r_mom,w_mom,Rik[0],Coord_atom[0][0],Coord_atom[2][0],BASIS_INFO[dm2[i].indexes[0]].Expon,BASIS_INFO[dm2[i].indexes[2]].Expon,nx_exp[0],nx_exp[2]);
      M2ik[1]=M2ab(Nroot_Lmax_plus_2,r_mom,w_mom,Rik[1],Coord_atom[0][1],Coord_atom[2][1],BASIS_INFO[dm2[i].indexes[0]].Expon,BASIS_INFO[dm2[i].indexes[2]].Expon,ny_exp[0],ny_exp[2]);
      M2ik[2]=M2ab(Nroot_Lmax_plus_2,r_mom,w_mom,Rik[2],Coord_atom[0][2],Coord_atom[2][2],BASIS_INFO[dm2[i].indexes[0]].Expon,BASIS_INFO[dm2[i].indexes[2]].Expon,nz_exp[0],nz_exp[2]);
      //r2
      Sjl[0]=Sab(Nroot_Lmax_plus_2,r_mom,w_mom,Rjl[0],Coord_atom[1][0],Coord_atom[3][0],BASIS_INFO[dm2[i].indexes[1]].Expon,BASIS_INFO[dm2[i].indexes[3]].Expon,nx_exp[1],nx_exp[3]);
      Sjl[1]=Sab(Nroot_Lmax_plus_2,r_mom,w_mom,Rjl[1],Coord_atom[1][1],Coord_atom[3][1],BASIS_INFO[dm2[i].indexes[1]].Expon,BASIS_INFO[dm2[i].indexes[3]].Expon,ny_exp[1],ny_exp[3]);
      Sjl[2]=Sab(Nroot_Lmax_plus_2,r_mom,w_mom,Rjl[2],Coord_atom[1][2],Coord_atom[3][2],BASIS_INFO[dm2[i].indexes[1]].Expon,BASIS_INFO[dm2[i].indexes[3]].Expon,nz_exp[1],nz_exp[3]);
      Mjl[0]=Mab(Nroot_Lmax_plus_2,r_mom,w_mom,Rjl[0],Coord_atom[1][0],Coord_atom[3][0],BASIS_INFO[dm2[i].indexes[1]].Expon,BASIS_INFO[dm2[i].indexes[3]].Expon,nx_exp[1],nx_exp[3]);
      Mjl[1]=Mab(Nroot_Lmax_plus_2,r_mom,w_mom,Rjl[1],Coord_atom[1][1],Coord_atom[3][1],BASIS_INFO[dm2[i].indexes[1]].Expon,BASIS_INFO[dm2[i].indexes[3]].Expon,ny_exp[1],ny_exp[3]);
      Mjl[2]=Mab(Nroot_Lmax_plus_2,r_mom,w_mom,Rjl[2],Coord_atom[1][2],Coord_atom[3][2],BASIS_INFO[dm2[i].indexes[1]].Expon,BASIS_INFO[dm2[i].indexes[3]].Expon,nz_exp[1],nz_exp[3]);
      //rho2
      res[0]=res[0]+Aijkl*Sik[0]*Sik[1]*Sik[2]*Sjl[0]*Sjl[1]*Sjl[2];
      //rho2 x1
      res[1]=res[1]+Aijkl*Mik[0]*Sik[1]*Sik[2]*Sjl[0]*Sjl[1]*Sjl[2];
      //rho2 y1
      res[2]=res[2]+Aijkl*Sik[0]*Mik[1]*Sik[2]*Sjl[0]*Sjl[1]*Sjl[2];
      //rho2 z1
      res[3]=res[3]+Aijkl*Sik[0]*Sik[1]*Mik[2]*Sjl[0]*Sjl[1]*Sjl[2];
      //rho2 x1x1
      res[4]=res[4]+Aijkl*M2ik[0]*Sik[1]*Sik[2]*Sjl[0]*Sjl[1]*Sjl[2];
      //rho2 x1y1
      res[5]=res[5]+Aijkl*Mik[0]*Mik[1]*Sik[2]*Sjl[0]*Sjl[1]*Sjl[2];
      //rho2 x1z1
      res[6]=res[6]+Aijkl*Mik[0]*Sik[1]*Mik[2]*Sjl[0]*Sjl[1]*Sjl[2];
      //rho2 y1y1
      res[7]=res[7]+Aijkl*Sik[0]*M2ik[1]*Sik[2]*Sjl[0]*Sjl[1]*Sjl[2];
      //rho2 y1z1
      res[8]=res[8]+Aijkl*Sik[0]*Mik[1]*Mik[2]*Sjl[0]*Sjl[1]*Sjl[2];
      //rho2 z1z1
      res[9]=res[9]+Aijkl*Sik[0]*Sik[1]*M2ik[2]*Sjl[0]*Sjl[1]*Sjl[2];
      //rho2 x1x2
      res[10]=res[10]+Aijkl*Mik[0]*Sik[1]*Sik[2]*Mjl[0]*Sjl[1]*Sjl[2];
      //rho2 x1y2
      res[11]=res[11]+Aijkl*Mik[0]*Sik[1]*Sik[2]*Sjl[0]*Mjl[1]*Sjl[2];
      //rho2 x1z2
      res[12]=res[12]+Aijkl*Mik[0]*Sik[1]*Sik[2]*Sjl[0]*Sjl[1]*Mjl[2];
      //rho2 y1y2
      res[13]=res[13]+Aijkl*Sik[0]*Mik[1]*Sik[2]*Sjl[0]*Mjl[1]*Sjl[2];
      //rho2 y1z2
      res[14]=res[14]+Aijkl*Sik[0]*Mik[1]*Sik[2]*Sjl[0]*Sjl[1]*Mjl[2];
      //rho2 z1z2
      res[15]=res[15]+Aijkl*Sik[0]*Sik[1]*Mik[2]*Sjl[0]*Sjl[1]*Mjl[2];
     }
     //Print results
     for(i=0;i<16;i++)
     {
      if(abs(res[i])<pow(TEN,-TEN))
      {res[i]=ZERO;}
     }
     Nelectrons=HALF*(ONE+sqrt(ONE+FOUR*res[0]));
     cout<<"#*****************************************************************************#";
     cout<<endl;
     cout<<setprecision(10)<<fixed<<scientific;
     cout<<"int rho2(1,2)    d1 d2    :"<<setw(18)<<res[0]<<endl;
     cout<<"int rho2(1,2) x1 d1 d2    :"<<setw(18)<<res[1]<<endl;
     cout<<"int rho2(1,2) y1 d1 d2    :"<<setw(18)<<res[2]<<endl;
     cout<<"int rho2(1,2) z1 d1 d2    :"<<setw(18)<<res[3]<<endl;
     cout<<"int rho2(1,2) x1x1 d1 d2  :"<<setw(18)<<res[4]<<endl;
     cout<<"int rho2(1,2) x1y1 d1 d2  :"<<setw(18)<<res[5]<<endl;
     cout<<"int rho2(1,2) x1z1 d1 d2  :"<<setw(18)<<res[6]<<endl;
     cout<<"int rho2(1,2) y1y1 d1 d2  :"<<setw(18)<<res[7]<<endl;
     cout<<"int rho2(1,2) y1z1 d1 d2  :"<<setw(18)<<res[8]<<endl;
     cout<<"int rho2(1,2) z1z1 d1 d2  :"<<setw(18)<<res[9]<<endl;
     cout<<"int rho2(1,2) x1x2 d1 d2  :"<<setw(18)<<res[10]<<endl;
     cout<<"int rho2(1,2) x1y2 d1 d2  :"<<setw(18)<<res[11]<<endl;
     cout<<"int rho2(1,2) x1z2 d1 d2  :"<<setw(18)<<res[12]<<endl;
     cout<<"int rho2(1,2) y1y2 d1 d2  :"<<setw(18)<<res[13]<<endl;
     cout<<"int rho2(1,2) y1z2 d1 d2  :"<<setw(18)<<res[14]<<endl;
     cout<<"int rho2(1,2) z1z2 d1 d2  :"<<setw(18)<<res[15]<<endl;
     cout<<"int rho(1) d1             :"<<setw(18)<<Nelectrons<<endl;
     cout<<endl;
     cout<<"Total Position Spread  :"<<endl;
     cout<<"[O. Brea et al, J. Chem. Theory Comput., 9, 5286 (2013).]"<<endl;
     cout<<endl;
     //xx
     TPS[0][0]=res[4]/(Nelectrons-ONE)+res[10]-res[1]*res[1]/pow(Nelectrons-ONE,TWO);
     //xy
     TPS[0][1]=res[5]/(Nelectrons-ONE)+res[11]-res[1]*res[2]/pow(Nelectrons-ONE,TWO);
     TPS[1][0]=res[5]/(Nelectrons-ONE)+res[11]-res[1]*res[2]/pow(Nelectrons-ONE,TWO);
     //xz
     TPS[0][2]=res[6]/(Nelectrons-ONE)+res[12]-res[1]*res[3]/pow(Nelectrons-ONE,TWO);
     TPS[2][0]=res[6]/(Nelectrons-ONE)+res[12]-res[1]*res[3]/pow(Nelectrons-ONE,TWO);
     //yy
     TPS[1][1]=res[7]/(Nelectrons-ONE)+res[13]-res[2]*res[2]/pow(Nelectrons-ONE,TWO);
     //yz
     TPS[1][2]=res[8]/(Nelectrons-ONE)+res[14]-res[2]*res[3]/pow(Nelectrons-ONE,TWO);
     TPS[2][1]=res[8]/(Nelectrons-ONE)+res[14]-res[2]*res[3]/pow(Nelectrons-ONE,TWO);
     //zz
     TPS[2][2]=res[9]/(Nelectrons-ONE)+res[15]-res[3]*res[3]/pow(Nelectrons-ONE,TWO);
     for(i=0;i<3;i++)
     {
      for(j=0;j<3;j++)
      {
       if(abs(TPS[i][j])<pow(TEN,-EIGHT)){TPS[i][j]=ZERO;}
      }
     }
     //Print TPS
     cout<<" ( "<<setw(18)<<TPS[0][0]<<" "<<setw(18)<<TPS[0][1]<<" ";
     cout<<setw(18)<<TPS[0][2]<<" )"<<endl;
     cout<<" ( "<<setw(18)<<TPS[1][0]<<" "<<setw(18)<<TPS[1][1]<<" ";
     cout<<setw(18)<<TPS[1][2]<<" )"<<endl;
     cout<<" ( "<<setw(18)<<TPS[2][0]<<" "<<setw(18)<<TPS[2][1]<<" ";
     cout<<setw(18)<<TPS[2][2]<<" )"<<endl;
     cout<<endl;
     cout<<endl;
     delete[] dm2;
    }
    else
    {
     cout<<endl;
     open_dm2.close();
     cout<<"Error! Could not open file: "<<name_dm2<<endl;
     cout<<endl;
    }
    delete[] BASIS_INFO;
    delete[] r_mom;
    delete[] w_mom;
   }
   else
   {
    cout<<endl;
    open_basis.close();
    cout<<"Error! Could not open file: "<<name_basis<<endl;
    cout<<endl;
   }
  }
  ////////////////////////////////////
  //Not Second moments nor Intracule//
  ////////////////////////////////////
  if(!Input_commands.second_moments && !Input_commands.intracule && !Input_commands.extracule)
  {
   cout<<"Please include either the $Intracule, $Extracule or $Moments keyword in the input!"<<endl;
  }
  ////////////////
  //Check times //
  ////////////////
  if(Input_commands.time_extra){cout<<"Time 4 (D H M S): "<<endl;system("date +%j' '%H' '%M' '%S ");cout<<endl;}
  ////////////////
  //End time    //
  ////////////////
 }
 else
 {
  cout<<endl;
  cout<<"Please include input commands file"<<endl;
  cout<<endl;
  cout<<endl;
 }
if(ID==0)
{
 /////////////////////
 // End clock       //
 /////////////////////
 system("date +%j' '%H' '%M' '%S>date_intrac.date");
 date_file.open("date_intrac.date");
 date_file>>DATE[1][0]>>DATE[1][1]>>DATE[1][2]>>DATE[1][3];
 date_file.close();
 system("rm date_intrac.date");
 calc_time(DATE);
 cout<<"#################################################"<<endl;
 cout<<endl;
 cout<<"                             __...------------._"<<endl;
 cout<<"                         ,-'                   `-.  "<<endl;
 cout<<"  It is a Trap!       ,-'                         `."<<endl;
 cout<<"                    ,'                            ,-`."<<endl;
 cout<<"                   ;                              `-' `."<<endl;
 cout<<"                  ;                                 .-. \\"<<endl;
 cout<<"                 ;                           .-.    `-'  \\"<<endl;
 cout<<"                ;                            `-'          \\"<<endl;
 cout<<"               ;                                          `."<<endl;
 cout<<"               ;                                           :"<<endl;
 cout<<"              ;                                            |"<<endl;
 cout<<"             ;                                             ;"<<endl;
 cout<<"            ;                            ___              ;"<<endl;
 cout<<"           ;                        ,-;-','.`.__          |"<<endl;
 cout<<"       _..;                      ,-' ;`,'.`,'.--`.        |"<<endl;
 cout<<"      ///;           ,-'   `. ,-'   ;` ;`,','_.--=:      /"<<endl;
 cout<<"     |'':          ,'        :     ;` ;,;,,-'_.-._`.   ,'"<<endl;
 cout<<"     '  :         ;_.-.      `.    :' ;;;'.ee.    \\|  /"<<endl;
 cout<<"      \\.'    _..-'/8o. `.     :    :! ' ':8888)   || /"<<endl;
 cout<<"       ||`-''    \\88o\\ :     :    :! :  :`""'    ;;/"<<endl;
 cout<<"       ||         \\\"88o\\;     `.    \\ `. `.      ;,'"<<endl;
 cout<<"       /)   ___    `.\"'/(--.._ `.    `.`.  `-..-' ;--."<<endl;
 cout<<"       \\(.=\"\"\"\"\"==.. `'-'     `.|      `-`-..__.-' `. `."<<endl;
 cout<<"        |          `\"==.__      )                    )  ;"<<endl;
 cout<<"        |   ||           `\"=== '                   .'  .'"<<endl;
 cout<<"        /\\,,||||  | |           \\                .'   .'"<<endl;
 cout<<"        | |||'|' |'|'           \\|             .'   _.' \\"<<endl;
 cout<<"        | |\' |  |           || ||           .'    .'    \\"<<endl;
 cout<<"        ' | \\ ' |'  .   ``-- `| ||         .'    .'       \\"<<endl;
 cout<<"          '  |  ' |  .    ``-.._ |  ;    .'    .'          `."<<endl;
 cout<<"       _.--,;`.       .  --  ...._,'   .'    .'              `.__"<<endl;
 cout<<"     ,'  ,';   `.     .   --..__..--'.'    .'                __/_\\"<<endl;
 cout<<"   ,'   ; ;     |    .   --..__.._.'     .'                ,'     `."<<endl;
 cout<<"  /    ; :     ;     .    -.. _.'     _.'                 /         `"<<endl;
 cout<<" /     :  `-._ |    .    _.--'     _.'                   |"<<endl;
 cout<<"/       `.    `--....--''       _.'                      |"<<endl;
 cout<<"          `._              _..-'                         |"<<endl;
 cout<<"             `-..____...-''                              |"<<endl;
 cout<<"                                                         |"<<endl;
 cout<<"                               mGk                       |"<<endl;
 cout<<endl;
 cout<<endl;
 cout<<"#################################################"<<endl;
 cout<<endl;
 cout<<setprecision(0)<<fixed;
 cout<<"  It took "<<setw(4)<<DATE[0][0]<<" days "<<setw(4)<<DATE[0][1]<<" hours "<<setw(4)<<DATE[0][2]<<" min. ";
 cout<<setw(4)<<DATE[0][3]<<" secs. "<<endl;
 cout<<endl;
 cout<<"#################################################"<<endl;
 cout<<endl;
 cout<<" Normal termination of RHO2_OPS code          "<<endl;
 cout<<endl;
 cout<<"#################################################"<<endl;
}
#ifdef HAVE_MPI
MPI::Finalize();
#endif
 return 0;
}

//Compute norm factor of the primitive [calls dfact() to do the double factorial]
double normGauss(double &expon, int &n, int &l,int &m)
{
 double norm;
 int n_df,m_df,l_df;
 norm=(double)n+l+m;
 n_df=dfact(2*n-1);
 l_df=dfact(2*l-1);
 m_df=dfact(2*m-1);
 norm=pow((TWO*expon/PI),(THREE/FOUR))*pow((FOUR*expon),(norm/TWO));
 norm=norm/pow(n_df*m_df*l_df,HALF);
 return norm;
}

//Check number of terms in the DM2 and store them
void terms_dm2(string name_file)
{
 int element[2],element_prime[2];
 double Dijkl;//,Trace=ZERO;
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
   if(abs(Dijkl)>=threshold)
   {
    nterms++;
//    if(element[0]==element_prime[0] && element[1]==element_prime[1])
//    {Trace=Trace+Dijkl;}
   }
  }
if(ID==0)
{
  cout<<"Number of terms in DM2 : "<<setw(17)<<nterms<<endl;
}
  input_data.close();
 }
}

void fill_in_dm2(string name_file)
{
 int element[2],element_prime[2];
 long int iterm=0,lnterms;
 double Dijkl;
 lnterms=ID;
 element[0]=10;element[1]=10;element_prime[0]=10;element_prime[1]=10;
 ifstream input_data(name_file.c_str(), ios::binary);
 if(input_data.good())
 {
  nterms=0;
  while((element[0]!=0 || element_prime[0]!=0) || (element[1]!=0 || element_prime[1]!=0))
  {
   input_data.seekg(RECORD_DELIMITER_LENGTH, ios::cur);
   input_data.read((char*) &element[0], sizeof(element[0]));
   input_data.read((char*) &element[1], sizeof(element[1]));
   input_data.read((char*) &element_prime[0], sizeof(element_prime[0]));
   input_data.read((char*) &element_prime[1], sizeof(element_prime[1]));
   input_data.read((char*) &Dijkl, sizeof(Dijkl));
   input_data.seekg(RECORD_DELIMITER_LENGTH, ios::cur);
   if(abs(Dijkl)>=threshold)
   {
    iterm++;
    if((iterm-1)==lnterms)
    {
     dm2[nterms].indexes[0]=element[0]-1;
     dm2[nterms].indexes[1]=element[1]-1;
     dm2[nterms].indexes[2]=element_prime[0]-1;
     dm2[nterms].indexes[3]=element_prime[1]-1;
     dm2[nterms].Dijkl=Dijkl;
     nterms++;
     lnterms=lnterms+nproc;
    }
   }
  }
  input_data.close();
 }
}

//Vijkl(S) is the exact integrated (in s coord.) value of the polinomyal (Eq. 16 in Cioslowski's paper)
/*double Vint(double &S, int &order, double *s, double *w,double &Sijkl, double &alpha_ijkl,double &zeta_ijkl,
double &Si,double &Sj,double &Sk, double &Sl, int &nxi,int &nxj,int &nxk,int &nxl)
{
 int i;
 double res=ZERO;
 for(i=0;i<order;i++)
 {
  res=res+w[i]
     *pow(pow(zeta_ijkl,-HALF)*s[i]+(alpha_ijkl-HALF)*S+(Sijkl-Si),(double)nxi)
     *pow(pow(zeta_ijkl,-HALF)*s[i]+(alpha_ijkl+HALF)*S+(Sijkl-Sj),(double)nxj)
     *pow(pow(zeta_ijkl,-HALF)*s[i]+(alpha_ijkl-HALF)*S+(Sijkl-Sk),(double)nxk)
     *pow(pow(zeta_ijkl,-HALF)*s[i]+(alpha_ijkl+HALF)*S+(Sijkl-Sl),(double)nxl);
 }
 return res;
}*/

//Vijkl(S) is the exact integrated (in s coord.) value of the polinomyal (Eq. 16 in Cioslowski's paper) 3D version.
double Vint3D(double S_intra[3], int &order, double *s, double *w,double Sijkl[3], double &alpha_ijkl,double &zeta_ijkl,
double Primitive_coords[4][3], int nsx_ijkl[4],int nsy_ijkl[4],int nsz_ijkl[4])
{
 int i,j;
 double res[3]={ZERO},zeta_ijkl_minus_half,zeta_ijkl_minus_half_s,alpha_ijkl_min_half,alpha_ijkl_plus_half;
 double alpha_ijkl_min_half_S[3],alpha_ijkl_plus_half_S[3];
 zeta_ijkl_minus_half=pow(zeta_ijkl,-HALF);
 alpha_ijkl_min_half=alpha_ijkl-HALF;
 alpha_ijkl_plus_half=alpha_ijkl+HALF;
 for(i=0;i<order;i++)
 {
  zeta_ijkl_minus_half_s=zeta_ijkl_minus_half*s[i];
  for(j=0;j<3;j++)
  {
   alpha_ijkl_min_half_S[j]=alpha_ijkl_min_half*S_intra[j];
   alpha_ijkl_plus_half_S[j]=alpha_ijkl_plus_half*S_intra[j];
  }
  //The point s where the quadrature is evaluated is fixed for all
  // This is Vijkl(X)
  res[0]=res[0]+w[i]
     *pow(zeta_ijkl_minus_half_s+alpha_ijkl_min_half_S[0] +(Sijkl[0]-Primitive_coords[0][0]),(double)nsx_ijkl[0])
     *pow(zeta_ijkl_minus_half_s+alpha_ijkl_plus_half_S[0]+(Sijkl[0]-Primitive_coords[1][0]),(double)nsx_ijkl[1])
     *pow(zeta_ijkl_minus_half_s+alpha_ijkl_min_half_S[0] +(Sijkl[0]-Primitive_coords[2][0]),(double)nsx_ijkl[2])
     *pow(zeta_ijkl_minus_half_s+alpha_ijkl_plus_half_S[0]+(Sijkl[0]-Primitive_coords[3][0]),(double)nsx_ijkl[3]);
  // This is  Vijkl(Y)
  res[1]=res[1]+w[i]
     *pow(zeta_ijkl_minus_half_s+alpha_ijkl_min_half_S[1] +(Sijkl[1]-Primitive_coords[0][1]),(double)nsy_ijkl[0])
     *pow(zeta_ijkl_minus_half_s+alpha_ijkl_plus_half_S[1]+(Sijkl[1]-Primitive_coords[1][1]),(double)nsy_ijkl[1])
     *pow(zeta_ijkl_minus_half_s+alpha_ijkl_min_half_S[1] +(Sijkl[1]-Primitive_coords[2][1]),(double)nsy_ijkl[2])
     *pow(zeta_ijkl_minus_half_s+alpha_ijkl_plus_half_S[1]+(Sijkl[1]-Primitive_coords[3][1]),(double)nsy_ijkl[3]);
  // This is  Vijkl(Z)
  res[2]=res[2]+w[i]
     *pow(zeta_ijkl_minus_half_s+alpha_ijkl_min_half_S[2] +(Sijkl[2]-Primitive_coords[0][2]),(double)nsz_ijkl[0])
     *pow(zeta_ijkl_minus_half_s+alpha_ijkl_plus_half_S[2]+(Sijkl[2]-Primitive_coords[1][2]),(double)nsz_ijkl[1])
     *pow(zeta_ijkl_minus_half_s+alpha_ijkl_min_half_S[2] +(Sijkl[2]-Primitive_coords[2][2]),(double)nsz_ijkl[2])
     *pow(zeta_ijkl_minus_half_s+alpha_ijkl_plus_half_S[2]+(Sijkl[2]-Primitive_coords[3][2]),(double)nsz_ijkl[3]);
 }
 return res[0]*res[1]*res[2];
}

//Vijkl2(S) is the exact integrated (in s coord.) value of the polinomyal (Eq. 16 in Cioslowski's paper) 3D version.
double Vint3D_2(double S_extra[3], int &order, double *s, double *w,double Sijkl[3], double &alpha_ijkl,double &zeta_ijkl,
double Primitive_coords[4][3], int nsx_ijkl[4],int nsy_ijkl[4],int nsz_ijkl[4])
{
 int i,j;
 double res[3]={ZERO},zeta_ijkl_minus_half,zeta_ijkl_minus_half_s,minus_two_alpha_ijkl_plus_one,two_alpha_ijkl_plus_one;
 double min_two_alpha_ijkl_plus_one_S[3],two_alpha_ijkl_plus_one_S[3];
 zeta_ijkl_minus_half=pow(zeta_ijkl,-HALF);
 minus_two_alpha_ijkl_plus_one=-TWO*alpha_ijkl+ONE;
 two_alpha_ijkl_plus_one=TWO*alpha_ijkl+ONE;
 for(i=0;i<order;i++)
 {
  zeta_ijkl_minus_half_s=zeta_ijkl_minus_half*s[i];
  for(j=0;j<3;j++)
  {
   min_two_alpha_ijkl_plus_one_S[j]=minus_two_alpha_ijkl_plus_one*S_extra[j];
   two_alpha_ijkl_plus_one_S[j]=two_alpha_ijkl_plus_one*S_extra[j];
  }
  //The point s where the quadrature is evaluated is fixed for all
  // This is Vijkl(X)
  res[0]=res[0]+w[i]
     *pow(-zeta_ijkl_minus_half_s+min_two_alpha_ijkl_plus_one_S[0]-(Sijkl[0]+Primitive_coords[0][0]),(double)nsx_ijkl[0])
     *pow( zeta_ijkl_minus_half_s+two_alpha_ijkl_plus_one_S[0]    +(Sijkl[0]-Primitive_coords[1][0]),(double)nsx_ijkl[1])
     *pow(-zeta_ijkl_minus_half_s+min_two_alpha_ijkl_plus_one_S[0]-(Sijkl[0]+Primitive_coords[2][0]),(double)nsx_ijkl[2])
     *pow( zeta_ijkl_minus_half_s+two_alpha_ijkl_plus_one_S[0]    +(Sijkl[0]-Primitive_coords[3][0]),(double)nsx_ijkl[3]);
  // This is  Vijkl(Y)
  res[1]=res[1]+w[i]
     *pow(-zeta_ijkl_minus_half_s+min_two_alpha_ijkl_plus_one_S[1]-(Sijkl[1]+Primitive_coords[0][1]),(double)nsy_ijkl[0])
     *pow( zeta_ijkl_minus_half_s+two_alpha_ijkl_plus_one_S[1]    +(Sijkl[1]-Primitive_coords[1][1]),(double)nsy_ijkl[1])
     *pow(-zeta_ijkl_minus_half_s+min_two_alpha_ijkl_plus_one_S[1]-(Sijkl[1]+Primitive_coords[2][1]),(double)nsy_ijkl[2])
     *pow( zeta_ijkl_minus_half_s+two_alpha_ijkl_plus_one_S[1]    +(Sijkl[1]-Primitive_coords[3][1]),(double)nsy_ijkl[3]);
  // This is  Vijkl(Z)
  res[2]=res[2]+w[i]
     *pow(-zeta_ijkl_minus_half_s+min_two_alpha_ijkl_plus_one_S[2]-(Sijkl[2]+Primitive_coords[0][2]),(double)nsz_ijkl[0])
     *pow( zeta_ijkl_minus_half_s+two_alpha_ijkl_plus_one_S[2]    +(Sijkl[2]-Primitive_coords[1][2]),(double)nsz_ijkl[1])
     *pow(-zeta_ijkl_minus_half_s+min_two_alpha_ijkl_plus_one_S[2]-(Sijkl[2]+Primitive_coords[2][2]),(double)nsz_ijkl[2])
     *pow( zeta_ijkl_minus_half_s+two_alpha_ijkl_plus_one_S[2]    +(Sijkl[2]-Primitive_coords[3][2]),(double)nsz_ijkl[3]);
 }
 return res[0]*res[1]*res[2];
}

//Jab is the upper bound function defined by A8 in Cioslowski's the paper.
double Jab(double &expa, double &expb,int &lxa,int &lya,int &lza,int &lxb,
int &lyb,int &lzb,double &Xa,double &Ya,double &Za,double &Xb,double &Yb,double &Zb)
{
 double res=ZERO,zeta_ab,e_ab,Rab[3],Xmax=ZERO,Ymax=ZERO,Zmax=ZERO,expa_sqrt_two_zab,expb_sqrt_two_zab,abs_Rb_Ra[3];
 abs_Rb_Ra[0]=abs(Xb-Xa);
 abs_Rb_Ra[1]=abs(Yb-Ya);
 abs_Rb_Ra[2]=abs(Zb-Za);
 zeta_ab=expa+expb;
 expa_sqrt_two_zab=expa*pow(TWO/zeta_ab,HALF);
 expb_sqrt_two_zab=expb*pow(TWO/zeta_ab,HALF);
 Xmax=Rmax[lxa][lxb];
 Ymax=Rmax[lya][lyb];
 Zmax=Rmax[lza][lzb];
 e_ab=expa*expb/zeta_ab;
 Rab[0]=Xa-Xb;
 Rab[1]=Ya-Yb;
 Rab[2]=Za-Zb;
 res=pow(PI,HALF*THREE)*pow(TWO*zeta_ab,-(double)(lxa+lxb+lya+lyb+lza+lzb)-HALF*THREE)
    *exp(-TWO*e_ab*scalarP(Rab))
    *pow(Xmax+expb_sqrt_two_zab*abs_Rb_Ra[0],TWO*(double)lxa)
    *pow(Xmax+expa_sqrt_two_zab*abs_Rb_Ra[0],TWO*(double)lxb)
    *pow(Ymax+expb_sqrt_two_zab*abs_Rb_Ra[1],TWO*(double)lya)
    *pow(Ymax+expa_sqrt_two_zab*abs_Rb_Ra[1],TWO*(double)lyb)
    *pow(Zmax+expb_sqrt_two_zab*abs_Rb_Ra[2],TWO*(double)lza)
    *pow(Zmax+expa_sqrt_two_zab*abs_Rb_Ra[2],TWO*(double)lzb);
 return res;
}
//Sab is the overlap between two primitive components
double Sab(int &order, double *s, double *w,double &S_ab, double &Sa,double &Sb,double &exp_a,double &exp_b, int &nxa,int &nxb)
{
 int i;
 double res=ZERO,z_ab;
 z_ab=exp_a+exp_b;
 for(i=0;i<order;i++)
 {
  res=res+w[i]*pow(pow(z_ab,-HALF)*s[i]+S_ab-Sa,(double)nxa)*pow(pow(z_ab,-HALF)*s[i]+S_ab-Sb,(double)nxb);
 }
 return res;
}
//Mab is the moment between two primitive components
double Mab(int &order, double *s, double *w,double &S_ab, double &Sa,double &Sb,double &exp_a,double &exp_b, int &nxa,int &nxb)
{
 int i;
 double res=ZERO,z_ab;
 z_ab=exp_a+exp_b;
 for(i=0;i<order;i++)
 {
  res=res+w[i]*pow(pow(z_ab,-HALF)*s[i]+S_ab-Sa,(double)nxa)*pow(pow(z_ab,-HALF)*s[i]+S_ab-Sb,(double)nxb)*(pow(z_ab,-HALF)*s[i]+S_ab);
 }
 return res;
}
//M2ab is the moment between two primitive components
double M2ab(int &order, double *s, double *w,double &S_ab, double &Sa,double &Sb,double &exp_a,double &exp_b, int &nxa,int &nxb)
{
 int i;
 double res=ZERO,z_ab;
 z_ab=exp_a+exp_b;
 for(i=0;i<order;i++)
 {
  res=res+w[i]*pow(pow(z_ab,-HALF)*s[i]+S_ab-Sa,(double)nxa)*pow(pow(z_ab,-HALF)*s[i]+S_ab-Sb,(double)nxb)*pow(pow(z_ab,-HALF)*s[i]+S_ab,TWO);
 }
 return res;
}
//Check for available order in Lebedev spherical quadrature
void grid_avail(int &Order)
{
  if(Order<=6)
  {
   Order =    6;
  }
  else if(Order<=14)
  {
   Order =   14;
  }
  else if(Order<=26)
  {
   Order =   26;
  }
  else if(Order<=38)
  {
   Order =   38;
  }
  else if(Order<=50)
  {
   Order =   50;
  }
  else if(Order<=74)
  {
   Order =   74;
  }
  else if(Order<=86)
  {
   Order =   86;
  }
  else if(Order<=110)
  {
   Order =  110;
  }
  else if(Order<=146)
  {
   Order =  146;
  }
  else if(Order<=170)
  {
   Order =  170;
  }
  else if(Order<=194)
  {
   Order =  194;
  }
  else if(Order<=230)
  {
   Order =  230;
  }
  else if(Order<=266)
  {
   Order =  266;
  }
  else if(Order<=302)
  {
   Order =  302;
  }
  else if(Order<=350)
  {
   Order =  350;
  }
  else if(Order<=434)
  {
   Order =  434;
  }
  else if(Order<=590)
  {
   Order =  590;
  }
  else if(Order<=770)
  {
   Order =  770;
  }
  else if(Order<=974)
  {
   Order =  974;
  }
  else if(Order<=1202)
  {
   Order = 1202;
  }
  else if(Order<=1454)
  {
   Order = 1454;
  }
  else if(Order<=1730)
  {
   Order = 1730;
  }
  else if(Order<=2030)
  {
   Order = 2030;
  }
  else if(Order<=2354)
  {
   Order = 2354;
  }
  else if(Order<=2702)
  {
   Order = 2702;
  }
  else if(Order<=3074)
  {
   Order = 3074;
  }
  else if(Order<=3470)
  {
   Order = 3470;
  }
  else if(Order<=3890)
  {
   Order = 3890;
  }
  else if(Order<=4334)
  {
   Order = 4334;
  }
  else if(Order<=4802)
  {
   Order = 4802;
  }
  else if(Order<=5294)
  {
   Order = 5294;
  }
  else
  {
   Order = 5810;
  }
if(ID==0)
{
  cout<<"Angular grid available :"<<setw(18)<<Order<<endl;
}
}

// User defined angular grid using Legendre Polinomials
void user_defined_angular(int &nang,int &nang2,double *x,double *y,double *z,double *w_theta_phi,string name_basis)
{
 int i,j,k;
 double *temp_w1,*temp_w2,*temp_theta,*temp_phi,zero=ZERO,one=ONE;
 temp_w1=new double[nang];
 temp_theta=new double[nang];
 temp_w2=new double[nang2];
 temp_phi=new double[nang2];
if(ID==0)
{
 cout<<"Angular grid set to    :"<<setw(18)<<nang*nang2<<endl;
 //////////////
 //Theta Grid//
 //////////////
 legendre_quadrature(name_basis.substr(0,name_basis.length()-6),nang,zero,one);
 //Read quadrature info
 ifstream read_ang_quad;
 // Read weights
 read_ang_quad.open((name_basis.substr(0,name_basis.length()-6)+"_w.txt").c_str());
 for(i=0;i<nang;i++)
 {
  read_ang_quad>>temp_w1[i];
 }
 read_ang_quad.close();
 // Read roots
 read_ang_quad.open((name_basis.substr(0,name_basis.length()-6)+"_x.txt").c_str());
 for(i=0;i<nang;i++)
 {
  read_ang_quad>>temp_theta[i];
  temp_theta[i]=acos(ONE-TWO*temp_theta[i]);
 }
 read_ang_quad.close();
 // Remove quadrature files
 system(("rm "+name_basis.substr(0,name_basis.length()-6)+"_r.txt").c_str());
 system(("rm "+name_basis.substr(0,name_basis.length()-6)+"_w.txt").c_str());
 system(("rm "+name_basis.substr(0,name_basis.length()-6)+"_x.txt").c_str());
 ////////////
 //Phi Grid//
 ////////////
 legendre_quadrature(name_basis.substr(0,name_basis.length()-6),nang2,zero,one);
 //Read quadrature info
 // Read weights
 read_ang_quad.open((name_basis.substr(0,name_basis.length()-6)+"_w.txt").c_str());
 for(i=0;i<nang2;i++)
 {
  read_ang_quad>>temp_w2[i];
 }
 read_ang_quad.close();
 // Read roots
 read_ang_quad.open((name_basis.substr(0,name_basis.length()-6)+"_x.txt").c_str());
 for(i=0;i<nang2;i++)
 {
  read_ang_quad>>temp_phi[i];
  temp_phi[i]=TWO*PI*temp_phi[i];
 }
 read_ang_quad.close();
 // Remove quadrature files
 system(("rm "+name_basis.substr(0,name_basis.length()-6)+"_r.txt").c_str());
 system(("rm "+name_basis.substr(0,name_basis.length()-6)+"_w.txt").c_str());
 system(("rm "+name_basis.substr(0,name_basis.length()-6)+"_x.txt").c_str());
}
#ifdef HAVE_MPI
MPI_Bcast(temp_w1,nang,MPI_DOUBLE,0,MPI_COMM_WORLD);
MPI_Bcast(temp_theta,nang,MPI_DOUBLE,0,MPI_COMM_WORLD);
MPI_Bcast(temp_w2,nang2,MPI_DOUBLE,0,MPI_COMM_WORLD);
MPI_Bcast(temp_phi,nang2,MPI_DOUBLE,0,MPI_COMM_WORLD);
#endif
 //////////////
 //Super Grid//
 //////////////
 k=0;
 for(i=0;i<nang;i++)
 {
  for(j=0;j<nang2;j++)
  {
   x[k]=sin(temp_theta[i])*cos(temp_phi[j]); 
   y[k]=sin(temp_theta[i])*sin(temp_phi[j]);
   z[k]=cos(temp_theta[i]);
   w_theta_phi[k]=temp_w1[i]*temp_w2[j];
   k++;
  }
 }
 //Delete temporal matrices 
 delete[] temp_w1;temp_w1=NULL;
 delete[] temp_w2;temp_w2=NULL;
 delete[] temp_theta;temp_theta=NULL;
 delete[] temp_phi;temp_phi=NULL;
}

void initialize(double Rmax[5][5])
{
 // To understand this part remember that for Jab the polynomial to integrate is of order
 // 2(la+lb) then the exact quadrature for Jab needs N_roots=[2(la+lb)+1]2.
 // for instance for la=0 and lb=0 N_roots=1/2 => we take ONE root which is ZERO
 // la=1 and lb = 0 => 0.7071067811865475
 // etc...
 //int order=2(nx,ny);
 //double    alpha=ZERO;
 //double    a=ZERO;
 //double    b=ONE;
 //gauss_hermite_rule("temp.txt",alpha,a,b,order);
 Rmax[0][0]=ZERO;              Rmax[0][1]=0.7071067811865475;Rmax[0][2]=1.2247448713915890;Rmax[0][3]=1.6506801238857840;Rmax[0][4]=2.0201828704560860;
 Rmax[1][0]=0.7071067811865475;Rmax[1][1]=1.2247448713915890;Rmax[1][2]=1.6506801238857840;Rmax[1][3]=2.0201828704560860;Rmax[1][4]=2.3506049736744930;
 Rmax[2][0]=1.2247448713915890;Rmax[2][1]=1.6506801238857840;Rmax[2][2]=2.0201828704560860;Rmax[2][3]=2.3506049736744930;Rmax[2][4]=2.6519613568352340;
 Rmax[3][0]=1.6506801238857840;Rmax[3][1]=2.0201828704560860;Rmax[3][2]=2.3506049736744930;Rmax[3][3]=2.6519613568352340;Rmax[3][4]=2.9306374202572430;
 Rmax[4][0]=2.0201828704560860;Rmax[4][1]=2.3506049736744930;Rmax[4][2]=2.6519613568352340;Rmax[4][3]=2.9306374202572430;Rmax[4][4]=3.1909932017815270;
}

//Calculate the time that the program took
void calc_time(double DATE[2][4])
{
 long double timeI,timeF,timeTOTAL;
 timeI=DATE[0][3]+DATE[0][2]*SIXTY+DATE[0][1]*pow(SIXTY,TWO)
      +DATE[0][0]*pow(SIXTY,TWO)*THREE*EIGHT;
 timeF=DATE[1][3]+DATE[1][2]*SIXTY+DATE[1][1]*pow(SIXTY,TWO)
      +DATE[1][0]*pow(SIXTY,TWO)*THREE*EIGHT;
 timeTOTAL=timeF-timeI;
 DATE[0][0]=(int)(timeTOTAL/(pow(SIXTY,TWO)*THREE*EIGHT)) ;
 timeTOTAL=timeTOTAL-DATE[0][0]*pow(SIXTY,TWO)*THREE*EIGHT;
 DATE[0][1]=(int)(timeTOTAL/(pow(SIXTY,TWO)));
 timeTOTAL=timeTOTAL-DATE[0][1]*pow(SIXTY,TWO);
 DATE[0][2]=(int)(timeTOTAL/SIXTY);
 timeTOTAL=timeTOTAL-DATE[0][2]*SIXTY;
 DATE[0][3]=timeTOTAL;
}


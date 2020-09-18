#include"Mathematical_Functions.h"
#define SQR(x)      ((x)*(x))                        // x^2
//////////////////////////
//Functions description //
//////////////////////////
//Euclidean norm for 3D
double scalarP(double vec[3])
{
return pow(vec[0],TWO)+pow(vec[1],TWO)+pow(vec[2],TWO);
}
double norm3D(double vec[3])
{
return sqrt(pow(vec[0],TWO)+pow(vec[1],TWO)+pow(vec[2],TWO));
}
//Double factorial (2k-1)!!
int dfact(int n)
{

 int i,ndfact;
 ndfact=1;
 for(i=n;i>0;i=i-2)
 {
  ndfact=ndfact*i;
 }
 return ndfact;
}
//Factorial
int factorial(int n)
{
 if(n>1)
 {n=n*factorial(n-1);}
 else{n=1;}
 return n;
}
//Permutes the indexes (Orders indexes first from lower to higher)
void permutation(int *indexes,int order,int &permutation_term)
{
 int i=0;
 sort(indexes,indexes+order);
 while(i<permutation_term)
 {
  next_permutation(indexes,indexes+order);
  i++;
 }
}
//Pochhammer Symbol (q)n rising
double PochhammerS(double n,double x)
{
 return Gamma(x+n)/Gamma(x);
}
//Second order polynom (ax^2+bx+c=0)
void second_order_polynom(double abc[3],double &root1,double &root2)
{
  double a=abc[0];
  double b=abc[1];
  double c=abc[2];
 if(pow(b,2)>4*a*c)
 {
 root1=(-b+sqrt(pow(b,2)-4*a*c))/2*a;
 root2=(-b-sqrt(pow(b,2)-4*a*c))/2*a;
 }
}
//Sign function
//The kind of the return value is that of A and B. If B\ge 0 then the result is ABS(A), else it is -ABS(A).
double sign(double a,double b)
{
 if(b>=0)
 {a=abs(a);}
 else
 {a=-abs(a);}
 return a;
}
//Print matrix
void print_mat(int &n,double **A)
{
 int i,j;
 for(i=0;i<n;i++)
 {
  for(j=0;j<n;j++)
  {
   cout<<setprecision(4)<<fixed<<A[i][j]<<"\t";
  }cout<<endl;
 }
}
//Matrix diagonalization and calculation of eigenvectors
void jacobi(int n, double **m, double **v)
{
 int i,j,k,ip,iq,maxiter=10000;
 double tol = pow(TEN,-(TWO*SIX)),apq,t,alpha,c,s,tau,d,delta,temp1,temp2;
 //Define v(matrix eigenvectors) as identity matrix
 for(i=0;i<n;i++)
 {
  for(j=0;j<n;j++)
  {v[i][j]=ZERO;}
   v[i][i] = ONE;
 }
 //Calculate Delta
 delta=ZERO;
 for(i=0;i<n;i++)
 {
  for(j=i+1;j<n;j++)
  {delta=delta+m[i][j]*m[i][j];}
 }
 //Iterations cycle
 for(i=0;i<maxiter;i++)
 {
  apq=ZERO;
  ip=0;
  iq=1;
  for(j=0;j<n;j++)
  {
   for(k=j+1;k<n;k++)
   {
    if(abs(m[j][k])>apq)
    {
     ip = j;
     iq = k;
     apq = abs(m[ip][iq]);
    }
   }
  }
  //Determine c, s, c, tau
  t=ONE;
  if(m[ip][ip]!=m[iq][iq])
  {
   alpha=(m[iq][iq]-m[ip][ip])/(TWO*m[ip][iq]);
   t=-alpha+alpha/abs(alpha)*sqrt(alpha*alpha+ONE);
  }
  c=ONE/sqrt(t*t+ONE);
  s=t*c;
  tau=s/(ONE+c);
  d=m[ip][iq];
  //Update matrix m en and the upper diagonal part
  m[ip][ip]=m[ip][ip]-t*m[ip][iq];
  m[iq][iq]=m[iq][iq]+t*m[ip][iq];
  for(j=0;j<ip;j++)
  {
   temp1=m[j][ip];
   temp2=m[j][iq];
   m[j][ip]=temp1-s*(temp2+tau*temp1);
   m[j][iq]=temp2+s*(temp1-tau*temp2);
  }
  for(j=ip+1;j<iq;j++)
  {
   temp1=m[ip][j];
   temp2=m[j][iq];
   m[j][iq]=temp2+s*(temp1-tau*temp2);
   m[ip][j]=temp1-s*(temp2+tau*temp1);
  }
  for(j=iq+1;j<n;j++)
  {
   temp1=m[ip][j];
   temp2=m[iq][j];
   m[iq][j]=temp2+s*(temp1-tau*temp2);
   m[ip][j]=temp1-s*(temp2+tau*temp1);
  }
  m[ip][iq]=ZERO;
  //Update v
  for(j=0;j<n;j++)
  {
   temp1=v[j][ip];
   temp2=v[j][iq];
   v[j][ip]=c*temp1-s*temp2;
   v[j][iq]=s*temp1+c*temp2;
  }
  //Update delta
  delta=delta-d*d;
  //If it has converge, update the lower diagonal part of m and finish jacobi()
  if(abs(delta)<=tol)
  {
   for(j=0;j<n;j++)
   {
    for(k=j+1;k<n;k++)
    {m[k][j] = m[j][k];}
   }
   return;
  }
 }
}

double dot(int &n,double *a,double *b)
{
 int i;
 double s=ZERO;
 for(i=0;i<n;i++)
 {s=s+a[i]*b[i];}
 return s;
}

void proyect(int &n,double *a,double *b,double *res)
{
 int i;
 for(i=0;i<n;i++)
 {
  res[i]=(dot(n,a,b)/dot(n,b,b))*b[i];
 }
}

void matmul(int n,double **A,double **B, double **RES)
{
 int i,j,k;
 double *row,*column;
 row=new double[n];
 column=new double[n];
 for(i=0;i<n;i++)
 {
  row[i]=ZERO;
  column[i]=ZERO;
 }
 for(i=0;i<n;i++)
 {
  for(j=0;j<n;j++)
  {
   for(k=0;k<n;k++)
   {
    row[k]=A[i][k];
    column[k]=B[k][j];
   }
   RES[i][j]=dot(n,row,column);
  }
 }
 delete[] row;
 delete[] column;
}
//A-A^T==RES
void mat_transpose(int n,double **A, double **RES)
{
 int i,j;
 for(i=0;i<n;i++)
 {
  for(j=0;j<n;j++)
  {
   RES[j][i]=A[i][j];
  }
 }
}
//RES==A
void mat_equal(int n,double **A, double **RES)
{
 int i,j;
 for(i=0;i<n;i++)
 {
  for(j=0;j<n;j++)
  {
   RES[i][j]=A[i][j];
  }
 }
}
//Matrix Inversion by Gauss-Jordan elimination
void mat_inverse(int n,double **A,double **Ainv)
{
 int i,j,k;
 bool noninvrt=false;
 double **temp,pivot;
 temp=new double*[n];
 for(i=0;i<n;i++)
 {
  temp[i]=new double[2*n];
 }
 for(i=0;i<n;i++)
 {
  for(j=0;j<2*n;j++)
  {
   if(j<n)
   {
    temp[i][j]=A[i][j];
   }
   else
   {
    temp[i][j]=ZERO;
    if(j==(i+n))
    {
     temp[i][j]=ONE;
    }
   }
  }
 }
 for(i=0;i<n;i++)
 {
  pivot=temp[i][i];
  k=i;
  for(j=i;j<n;j++)
  {
   if(temp[j][i]>pivot)
   {
    k=j;
    pivot=temp[j][i];
   }
  }
  if(pivot!=ZERO)
  {
   for(j=0;j<2*n;j++)
   {
    pivot=temp[k][j];
    temp[k][j]=temp[i][j];
    temp[i][j]=pivot;
   }
   pivot=temp[i][i];
   for(j=0;j<2*n;j++)
   {
    temp[i][j]=temp[i][j]/pivot;
   }
   for(j=i+1;j<n;j++)
   {
    pivot=-temp[j][i];
    for(k=i;k<2*n;k++)
    {
     temp[j][k]=temp[j][k]+pivot*temp[i][k];
    }
   }
  }
  else
  {
   cout<<"Warning! Not an invertible matrix!"<<endl;
   noninvrt=true;
   break;
  }
 }
 if(!noninvrt)
 {
  for(i=0;i<n;i++)
  {
   for(j=0;j<i;j++)
   {
    pivot=-temp[n-1-i][n-1-j];
    for(k=0;k<2*n;k++)
    {
     temp[n-1-i][k]=temp[n-1-i][k]+pivot*temp[n-1-j][k];
    }
   }
  }
 }
 for(i=0;i<n;i++)
 {
  for(j=n;j<2*n;j++)
  {
   Ainv[i][j-n]=temp[i][j];
  }
 }
 for(i=0;i<n;i++)
 {delete[] temp[i];}
 delete[] temp;
}
//Spherical harmonics
double Y_lm(double l,double m,double teta, double phi)
{
 double eval;
 eval=pow((TWO*l+ONE)*(double)(factorial((int)l-(int)abs(m)))/(FOUR*PI*(double)(factorial((int)l+(int)abs(m)))),HALF);
 eval=eval*P_lm(l,abs(m),cos(teta));
 if(m>ZERO)
 {
  eval=pow(TWO,HALF)*eval*cos(m*phi);
 }
 else if(m==ZERO)
 {
  eval=eval;
 }
 else
 {
  eval=pow(TWO,HALF)*eval*sin(abs(m)*phi);
 }
 return eval;
}
//Asociated legendre polynomials
double P_lm(double l, double m, double x)
{
 double eval;
 if(l==ZERO && m==ZERO)
 {eval=ONE;}
 else if(l==ONE && m==ZERO)
 {eval=x;}
 else if(l>ONE && m==ZERO)
 {eval=((TWO*l-ONE)/l)*x*P_lm(l-ONE,m,x)-(l-ONE)*P_lm(l-TWO,m,x)/l;}
 else if(m>ZERO)
 {
  if(x==ONE || x==-ONE)
  {
   if(x>ZERO)
   {x=x-pow(TEN,-TEN);}
   else
   {x=x+pow(TEN,-TEN);}
  }
  eval=((l-(m-ONE))*x*P_lm(l,m-ONE,x)-(l+m-ONE)*P_lm(l-ONE,m-ONE,x))/pow(ONE-x*x,HALF);
 }
 else if(m<ZERO)
 {eval=pow(-ONE,m)*((double)(factorial((int)l-(int)m))/(double)(factorial((int)l+(int)m)))*P_lm(l,m,x);}
 else
 {eval=ZERO;}
 return eval;
}
//Hermite polynomials
double H_n(double n,double x)
{
 double eval;
 if(n==ZERO)
 {eval=ONE;}
 else if(n==ONE)
 {eval=TWO*x;}
 else
 {eval=TWO*x*H_n(n-ONE,x)-TWO*(n-ONE)*H_n(n-TWO,x);}
 return eval;
}
//Generalized Laguerre polynomials
double L_n_alpha(double n,double alpha,double x)
{
 double eval;
 if(n==ZERO)
 {eval=ONE;}
 else
 {eval=((n+alpha)*L_n_alpha(n-ONE,alpha,x)-x*L_n_alpha(n-ONE,alpha+ONE,x))/n;}
 return eval;
}
//Chebyshev of first kind polynomials
double T_n(double n,double x)
{
 double eval;
 if(n==ZERO)
 {eval=ONE;}
 else if(n==ONE)
 {eval=x;}
 else
 {eval=TWO*x*T_n(n-ONE,x)-T_n(n-TWO,x);}
 return eval;
}
//Chebyshev of second kind polynomials
double U_n(double n,double x)
{
 double eval;
 if(n==ZERO)
 {eval=ONE;}
 else if(n==ONE)
 {eval=TWO*x;}
 else
 {eval=TWO*x*U_n(n-ONE,x)-U_n(n-TWO,x);}
 return eval;
}
//Gamma function
double Gamma(double x)
{
 int i;
 double eval,t;
 if(x>=ONE)
 {
// Function integrate 0 to 1 [pow(t/(ONE-t),x-ONE)*exp(-(t/(ONE-t)))*(ONE/pow(ONE-t,TWO))] dt
  eval=ZERO;
  for(i=1;i<10000000;i++)
  {
   t=((double)i)/pow(TEN,SEVEN);
   eval=eval+pow(t/(ONE-t),x-ONE)*exp(-(t/(ONE-t)))*(ONE/pow(ONE-t,TWO));
  }
  eval=eval+pow(ONE/pow(ONE,-EIGHT),x-ONE)*exp(-(ONE/pow(ONE,-EIGHT)))*(ONE/pow(pow(ONE,-EIGHT),TWO))/TWO;
  eval=(ONE/pow(TEN,SEVEN))*eval;
 }
 else
 {eval=Gamma(x+ONE)/x;}
 return eval;
}
//Bessel first kind Jv(x) function
double J_v(double v, double x)
{
 int i;
 double eval,eval2,t;
// Function integrate 0 to pi [ cos(nt-x*sin(t)) ] dt
  eval=ZERO;
  for(i=1;i<1000000;i++)
  {
   t=((double)i)*PI/pow(TEN,SIX);
   eval=eval+cos(x*sin(t)-v*t);
  }
  eval=eval+(cos(-v*PI)/TWO)+HALF;
  eval=(eval/pow(TEN,SIX));
// Function: sin(v*pi)/pi integrate 0 to inf [ exp(-x*sinh(t)-vt) ] dt
 if((v-(int)v)!=ZERO)
 {
  eval2=ZERO;
  for(i=1;i<1000000;i++)
  {
   t=((double)i)/pow(TEN,SIX);
   eval2=eval2+exp(-x*sinh(t/(ONE-t))-v*t/(ONE-t))/(pow(ONE-t,TWO));
  }
  t=ONE-pow(TEN,-EIGHT);
  eval2=eval2+(HALF+exp(-x*sinh(t/(ONE-t))-v*t/(ONE-t))/(TWO*pow(ONE-t,TWO)));
  eval2=-sin(PI*v)*eval2/(PI*pow(TEN,SIX));
  eval=eval+eval2;
 }
 return eval;
}
//Bessel second kind Yv(x) function
double Y_v(double v, double x)
{
 double eval,diff;
 diff=(int)v;
 if((v-diff)==ZERO){v=v+pow(TEN,-EIGHT);}
 eval=(J_v(v,x)*cos(v*PI)-J_v(-v,x))/sin(v*PI);
 return eval;
}
//Spherical Bessel jn(x)
double j_n(double n,double x)
{
 return pow(PI/(TWO*x),HALF)*J_v(n+HALF,x);
}
//Spherical Bessel yn(x)
double y_n(double n,double x)
{
 return pow(-ONE,n+ONE)*pow(PI/(TWO*x),HALF)*J_v(-n-HALF,x);
}
/////////////////////////////////////////////////////////////////////////
// Analytic Solution for 3D Quantum Harmonic Oscillator                //
/////////////////////////////////////////////////////////////////////////
double Q_HARM_OSC_3D(int nx, int ny, int nz, double omega, double xyz[3])
{
 double eval,phi_nx,phi_ny,phi_nz;
 phi_nx=pow(pow(TWO,nx)*factorial(nx),-HALF)*pow(omega/PI,ONE/FOUR);
 phi_nx=phi_nx*H_n(nx,xyz[0]*pow(omega,HALF))*exp(-HALF*omega*xyz[0]*xyz[0]);
 phi_ny=pow(pow(TWO,ny)*factorial(ny),-HALF)*pow(omega/PI,ONE/FOUR);
 phi_ny=phi_ny*H_n(ny,xyz[1]*pow(omega,HALF))*exp(-HALF*omega*xyz[1]*xyz[1]);
 phi_nz=pow(pow(TWO,nz)*factorial(nz),-HALF)*pow(omega/PI,ONE/FOUR);
 phi_nz=phi_nz*H_n(nz,xyz[2]*pow(omega,HALF))*exp(-HALF*omega*xyz[2]*xyz[2]);
 eval=phi_nx*phi_ny*phi_nz;
 return eval;
}
//Check if is NAN for log
void nan(string in,bool &checked)
{
 if(in[0]=='a')
 {checked=false;}
 else if(in[0]=='b')
 {checked=false;}
 else if(in[0]=='c')
 {checked=false;}
 else if(in[0]=='d')
 {checked=false;}
 else if(in[0]=='e')
 {checked=false;}
 else if(in[0]=='f')
 {checked=false;}
 else if(in[0]=='g')
 {checked=false;}
 else if(in[0]=='h')
 {checked=false;}
 else if(in[0]=='i')
 {checked=false;}
 else if(in[0]=='j')
 {checked=false;}
 else if(in[0]=='k')
 {checked=false;}
 else if(in[0]=='l')
 {checked=false;}
 else if(in[0]=='m')
 {checked=false;}
 else if(in[0]=='n')
 {checked=false;}
 else if(in[0]=='o')
 {checked=false;}
 else if(in[0]=='p')
 {checked=false;}
 else if(in[0]=='q')
 {checked=false;}
 else if(in[0]=='r')
 {checked=false;}
 else if(in[0]=='s')
 {checked=false;}
 else if(in[0]=='t')
 {checked=false;}
 else if(in[0]=='u')
 {checked=false;}
 else if(in[0]=='w')
 {checked=false;}
 else if(in[0]=='x')
 {checked=false;}
 else if(in[0]=='y')
 {checked=false;}
 else if(in[0]=='z')
 {checked=false;}
else{checked=true;}
}

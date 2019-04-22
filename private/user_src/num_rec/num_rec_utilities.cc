#include "num_rec_utilities.h"


namespace oomph
{



//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////


//=================================================================
/// Evaluate the fitting function and its derivatives
/// w.r.t. fitting parameters (done by FD by default; can be
/// overloaded)
//=================================================================
double LevenbergMarquardtFittingFunctionObject::fitting_function(
 const double& x,
 Vector<double>& dfit_dparam)
{
 // Get reference value
 double fct=fitting_function(x);
 
 // FD step
 double eps_fd=1.0e-8;
 
 // Do FD loop
 unsigned n_param=Parameter.size();
 for (unsigned i=0;i<n_param;i++)
  {
   double backup=Parameter[i];
   Parameter[i]+=eps_fd;
   double fct_pls=fitting_function(x);
   dfit_dparam[i]=(fct_pls-fct)/eps_fd;
   Parameter[i]=backup;
  }
 return fct;
}




//=====================================================================
/// Fit the parameters to the pairs of (x,y) data specified, 
/// using max_iter Levenberg Marquardt iterations
//=====================================================================
void LevenbergMarquardtFitter::fit_it(
 const Vector<std::pair<double,double> >& fitting_data,
 const unsigned& max_iter, const bool& quiet)
{
 if (Fitting_function_object_pt==0)
  {
   throw OomphLibError("Fitting_function_object_pt==0",
                       "Problem::distribute()",
                       OOMPH_EXCEPTION_LOCATION);    
  }

 // Number of parameters:
 unsigned nparam=Fitting_function_object_pt->nparameter();

 // By default regard all parameters as fittable -- can generalise
 // this at some point.
 std::vector<bool> ia(nparam,true);
   
 // Chi squared
 double chisq=0.0;

 // Number of data pairs
 unsigned ndata=fitting_data.size();

 // Vector of standard deviations -- just set to one
 Vector<double> sig(ndata,1.0);

 // Move to vectors (as required by Num Rec interface
 Vector<double> x(ndata), y(ndata);
 for (unsigned i=0;i<ndata;i++)
  {
   x[i]=fitting_data[i].first;
   y[i]=fitting_data[i].second;
  }

 // "Workspace" for numerical recipes
 int ma = static_cast<int>(nparam); 
 DenseDoubleMatrix covar(ma,ma);
 DenseDoubleMatrix alpha(ma,ma);


 if (!quiet)
  {
   oomph_info << "Chi_squared" << " ";
   for (unsigned i=0;i<unsigned(ma);i++)
    { 
     oomph_info << " parameter " << i << " ";
    }       
   oomph_info << std::endl;
  }

 // Start iteration with almda negative for setup
 double alamda=-0.1;
 for (unsigned iter=0;iter<max_iter;iter++)
  {
   // This is where Num Rec code starts so now it gets really ugly
   static int mfit;
   static double ochisq;
   int j,k,l;
          
   static Vector<double> oneda(ma);
   static Vector<double> atry(ma), beta(ma), da(ma);

   // Initialisation
   if (alamda < 0.0)
    {
     mfit=0;
     for (j=0;j<ma;j++)
      {
       if (ia[j]) mfit++;
      }
     alamda=0.001;
     mrqcof(x,y,sig,Fitting_function_object_pt->parameter(),ia,
            alpha,beta,chisq);
     ochisq=chisq;
     for (j=0;j<ma;j++)
      {
       atry[j]=Fitting_function_object_pt->parameter(j); 
      }
    } 
     
   DenseDoubleMatrix temp(mfit,mfit);
   for (j=0;j<mfit;j++)
    {
     for (k=0;k<mfit;k++)
      {
       covar(j,k)=alpha(j,k);
      }
     covar(j,j)=alpha(j,j)*(1.0+alamda);
     for (k=0;k<mfit;k++)
      {
       temp(j,k)=covar(j,k);
      }
     oneda[j]=beta[j];
    }
     
   // Linear solver
   temp.solve(oneda); 
     
   for (j=0;j<mfit;j++)
    {
     for (k=0;k<mfit;k++)
      {
       covar(j,k)=temp(j,k);
      }
     da[j]=oneda[j];
    }
     
   // Converged
   if (alamda == 0.0)
    {
     return;
    }
     
          
   for (j=0,l=0;l<ma;l++)
    {
     if (ia[l]) atry[l]=Fitting_function_object_pt->parameter(l)+da[j++];
    }
   mrqcof(x,y,sig,atry,ia,covar,da,chisq);
   if (chisq < ochisq)
    {
     alamda *= 0.1;
     ochisq=chisq;
     for (j=0;j<mfit;j++)
      {
       for (k=0;k<mfit;k++) 
        {
         alpha(j,k)=covar(j,k);
        }
       beta[j]=da[j];
      }
     
     for (l=0;l<ma;l++)
      { 
       Fitting_function_object_pt->parameter(l)=atry[l];
      }
       
     if (!quiet)
      {
       //Store the current output flags
       std::ios_base::fmtflags ff = std::cout.flags();
       // Output with fixed width
       std::cout.setf(std::ios_base::scientific,std::ios_base::floatfield);
       std::cout.width(15);
       std::cout << chisq << " ";
       for (l=0;l<ma;l++)
        { 
         std::cout << atry[l] << " ";
        }       
       std::cout << std::endl;
       // Reset
       std::cout.setf(ff, std::ios_base::floatfield);
       std::cout.width(0);
      }
    }
   else
    {
     alamda *= 10.0;
     chisq=ochisq;
    }
     
  }

}

//==================================================================
/// Private helper function -- don't look into it...
//==================================================================
void LevenbergMarquardtFitter::mrqcof(Vector<double>& x, 
                                      Vector<double>& y, 
                                      Vector<double>& sig, 
                                      Vector<double>& a,
                                      std::vector<bool>& ia, 
                                      DenseDoubleMatrix& alpha, 
                                      Vector<double>& beta, 
                                      double& chisq)
{
 int i,j,k,l,m,mfit=0;
 double ymod,wt,sig2i,dy;
 
 int ndata=x.size();
 int ma=a.size();
 Vector<double> dyda(ma);
 for (j=0;j<ma;j++)
  {
   if (ia[j]) mfit++;
  }
 
 for (j=0;j<mfit;j++)
  {
   for (k=0;k<=j;k++)
    {
     alpha(j,k)=0.0;
    }
   beta[j]=0.0;
  }
 
 chisq=0.0;
 for (i=0;i<ndata;i++) 
  {
   Vector<double> backup=Fitting_function_object_pt->parameter();
   Fitting_function_object_pt->parameter()=a;
   ymod=Fitting_function_object_pt->fitting_function(x[i],dyda);
   Fitting_function_object_pt->parameter()=backup;
   sig2i=1.0/(sig[i]*sig[i]);
   dy=y[i]-ymod;
   for (j=0,l=0;l<ma;l++)
    {
     if (ia[l])
      {
       wt=dyda[l]*sig2i;
       for (k=0,m=0;m<l+1;m++)
        {
         if (ia[m]) alpha(j,k++) += wt*dyda[m];
        }
       beta[j++] += dy*wt;
      }
    }
   chisq += dy*dy*sig2i;
  }
 
 
 for (j=1;j<mfit;j++)
  {
   for (k=0;k<j;k++)
    {
     alpha(k,j)=alpha(j,k);
    }
  }
}


//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////


//=========================================================================
/// Namespaces for (Numerical-Recipes-based) functions for
/// eigensolver, based on Jacobi rotations.
//=========================================================================
namespace JacobiEigenSolver
{

 /// Perform one Jacobi rotation on matrix a
 inline void rot(DenseDoubleMatrix&a, const double s, const double tau, 
                 const unsigned long i, const unsigned long j, 
                 const unsigned long k, const unsigned long l)
 {
  double g,h;
  
  g=a(i,j);
  h=a(k,l);
  a(i,j)=g-s*(h+g*tau);
  a(k,l)=h+s*(g-h*tau);

 }


/// \short Use Jacobi rotations to determine eigenvalues and eigenvectors of 
/// matrix a. d[i]=i-th eigenvalue; v(i,j)=i-th component of j-th eigenvector
/// (note that this is the transpose of what we'd like to have...);
/// nrot=number of rotations used. 
 void jacobi(DenseDoubleMatrix& a, Vector<double>& d, 
             DenseDoubleMatrix& v, unsigned long& nrot)
 {
#ifdef PARANOID
  // Check Matrix a is square
  if (a.ncol()!=a.nrow())
   {
    throw OomphLibError(
     "This matrix is not square, the matrix MUST be square!",
     "JacobiEigenSolver::jacobi()",
     OOMPH_EXCEPTION_LOCATION);
   }
#endif
 
  // If matrix v is wrong size, correct it!
  if (v.ncol()!=a.ncol() || v.nrow()!=a.nrow())
   {
    v.resize(a.nrow(),a.nrow(),0.0);
   }

  unsigned long i,j,ip,iq;
  double tresh,theta,tau,t,sm,s,h,g,c;
 
  unsigned long n=d.size();
  Vector<double> b(n);
  Vector<double> z(n);
  for (ip=0;ip<n;ip++) {
   for (iq=0;iq<n;iq++) v(ip,iq)=0.0;
   v(ip,ip)=1.0;
  }
  for (ip=0;ip<n;ip++) {
   b[ip]=d[ip]=a(ip,ip);
   z[ip]=0.0;
  }
  nrot=0;
  for (i=1;i<=50;i++) {
   sm=0.0;
   for (ip=0;ip<n-1;ip++) {
    for (iq=ip+1;iq<n;iq++)
     sm += std::fabs(a(ip,iq));
   }
   if (sm == 0.0)
    return;
   if (i < 4)
    tresh=0.2*sm/(n*n);
   else
    tresh=0.0;
   for (ip=0;ip<n-1;ip++) {
    for (iq=ip+1;iq<n;iq++) {
     g=100.0*std::fabs(a(ip,iq));
     if (i > 4 && (std::fabs(d[ip])+g) == std::fabs(d[ip])
         && (std::fabs(d[iq])+g) == std::fabs(d[iq]))
      a(ip,iq)=0.0;
     else if (std::fabs(a(ip,iq)) > tresh) {
      h=d[iq]-d[ip];
      if ((std::fabs(h)+g) == std::fabs(h))
       t=(a(ip,iq))/h;
      else {
       theta=0.5*h/(a(ip,iq));
       t=1.0/(std::fabs(theta)+std::sqrt(1.0+theta*theta));
       if (theta < 0.0) t = -t;
      }
      c=1.0/std::sqrt(1+t*t);
      s=t*c;
      tau=s/(1.0+c);
      h=t*a(ip,iq);
      z[ip] -= h;
      z[iq] += h;
      d[ip] -= h;
      d[iq] += h;
      a(ip,iq)=0.0;
      for (j=0;j<ip;j++)
       rot(a,s,tau,j,ip,j,iq);
      for (j=ip+1;j<iq;j++)
       rot(a,s,tau,ip,j,j,iq);
      for (j=iq+1;j<n;j++)
       rot(a,s,tau,ip,j,iq,j);
      for (j=0;j<n;j++)
       rot(v,s,tau,j,ip,j,iq);
      ++nrot;
     }
    }
   }
   for (ip=0;ip<n;ip++) {
    b[ip] += z[ip];
    d[ip]=b[ip];
    z[ip]=0.0;
   }
  }
  throw OomphLibError(
   "Too many iterations in routine jacobi",
   "JacobiEigenSolver::jacobi()",
   OOMPH_EXCEPTION_LOCATION);
 }

}

/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////



}


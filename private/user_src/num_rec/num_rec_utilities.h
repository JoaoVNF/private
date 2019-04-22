//Header file for num rec inspired code

#ifndef OOMPH_NR_HEADER
#define OOMPH_NR_HEADER

									       
// Config header generated by autoconfig
#ifdef HAVE_CONFIG_H
 #include <oomph-lib-config.h>
#endif


//OOMPH-LIB headers
#include "../generic/nodes.h"
#include "../generic/Qelements.h"
#include "../generic/oomph_utilities.h"



namespace oomph
{



//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////


//============================================================
/// Spline class -- padded to remain constant outside range
//============================================================
class Spline
{
 
public:
 
 /// \short Constructor: Pass x and y coordinates of function to be
 /// fitted (natural spline -- extension to prescribed slope is easy...)
 Spline(const Vector<double> &x_sample, const Vector<double> &y_sample,
        bool zero_slope_at_ends=true) :
  X_sample(x_sample), Y_sample(y_sample)
  {

   unsigned n=x_sample.size();

   // If the first entry in the X sample vector is larger than
   // the last one we need to sort
   if(X_sample[0] > X_sample[n-1])
    {
     // Vector of pairs
     Vector<std::pair<double, double> > vector_of_pairs(n);
     std::pair<double, double> tmp;
     
     // Set up the pairs
     for(unsigned i=0;i<n;i++)
      {
       tmp = std::make_pair(X_sample[i], Y_sample[i]);
       vector_of_pairs[i]=tmp;
      }
     
     // Sort the vector of pairs using standard sort algorithm
     // (based on first entry in pair (x-value))
     sort(vector_of_pairs.begin(), vector_of_pairs.end());
     
     // Update the vectors accordingly
     for(unsigned i=0;i<n;i++)
      {
       X_sample[i] = vector_of_pairs[i].first;
       Y_sample[i] = vector_of_pairs[i].second;
      }

    }

   Y_at_x_max=Y_sample[0];
   Y_at_x_min=Y_sample[0];
   X_min=X_sample[0];
   X_max=X_sample[0];

   int k=0;
   double p=0.0;
   double qn=0.0;
   double sig=0.0;
   double un=0.0;;
   
   Spline_coeff.resize(n);
   Vector<double> u(n-1);
   
   // Natural spline
   if (!zero_slope_at_ends)
    Spline_coeff[0]=u[0]=0.0;
   else 
    {
     /// set slope at end to be zero
     double yp1=0.0;
     Spline_coeff[0] = -0.5;
     u[0]=(3.0/(X_sample[1]-X_sample[0]))*
      ((Y_sample[1]-Y_sample[0])/(X_sample[1]-X_sample[0])-yp1);
    }
   for (unsigned i=1;i<n-1;i++) 
    {

     if (X_sample[i]<X_min)
      {
       X_min=X_sample[i];
       Y_at_x_min=Y_sample[i];
      }
     if (X_sample[i]>X_max) 
      {
       X_max=X_sample[i];
       Y_at_x_max=Y_sample[i];
      }

     if ((X_sample[i+1]-X_sample[i])==0.0)
      {
       std::cout<<X_sample[i]<<" "<<Y_sample[i]<<std::endl;
       std::cout<<X_sample[i+1]<<" "<<Y_sample[i+1]<<std::endl;
       std::ostringstream error_stream;
       error_stream<< "X_sample[i+1] = X_sample[i] for i = " 
                   << i << std::endl;
       throw OomphLibError(
        error_stream.str(),
        OOMPH_CURRENT_FUNCTION,
        OOMPH_EXCEPTION_LOCATION);
      }
     if ((X_sample[i+1]-X_sample[i-1])==0.0)
      {
       std::ostringstream error_stream;
       error_stream<< "X_sample[i+1] = X_sample[i-1] for i = " 
                   << i << std::endl;
       throw OomphLibError(
        error_stream.str(),
        OOMPH_CURRENT_FUNCTION,
        OOMPH_EXCEPTION_LOCATION);
      }
     sig=(X_sample[i]-X_sample[i-1])/(X_sample[i+1]-X_sample[i-1]);
     p=sig*Spline_coeff[i-1]+2.0;
     Spline_coeff[i]=(sig-1.0)/p;
     u[i]=(Y_sample[i+1]-Y_sample[i])/(X_sample[i+1]-X_sample[i]) - 
      (Y_sample[i]-Y_sample[i-1])/(X_sample[i]-X_sample[i-1]);
     u[i]=(6.0*u[i]/(X_sample[i+1]-X_sample[i-1])-sig*u[i-1])/p;
    }
   // Natural spline
   if (!zero_slope_at_ends)
    qn=un=0.0;
   else 
    {
     double ypn=0.0;
     qn=0.5;
     un=(3.0/(X_sample[n-1]-X_sample[n-2]))*
      (ypn-(Y_sample[n-1]-Y_sample[n-2])/(X_sample[n-1]-X_sample[n-2]));
    }

   {
    unsigned i=n-1;
    if (X_sample[i]<X_min)
     {
      X_min=X_sample[i];
      Y_at_x_min=Y_sample[i];
     }
    if (X_sample[i]>X_max) 
     {
      X_max=X_sample[i];
      Y_at_x_max=Y_sample[i];
     }
   }

   Spline_coeff[n-1]=(un-qn*u[n-2])/(qn*Spline_coeff[n-2]+1.0);
   for (k=n-2;k>=0;k--)
    {
     Spline_coeff[k]=Spline_coeff[k]*Spline_coeff[k+1]+u[k];
    }
  }
 
 
 /// Spline evalution -- padded to remain constant outside range
 double spline(const double& x)
  {

   if (x>X_max) return Y_at_x_max;
   if (x<X_min) return Y_at_x_min;

   double y=0.0;
   int k=0;
   double h=0.0;
   double b=0.0;
   double a=0.0;
   
   int n=X_sample.size();
   int klo=0;
   int khi=n-1;
   while (khi-klo > 1)
    {
     k=(khi+klo) >> 1;
     if (X_sample[k] > x) 
      {
       khi=k;
      }
     else
      {
       klo=k;
      }
    }
   h=X_sample[khi]-X_sample[klo];
   // never get here... now caught in setup (constructor)
   //if (h == 0.0)
   // {
   //  oomph_info << "Bad X_sample input to routine splint\n";
   //  abort();
   // }
   a=(X_sample[khi]-x)/h;
   b=(x-X_sample[klo])/h;
   y=a*Y_sample[klo]+b*Y_sample[khi]+((a*a*a-a)*Spline_coeff[klo]
                                      +(b*b*b-b)*Spline_coeff[khi])*(h*h)/6.0;
   
   return y;
  }
 
 
private:
 
 /// x samples
 Vector<double> X_sample;
 
 /// y samples
 Vector<double> Y_sample;
 
 /// Spline coefficients
 Vector<double> Spline_coeff;
 
 /// Value for padding to the right of data
 double Y_at_x_max;

 /// Value for padding to the left of data
 double Y_at_x_min;

 /// Min sample value
 double X_min;

 /// Max sample value
 double X_max;

};


//============================================================
/// Parametric spline class
//============================================================
class ParametricSpline
{
 
public:
 
 /// \short Constructor: Pass x and y coordinates of function to be
 /// fitted (natural spline -- extension to prescribed slope is easy...)
 ParametricSpline(const Vector<Vector<double> > &vertices)
  {
   /// Number of vertices
   unsigned n_sample = vertices.size();

   /// Initialise arc length to zero
   double zeta_total = 0.0;

   X_sample.resize(n_sample);
   Y_sample.resize(n_sample);
   Zeta_sample.resize(n_sample);

   for(unsigned i=0; i<n_sample; i++)
    {
     X_sample[i] = vertices[i][0];
     Y_sample[i] = vertices[i][1];

     if(i>0)
      {
       double x0=vertices[i-1][0];
       double y0=vertices[i-1][1];
       double x1=vertices[i][0];
       double y1=vertices[i][1];
       
       zeta_total += sqrt((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0));
      }
     Zeta_sample[i] = zeta_total;

    }

   Zeta_max = zeta_total;

   /// Initialise the XZetaSpline and YZetaSpline
   X_zeta_spline = new Spline(Zeta_sample, X_sample);
   Y_zeta_spline = new Spline(Zeta_sample, Y_sample);
  }
 
 
 /// Spline evalution -- padded to remain constant outside range
 Vector<double> position(const double& zeta)
  {
   double x = X_zeta_spline->spline(zeta*Zeta_max);
   double y = Y_zeta_spline->spline(zeta*Zeta_max);

   Vector<double> pos(2, 0.0);

   pos[0]=x;
   pos[1]=y;

   return pos;
   
  }
 
 
private:
 
 /// x samples
 Vector<double> X_sample;
 
 /// y samples
 Vector<double> Y_sample;
 
 /// Zeta samples
 Vector<double> Zeta_sample;

 /// Spline representation for x as a fct of zeta
 Spline* X_zeta_spline;

 /// Spline representation for y as a fct of zeta
 Spline* Y_zeta_spline;

 /// Maximum arc length
 double Zeta_max;

};


//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////



//=====================================================================
/// Base class for functions whose parameters can be fitted
/// by Levenberg Marquardt technique.
//=====================================================================
class LevenbergMarquardtFittingFunctionObject
{
 
public:


 /// Constructor: Specify number of fitting parameters
 LevenbergMarquardtFittingFunctionObject(const unsigned& n_param)
  {
   Parameter.resize(n_param);
   for (unsigned i=0;i<n_param;i++)
    {
     Parameter[i]=0.0;
    }
  }

 /// Empty destructor
 virtual ~LevenbergMarquardtFittingFunctionObject() {}

 /// \short Evaluate the fitting function for the current set
 /// of parameters: Pure virtual, must be implemented.
 virtual double fitting_function(const double& x)=0;

 /// \short Evaluate the fitting function and its derivatives
 /// w.r.t. fitting parameters (done by FD by default; can be
 /// overloaded)
 virtual double fitting_function(const double& x,
                                 Vector<double>& dfit_dparam);
                   
 /// \short Number of parameters in fitting function. Pure virtual, must be
 /// implemented.
 virtual unsigned nparameter()=0;

 /// Access to i-th fitting parameter
 double& parameter(const unsigned& i)
  {
   return Parameter[i];
  }

 /// Access to vector of fitting parameters
 Vector<double>& parameter()
  {
   return Parameter;
  }


protected:

 /// Vector of fitting parameters
 Vector<double> Parameter;

};





//=====================================================================
/// Damped oscillatory function whose parameters can be
/// fitted with Levenberg Marquardt.
//=====================================================================
class DampedOscillatoryFittingFunctionObject : 
 virtual public LevenbergMarquardtFittingFunctionObject
{
 
public:
 
 /// Constructor: Number of fitting parameters is five. 
 DampedOscillatoryFittingFunctionObject() :  
  LevenbergMarquardtFittingFunctionObject(5)
  {}


 /// \short Evaluate the fitting function for the current set
 /// of parameters
 double fitting_function(const double& x)
  {
   return Parameter[0]+
    exp(Parameter[1]*x)*Parameter[2]*sin(Parameter[3]*x+Parameter[4]);
  }
                                 
 /// \short Overload all interfaces of the fitting function, call the default
 /// finite difference version
 double fitting_function(const double &x,
                         Vector<double> &dfit_dparam)
  {
   return 
    LevenbergMarquardtFittingFunctionObject::fitting_function(x,dfit_dparam);
  }

 /// Number of parameters in fitting function
 virtual unsigned nparameter()
  {
   return 5;
  }


};


//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////


//=====================================================================
/// Class that allows fitting of free parameters in function
/// (represented by a LevenbergMarquardtFittingFunctionObject)
/// to given (x,y) data.
//=====================================================================
class LevenbergMarquardtFitter
{
 
public:

 /// Empty constructor 
 LevenbergMarquardtFitter(): Fitting_function_object_pt(0)
  {}
 
 /// \short Access to pointer to LevenbergMarquardtFittingFunctionObject
 /// whose parameters we want to determine by fit to data.
 LevenbergMarquardtFittingFunctionObject*& fitting_function_object_pt()
  {
   return Fitting_function_object_pt;
  }

 /// \short Fit the parameters to the pairs of (x,y) data specified, 
 /// using max_iter Levenberg Marquardt iterations
 void fit_it(const Vector<std::pair<double,double> >& fitting_data,
             const unsigned& max_iter,
             const bool& quiet=true);

private:

 /// Pointer to LevenbergMarquardtFittingFunctionObject
 LevenbergMarquardtFittingFunctionObject* Fitting_function_object_pt;


 /// Private helper function -- don't look into it...
 void mrqcof(Vector<double>& x, 
             Vector<double>& y, 
             Vector<double>& sig, 
             Vector<double>& a,
             std::vector<bool>& ia, 
             DenseDoubleMatrix& alpha, 
             Vector<double>& beta, 
             double& chisq);

 
};




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
 extern inline void rot(DenseDoubleMatrix&a, const double s, const double tau, 
                 const unsigned long i, const unsigned long j, 
                        const unsigned long k, const unsigned long l);

/// \short Use Jacobi rotations to determine eigenvalues and eigenvectors of 
/// matrix a. d[i]=i-th eigenvalue; v(i,j)=i-th component of j-th eigenvector
/// (note that this is the transpose of what we'd like to have...);
/// nrot=number of rotations used. 
 extern void jacobi(DenseDoubleMatrix& a, Vector<double>& d, 
             DenseDoubleMatrix& v, unsigned long& nrot);

}

/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////



}


#endif
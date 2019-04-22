//Header file for Spectral Helmholtz elements
#ifndef OOMPH_SPECTRAL_HELMHOLTZ_ELEMENTS_HEADER
#define OOMPH_SPECTRAL_HELMHOLTZ_ELEMENTS_HEADER

// Config header generated by autoconfig
#ifdef HAVE_CONFIG_H
  #include <oomph-lib-config.h>
#endif

//OOMPH-LIB headers
#include "helmholtz_elements.h"
#include "../generic/Qspectral_elements.h"

namespace oomph
{

//======================================================================
/// QSpectralHelmholtzElement elements are linear/quadrilateral/brick-shaped 
/// Helmholtz elements with isoparametric spectral interpolation for the 
/// function. Note that the implementation is HelmholtzEquations<DIM> does
/// not use sum factorisation for the evaluation of the residuals and is,
/// therefore, not optimal for higher dimensions.
//======================================================================
template <unsigned DIM, unsigned NNODE_1D>
 class QSpectralHelmholtzElement : public virtual QSpectralElement<DIM,NNODE_1D>,
 public virtual HelmholtzEquations<DIM>
{

  private:
 
 /// \short Static array of ints to hold number of variables at 
 /// nodes: Initial_Nvalue[n]
 static const unsigned Initial_Nvalue;
 
  public:


 ///\short  Constructor: Call constructors for QSpectralElement and 
 /// Helmholtz equations
 QSpectralHelmholtzElement() : QSpectralElement<DIM,NNODE_1D>(), 
  HelmholtzEquations<DIM>()
  {}
 
 /// Broken copy constructor
 QSpectralHelmholtzElement(const QSpectralHelmholtzElement<DIM,NNODE_1D>& dummy) 
  { 
   BrokenCopy::broken_copy("QSpectralHelmholtzElement");
  } 
 
 /// Broken assignment operator
 void operator=(const QSpectralHelmholtzElement<DIM,NNODE_1D>&) 
  {
   BrokenCopy::broken_assign("QSpectralHelmholtzElement");
  }

/*  /// \short Overload assign_all_generic_local_eqn_numbers since it is */
/*  ///        already overloaded in both HelmholtzEquations AND */
/*  ///        QSpectralElement */
/*  virtual inline void assign_all_generic_local_eqn_numbers() */
/*   { */
/*    //Assign unique external data */
/*    assign_unique_external_data_helper(); */
/*    //Call the QSpectralElement assign */
/*    QSpectralElement<DIM,NNODE_1D>::assign_all_generic_local_eqn_numbers(); */
/*   } */

 /// \short  Required  # of `values' (pinned or dofs) 
 /// at node n
 inline unsigned required_nvalue(const unsigned &n) const 
  {return Initial_Nvalue;}

 /// \short Output function:  
 ///  x,y,u   or    x,y,z,u
 void output(std::ostream &outfile)
  {HelmholtzEquations<DIM>::output(outfile);}

 ///  \short Output function:  
 ///   x,y,u   or    x,y,z,u at n_plot^DIM plot points
 void output(std::ostream &outfile, const unsigned &n_plot)
  {HelmholtzEquations<DIM>::output(outfile,n_plot);}


 /// \short C-style output function:  
 ///  x,y,u   or    x,y,z,u
 void output(FILE* file_pt)
  {HelmholtzEquations<DIM>::output(file_pt);}


 ///  \short C-style output function:  
 ///   x,y,u   or    x,y,z,u at n_plot^DIM plot points
 void output(FILE* file_pt, const unsigned &n_plot)
  {HelmholtzEquations<DIM>::output(file_pt,n_plot);}


 /// \short Output function for an exact solution:
 ///  x,y,u_exact   or    x,y,z,u_exact at n_plot^DIM plot points
 void output_fct(std::ostream &outfile, const unsigned &n_plot,
                 FiniteElement::SteadyExactSolutionFctPt exact_soln_pt)
  {HelmholtzEquations<DIM>::output_fct(outfile,n_plot,exact_soln_pt);}



 /// \short Output function for a time-dependent exact solution.
 ///  x,y,u_exact   or    x,y,z,u_exact at n_plot^DIM plot points
 /// (Calls the steady version)
 void output_fct(std::ostream &outfile, const unsigned &n_plot,
                 const double& time,
                 FiniteElement::UnsteadyExactSolutionFctPt exact_soln_pt)
  {HelmholtzEquations<DIM>::output_fct(outfile,n_plot,time,exact_soln_pt);}


protected:

/// Shape, test functions & derivs. w.r.t. to global coords. Return Jacobian.
 inline double dshape_and_dtest_eulerian_helmholtz(
  const Vector<double> &s, Shape &psi, DShape &dpsidx, 
  Shape &test, DShape &dtestdx) const;
 

 /// \short Shape, test functions & derivs. w.r.t. to global coords. at
 /// integration point ipt. Return Jacobian.
 inline double dshape_and_dtest_eulerian_at_knot_helmholtz(
  const unsigned& ipt, Shape &psi, DShape &dpsidx, 
  Shape &test, DShape &dtestdx) 
  const;

};


//Inline functions:


//======================================================================
/// Define the shape functions and test functions and derivatives
/// w.r.t. global coordinates and return Jacobian of mapping.
///
/// Galerkin: Test functions = shape functions
//======================================================================
template<unsigned DIM, unsigned NNODE_1D>
double QSpectralHelmholtzElement<DIM,NNODE_1D>::
 dshape_and_dtest_eulerian_helmholtz(
  const Vector<double> &s,
  Shape &psi, 
  DShape &dpsidx,
  Shape &test, 
  DShape &dtestdx) const
{
 //Call the geometrical shape functions and derivatives  
 double J = this->dshape_eulerian(s,psi,dpsidx);

 //Loop over the test functions and derivatives and set them equal to the
 //shape functions
 for(unsigned i=0;i<NNODE_1D;i++)
  {
   test[i] = psi[i]; 
   for(unsigned j=0;j<DIM;j++)
    {
     dtestdx(i,j) = dpsidx(i,j);
    }
  }
 
 //Return the jacobian
 return J;
}

//======================================================================
/// Define the shape functions and test functions and derivatives
/// w.r.t. global coordinates and return Jacobian of mapping.
///
/// Galerkin: Test functions = shape functions
//======================================================================
template<unsigned DIM, unsigned NNODE_1D>
double QSpectralHelmholtzElement<DIM,NNODE_1D>::
 dshape_and_dtest_eulerian_at_knot_helmholtz(
  const unsigned &ipt,
  Shape &psi, 
  DShape &dpsidx,
  Shape &test, 
  DShape &dtestdx) const
{
 //Call the geometrical shape functions and derivatives  
 double J = this->dshape_eulerian_at_knot(ipt,psi,dpsidx);

 //Set the pointers of the test functions
 test = psi;
 dtestdx = dpsidx;

 //Return the jacobian
 return J;
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////



//=======================================================================
/// Face geometry for the QSpectralHelmholtzElement elements: The spatial 
/// dimension of the face elements is one lower than that of the
/// bulk element but they have the same number of points
/// along their 1D edges.
//=======================================================================
template<unsigned DIM, unsigned NNODE_1D>
class FaceGeometry<QSpectralHelmholtzElement<DIM,NNODE_1D> >: 
 public virtual QSpectralElement<DIM-1,NNODE_1D>
{

  public:
 
 /// \short Constructor: Call the constructor for the
 /// appropriate lower-dimensional QElement
 FaceGeometry() : QSpectralElement<DIM-1,NNODE_1D>() {}

};


//=======================================================================
/// Face geometry for the 1D QHelmholtzElement elements: Point elements
//=======================================================================
template<unsigned NNODE_1D>
class FaceGeometry<QSpectralHelmholtzElement<1,NNODE_1D> >: 
 public virtual PointElement
{

  public:
 
 /// \short Constructor: Call the constructor for the
 /// appropriate lower-dimensional QElement
 FaceGeometry() : PointElement() {}

};

}

#endif

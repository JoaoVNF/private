//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC//           Version 0.90. August 3, 2009.
//LIC// 
//LIC// Copyright (C) 2006-2009 Matthias Heil and Andrew Hazel
//LIC// 
//LIC// This library is free software; you can redistribute it and/or
//LIC// modify it under the terms of the GNU Lesser General Public
//LIC// License as published by the Free Software Foundation; either
//LIC// version 2.1 of the License, or (at your option) any later version.
//LIC// 
//LIC// This library is distributed in the hope that it will be useful,
//LIC// but WITHOUT ANY WARRANTY; without even the implied warranty of
//LIC// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//LIC// Lesser General Public License for more details.
//LIC// 
//LIC// You should have received a copy of the GNU Lesser General Public
//LIC// License along with this library; if not, write to the Free Software
//LIC// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
//LIC// 02110-1301  USA.
//LIC// 
//LIC// The authors may be contacted at oomph-lib@maths.man.ac.uk.
//LIC// 
//LIC//====================================================================
//Header file for UnsteadyHeat elements
#ifndef OOMPH_UNSTEADY_HEAT_ELEMENTS_HEADER
#define OOMPH_UNSTEADY_HEAT_ELEMENTS_HEADER
				       				       
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

//=============================================================
/// A class for all isoparametric elements that solve the 
/// UnsteadyHeat equations.
/// \f[ 
/// \frac{\partial^2 u}{\partial x_i^2}=\frac{\partial u}{\partial t}+f(t,x_j)
/// \f] 
/// This contains the generic maths. Shape functions, geometric
/// mapping etc. must get implemented in derived class.
/// Note that this class assumes an isoparametric formulation, i.e. that
/// the scalar unknown is interpolated using the same shape funcitons
/// as the position.
//=============================================================
template <unsigned DIM>
class DummyUnsteadyHeatEquations : public virtual FiniteElement
{

public:

 /// \short Function pointer to coefficient function fct(t,x,f(x,t)) -- 
 /// x is a Vector! 
 typedef void (*UnsteadyHeatCoefficientFctPt)(const double& time,
                                         const Vector<double>& x,
                                         double& u);


 /// \short Constructor: Initialises the Coefficient_fct_pt to null and 
 /// sets flag to use ALE formulation of the equations.
 DummyUnsteadyHeatEquations() : Coefficient_fct_pt(0), ALE_is_disabled(false) {}
 

 /// Broken copy constructor
 DummyUnsteadyHeatEquations(const DummyUnsteadyHeatEquations& dummy) 
  { 
   BrokenCopy::broken_copy("DummyUnsteadyHeatEquations");
  } 
 
 /// Broken assignment operator
 void operator=(const DummyUnsteadyHeatEquations&) 
  {
   BrokenCopy::broken_assign("DummyUnsteadyHeatEquations");
  }

  /// \short Return the index at which the unknown value
 /// is stored. The default value, 0, is appropriate for single-physics
 /// problems, when there is only one variable, the value that satisfies the
 /// unsteady heat equation. 
 /// In derived multi-physics elements, this function should be overloaded
 /// to reflect the chosen storage scheme. Note that these equations require
 /// that the unknown is always stored at the same index at each node.
 virtual inline unsigned u_index_ust_heat() const {return 0;}
 
 /// \short du/dt at local node n. 
 /// Uses suitably interpolated value for hanging nodes.
 double du_dt_ust_heat(const unsigned &n) const
  {
   // Get the data's timestepper
   TimeStepper* time_stepper_pt= this->node_pt(n)->time_stepper_pt();

   //Initialise dudt
   double dudt=0.0;
   
   //Loop over the timesteps, if there is a non Steady timestepper
   if (!time_stepper_pt->is_steady())
     {
     //Find the index at which the variable is stored
     const unsigned u_nodal_index = u_index_ust_heat();
     
     // Number of timsteps (past & present)
     const unsigned n_time = time_stepper_pt->ntstorage();
     
     //Add the contributions to the time derivative
     for(unsigned t=0;t<n_time;t++)
      {
       dudt += time_stepper_pt->weight(1,t)*nodal_value(t,n,u_nodal_index);
      }
    }
   return dudt;
  }

 /// \short Disable ALE, i.e. assert the mesh is not moving -- you do this
 /// at your own risk!
 void disable_ALE()
  {
   ALE_is_disabled=true;
  }


 /// \short (Re-)enable ALE, i.e. take possible mesh motion into account
 /// when evaluating the time-derivative. Note: By default, ALE is 
 /// enabled, at the expense of possibly creating unnecessary work 
 /// in problems where the mesh is, in fact, stationary. 
 void enable_ALE()
  {
   ALE_is_disabled=false;
  }


 /// Output with default number of plot points
 void output(std::ostream &outfile) 
  {
   unsigned nplot=5;
   output(outfile,nplot);
  }


 /// \short Output FE representation of soln: x,y,u or x,y,z,u at 
 /// n_plot^DIM plot points
 void output(std::ostream &outfile, const unsigned &nplot);

 /// C_style output with default number of plot points
 void output(FILE* file_pt)
  {
   unsigned n_plot=5;
   output(file_pt,n_plot);
  }


 /// \short C-style output FE representation of soln: x,y,u or x,y,z,u at 
 /// n_plot^DIM plot points
 void output(FILE* file_pt, const unsigned &n_plot);


 /// Output exact soln: x,y,u_exact or x,y,z,u_exact at nplot^DIM plot points
 void output_fct(std::ostream &outfile, const unsigned &nplot, 
                 FiniteElement::SteadyExactSolutionFctPt 
                 exact_soln_pt);


 /// \short Output exact soln: x,y,u_exact or x,y,z,u_exact at 
 /// nplot^DIM plot points (time-dependent version)
 virtual 
  void output_fct(std::ostream &outfile, const unsigned &nplot,
                  const double& time, 
                  FiniteElement::UnsteadyExactSolutionFctPt 
                  exact_soln_pt);


 /// Get error against and norm of exact solution
 void compute_error(std::ostream &outfile, 
                    FiniteElement::SteadyExactSolutionFctPt 
                    exact_soln_pt,
                    double& error, double& norm);


 /// Get error against and norm of exact solution
 void compute_error(std::ostream &outfile, 
                    FiniteElement::UnsteadyExactSolutionFctPt 
                    exact_soln_pt,
                    const double& time, double& error, double& norm);


 /// Access function: Pointer to coefficient function
 UnsteadyHeatCoefficientFctPt& coefficient_fct_pt() {return Coefficient_fct_pt;}


 /// Access function: Pointer to coefficient function. Const version
 UnsteadyHeatCoefficientFctPt coefficient_fct_pt() const {return Coefficient_fct_pt;}


 /// \short Get coefficient term at continous time t and (Eulerian) position x.
 /// Virtual so it can be overloaded in derived multiphysics elements. 
 virtual inline void get_coefficient_ust_heat(const double& t,
                                         const unsigned& ipt,
                                         const Vector<double>& x,
                                         double& coefficient) const
  {
   //If no coefficient function has been set, return zero
   if(Coefficient_fct_pt==0) {coefficient = 0.0;}
   else
    {
     // Get coefficient strength
     (*Coefficient_fct_pt)(t,x,coefficient);
    }
  }

 /// Get flux: flux[i] = du/dx_i
 void get_flux(const Vector<double>& s, Vector<double>& flux) const
  {
   //Find out how many nodes there are in the element
   unsigned n_node = nnode();

   //Find the index at which the variable is stored
   unsigned u_nodal_index = u_index_ust_heat();

   //Set up memory for the shape and test functions
   Shape psi(n_node);
   DShape dpsidx(n_node,DIM);
 
   //Call the derivatives of the shape and test functions
   dshape_eulerian(s,psi,dpsidx);
     
   //Initialise to zero
   for(unsigned j=0;j<DIM;j++) {flux[j] = 0.0;}
   
   // Loop over nodes
   for(unsigned l=0;l<n_node;l++) 
    {
     //Loop over derivative directions
     for(unsigned j=0;j<DIM;j++)
      {                               
       flux[j] += nodal_value(l,u_nodal_index)*dpsidx(l,j);
      }
    }
  }


 /// Compute element residual Vector (wrapper)
 void fill_in_contribution_to_residuals(Vector<double> &residuals)
  {
   //Call the generic residuals function with flag set to 0
   //using a dummy matrix argument
   fill_in_generic_residual_contribution_ust_heat(
    residuals,GeneralisedElement::Dummy_matrix,0);
  }


 /// Compute element residual Vector and element Jacobian matrix (wrapper)
 void fill_in_contribution_to_jacobian(Vector<double> &residuals,
                                   DenseMatrix<double> &jacobian)
  {
   //Call the generic routine with the flag set to 1
   fill_in_generic_residual_contribution_ust_heat(residuals,jacobian,1);
  }
 

 /// Return FE representation of function value u(s) at local coordinate s
 inline double interpolated_u_ust_heat(const Vector<double> &s) const
  {
   //Find number of nodes
   unsigned n_node = nnode();

   //Find the index at which the variable is stored
   unsigned u_nodal_index = u_index_ust_heat();

   //Local shape function
   Shape psi(n_node);

   //Find values of shape function
   shape(s,psi);

   //Initialise value of u
   double interpolated_u = 0.0;

   //Loop over the local nodes and sum
   for(unsigned l=0;l<n_node;l++) 
    {
     interpolated_u += nodal_value(l,u_nodal_index)*psi[l];
    }

   return(interpolated_u);
  }

 /// \short Self-test: Return 0 for OK
 unsigned self_test();


protected:

 /// \short Shape/test functions and derivs w.r.t. to global coords at 
 /// local coord. s; return  Jacobian of mapping
 virtual double dshape_and_dtest_eulerian_ust_heat(const Vector<double> &s, 
                                                   Shape &psi, 
                                                   DShape &dpsidx, 
                                                   Shape &test, 
                                                   DShape &dtestdx) const=0;


 /// \short Shape/test functions and derivs w.r.t. to global coords at 
 /// integration point ipt; return  Jacobian of mapping
 virtual double dshape_and_dtest_eulerian_at_knot_ust_heat(
  const unsigned &ipt, 
  Shape &psi, 
  DShape &dpsidx,
  Shape &test, 
  DShape &dtestdx)
  const=0;

 /// \short Compute element residual Vector only (if flag=and/or element 
 /// Jacobian matrix 
 virtual void fill_in_generic_residual_contribution_ust_heat(
  Vector<double> &residuals, DenseMatrix<double> &jacobian, 
  unsigned flag); 

 /// Pointer to coefficient function:
 UnsteadyHeatCoefficientFctPt Coefficient_fct_pt;

 /// \short Boolean flag to indicate if ALE formulation is disabled when 
 /// time-derivatives are computed. Only set to true if you're sure
 /// that the mesh is stationary.
 bool ALE_is_disabled;

};






///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////



//======================================================================
/// DummyQUnsteadyHeatElement elements are linear/quadrilateral/brick-shaped 
/// UnsteadyHeat elements with isoparametric interpolation for the function.
//======================================================================
template <unsigned DIM, unsigned NNODE_1D>
 class DummyQUnsteadyHeatElement : public virtual QElement<DIM,NNODE_1D>,
 public virtual DummyUnsteadyHeatEquations<DIM>
{
  private:

 /// \short Static array of ints to hold number of variables at 
 /// nodes: Initial_Nvalue[n]
 static const unsigned Initial_Nvalue;
 
  public:
 
 ///\short  Constructor: Call constructors for QElement and 
 /// UnsteadyHeat equations
 DummyQUnsteadyHeatElement() : QElement<DIM,NNODE_1D>(), 
  DummyUnsteadyHeatEquations<DIM>()
  { }

 /// Broken copy constructor
 DummyQUnsteadyHeatElement(const DummyQUnsteadyHeatElement<DIM,NNODE_1D>& dummy) 
  { 
   BrokenCopy::broken_copy("DummyQUnsteadyHeatElement");
  } 
 
 /// Broken assignment operator
 void operator=(const DummyQUnsteadyHeatElement<DIM,NNODE_1D>&) 
  {
   BrokenCopy::broken_assign("DummyQUnsteadyHeatElement");
  }

 /// \short  Required  # of `values' (pinned or dofs) 
 /// at node n
 inline unsigned required_nvalue(const unsigned &n) const 
  {return Initial_Nvalue;}

 /// \short Output function:  
 ///  x,y,u   or    x,y,z,u
 void output(std::ostream &outfile)
  {DummyUnsteadyHeatEquations<DIM>::output(outfile);}


 ///  \short Output function:  
 ///   x,y,u   or    x,y,z,u at n_plot^DIM plot points
 void output(std::ostream &outfile, const unsigned &n_plot)
  {DummyUnsteadyHeatEquations<DIM>::output(outfile,n_plot);}


 /// \short C-style output function:  
 ///  x,y,u   or    x,y,z,u
 void output(FILE* file_pt)
  {DummyUnsteadyHeatEquations<DIM>::output(file_pt);}


 ///  \short C-style output function:  
 ///   x,y,u   or    x,y,z,u at n_plot^DIM plot points
 void output(FILE* file_pt, const unsigned &n_plot)
  {DummyUnsteadyHeatEquations<DIM>::output(file_pt,n_plot);}


 /// \short Output function for an exact solution:
 ///  x,y,u_exact   or    x,y,z,u_exact at n_plot^DIM plot points
 void output_fct(std::ostream &outfile, const unsigned &n_plot,
                 FiniteElement::SteadyExactSolutionFctPt 
                 exact_soln_pt)
  {DummyUnsteadyHeatEquations<DIM>::output_fct(outfile,n_plot,exact_soln_pt);}



 /// \short Output function for a time-dependent exact solution.
 ///  x,y,u_exact   or    x,y,z,u_exact at n_plot^DIM plot points
 /// (Calls the steady version)
 void output_fct(std::ostream &outfile, const unsigned &n_plot,
                 const double& time,
                 FiniteElement::UnsteadyExactSolutionFctPt exact_soln_pt)
  {DummyUnsteadyHeatEquations<DIM>::output_fct(outfile,n_plot,time,exact_soln_pt);}



protected:

 /// Shape, test functions & derivs. w.r.t. to global coords. Return Jacobian.
 inline double dshape_and_dtest_eulerian_ust_heat(const Vector<double> &s, 
                                                  Shape &psi, 
                                                  DShape &dpsidx, 
                                                  Shape &test, 
                                                  DShape &dtestdx) const;
 

 /// \short Shape/test functions and derivs w.r.t. to global coords at 
 /// integration point ipt; return  Jacobian of mapping
 inline double dshape_and_dtest_eulerian_at_knot_ust_heat(const unsigned &ipt, 
                                                          Shape &psi, 
                                                          DShape &dpsidx,
                                                          Shape &test, 
                                                          DShape &dtestdx)
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
double DummyQUnsteadyHeatElement<DIM,NNODE_1D>::
 dshape_and_dtest_eulerian_ust_heat(const Vector<double> &s,
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
template<unsigned DIM,unsigned NNODE_1D>
double DummyQUnsteadyHeatElement<DIM,NNODE_1D>::
 dshape_and_dtest_eulerian_at_knot_ust_heat(
 const unsigned &ipt,
 Shape &psi, 
 DShape &dpsidx,
 Shape &test, 
 DShape &dtestdx) const
{
 //Call the geometrical shape functions and derivatives  
 double J = this->dshape_eulerian_at_knot(ipt,psi,dpsidx);

 //Set the test functions equal to the shape functions 
 //(sets internal pointers)
 test = psi;
 dtestdx = dpsidx;

 //Return the jacobian
 return J;
}


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////



//=======================================================================
/// Face geometry for the DummyQUnsteadyHeatElement elements: The spatial 
/// dimension of the face elements is one lower than that of the
/// bulk element but they have the same number of points
/// along their 1D edges.
//=======================================================================
template<unsigned DIM, unsigned NNODE_1D>
class FaceGeometry<DummyQUnsteadyHeatElement<DIM,NNODE_1D> >: 
 public virtual QElement<DIM-1,NNODE_1D>
{

  public:
 
 /// \short Constructor: Call the constructor for the
 /// appropriate lower-dimensional QElement
 FaceGeometry() : QElement<DIM-1,NNODE_1D>() {}

};

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////


//=======================================================================
/// Face geometry for the 1D DummyQUnsteadyHeatElement elements: Point elements
//=======================================================================
template<unsigned NNODE_1D>
class FaceGeometry<DummyQUnsteadyHeatElement<1,NNODE_1D> >: 
 public virtual PointElement
{

  public:
 
 /// \short Constructor: Call the constructor for the
 /// appropriate lower-dimensional QElement
 FaceGeometry() : PointElement() {}

};


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

}

#endif

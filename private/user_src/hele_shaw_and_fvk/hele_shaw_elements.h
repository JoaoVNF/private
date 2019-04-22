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
//Header file for HeleShaw elements
#ifndef OOMPH_HELE_SHAW_ELEMENTS_HEADER
#define OOMPH_HELE_SHAW_ELEMENTS_HEADER


// Config header generated by autoconfig
#ifdef HAVE_CONFIG_H
#include <oomph-lib-config.h>
#endif


//OOMPH-LIB headers

// hierher uncomment these
//#include "../generic/nodes.h"
//#include "../generic/Qelements.h"
//#include "../generic/oomph_utilities.h"

#include "../../../src/generic/projection.h"
#include "../../../src/generic/nodes.h"
#include "../../../src/generic/Qelements.h"
#include "../../../src/generic/oomph_utilities.h"

namespace oomph
{
 
//=============================================================
/// A class for all isoparametric elements that solve the
/// HeleShaw equations.
/// \f[
/// dh_dt + div ( b^3 grad p)=0
/// \f]
/// This contains the generic maths. Shape functions, geometric
/// mapping etc. must get implemented in derived class.
//=============================================================
 class HeleShawEquations : public virtual FiniteElement
 {
  
   public:
  
  /// \short Function pointer to function which provides h(x,t) and dh/dt,
  /// for vector x.
  typedef void (*UpperWallFctPt)(const Vector<double>& x, double& h,
                                 double& dhdt);
  
  typedef void (*UpperWallFluxFctPt)(const Vector<double>&x,
                                     double&h,
                                     double& dhdt,
                                     Vector<double>& dhdx,
                                     Vector<double>& d_dhdt_dx );
  
 /* void hierher_shite() */
 /* { */
 /*  oomph_info << "hierher start shite: " << std::endl; */
 /*  // Loop over nodes */
 /*  unsigned n_node=this->nnode(); */
 /*  for(unsigned l=0; l<n_node; l++) */
 /*   { */
 /*    // Loop over directions */
 /*    for(unsigned j=0; j<2; j++) */
 /*     { */
 /*      oomph_info << "l, j " << l << " " << j << std::endl; */
 /*      double tmp=this->raw_nodal_position(l,j);  */
 /*      oomph_info << "tmp " << tmp << std::endl; */
 /*     } */
 /*   } */
 /*  oomph_info << "hierher end shite "<< std::endl; */
 /* }  */

  /// Constructor
  HeleShawEquations() : Upper_wall_fct_pt(0),
   Hele_shaw_disabled(false)
   {}
   
   /// Broken copy constructor
   HeleShawEquations(const HeleShawEquations& dummy)
    {
     BrokenCopy::broken_copy("HeleShawEquations");
    }
      
   /// Broken assignment operator
   void operator=(const HeleShawEquations&)
    {
     BrokenCopy::broken_assign("HeleShawEquations");
    }
   
   /// \short Return the index at which the unknown value
   /// is stored. The default value, 0, is appropriate for single-physics
   /// problems, when there is only one variable, the value that satisfies
   /// the hele_shaw equation.
   /// In derived multi-physics elements, this function should be overloaded
   /// to reflect the chosen storage scheme. Note that these equations require
   /// that the unknown is always stored at the same index at each node.
   virtual inline unsigned p_index_hele_shaw() const {return 0;}
   

   /// \short Switch off Hele Shaw equations (needed/useful in FSI version
   /// of this problem where we're not actually solving the HS equations
   /// in the regions occupied by the air bubble
   void disable_hele_shaw()
   {
    Hele_shaw_disabled=true;
   }


   /// Output with default number of plot points
   void output(std::ostream &outfile)
   {
    const unsigned n_plot = 3;
    output(outfile,n_plot);
   }
   
   /// \short Output FE representation of soln: x,y,u or x,y,z,u at
   /// n_plot^2 plot points
   void output(std::ostream &outfile, const unsigned &n_plot);
   
   /// C_style output with default number of plot points
   void output(FILE* file_pt)
   {
    const unsigned n_plot = 3;
    output(file_pt,n_plot);
   }
   
   /// \short C-style output FE representation of soln: x,y,u or x,y,z,u at
   /// n_plot^2 plot points
   void output(FILE* file_pt, const unsigned &n_plot);
   
   /// Output exact soln: x,y,u_exact or x,y,z,u_exact at n_plot^2 plot points
   void output_fct(std::ostream &outfile, const unsigned &n_plot,
                   FiniteElement::SteadyExactSolutionFctPt exact_soln_pt);
   
   /// \short Output exact soln: x,y,u_exact or x,y,z,u_exact at
   /// n_plot^2 plot points (dummy time-dependent version to
   /// keep intel compiler happy)
   virtual void output_fct(std::ostream &outfile, const unsigned &n_plot,
                           const double& time,
                           FiniteElement::UnsteadyExactSolutionFctPt
                           exact_soln_pt)
   {
    throw OomphLibError(
     "There is no time-dependent output_fct() for HeleShaw elements ",
     "HeleShawEquations::output_fct()",
     OOMPH_EXCEPTION_LOCATION);
   }
      
   /// Get error against and norm of exact solution
   void compute_error(std::ostream &outfile,
                      FiniteElement::SteadyExactSolutionFctPt exact_soln_pt,
                      double& error, double& norm);
   
   /// Dummy, time dependent error checker
   void compute_error(std::ostream &outfile,
                      FiniteElement::UnsteadyExactSolutionFctPt exact_soln_pt,
                      const double& time, double& error, double& norm)
   {
    throw OomphLibError(
     "There is no time-dependent compute_error() for HeleShaw elements",
     "HeleShawEquations::compute_error()",
     OOMPH_EXCEPTION_LOCATION);
   }

   /// Access function: Pointer to source function
   UpperWallFctPt& upper_wall_fct_pt() {return Upper_wall_fct_pt;}
   
   /// Access function: Pointer to source function. Const version
   UpperWallFctPt upper_wall_fct_pt() const {return Upper_wall_fct_pt;}
       
   /// \short Compute gap width b and its Eulerian (!) derivative db/dt.
   /// Shape function and deriv w.r.t. to global coords is passed in for 
   /// efficiency when it comes to overloading this function in an 
   /// FSI context. 
   inline virtual void get_upper_wall_data(const Vector<double>& s,
                                           const Vector<double>& x,
                                           const Shape& psi,
                                           const DShape& dpsidx,
                                           double& b,
                                           double& dbdt) const
    {
     // If no function has been set, assume constant thickness
     if (Upper_wall_fct_pt==0)
      {
       // Unit gapwidth
       b = 1.0;
       
       // Zero velocity
       dbdt = 0.0;
      }
     else
      {
       // Get data from function
       (*Upper_wall_fct_pt)(x,b,dbdt);
      }
    }
   
    
   /// Get pressure flux: gradient[i] = dp/dx_i
   /// This is useful to compute the velocity components, and can also be used
   /// as a flux vector for the Z2 error estimator(see eg Thele_shaw_elements).
   /// We could also use velocity as the flux vector for the Z2 error estimator
   void get_pressure_gradient(const Vector<double>& s,
                              Vector<double>& gradient) const
   {
    // Find out how many nodes there are in the element
    const unsigned n_node = nnode();
     
    // Get the index at which the unknown is stored
    const unsigned p_nodal_index = p_index_hele_shaw();;
     
    // Set up memory for the shape and test functions
    Shape psi(n_node);
    DShape dpsidx(n_node,2);
     
    // Call the derivatives of the shape and test functions
    dshape_eulerian(s,psi,dpsidx);
     
    // Initialise to zero
    for(unsigned j=0; j<2; j++)
     {
      gradient[j] = 0.0;
     }
     
    // Loop over nodes
    for(unsigned l=0; l<n_node; l++)
     {
      // Loop over derivative directions
      for(unsigned j=0; j<2; j++)
       {
        gradient[j] += this->nodal_value(l,p_nodal_index)*dpsidx(l,j);
       }
     }
   }
    
    
   /// The current nondimensionalisation has velocity[i] = -h^2 *dp/dx_i
   void get_velocity(const Vector<double>& s, Vector<double>& velocity) const
   {
     
    /// To find the velocity, we multiply the pressure gradient by b^2. We need
    /// to interpolate to find x(s) in order to call b(x,t) via 
    /// get_upper_wall_data.
     
    // Find out how many nodes there are in the element
    const unsigned n_node = nnode();
     
    // Set up memory for the shape and test functions
    Shape psi(n_node);
    DShape dpsidx(n_node,2);
     
    // Call the derivatives of the shape and test functions
    dshape_eulerian(s,psi,dpsidx);
     
    // Calculate local values of unknown
    // Allocate and initialise to zero
    Vector<double> interpolated_x(2,0.0);
    Vector<double> pressure_gradient(2,0.0);
     
    // Initialise
    double dbdt = 0.0;
    double b = 0.0;
     
    // Loop over nodes to assemble the coordinate
    for(unsigned l=0; l<n_node; l++)
     {
      // Loop over coordinate directions
      for(unsigned j=0; j<2; j++)
       {
        interpolated_x[j] += raw_nodal_position(l,j)*psi(l);
       }
     }
     
    // Now get the gap width; local coordinate of integration point
    // is only used in multiphysics
    get_upper_wall_data(s,interpolated_x,psi,dpsidx,b,dbdt);
    get_pressure_gradient(s,pressure_gradient);
     
    /// Now assemble the velocity components.
    for(unsigned j=0; j<2; j++)
     {
      velocity[j] = -b*b*pressure_gradient[j];
     }
     
   }    
    
   /// Add the element's contribution to its residual vector (wrapper)
   void fill_in_contribution_to_residuals(Vector<double> &residuals)
   {
    // Call the generic residuals function with flag set to 0
    // using a dummy matrix argument
    fill_in_generic_residual_contribution_hele_shaw(
     residuals,GeneralisedElement::Dummy_matrix,0);
   }    
    
   /// \short Add the element's contribution to its residual vector and
   /// element Jacobian matrix (wrapper)
   void fill_in_contribution_to_jacobian(Vector<double> &residuals,
                                         DenseMatrix<double> &jacobian)
   {
    // Call the generic routine with the flag set to 1
    fill_in_generic_residual_contribution_hele_shaw(residuals,jacobian,1);
   }
    
   /// \short Return FE representation of function value u_hele_shaw(s)
   /// at local coordinate s
   inline double interpolated_p_hele_shaw(const Vector<double> &s) const
    {
     // Find number of nodes
     const unsigned n_node = nnode();
      
     // Get the index at which the hele_shaw unknown is stored
     const unsigned p_nodal_index = p_index_hele_shaw();
      
     // Local shape function
     Shape psi(n_node);
      
     // Find values of shape function
     shape(s,psi);
      
     // Initialise value of u
     double interpolated_p = 0.0;
      
     // Loop over the local nodes and sum
     for(unsigned l=0; l<n_node; l++)
      {
       interpolated_p += this->nodal_value(l,p_nodal_index)*psi[l];
      }
      
     return(interpolated_p);
    }
    

   /// \short Self-test: Return 0 for OK
   unsigned self_test();
    
    
   protected:
    
   /// \short Shape/test functions and derivs w.r.t. to global coords at
   /// local coord. s; return  Jacobian of mapping
   virtual double dshape_and_dtest_eulerian_hele_shaw(
    const Vector<double> &s,
    Shape &psi,
    DShape &dpsidx, 
    Shape &test,
    DShape &dtestdx) const = 0;
    
   /// \short Shape/test functions and derivs w.r.t. to global coords at
   /// integration point ipt; return  Jacobian of mapping
   virtual double dshape_and_dtest_eulerian_at_knot_hele_shaw(
    const unsigned &ipt,
    Shape &psi,
    DShape &dpsidx,
    Shape &test,
    DShape &dtestdx) const = 0;
   
   /// \short Shape/test functions and derivs w.r.t. to global coords at
   /// integration point ipt; return Jacobian of mapping (J). Also compute
   /// derivatives of dpsidx, dtestdx and J w.r.t. nodal coordinates.
   virtual double dshape_and_dtest_eulerian_at_knot_hele_shaw(
    const unsigned &ipt,
    Shape &psi,
    DShape &dpsidx,
    RankFourTensor<double> &d_dpsidx_dX,
    Shape &test,
    DShape &dtestdx,
    RankFourTensor<double> &d_dtestdx_dX,
    DenseMatrix<double> &djacobian_dX) const = 0;
    
   /// \short Compute element residual Vector only (if flag=and/or element
   /// Jacobian matrix
   virtual void fill_in_generic_residual_contribution_hele_shaw(
    Vector<double> &residuals, DenseMatrix<double> &jacobian,
    const unsigned& flag);
    
   /// Pointer to function that specifies the gap width and wall velocity
   UpperWallFctPt Upper_wall_fct_pt;
    
   /// \short Boolean to  switch off Hele Shaw equations (needed/useful
   /// in FSI version
   /// of this problem where we're not actually solving the HS equations
   /// in the regions occupied by the air bubble. Defaults to false.
   bool Hele_shaw_disabled;

    
 };
 
 
 
 
 
 
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
 
 
 
//======================================================================
/// QHeleShawElement elements are linear/quadrilateral/brick-shaped
/// HeleShaw elements with isoparametric interpolation for the function.
//======================================================================
 template <unsigned NNODE_1D>
  class QHeleShawElement : public virtual QElement<2,NNODE_1D>,
  public virtual HeleShawEquations
  {
   
    private:
   
   /// \short Static int that holds the number of variables at
   /// nodes: always the same
   static const unsigned Initial_Nvalue;
   

    public:  
   
   /// \short  Constructor: Call constructors for QElement and
   /// HeleShaw equations
   QHeleShawElement() : QElement<2,NNODE_1D>(), HeleShawEquations()
    {}
    
    /// Broken copy constructor
    QHeleShawElement(const QHeleShawElement<NNODE_1D>& dummy)
     {
      BrokenCopy::broken_copy("QHeleShawElement");
     }
    
    /// Broken assignment operator
    void operator = (const QHeleShawElement<NNODE_1D>&)
     {
      BrokenCopy::broken_assign("QHeleShawElement");
     } 
   
    /// \short  Required  # of `values' (pinned or dofs)
    /// at node n
    inline unsigned required_nvalue(const unsigned &n) const
     {return Initial_Nvalue;}
    
    /// \short Output function:
    ///  x,y,u   or    x,y,z,u
    void output(std::ostream &outfile)
    {HeleShawEquations::output(outfile);}
    
    ///  \short Output function:
    ///   x,y,u   or    x,y,z,u at n_plot^2 plot points
    void output(std::ostream &outfile, const unsigned &n_plot)
    {HeleShawEquations::output(outfile,n_plot);}
        
    /// \short C-style output function:
    ///  x,y,u   or    x,y,z,u
    void output(FILE* file_pt)
    {HeleShawEquations::output(file_pt);}
        
    ///  \short C-style output function:
    ///   x,y,u   or    x,y,z,u at n_plot^2 plot points
    void output(FILE* file_pt, const unsigned &n_plot)
    {HeleShawEquations::output(file_pt,n_plot);}
        
    /// \short Output function for an exact solution:
    ///  x,y,u_exact   or    x,y,z,u_exact at n_plot^2 plot points
    void output_fct(std::ostream &outfile, const unsigned &n_plot,
                    FiniteElement::SteadyExactSolutionFctPt exact_soln_pt)
    {HeleShawEquations::output_fct(outfile,n_plot,exact_soln_pt);}
        
    /// \short Output function for a time-dependent exact solution.
    ///  x,y,u_exact   or    x,y,z,u_exact at n_plot^2 plot points
    /// (Calls the steady version)
    void output_fct(std::ostream &outfile, const unsigned &n_plot,
                    const double& time,
                    FiniteElement::UnsteadyExactSolutionFctPt exact_soln_pt)
    {HeleShawEquations::output_fct(outfile,n_plot,time,exact_soln_pt);}
    
        
    protected:
    
    /// \short Shape, test functions & derivs. w.r.t. to global coords. 
    /// Return Jacobian
    inline double dshape_and_dtest_eulerian_hele_shaw(
     const Vector<double> &s, Shape &psi, DShape &dpsidx,
     Shape &test, DShape &dtestdx) const;
    
    /// \short Shape, test functions & derivs. w.r.t. to global coords. at
    /// integration point ipt. Return Jacobian.
    inline double dshape_and_dtest_eulerian_at_knot_hele_shaw(
     const unsigned& ipt,
     Shape &psi,
     DShape &dpsidx,
     Shape &test,
     DShape &dtestdx) const;
        
    /// \short Shape/test functions and derivs w.r.t. to global coords at
    /// integration point ipt; return Jacobian of mapping (J). Also compute
    /// derivatives of dpsidx, dtestdx and J w.r.t. nodal coordinates.
    inline double dshape_and_dtest_eulerian_at_knot_hele_shaw(
     const unsigned &ipt,
     Shape &psi,
     DShape &dpsidx,
     RankFourTensor<double> &d_dpsidx_dX,
     Shape &test,
     DShape &dtestdx,
     RankFourTensor<double> &d_dtestdx_dX,
     DenseMatrix<double> &djacobian_dX) const;
  };
 
 
 
 
//Inline functions:
 
 
//======================================================================
/// Define the shape functions and test functions and derivatives
/// w.r.t. global coordinates and return Jacobian of mapping.
///
/// Galerkin: Test functions = shape functions
//======================================================================
 template<unsigned NNODE_1D>
  double QHeleShawElement<NNODE_1D>::dshape_and_dtest_eulerian_hele_shaw(
   const Vector<double> &s,
   Shape &psi,
   DShape &dpsidx,
   Shape &test,
   DShape &dtestdx) const
  {
   // Call the geometrical shape functions and derivatives
   const double J = this->dshape_eulerian(s,psi,dpsidx);
   
   // Set the test functions equal to the shape functions
   test = psi;
   dtestdx = dpsidx;
   
   // Return the jacobian
   return J;
  }
 
 
 
 
//======================================================================
/// Define the shape functions and test functions and derivatives
/// w.r.t. global coordinates and return Jacobian of mapping.
///
/// Galerkin: Test functions = shape functions
//======================================================================
 template<unsigned NNODE_1D>
  double QHeleShawElement<NNODE_1D>::
  dshape_and_dtest_eulerian_at_knot_hele_shaw(
   const unsigned &ipt,
   Shape &psi,
   DShape &dpsidx,
   Shape &test,
   DShape &dtestdx) const
  {
   // Call the geometrical shape functions and derivatives
   const double J = this->dshape_eulerian_at_knot(ipt,psi,dpsidx);
   
   // Set the pointers of the test functions
   test = psi;
   dtestdx = dpsidx;
   
   // Return the jacobian
   return J;
  }



 
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
 
 
/// Define the shape functions (psi) and test functions (test) and
/// their derivatives w.r.t. global coordinates (dpsidx and dtestdx)
/// and return Jacobian of mapping (J). Additionally compute the
/// derivatives of dpsidx, dtestdx and J w.r.t. nodal coordinates.
///
/// Galerkin: Test functions = shape functions
//======================================================================
 template<unsigned NNODE_1D>
  double QHeleShawElement<NNODE_1D>::
  dshape_and_dtest_eulerian_at_knot_hele_shaw(
   const unsigned &ipt,
   Shape &psi,
   DShape &dpsidx,
   RankFourTensor<double> &d_dpsidx_dX,
   Shape &test,
   DShape &dtestdx,
   RankFourTensor<double> &d_dtestdx_dX,
   DenseMatrix<double> &djacobian_dX) const
  {
   // Call the geometrical shape functions and derivatives
   const double J = this->dshape_eulerian_at_knot(ipt,psi,dpsidx,
                                                  djacobian_dX,d_dpsidx_dX);
   
   // Set the pointers of the test functions
   test = psi;
   dtestdx = dpsidx;
   d_dtestdx_dX = d_dpsidx_dX;
   
   // Return the jacobian
   return J;
  }
 
 
 
//=======================================================================
/// Face geometry for the QHeleShawElement elements: The spatial
/// dimension of the face elements is one lower than that of the
/// bulk element but they have the same number of points
/// along their 1D edges.
//=======================================================================
 template<unsigned NNODE_1D>
  class FaceGeometry<QHeleShawElement<NNODE_1D> >:
 public virtual QElement<1,NNODE_1D>
  {   
    public:
   
   /// \short Constructor: Call the constructor for the
   /// appropriate lower-dimensional QElement
   FaceGeometry() : QElement<1,NNODE_1D>() {}    
  };
 
 
 
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
 
 
//==========================================================
/// HeleShaw upgraded to become projectable
//==========================================================
 template<class HELE_SHAW_ELEMENT>
  class ProjectableHeleShawElement :
 public virtual ProjectableElement<HELE_SHAW_ELEMENT>
 {  
   public:
  
  /// \short Specify the values associated with field fld.
  /// The information is returned in a vector of pairs which comprise
  /// the Data object and the value within it, that correspond to field fld.
  Vector<std::pair<Data*,unsigned> > data_values_of_field(const unsigned& fld)
   {
#ifdef PARANOID
    if (fld!=0)
     {
      std::stringstream error_stream;
      error_stream
       << "HeleShaw elements only store a single field so fld must be 0 rather"
       << " than " << fld << std::endl;
      throw OomphLibError(
       error_stream.str(),
       "ProjectableHeleShawElement::data_values_of_field()",
       OOMPH_EXCEPTION_LOCATION);
     }
#endif
    
    // Create the vector
    Vector<std::pair<Data*,unsigned> > data_values;
    
    // Loop over all nodes
    unsigned nnod = this->nnode();
    for (unsigned j=0; j<nnod; j++)
     {
      // Add the data value associated field: The node itself
      data_values.push_back(std::make_pair(this->node_pt(j),fld));
     }
    
    // Return the vector
    return data_values;
   }
  
  /// \short Number of fields to be projected: Just one
  unsigned nfields_for_projection()
  {
   return 1;
  }
  
  /// \short Number of history values to be stored for fld-th field.
  unsigned nhistory_values_for_projection(const unsigned &fld)
  {
#ifdef PARANOID
   if (fld!=0)
    {
     std::stringstream error_stream;
     error_stream
      << "HeleShaw elements only store a single field so fld must be 0 rather"
      << " than " << fld << std::endl;
     throw OomphLibError(
      error_stream.str(),
      "ProjectableHeleShawElement::nhistory_values_for_projection()",
      OOMPH_EXCEPTION_LOCATION);
    }
#endif
   return this->node_pt(0)->ntstorage();
  }
  
  /// \short Number of positional history values
  unsigned nhistory_values_for_coordinate_projection()
  {
   return this->node_pt(0)->position_time_stepper_pt()->ntstorage();
  }
  
  /// \short Return Jacobian of mapping and shape functions of field fld
  /// at local coordinate s
  double jacobian_and_shape_of_field(const unsigned &fld,
                                     const Vector<double> &s,
                                     Shape &psi)
  {
#ifdef PARANOID
   if (fld!=0)
    {
     std::stringstream error_stream;
     error_stream
      << "HeleShaw elements only store a single field so fld must be 0 rather"
      << " than " << fld << std::endl;
     throw OomphLibError(
      error_stream.str(),
      "ProjectableHeleShawElement::jacobian_and_shape_of_field()",
      OOMPH_EXCEPTION_LOCATION);
    }
#endif
   unsigned n_dim = this->dim();
   unsigned n_node = this->nnode();
   Shape test(n_node);
   DShape dpsidx(n_node,n_dim), dtestdx(n_node,n_dim);
   double J = this->dshape_and_dtest_eulerian_hele_shaw(s,psi,dpsidx,
                                                        test,dtestdx);
   return J;
  }
    
  /// \short Return interpolated field fld at local coordinate s, at time level
  /// t (t=0: present; t>0: history values)
  double get_field(const unsigned &t,
                   const unsigned &fld,
                   const Vector<double>& s)
  {
#ifdef PARANOID
   if (fld!=0)
    {
     std::stringstream error_stream;
     error_stream
      << "HeleShaw elements only store a single field so fld must be 0 rather"
      << " than " << fld << std::endl;
     throw OomphLibError(
      error_stream.str(),
      "ProjectableHeleShawElement::jget_field()",
      OOMPH_EXCEPTION_LOCATION);
    }
#endif
   // Find the index at which the variable is stored
   unsigned p_nodal_index = this->p_index_hele_shaw();
   
   // Local shape function
   unsigned n_node=this->nnode();
   Shape psi(n_node);
   
   // Find values of shape function
   this->shape(s,psi);
   
   // Initialise value of u
   double interpolated_p = 0.0;
   
   // Sum over the local nodes
   for(unsigned l=0; l<n_node; l++)
    {
     interpolated_p += this->nodal_value(l,p_nodal_index)*psi[l];
    }
   return interpolated_p;
  }
    
  
  /// Return number of values in field fld: One per node
  unsigned nvalue_of_field(const unsigned &fld)
  {
#ifdef PARANOID
   if (fld!=0)
    {
     std::stringstream error_stream;
     error_stream
      << "HeleShaw elements only store a single field so fld must be 0 rather"
      << " than " << fld << std::endl;
     throw OomphLibError(
      error_stream.str(),
      "ProjectableHeleShawElement::nvalue_of_field()",
      OOMPH_EXCEPTION_LOCATION);
    }
#endif
   return this->nnode();
  }
  
  
  /// Return local equation number of value j in field fld.
  int local_equation(const unsigned &fld,
                     const unsigned &j)
  {
#ifdef PARANOID
   if (fld!=0)
    {
     std::stringstream error_stream;
     error_stream
      << "HeleShaw elements only store a single field so fld must be 0 rather"
      << " than " << fld << std::endl;
     throw OomphLibError(
      error_stream.str(),
      "ProjectableHeleShawElement::local_equation()",
      OOMPH_EXCEPTION_LOCATION);
    }
#endif
   const unsigned p_nodal_index = this->p_index_hele_shaw();
   return this->nodal_local_eqn(j,p_nodal_index);
  }
  
 };
 
 
//=======================================================================
/// Face geometry for element is the same as that for the underlying
/// wrapped element
//=======================================================================
 template<class ELEMENT>
  class FaceGeometry<ProjectableHeleShawElement<ELEMENT> >
  : public virtual FaceGeometry<ELEMENT>
 {
   public:
  FaceGeometry() : FaceGeometry<ELEMENT>() {}
 };
 
 
//=======================================================================
/// Face geometry of the Face Geometry for element is the same as
/// that for the underlying wrapped element
//=======================================================================
 template<class ELEMENT>
  class FaceGeometry<FaceGeometry<ProjectableHeleShawElement<ELEMENT> > >
  : public virtual FaceGeometry<FaceGeometry<ELEMENT> >
 {
   public:
  FaceGeometry() : FaceGeometry<FaceGeometry<ELEMENT> >() {}
 };
 
 
 
}






#endif


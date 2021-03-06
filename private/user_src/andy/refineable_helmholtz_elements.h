//Header file for refineable QHelmholtzElement elements

#ifndef OOMPH_REFINEABLE_HELMHOLTZ_ELEMENTS_HEADER
#define OOMPH_REFINEABLE_HELMHOLTZ_ELEMENTS_HEADER

// Config header generated by autoconfig
#ifdef HAVE_CONFIG_H
  #include <oomph-lib-config.h>
#endif


//oomph-lib headers
#include "../generic/refineable_quad_element.h"
#include "../generic/refineable_brick_element.h"
#include "../generic/error_estimator.h"
#include "helmholtz_elements.h"

namespace oomph
{

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////



//======================================================================
/// Refineable version of Helmholtz equations
///
///
//======================================================================
template <unsigned DIM>
class RefineableHelmholtzEquations : public virtual HelmholtzEquations<DIM>,
                                   public virtual RefineableElement,
                                   public virtual ElementWithZ2ErrorEstimator
{
  public:

 /// \short Constructor, simply call other constructors
 RefineableHelmholtzEquations() : HelmholtzEquations<DIM>(),
  RefineableElement(), ElementWithZ2ErrorEstimator() 
  { } 

 /// Broken copy constructor
 RefineableHelmholtzEquations(const RefineableHelmholtzEquations<DIM>& dummy) 
  { 
   BrokenCopy::broken_copy("RefineableHelmholtzEquations");
  } 
 
 /// Broken assignment operator
 void operator=(const RefineableHelmholtzEquations<DIM>&) 
  {
   BrokenCopy::broken_assign("RefineableHelmholtzEquations");
  }
 
 /// Number of 'flux' terms for Z2 error estimation 
 unsigned num_Z2_flux_terms() {return DIM;}

 /// Get 'flux' for Z2 error recovery:  Standard flux.from Helmholtz equations
 void get_Z2_flux(const Vector<double>& s, Vector<double>& flux)
  {this->get_flux(s,flux);}

/*  /// \short Overload assign_all_generic_local_eqn_numbers, since this */
/*  ///        function is overloaded in both RefineableElement AND */
/*  ///        HelmholtzEquations */
/*  virtual inline void assign_all_generic_local_eqn_numbers() */
/*   { */
/*    //Assign unique external data */
/*    assign_unique_external_data_helper(); */
/*    //Call the RefineableElement version (? - but this will assign nodals?) */
/*    RefineableElement::assign_all_generic_local_eqn_numbers(); */
/*   } */

/// \short Get the function value u in Vector.
/// Note: Given the generality of the interface (this function
/// is usually called from black-box documentation or interpolation routines),
/// the values Vector sets its own size in here.
void get_interpolated_values(const Vector<double>&s,  Vector<double>& values)
 {
  // Set size of Vector: u
  values.resize(1);
  
  //Find number of nodes
  unsigned n_node = nnode();
  
  //Local shape function
  Shape psi(n_node);
  
  //Find values of shape function
  shape(s,psi);
  
  //Initialise value of u
  values[0] = 0.0;

  //Find the index at which the helmholtz unknown is stored
  unsigned u_nodal_index = this->u_index_helmholtz();
  
  //Loop over the local nodes and sum up the values
  for(unsigned l=0;l<n_node;l++)
   {
    values[0] += this->nodal_value(l,u_nodal_index)*psi[l];
   }
 }


 /// \short Get the function value u in Vector.
 /// Note: Given the generality of the interface (this function
 /// is usually called from black-box documentation or interpolation routines),
 /// the values Vector sets its own size in here.
 void get_interpolated_values(const unsigned& t, const Vector<double>&s, 
                              Vector<double>& values)
  {
   if (t!=0)
    {
     std::string error_message =
      "Time-dependent version of get_interpolated_values() ";
     error_message += "not implemented for this element \n";
     throw 
      OomphLibError(error_message,
                    "RefineableHelmholtzEquations::get_interpolated_values()",
                    OOMPH_EXCEPTION_LOCATION);
    }
   else
    {
     //Make sure that we call this particular object's steady 
     //get_interpolated_values (it could get overloaded lower down)
     RefineableHelmholtzEquations<DIM>::get_interpolated_values(s,values);
    }
  }

 
 ///  Further build: Copy source function pointer from father element
 void further_build()
  {
   this->Source_fct_pt=dynamic_cast<RefineableHelmholtzEquations<DIM>*>(
    this->father_element_pt())->source_fct_pt();
  }


  private:


/// \short Add element's contribution to elemental residual vector and/or 
/// Jacobian matrix 
/// flag=1: compute both
/// flag=0: compute only residual vector
 void fill_in_generic_residual_contribution_helmholtz(
  Vector<double> &residuals, DenseMatrix<double> &jacobian, 
  unsigned flag); 
 
};


//======================================================================
/// Refineable version of 2D QHelmholtzElement elements
///
///
//======================================================================
template <unsigned DIM, unsigned NNODE_1D>
 class RefineableQHelmholtzElement : 
 public QHelmholtzElement<DIM,NNODE_1D>,
 public virtual RefineableHelmholtzEquations<DIM>,
 public virtual RefineableQElement<DIM>
{
  public:

 /// \short Constructor, simply call the other constructors 
 RefineableQHelmholtzElement() : 
  RefineableElement(),
  RefineableHelmholtzEquations<DIM>(),
  RefineableQElement<DIM>(),
  QHelmholtzElement<DIM,NNODE_1D>()
   {} 


 /// Broken copy constructor
 RefineableQHelmholtzElement(const RefineableQHelmholtzElement<DIM,NNODE_1D>& 
                           dummy) 
  { 
   BrokenCopy::broken_copy("RefineableQuadHelmholtzElement");
  } 
 
 /// Broken assignment operator
 void operator=(const RefineableQHelmholtzElement<DIM,NNODE_1D>&) 
  {
   BrokenCopy::broken_assign("RefineableQuadHelmholtzElement");
  }
 
 /// Number of continuously interpolated values: 1
 unsigned ncont_interpolated_values() const {return 1;}

 /// \short Number of vertex nodes in the element
 unsigned nvertex_node() const
  {return QHelmholtzElement<DIM,NNODE_1D>::nvertex_node();}

 /// \short Pointer to the j-th vertex node in the element
 Node* vertex_node_pt(const unsigned& j) const
  {return QHelmholtzElement<DIM,NNODE_1D>::vertex_node_pt(j);}

 /// Rebuild from sons: empty
 void rebuild_from_sons(Mesh* &mesh_pt) {}

 /// \short Order of recovery shape functions for Z2 error estimation:
 /// Same order as shape functions.
 unsigned nrecovery_order() {return (NNODE_1D-1);}

 ///  \short Perform additional hanging node procedures for variables
 /// that are not interpolated by all nodes. Empty.
 void further_setup_hanging_nodes(){}

};
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////



//=======================================================================
/// Face geometry for the RefineableQuadHelmholtzElement elements: The spatial 
/// dimension of the face elements is one lower than that of the
/// bulk element but they have the same number of points
/// along their 1D edges.
//=======================================================================
template<unsigned DIM, unsigned NNODE_1D>
class FaceGeometry<RefineableQHelmholtzElement<DIM,NNODE_1D> >: 
 public virtual QElement<DIM-1,NNODE_1D>
{

  public:
 
 /// \short Constructor: Call the constructor for the
 /// appropriate lower-dimensional QElement
 FaceGeometry() : QElement<DIM-1,NNODE_1D>() {}

};

}

#endif


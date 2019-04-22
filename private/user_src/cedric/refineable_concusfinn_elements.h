
#ifndef OOMPH_REFINEABLE_YOUNGLAPLACE_ELEMENTS_HEADER
#define OOMPH_REFINEABLE_YOUNGLAPLACE_ELEMENTS_HEADER

// Config header generated by autoconfig
#ifdef HAVE_CONFIG_H
  #include <oomph-lib-config.h>
#endif


//oomph-lib headers
#include "generic/refineable_quad_element.h"
#include "generic/refineable_brick_element.h"
#include "generic/error_estimator.h"
#include "concusfinn_elements.h"


namespace oomph
{


///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////



//======================================================================
/// Refineable version of YoungLaplace equations
///
///
//======================================================================
class RefineableYoungLaplaceEquations : public virtual YoungLaplaceEquations,
                                    public virtual RefineableElement,
                                    public virtual ElementWithZ2ErrorEstimator
{
  public:

 /// \short Constructor: Pass refinement level to refineable element
 /// (default 0 = root)
 RefineableYoungLaplaceEquations(const int& refine_level=0) : 
  YoungLaplaceEquations(),
  RefineableElement(),
  ElementWithZ2ErrorEstimator()
  {} 


 /// Broken copy constructor
 RefineableYoungLaplaceEquations(const RefineableYoungLaplaceEquations& dummy)
  { 
   BrokenCopy::broken_copy("RefineableYoungLaplaceEquations");
  } 
 
 /// Broken assignment operator
 void operator=(const RefineableYoungLaplaceEquations&) 
  {
   BrokenCopy::broken_assign("RefineableYoungLaplaceEquations");
  }
 
 /// Compute element residual vector taking hanging nodes into account
 void fill_in_contribution_to_residuals(Vector<double> &residuals);

 /// Number of 'flux' terms for Z2 error estimation 
 unsigned num_Z2_flux_terms() {return 2;}

 /// Get 'flux' for Z2 error recovery:  Standard flux.from YoungLaplace equations
 void get_Z2_flux(const Vector<double>& s, Vector<double>& flux)
  {this->get_flux(s,flux);}

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
  
  //Loop over the local nodes and sum
  for(unsigned l=0;l<n_node;l++)
   {
    values[0] += this->u(l)*psi[l];
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
     throw OomphLibError(
      "These equations are steady => No time dependence",
      "RefineableYoungLaplaceEquations::get_interpolated_values()",
      OOMPH_EXCEPTION_LOCATION);
    }
   else
    {
     get_interpolated_values(s,values);
    }
  }

 
 /// \short Further build: Copy function pointers to spine and spine base
 /// functions from father. Kappa is passed across with set_kappa(...)
 /// to ensure that it's added to the element's external Data.
 void further_build()
  {
   // Copy kappa with set_kappa() to ensure that it's added to the
   // element's external Data.
   this->set_kappa(dynamic_cast<RefineableYoungLaplaceEquations*>(
    this->father_element_pt())->kappa_pt());
   
   // Copy spine functions
   this->Spine_fct_pt=dynamic_cast<RefineableYoungLaplaceEquations*>(
    this->father_element_pt())->spine_fct_pt();
   
   this->Spine_base_fct_pt=dynamic_cast<RefineableYoungLaplaceEquations*>(
   this->father_element_pt())->spine_base_fct_pt();
  }

};


//======================================================================
/// Refineable version of 2D QYoungLaplaceElement elements
///
///
//======================================================================
template <unsigned NNODE_1D>
class RefineableQYoungLaplaceElement : 
public QYoungLaplaceElement<NNODE_1D>,
           public virtual RefineableYoungLaplaceEquations,
           public virtual RefineableQElement<2>
{
  public:

 /// \short Constructor: Pass refinement level to refineable quad element
 /// (default 0 = root)
  RefineableQYoungLaplaceElement() : 
   RefineableElement(),
   RefineableYoungLaplaceEquations(),
   RefineableQElement<2>(),
   QYoungLaplaceElement<NNODE_1D>()
   {} 

 /// Broken copy constructor
 RefineableQYoungLaplaceElement(const RefineableQYoungLaplaceElement<NNODE_1D>&
                              dummy) 
  { 
   BrokenCopy::broken_copy("RefineableQuadYoungLaplaceElement");
  } 
 
 /// Broken assignment operator
 void operator=(const RefineableQYoungLaplaceElement<NNODE_1D>&) 
  {
   BrokenCopy::broken_assign("RefineableQuadYoungLaplaceElement");
  }
 
 /// Number of continuously interpolated values: 1
 unsigned ncont_interpolated_values() const {return 1;}

 /// \short Number of vertex nodes in the element
 unsigned nvertex_node() const
  {return QYoungLaplaceElement<NNODE_1D>::nvertex_node();}

 /// \short Pointer to the j-th vertex node in the element
 Node* vertex_node_pt(const unsigned& j) const
  {return QYoungLaplaceElement<NNODE_1D>::vertex_node_pt(j);}

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
/// Face geometry for the RefineableQuadYoungLaplaceElement elements: The spatial
/// dimension of the face elements is one lower than that of the
/// bulk element but they have the same number of points
/// along their 1D edges.
//=======================================================================
template<unsigned NNODE_1D>
class FaceGeometry<RefineableQYoungLaplaceElement<NNODE_1D> >: 
 public virtual QElement<1,NNODE_1D>
{

  public:
 
 /// \short Constructor: Call the constructor for the
 /// appropriate lower-dimensional QElement
 FaceGeometry() : QElement<1,NNODE_1D>() {}

};


}

#endif


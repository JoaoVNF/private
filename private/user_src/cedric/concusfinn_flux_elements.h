// Header file for elements that are used to apply contact angle
// BC for Young Laplace eqns

#ifndef OOMPH_YOUNGLAPLACE_FLUX_ELEMENTS_HEADER
#define OOMPH_YOUNGLAPLACE_FLUX_ELEMENTS_HEADER


// Config header generated by autoconfig
#ifdef HAVE_CONFIG_H
  #include <oomph-lib-config.h>
#endif


// oomph-lib ncludes
#include "generic/Qelements.h"

namespace oomph
{

//======================================================================
/// \short A class for elements that allow the imposition of an 
/// contact angle bcs for Young Laplace elements.
/// The element geometry is obtained from the  FaceGeometry<ELEMENT> 
/// policy class.
//======================================================================
template <class ELEMENT>
class YoungLaplaceFluxElement : public virtual FaceGeometry<ELEMENT>, 
public virtual FaceElement
{
 
public:

 /// \short Constructor, takes the pointer to the "bulk" element, the 
 /// index of the fixed local coordinate and its value represented
 /// by an integer (+/- 1), indicating that the face is located
 /// at the max. or min. value of the "fixed" local coordinate
 /// in the bulk element.
 YoungLaplaceFluxElement(FiniteElement* const &bulk_el_pt, 
                         const int& face_index);
 

 ///\short  Broken empty constructor
 YoungLaplaceFluxElement()
  {
   throw OomphLibError(
    "Don't call empty constructor for YoungLaplaceFluxElement",
    "YoungLaplaceFluxElement::YoungLaplaceFluxElement",
    OOMPH_EXCEPTION_LOCATION);
  }

 /// Broken copy constructor
 YoungLaplaceFluxElement(const YoungLaplaceFluxElement& dummy) 
  { 
   BrokenCopy::broken_copy("YoungLaplaceFluxElement");
  } 
 
 /// Broken assignment operator
 void operator=(const YoungLaplaceFluxElement&) 
  {
   BrokenCopy::broken_assign("YoungLaplaceFluxElement");
  }


 /// \short Access function for the pointer to the prescribed contact angle
 /// (const version)
 double* prescribed_cos_gamma_pt() const
  {
   return Prescribed_cos_gamma_pt;
  }


 /// Access function for the pointer to the prescribed contact angle
 double*& prescribed_cos_gamma_pt() 
  {
   return Prescribed_cos_gamma_pt;
  }


 /// Add the element's contribution to its residual vector
 void fill_in_contribution_to_residuals(Vector<double> &residuals)
  {
    //old_fill_in_contribution_to_residuals(residuals);
    new_fill_in_contribution_to_residuals(residuals);
  }

 /// \short Get the local equation number of the (one and only) unknown
 /// stored at local node n (returns -1 if value is pinned).
 /// Can be overloaded in derived multiphysics elements.
 inline int u_local_eqn(const unsigned& n)
  {
   //Local equation number is the first value stored at the node
   return nodal_local_eqn(n,0);
  }

 /// Output function -- forward to broken version in FiniteElement
 /// until somebody decides what exactly they want to plot here...
 void output(std::ostream &outfile) {FiniteElement::output(outfile);}

 /// \short Output function -- forward to broken version in FiniteElement
 /// until somebody decides what exactly they want to plot here...
 void output(std::ostream &outfile, const unsigned &n_plot)
  {FiniteElement::output(outfile,n_plot);}


 /// C-style output function -- forward to broken version in FiniteElement
 /// until somebody decides what exactly they want to plot here...
 void output(FILE* file_pt) {FiniteElement::output(file_pt);}

 /// \short C-style output function -- forward to broken version in 
 /// FiniteElement until somebody decides what exactly they want to plot 
 /// here...
 void output(FILE* file_pt, const unsigned &n_plot)
  {FiniteElement::output(file_pt,n_plot);}

 /// Compute cosinus of actual contact angle
 double actual_cos_contact_angle(const Vector<double>& s);

 /// Calculate the length of the face element
 double get_projected_length();

 /// Get tangent and normal vectors to contact line
 void contact_line_vectors(const Vector<double>& s,
                           Vector<double>& tangent,
                           Vector<double>& normal)
  {
   // Dummy
   double norm_of_drds;
   Vector<double> spine(3);
   contact_line_vectors(s,tangent,normal,spine,norm_of_drds);
  }

 /// \short Get tangent and normal vectors to contact line. 
 /// Final argument gives the norm of the tangent
 void contact_line_vectors(const Vector<double>& s,
                           Vector<double>& tangent,
                           Vector<double>& normal,
			   double& norm_of_drds)
  {
   // Dummy
   Vector<double> spine(3);
   contact_line_vectors(s,tangent,normal,spine,norm_of_drds);
  }


 /// \short Get tangent and normal vectors to contact line and spine vector 
 /// (wall normal can then be obtained by cross-product). Final
 /// argument gives the norm of dR/ds where R is the vector to the
 /// contact line and s the local coordinate in the element
 void contact_line_vectors(const Vector<double>& s,
                           Vector<double>& tangent,
                           Vector<double>& normal,
                           Vector<double>& spine,
                           double& norm_of_drds);

protected:

 /// \short Define an access function to the first data value stored 
 /// at each node. In a more general "Equation" element, such abstraction
 /// is essential, because different Elements will store the same variables
 /// in different locations.
 double &u(const unsigned int &n) {return *this->node_pt(n)->value_pt(0);}

 /// Function to calculate the cos of the prescribed contact angle
 double prescribed_cos_gamma()
  {
   //If the function pointer is zero return zero
   if(Prescribed_cos_gamma_pt == 0)
    {
     return 0.0;
    }
   //Otherwise de-reference pointer
   else
    {
     return *Prescribed_cos_gamma_pt;
    }
  }


private:


 /// \short Add the element's contribution to its residual vector.
 /// old version
 //void old_fill_in_contribution_to_residuals(Vector<double> &residuals);


 /// \short Add the element's contribution to its residual vector.
 /// new version
 void new_fill_in_contribution_to_residuals(Vector<double> &residuals);

 
 /// Pointer to prescribed cos gamma
 double* Prescribed_cos_gamma_pt;


}; 


}

#endif

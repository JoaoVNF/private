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
// Header file for elements that are used to apply prescribed flux
// boundary conditions to the Poisson equations
#ifndef OOMPH_HELE_SHAW_FLUX_ELEMENTS_HEADER
#define OOMPH_HELE_SHAW_FLUX_ELEMENTS_HEADER


// Config header generated by autoconfig
#ifdef HAVE_CONFIG_H
  #include <oomph-lib-config.h>
#endif

// oomph-lib includes
#include "../../src/generic/Qelements.h"

namespace oomph
{

//======================================================================
/// \short A class for elements that allow the imposition of an
/// applied flux on the boundaries of Poisson elements.
/// The element geometry is obtained from the  FaceGeometry<ELEMENT>
/// policy class.
//======================================================================
template <class ELEMENT>
class HeleShawFluxElement : public virtual FaceGeometry<ELEMENT>,
public virtual FaceElement
{

public:

 /// \short Function pointer to the prescribed-flux function fct(w,f(x)) --
 /// w is the deflection of the sheet!
 typedef void (*PrescribedFluxFctPt)(const double& w, 
                                     const Vector<double>& dwdx, 
                                     double& flux);
   
 /// \short Function pointer to function which provides wall speed as a 
 /// function of x. This allows us to solve for bubble motion in a moving 
 /// frame. A constant wall speed does not affect the mass conservation 
 /// equations, but does feature in the kinematic equation for interface 
 /// motion.
 typedef void (*WallSpeedFctPt)(const Vector<double>& x, 
                                Vector<double>& U_wall);

 /// \short Constructor, takes the pointer to the "bulk" element and the
 /// index of the face to which the element is attached.
 HeleShawFluxElement(FiniteElement* const &bulk_el_pt,
                     const int& face_index);

 ///\short  Broken empty constructor
 HeleShawFluxElement()
  {
   throw OomphLibError(
    "Don't call empty constructor for HeleShawFluxElement",
    "HeleShawFluxElement::HeleShawFluxElement()",
    OOMPH_EXCEPTION_LOCATION);
  }

 /// Broken copy constructor
 HeleShawFluxElement(const HeleShawFluxElement& dummy)
  {
   BrokenCopy::broken_copy("HeleShawFluxElement");
  }

 /// Broken assignment operator
 void operator=(const HeleShawFluxElement&)
  {
   BrokenCopy::broken_assign("HeleShawFluxElement");
  }
    
 /// function to enable thin film effects
 void enable_thin_film_effects()
 {
  Apply_thin_film_effects=true;
 }
 
 /// function to enable thin film effects
 void disable_thin_film_effects()
 {
  Apply_thin_film_effects=false;
 }

 /// \short Specify the value of nodal zeta from the face geometry
 /// The "global" intrinsic coordinate of the element when
 /// viewed as part of a geometric object should be given by
 /// the FaceElement representation, by default (needed to break
 /// indeterminacy if bulk element is SolidElement)
 double zeta_nodal(const unsigned &n, const unsigned &k,
                          const unsigned &i) const
  {return FaceElement::zeta_nodal(n,k,i);}

 /// Access function for the prescribed-flux function pointer
 PrescribedFluxFctPt& flux_fct_pt() {return Flux_fct_pt;}


 /// Add the element's contribution to its residual vector
 inline void fill_in_contribution_to_residuals(Vector<double> &residuals)
  {
   //Call the generic residuals function with flag set to 0
   //using a dummy matrix argument
   fill_in_generic_residual_contribution_poisson_flux(
    residuals,GeneralisedElement::Dummy_matrix,0);
  }


 /* /// \short Add the element's contribution to its residual vector and its */
 /* /// Jacobian matrix */
 /* inline void fill_in_contribution_to_jacobian(Vector<double> &residuals, */
 /*                                          DenseMatrix<double> &jacobian) */
 /*  { */
 /*   //Call the generic routine with the flag set to 1 */
 /*   fill_in_generic_residual_contribution_poisson_flux(residuals,jacobian,1); */
 /*   //oomph_info<<"Contribution from FluxElement"<<std::endl; */
 /*   GeneralisedElement::fill_in_jacobian_from_external_by_fd(residuals, */
 /*                                                            jacobian,true); */
 /*  } */
 
 ///\short Compute the element's residual vector and the Jacobian matrix.
 /// Jacobian is computed by finite-differencing.
 void fill_in_contribution_to_jacobian(Vector<double> &residuals,
                                       DenseMatrix<double> &jacobian)
 {
  //Add the contribution to the residuals
  this->fill_in_contribution_to_residuals(residuals);
  
  //Allocate storage for the full residuals (residuals of entire element)
  unsigned n_dof = this->ndof();
  Vector<double> full_residuals(n_dof);
  
  //Get the residuals for the entire element
  this->get_residuals(full_residuals);
  
  //Get the solid entries in the jacobian using finite differences
  this->fill_in_jacobian_from_solid_position_by_fd(full_residuals,jacobian);
  
  //There could be internal data
  //(finite-difference the lot by default)
  this->fill_in_jacobian_from_internal_by_fd(full_residuals,jacobian,true);
  
  //There could also be external data
  //(finite-difference the lot by default)
  //oomph_info<<"Contribution from BulkElement"<<std::endl;
  this->fill_in_jacobian_from_external_by_fd(full_residuals,jacobian,true);
  
  //There could also be nodal data
  this->fill_in_jacobian_from_nodal_by_fd(full_residuals,jacobian);
 }

 void fill_in_contribution_to_jacobian_and_mass_matrix(
  Vector<double> &residuals,
  DenseMatrix<double> &jacobian,
  DenseMatrix<double> &mass_matrix
  )
 {
  //std::cout << "Calling jacobian and mass matrix in flux " << std::endl;
  fill_in_contribution_to_jacobian(residuals,jacobian);
 }
    
 /// The value of the inverse Capillary number
 const double& ca() const
 {
  return *Ca_pt;
 }
 
 /// Pointer to the inverse Capillary number
 double*& ca_pt() {return Ca_pt;}
    
 /// Aspect ratio: Reference gap width / in-plane lengthscale
 double aspect_ratio() const
 {
#ifdef PARANOID
  if (Aspect_ratio_pt==0)
   {
    throw OomphLibError(
     "Aspect_ratio_pt has not been set yet for HeleShaw elements",
     "HeleShawInterfaceElement::aspect_ratio()",
     OOMPH_EXCEPTION_LOCATION);
   }
#endif
  return *Aspect_ratio_pt;
 }
 
 /// \short Pointer to aspect ratio
 double*& aspect_ratio_pt()
 {
  return Aspect_ratio_pt;
 }    
 
 /// \short Pointer to aspect ratio. Const version.
 double* aspect_ratio_pt() const
 {
  return Aspect_ratio_pt;
 }
   
 /// Access function: Pointer to wall speed function
 WallSpeedFctPt& wall_speed_fct_pt() {return Wall_speed_fct_pt;}
   
 /// Access function: Pointer to wall speed function. Const version
 WallSpeedFctPt wall_speed_fct_pt() const {return Wall_speed_fct_pt;}
    
 /// Return the value of the in-plane speed of the plates
 void get_wall_velocity(Vector<double>&x, Vector<double>& U)
 {
  if (Wall_speed_fct_pt==0)
   {
    U[0]=0.0;
    U[1]=0.0;
   }
  else
   {
    (*Wall_speed_fct_pt)(x,U);
   }
 }
    
 /// Thin film homotopy parameter
 void thin_film_homotopy(double& parameter) const
 {
  if (Thin_film_homotopy_pt==0)
   {
    parameter = 1.0;
   }
  else
   {
    parameter = *Thin_film_homotopy_pt;
   }
 }
 
 /// \short Pointer to thin film homotopy parameter
 double*& thin_film_homotopy_pt()
 {
  return Thin_film_homotopy_pt;
 }    
 
 /// \short Pointer to thin film homotopy parameter. Const version.
 double* thin_film_homotopy_pt() const
 {
  return Thin_film_homotopy_pt;
 }
 
 /// Function to calculate the volume flux across the boundary
 /// Q = \int -h^3 dp/dx - h dS
 double compute_volume_flux()
 {  
  //Find out how many nodes there are
  const unsigned n_node = nnode();
  
  //Set up memory for the shape and test functions
  Shape psif(n_node), testf(n_node);
  
  //Set the value of Nintpt
  const unsigned n_intpt = integral_pt()->nweight();
  
  //Set the Vector to hold local coordinates
  Vector<double> s(Dim-1);

  // Vector for local coordinates in bulk element
  Vector<double> s_bulk(Dim);
   
  // Get pointer to assocated bulk element
  ELEMENT* bulk_el_pt=dynamic_cast<ELEMENT*>(bulk_element_pt());

  // Initialise volume flux
  double volume_flux = 0.0;
  
  //Loop over the integration points
  //--------------------------------
  for(unsigned ipt=0;ipt<n_intpt;ipt++)
   {
    //Assign values of s
    for(unsigned i=0;i<(Dim-1);i++) {s[i] = integral_pt()->knot(ipt,i);}

    //Get the bulk coordinates
    this->get_local_coordinate_in_bulk(s,s_bulk);
    
    //Get the integral weight
    double w = integral_pt()->weight(ipt);
    
    //Find the shape and test functions and return the Jacobian
    //of the mapping
    double J = shape_and_test(s,psif,testf);
    
    //Premultiply the weights and the Jacobian
    double W = w*J;
    
    Vector<double> interpolated_x(2,0.0);
    
    //Calculate velocities and derivatives
    for(unsigned l=0;l<n_node;l++)
     {
      // Loop over directional components
      for(unsigned i=0; i<2; i++)
       {
        // Spatial bits
        interpolated_x[i] += this->raw_nodal_position(l,i)*psif[l];
       }
     }
    
    // Get the pressure gradient
    Vector<double> pressure_gradient(2,0.0);
    bulk_el_pt->get_pressure_gradient(s_bulk,pressure_gradient);
     
    // Get the velocity of the plates
    Vector<double> U_wall(2,0.0);
    get_wall_velocity(interpolated_x,U_wall);
     
    // Now calculate the unit normal vector
    Vector<double> interpolated_n(2);
    outer_unit_normal(ipt,interpolated_n);
    
    // Get h from bulk! 
    double b = 0.0;
    double dbdt = 0.0;
    Vector<double> dbdx(2,0.0);
    
    //---------------------------------------------------------
    // NOTE: This is dangerous/wrong. We're really supposed to
    // pass in dpsi/dx (2D!) which we don't have/need here.
    // However it's only used for the computation of the ALE bits
    // in dh/dt and dh/dt isn't used here! Phew...
    //---------------------------------------------------------
    // Number of nodes in bulk element
    unsigned n_nod_bulk=bulk_el_pt->nnode();
    DShape dpsidx_bulk(n_nod_bulk,2);
    Shape psi_bulk(n_nod_bulk);
    bulk_el_pt->dshape_eulerian(s_bulk,psi_bulk,dpsidx_bulk);
    bulk_el_pt->get_upper_wall_data(s_bulk,interpolated_x,
                                    psi_bulk,dpsidx_bulk,b,dbdt,dbdx);

    //Now add to the contributions to the volume flux
    
    //Loop over the test functions
    for(unsigned l=0;l<n_node;l++)
     {
      // Loop over directional components
      for(unsigned i=0; i<2; i++)
       {
        //Add the contribution
        volume_flux -= pow(b,3.0)*pressure_gradient[i]*
         interpolated_n[i]*testf[l]*W;

        // Moving frame terms
        volume_flux -= U_wall[i]*b*interpolated_n[i]*testf[l]*W;
       }
     }
   }

  return volume_flux;
 }

 /// Function to calculate the cross section of the boundary
 /// A = \int b dS
 double compute_cross_section()
 {  
  //Find out how many nodes there are
  const unsigned n_node = nnode();
  
  //Set up memory for the shape and test functions
  Shape psif(n_node), testf(n_node);
  
  //Set the value of Nintpt
  const unsigned n_intpt = integral_pt()->nweight();
  
  //Set the Vector to hold local coordinates
  Vector<double> s(Dim-1);

  // Vector for local coordinates in bulk element
  Vector<double> s_bulk(Dim);
   
  // Get pointer to assocated bulk element
  ELEMENT* bulk_el_pt=dynamic_cast<ELEMENT*>(bulk_element_pt());

  // Initialise cross section
  double cross_section = 0.0;
  
  //Loop over the integration points
  //--------------------------------
  for(unsigned ipt=0;ipt<n_intpt;ipt++)
   {
    //Assign values of s
    for(unsigned i=0;i<(Dim-1);i++) {s[i] = integral_pt()->knot(ipt,i);}

    //Get the bulk coordinates
    this->get_local_coordinate_in_bulk(s,s_bulk);
    
    //Get the integral weight
    double w = integral_pt()->weight(ipt);
    
    //Find the shape and test functions and return the Jacobian
    //of the mapping
    double J = shape_and_test(s,psif,testf);
    
    //Premultiply the weights and the Jacobian
    double W = w*J;
    
    Vector<double> interpolated_x(2,0.0);
    
    //Calculate velocities and derivatives
    for(unsigned l=0;l<n_node;l++)
     {
      // Loop over directional components
      for(unsigned i=0; i<2; i++)
       {
        // Spatial bits
        interpolated_x[i] += this->nodal_position(l,i)*psif[l];
       }
     }
    
    // Get b from bulk!
    double b = 0.0;
    double dbdt = 0.0;
    Vector<double> dbdx(2,0.0);
    
    //---------------------------------------------------------
    // NOTE: This is dangerous/wrong. We're really supposed to
    // pass in dpsi/dx (2D!) which we don't have/need here.
    // However it's only used for the computation of the ALE bits
    // in dh/dt and dh/dt isn't used here! Phew...
    //---------------------------------------------------------
    // Number of nodes in bulk element
    unsigned n_nod_bulk=bulk_el_pt->nnode();
    DShape dpsidx_bulk(n_nod_bulk,2);
    Shape psi_bulk(n_nod_bulk);
    bulk_el_pt->dshape_eulerian(s_bulk,psi_bulk,dpsidx_bulk);
    bulk_el_pt->get_upper_wall_data(s_bulk,interpolated_x,
                                    psi_bulk,dpsidx_bulk,b,dbdt,dbdx);
    
    //Now add to the contributions to the volume flux

    if(Apply_thin_film_effects)
     {
      // Get the value of the Capillary number
      Vector<double> x(2,0.0);                                                    /////////////////////////////Joao changing the thin film dependency from ca() to ca()*U_frame_modulus
      Vector<double> U_wall(2,0.0);                                               /////////////////////////////Joao changing the thin film dependency from ca() to ca()*U_frame_modulus
      get_wall_velocity(x,U_wall);                                                /////////////////////////////Joao changing the thin film dependency from ca() to ca()*U_frame_modulus 
      double U_frame_modulus = sqrt(U_wall[0]*U_wall[0] + U_wall[1]*U_wall[1]);   /////////////////////////////Joao changing the thin film dependency from ca() to ca()*U_frame_modulus
      double Ca = ca()*U_frame_modulus;                                           /////////////////////////////Joao changing the thin film dependency from ca() to ca()*U_frame_modulus
//      double Ca = ca();
      b -= thin_film_effect(Ca)*b;
     }

    //Loop over the test functions
    for(unsigned l=0;l<n_node;l++)
     {
      //Add the contribution
      cross_section += b*testf[l]*W;
     }
   }

  return cross_section;
 }

 /// Function to calculate the mass flow due to thin films
 double compute_thin_film_mass_flow()
 {  
  //Find out how many nodes there are
  const unsigned n_node = nnode();
  
  //Set up memory for the shape and test functions
  Shape psif(n_node), testf(n_node);
  
  //Set the value of Nintpt
  const unsigned n_intpt = integral_pt()->nweight();
  
  //Set the Vector to hold local coordinates
  Vector<double> s(Dim-1);

  // Vector for local coordinates in bulk element
  Vector<double> s_bulk(Dim);
   
  // Get pointer to assocated bulk element
  ELEMENT* bulk_el_pt=dynamic_cast<ELEMENT*>(bulk_element_pt());

  // Get the value of the Capillary number
  double Ca = ca();

  // Initialise mass flow
  double mass_flow = 0.0;

  //Loop over the integration points
  //--------------------------------
  for(unsigned ipt=0;ipt<n_intpt;ipt++)
   {
    //Assign values of s
    for(unsigned i=0;i<(Dim-1);i++) {s[i] = integral_pt()->knot(ipt,i);}

    //Get the bulk coordinates
    this->get_local_coordinate_in_bulk(s,s_bulk);
    
    //Get the integral weight
    double w = integral_pt()->weight(ipt);
    
    //Find the shape and test functions and return the Jacobian
    //of the mapping
    double J = shape_and_test(s,psif,testf);
    
    //Premultiply the weights and the Jacobian
    double W = w*J;
    
    Vector<double> interpolated_x(2,0.0);
    
    //Calculate velocities and derivatives
    for(unsigned l=0;l<n_node;l++)
     {
      // Loop over directional components
      for(unsigned i=0; i<2; i++)
       {
        // Spatial bits
        interpolated_x[i] += this->nodal_position(l,i)*psif[l];
       }
     }
    
    // Get the pressure gradient
    Vector<double> pressure_gradient(2,0.0);
    bulk_el_pt->get_pressure_gradient(s_bulk,pressure_gradient);
     
    // Get the velocity of the plates
    Vector<double> U_wall(2,0.0);
    get_wall_velocity(interpolated_x,U_wall);
     
    // Now calculate the unit normal vector
    Vector<double> interpolated_n(2);
    outer_unit_normal(ipt,interpolated_n);
    
    // Get h from bulk! 
    double b = 0.0;
    double dbdt = 0.0;
    Vector<double> dbdx(2,0.0);
    
    //---------------------------------------------------------
    // NOTE: This is dangerous/wrong. We're really supposed to
    // pass in dpsi/dx (2D!) which we don't have/need here.
    // However it's only used for the computation of the ALE bits
    // in dh/dt and dh/dt isn't used here! Phew...
    //---------------------------------------------------------
    // Number of nodes in bulk element
    unsigned n_nod_bulk=bulk_el_pt->nnode();
    DShape dpsidx_bulk(n_nod_bulk,2);
    Shape psi_bulk(n_nod_bulk);
    //bulk_el_pt->shape(s_bulk,psi_bulk);
    // obacht correct?
    bulk_el_pt->dshape_eulerian(s_bulk,psi_bulk,dpsidx_bulk);
    bulk_el_pt->get_upper_wall_data(s_bulk,interpolated_x,
                                    psi_bulk,dpsidx_bulk,b,dbdt,dbdx);
    
    double combined_gradient = dbdx[0]*dbdx[1];
    double sum_square_gradients = dbdx[0]*dbdx[0]+dbdx[1]*dbdx[1];
    
    //Now add to the contributions to the volume flux

    // Modulus of the frame velocity                                                                  //////Joao
    double U_frame_modulus = sqrt(U_wall[0]*U_wall[0] + U_wall[1]*U_wall[1]);  /////////////////////////////Joao changing the thin film dependency from Ca to Ca*U_frame_modulus

    //Loop over the test functions
    for(unsigned l=0;l<n_node;l++)
     {
      // Loop over directional components
      for(unsigned i=0; i<2; i++)
       {
        //Add the contribution
        mass_flow -= pow(b,3.0)*pow(thin_film_effect(Ca*U_frame_modulus),3.0)*            //////Joao
         pressure_gradient[i]*interpolated_n[i]*testf[l]*W;
        mass_flow -= U_wall[i]*b*thin_film_effect(Ca*U_frame_modulus)*                    //////Joao
         interpolated_n[i]*testf[l]*W;
       }
     }
   }

  return mass_flow;
 }

 /// Output function -- forward to broken version in FiniteElement
 /// until somebody decides what exactly they want to plot here...
 void output(std::ostream &outfile)// {FiniteElement::output(outfile);}
 {
  //Find out how many nodes there are
  const unsigned n_node = nnode();
  
  //Set up memory for the shape and test functions
  Shape psif(n_node);
  
  //Set the value of Nintpt
  const unsigned n_intpt = integral_pt()->nweight();
  
  //Set the Vector to hold local coordinates
  Vector<double> s(Dim-1);

  // Vector for local coordinates in bulk element
  Vector<double> s_bulk(Dim);
   
  // Get pointer to assocated bulk element
  ELEMENT* bulk_el_pt=dynamic_cast<ELEMENT*>(bulk_element_pt());
  
  //Loop over the integration points
  //--------------------------------
  for(unsigned ipt=0;ipt<n_intpt;ipt++)
   {
    //Assign values of s
    for(unsigned i=0;i<(Dim-1);i++) {s[i] = integral_pt()->knot(ipt,i);}

    //Get the bulk coordinates
    this->get_local_coordinate_in_bulk(s,s_bulk);
    
    //Find the shape and test functions and return the Jacobian
    //of the mapping
    shape(s,psif);
    
    Vector<double> interpolated_x(2,0.0);
    
    //Calculate velocities and derivatives
    for(unsigned l=0;l<n_node;l++)
     {
      // Loop over directional components
      for(unsigned i=0; i<2; i++)
       {
        // Spatial bits
        interpolated_x[i] += this->raw_nodal_position(l,i)*psif[l];
       }
     }
    
    // Get the pressure gradient
    Vector<double> pressure_gradient(2,0.0);
    bulk_el_pt->get_pressure_gradient(s_bulk,pressure_gradient);
     
    // Get the velocity of the plates
    Vector<double> U_wall(2,0.0);
    get_wall_velocity(interpolated_x,U_wall);
     
    // Now calculate the unit normal vector
    Vector<double> interpolated_n(2);
    outer_unit_normal(ipt,interpolated_n);

    double b = 0.0;
    double dbdt = 0.0;
    Vector<double> dbdx(2,0.0);

    // Number of nodes in bulk element
    unsigned n_nod_bulk=bulk_el_pt->nnode();
    DShape dpsidx_bulk(n_nod_bulk,2);
    Shape psi_bulk(n_nod_bulk);
    bulk_el_pt->dshape_eulerian(s_bulk,psi_bulk,dpsidx_bulk);
    bulk_el_pt->get_upper_wall_data(s_bulk,interpolated_x,
                                    psi_bulk,dpsidx_bulk,b,dbdt,dbdx);

    outfile << interpolated_x[0] << " "
            << interpolated_x[1] << " "
            << pressure_gradient[0] << " "
            << pressure_gradient[1] << " "
            << b << " "
            << dbdx[0] << " "
            << dbdx[1] << " "
            << interpolated_n[0] << " "
            << interpolated_n[1] << std::endl;
   }
 }

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


protected:

 /// \short Function to compute the shape and test functions and to return
 /// the Jacobian of mapping between local and global (Eulerian)
 /// coordinates
 inline double shape_and_test(const Vector<double> &s, Shape &psi, Shape &test)
  const
  {
   //Find number of nodes
   unsigned n_node = nnode();

   //Get the shape functions
   shape(s,psi);

   //Set the test functions to be the same as the shape functions
   for(unsigned i=0;i<n_node;i++) {test[i] = psi[i];}

   //Return the value of the jacobian
   return J_eulerian(s);
  }


 /// \short Function to compute the shape and test functions and to return
 /// the Jacobian of mapping between local and global (Eulerian)
 /// coordinates
 inline double shape_and_test_at_knot(const unsigned &ipt,
                                      Shape &psi, Shape &test)
  const
  {
   //Find number of nodes
   unsigned n_node = nnode();

   //Get the shape functions
   shape_at_knot(ipt,psi);

   //Set the test functions to be the same as the shape functions
   for(unsigned i=0;i<n_node;i++) {test[i] = psi[i];}

   //Return the value of the jacobian
   return J_eulerian_at_knot(ipt);
  }


 /// Function to calculate the prescribed flux at a given spatial
 /// position
 void get_flux(const double& w, const Vector<double>& dwdx, 
               double& flux)
  {
   //If the function pointer is zero return zero
   if(Flux_fct_pt == 0)
    {
     //oomph_info << "obacht what is happening here!?"<<std::endl;
     flux=0.0;
    }
   //Otherwise call the function
   else
    {
     (*Flux_fct_pt)(w,dwdx,flux);
    }
  }

private:

 /// Flag to indicate whether we apply thin-film effects or not
 bool Apply_thin_film_effects;
 
 /// Thin film effect contribution when computing the effective film thickness
 double thin_film_effect(double Ca)
 {
  if(Apply_thin_film_effects)
   {
    double homotopy_parameter=0.0;
    thin_film_homotopy(homotopy_parameter);
    return homotopy_parameter*pow(Ca,2.0/3.0)/(0.76+2.16*pow(Ca,2.0/3.0));
   }
  else
   {
    return 0.0;
   }
 }

 /// \short Add the element's contribution to its residual vector.
 /// flag=1(or 0): do (or don't) compute the contribution to the
 /// Jacobian as well.
 void fill_in_generic_residual_contribution_poisson_flux(
  Vector<double> &residuals, DenseMatrix<double> &jacobian,
  const unsigned& flag);


 /// Function pointer to the (global) prescribed-flux function
 PrescribedFluxFctPt Flux_fct_pt;
   
 /// Function pointer to the (global) wall speed function
 WallSpeedFctPt Wall_speed_fct_pt;

 ///The spatial dimension of the problem
 unsigned Dim;

 ///The index at which the unknown is stored at the nodes
 unsigned U_index_poisson;
   
 /// \short Pointer to the aspect_ratio ratio: reference gap width / 
 /// in-plane lengthscale
 double* Aspect_ratio_pt;
   
 /// Pointer to the Capillary number
 double *Ca_pt;

 /// Pointer to a homotopy parameter to control thin film effects
 double* Thin_film_homotopy_pt;


};

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////



//===========================================================================
/// Constructor, takes the pointer to the "bulk" element, the
/// index of the fixed local coordinate and its value represented
/// by an integer (+/- 1), indicating that the face is located
/// at the max. or min. value of the "fixed" local coordinate
/// in the bulk element.
//===========================================================================
template<class ELEMENT>
HeleShawFluxElement<ELEMENT>::
HeleShawFluxElement(FiniteElement* const &bulk_el_pt,
                   const int &face_index) :
  FaceGeometry<ELEMENT>(), FaceElement()
  {
#ifdef PARANOID
   {
    //Check that the element is not a refineable 3d element
    ELEMENT* elem_pt = new ELEMENT;
    //If it's three-d
    if(elem_pt->dim()==3)
     {
      //Is it refineable
      if(dynamic_cast<RefineableElement*>(elem_pt))
       {
        //Issue a warning
        OomphLibWarning(
         "This flux element will not work correctly if nodes are hanging\n",
         "HeleShawFluxElement::Constructor",
         OOMPH_EXCEPTION_LOCATION);
       }
     }
   }
#endif

   // Let the bulk element build the FaceElement, i.e. setup the pointers
   // to its nodes (by referring to the appropriate nodes in the bulk
   // element), etc.
   bulk_el_pt->build_face_element(face_index,this);

   // Initialise the prescribed-flux function pointer to zero
   Flux_fct_pt = 0;

   // Initialise the wall speed function pointer to zero
   Wall_speed_fct_pt = 0;

   // Extract the dimension of the problem from the dimension of
   // the first node
   Dim = this->node_pt(0)->ndim();

   //Set up U_index_poisson. Initialise to zero, which probably won't change
   //in most cases, oh well, the price we pay for generality
   U_index_poisson = 0;

   //Cast to the appropriate PoissonEquation so that we can
   //find the index at which the variable is stored
   //We assume that the dimension of the full problem is the same
   //as the dimension of the node, if this is not the case you will have
   //to write custom elements, sorry
   switch(Dim)
    {
//     //One dimensional problem
//    case 1:
//    {
//     PoissonEquations<1>* eqn_pt =
//      dynamic_cast<PoissonEquations<1>*>(bulk_el_pt);
//     //If the cast has failed die
//     if(eqn_pt==0)
//      {
//       std::string error_string =
//        "Bulk element must inherit from PoissonEquations.";
//       error_string +=
//        "Nodes are one dimensional, but cannot cast the bulk element to\n";
//       error_string += "PoissonEquations<1>\n.";
//       error_string +=
//        "If you desire this functionality, you must implement it yourself\n";
//
//       throw OomphLibError(error_string,
//                           "HeleShawFluxElement::HeleShawFluxElement()",
//                           OOMPH_EXCEPTION_LOCATION);
//      }
//     //Otherwise read out the value
//     else
//      {
//       //Read the index from the (cast) bulk element
//       U_index_poisson = eqn_pt->u_index_poisson();
//      }
//    }
//    break;

    //Two dimensional problem
    case 2:
    {
     //Read the index from the (cast) bulk element.
     U_index_poisson = 
      dynamic_cast<ELEMENT*>(this->bulk_element_pt())->p_index_hele_shaw();
    }
    break;

//    //Three dimensional problem
//    case 3:
//    {
//     PoissonEquations<3>* eqn_pt =
//      dynamic_cast<PoissonEquations<3>*>(bulk_el_pt);
//     //If the cast has failed die
//     if(eqn_pt==0)
//      {
//       std::string error_string =
//        "Bulk element must inherit from PoissonEquations.";
//       error_string +=
//        "Nodes are three dimensional, but cannot cast the bulk element to\n";
//       error_string += "PoissonEquations<3>\n.";
//       error_string +=
//        "If you desire this functionality, you must implement it yourself\n";
//
//       throw OomphLibError(error_string,
//                        "HeleShawFluxElement::HeleShawFluxElement()",
//                        OOMPH_EXCEPTION_LOCATION);
//
//      }
//     else
//      {
//       //Read the index from the (cast) bulk element.
//       U_index_poisson = eqn_pt->u_index_poisson();
//      }
//    }
//    break;

    //Any other case is an error
    default:
     std::ostringstream error_stream;
     error_stream <<  "Dimension of node is " << Dim
                  << ". It should be 1,2, or 3!" << std::endl;

     throw OomphLibError(error_stream.str(),
                         "HeleShawFluxElement::HeleShawFluxElement()",
                         OOMPH_EXCEPTION_LOCATION);
     break;
    }

   /// Initialise the bubble tail node pointers to NULL
   //Top_tail_node_pt = 0;
   //Bottom_tail_node_pt = 0;

   Ca_pt = 0;

   // By default we don't apply any thin film effects
   Apply_thin_film_effects=false;   /////Joao====>>>> Before false
  }


//===========================================================================
/// Compute the element's residual vector and the (zero) Jacobian matrix.
//===========================================================================
template<class ELEMENT>
void HeleShawFluxElement<ELEMENT>::
fill_in_generic_residual_contribution_poisson_flux(
 Vector<double> &residuals, DenseMatrix<double> &jacobian,
 const unsigned& flag)
{
 //Find out how many nodes there are
 const unsigned n_node = nnode();

 //Set up memory for the shape and test functions
 Shape psif(n_node), testf(n_node);

 //Set the value of Nintpt
 const unsigned n_intpt = integral_pt()->nweight();

 //Set the Vector to hold local coordinates
 Vector<double> s(Dim-1);

 //Integers to hold the local equation and unknown numbers
 int local_eqn=0;

 // Locally cache the index at which the variable is stored
 const unsigned u_index_poisson = U_index_poisson;
     
 // Get pointer to bulk element
 ELEMENT* bulk_el_pt = dynamic_cast<ELEMENT*>(this->bulk_element_pt());

 //Loop over the integration points
 //--------------------------------
 for(unsigned ipt=0;ipt<n_intpt;ipt++)
  {

   //Assign values of s
   for(unsigned i=0;i<(Dim-1);i++) {s[i] = integral_pt()->knot(ipt,i);}

   //Get the integral weight
   double w = integral_pt()->weight(ipt);

   //Find the shape and test functions and return the Jacobian
   //of the mapping
   double J = shape_and_test(s,psif,testf);

   //Premultiply the weights and the Jacobian
   double W = w*J;
   
   Vector<double> interpolated_x(2,0.0);
       
   // Loop over the shape functions
   for(unsigned l=0; l<n_node; l++)
    {
     // Loop over directional components
     for(unsigned i=0; i<2; i++)
      {         
       // Spatial bits
       interpolated_x[i] += this->nodal_position(l,i)*psif(l);
      }
    }
   
   // Now calculate the unit normal vector
   Vector<double> interpolated_n(2);
   outer_unit_normal(ipt,interpolated_n);
   
   // Get b from bulk!
   Vector<double> s_bulk(2);
   get_local_coordinate_in_bulk(s,s_bulk);
   double b = 0.0;
   double dbdt = 0.0;
   Vector<double> dbdx(2,0.0);
   
   //---------------------------------------------------------
   // NOTE: This is dangerous/wrong. We're really supposed to
   // pass in dpsi/dx (2D!) which we don't have/need here.
   // However it's only used for the computation of the ALE bits
   // in dh/dt and dh/dt isn't used here! Phew...
   //---------------------------------------------------------
   // Number of nodes in bulk element
   unsigned n_nod_bulk=bulk_el_pt->nnode();
   DShape dpsidx_bulk(n_nod_bulk,2);
   Shape psi_bulk(n_nod_bulk);
   bulk_el_pt->dshape_eulerian(s_bulk,psi_bulk,dpsidx_bulk);
   bulk_el_pt->get_upper_wall_data(s_bulk,interpolated_x,
                                   psi_bulk,dpsidx_bulk,b,dbdt,dbdx);

   //Get the imposed flux
   double flux=0.0;
   get_flux(b,dbdx,flux);

   //Now add to the appropriate equations

   //Loop over the test functions
   for(unsigned l=0;l<n_node;l++)
    {
     local_eqn = nodal_local_eqn(l,u_index_poisson);
     /*IF it's not a boundary condition*/
     if(local_eqn >= 0)
      {
       // normal flux, Q = Q_x n_x + Q_y n_y
       residuals[local_eqn] -= flux*testf[l]*W;

       // Imposed traction doesn't depend upon the solution,
       // --> the Jacobian is always zero, so no Jacobian
       // terms are required
      }
    }
  }
}


} // end of namespace

#endif

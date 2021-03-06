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
#ifndef OOMPH_HELE_SHAW_FOEPPL_VON_KARMAN_ELEMENTS_HEADER
#define OOMPH_HELE_SHAW_FOEPPL_VON_KARMAN_ELEMENTS_HEADER


// Config header generated by autoconfig
#ifdef HAVE_CONFIG_H
#include <oomph-lib-config.h>
#endif

//#include "foeppl_von_karman.h"
#include "nico_elastic_hele_shaw.h"
#include "Thele_shaw_elements.h"

namespace oomph
{


namespace Junk
{ 
 bool Do_junk_string=false;
 std::string Junk_string;
}


//==============================================================
/// Combined triangular Hele Shaw and Foeppl von Karman elements
/// for elastic-walled Hele Shaw problem
//==============================================================
template<unsigned NNODE_1D>
class THeleShawFoepplvonKarmanDisplacementElement
 : public virtual TFoepplvonKarmanDisplacementElement<NNODE_1D>,
   public virtual THeleShawElement<NNODE_1D>
{

public:

 /// Constructor
 THeleShawFoepplvonKarmanDisplacementElement() : 
 TFoepplvonKarmanDisplacementElement<NNODE_1D>(),
  THeleShawElement<NNODE_1D>()
   {
    Q_pt = &Default_Physical_Constant_Value;
    Aspect_ratio_pt=0;
    Ca_pt=0;
    Apply_thin_film_effects=false;     ///////Joao =====>>>>  Before false
    Thin_film_homotopy_pt = 0;
   }
  

  /// Access to FSI parameter Q
  const double &q() const {return *Q_pt;}
  
  /// Pointer to FSI parameter Q
  double* &q_pt() {return Q_pt;}
  

  /// \short Overloaded function to return the "bubble volume".
  /// Overloads function in FvK elements and incorporates aspect
  /// ratio and offset.
  virtual double get_bounded_volume() const
  {

   //Number of nodes and integration points for the current element
   const unsigned n_node = this->nnode();
   const unsigned n_intpt = this->integral_pt()->nweight();
   
   //Shape functions and their derivatives
   Shape psi(n_node);
   DShape dpsidx(n_node,2);
   
   //The nodal index at which the displacement is stored
   const unsigned w_nodal_index = this->nodal_index_fvk();
   
   // obacht
   double Ca_local = 0.0;
   if(Apply_thin_film_effects)
    {
     double Bubble_tip_velocity = *Bubble_tip_velocity_pt;
     Ca_local = Bubble_tip_velocity*ca();
    }

   //Initalise the integral
   double integral = 0;
   
   //Loop over the integration points
   for(unsigned ipt=0;ipt<n_intpt;ipt++)
    {
     //Get the integral weight
     double w = this->integral_pt()->weight(ipt);
     
     //Get determinant of the Jacobian of the mapping
     double J = this->dshape_eulerian_at_knot(ipt,psi,dpsidx);
     
     //Premultiply the weight and Jacobian
     double W = w*J;
     
     //Initialise storage for the w value and nodal value
     double interpolated_w = 0;
     double w_nodal_value;
     
     //Loop over the shape functions/nodes
     for(unsigned l=0;l<n_node;l++)
      {
       //Get the current nodal value
       w_nodal_value = this->raw_nodal_value(l,w_nodal_index);

       //Add the contribution to interpolated w
       interpolated_w += w_nodal_value*psi(l);
      }
     
     //Add the contribution from the current integration point
     integral += (1.0+interpolated_w/aspect_ratio())*
      (1.0-this->thin_film_effect(Ca_local))*W;
    }
   return integral;
  }

  /// \short Output all quantities at Gauss points
  void output_at_gauss_points(std::ostream &outfile) const
   {
    // Find out how many nodes there are
    const unsigned n_node = this->nnode();
    
    // Set up memory for the shape and test functions
    Shape psi(n_node), test(n_node);
    DShape dpsidx(n_node,2), dtestdx(n_node,2);
    
    // Index at which the hele_shaw unknown is stored
    //const unsigned p_nodal_index = p_index_hele_shaw();

    // Vector of local coordinates and velocity
    Vector<double> s(2);
    Vector<double> velocity(2);
    
    // Loop over integration points
    const unsigned n_intpt = this->integral_pt()->nweight();

    //Loop over the integration points
    for(unsigned ipt=0;ipt<n_intpt;ipt++)
     {

      // Get local coords of int point
      Vector<double> s(2);
      for(unsigned i=0;i<2;i++)
       {
        s[i] = this->integral_pt()->knot(ipt,i);
       }
      
      // Get velocity
      this->get_velocity(s,velocity);

      // Get coords
      Vector<double> x(2);
      for(unsigned i=0;i<2;i++)
       {
        x[i]=this->interpolated_x(s,i);
       }

      // Get pressure
      double pressure=0.0;
      double dummy=0.0;
      this->get_pressure_fvk(ipt,x,dummy,pressure);
      

      // Call the derivatives of the shape and test functions
      //double J = 
      this->dshape_and_dtest_eulerian_at_knot_hele_shaw(
       ipt,psi,dpsidx,
       test,dtestdx);
      
      
      // Calculate local values of unknown
      // Allocate and initialise to zero
      Vector<double> interpolated_x(2,0.0);
      Vector<double> interpolated_dpdx(2,0.0);
      
      // Loop over nodes
      for(unsigned l=0; l<n_node; l++)
       {
        // Loop over directions
        for(unsigned j=0; j<2; j++)
         {
          interpolated_x[j] += this->raw_nodal_position(l,j)*psi(l);
         }
       }
      

      // Get gap width and (Eulerian!) wall velocity
      // (Pass in psi and dpsidx to avoid their recomputation
      // in fsi context)
      double h = 1.0;
      double dhdt = 0.0;
      Vector<double> dhdx(2,0.0);
      Junk::Do_junk_string=true;
      get_upper_wall_data(s,interpolated_x,psi,dpsidx,h,dhdt,dhdx);
      Junk::Do_junk_string=false;

      // Output
      for(unsigned i=0;i<2;i++)
       {
        outfile << x[i] << " ";
       }
      outfile << velocity[0] << " "
              << velocity[1] << " "
              << pressure << " "
              << this->interpolated_w_fvk(s) << " "
              << h << " " 
              << dhdt << " " 
              << Junk::Junk_string
              << "\n";
     }
   }

  /// \short Overloaded version of function that determines the pressure
  /// on the elastic membrane at Gauss point ipt, with spatial
  /// coordinate x.
  inline virtual void get_pressure_fvk(const unsigned& ipt,
                                       const Vector<double>& x,
                                       const double& w,
                                       double& pressure) const
   {
    // Get contribution from Foeppl von Karmann -- essentially
    // the external pressure
    double pressure_fvk=0.0;
    double dummy=0.0;
    FoepplvonKarmanDisplacementEquations::get_pressure_fvk(ipt,x,dummy,
                                                           pressure_fvk);
    
    // Which nodal value represents the Hele Shaw pressure?
    unsigned p_nodal_index = this->p_index_hele_shaw();

    // Hele Shaw pressure
    double pressure_hs = 0.0;

    // We're in the bubble, so get the pressure from the external data
    // obacht assumes that first external data that is added is the bubble
    // pressure
    if(this->Hele_shaw_disabled)
     {
      if (this->nexternal_data()!=0)
       {
        pressure_hs = this->external_data_pt(0)->value(0);
       }
      else
       {
        pressure_hs=0.0;
       }
     }
    // We're outside the bubble, so get the pressure from Hele-Shaw
    else 
     {
      //Set up memory for the shape and test functions
      const unsigned n_node = this->nnode();
      Shape psi(n_node);      
      Vector<double> s(2);
      for(unsigned i=0;i<2;i++)
       {
        s[i] = this->integral_pt()->knot(ipt,i);
       }
      
      this->shape(s,psi);
      
      //Initialise storage for the w value and nodal value
      double p_nodal_value = 0;
      
      for(unsigned l=0;l<n_node;l++)
       {
        //Get the current nodal value
        p_nodal_value = this->raw_nodal_value(l,p_nodal_index);
        //Add the contribution to interpolated p
        pressure_hs += p_nodal_value*psi(l);
       }
     }
    
    // Return the sum of the single-physics FvK pressure and the multi-physics
    // Hele-Shaw pressure, weighted by the FSI parameter, q
    pressure = pressure_fvk + this->q()*pressure_hs;
   }

  /// \short Get wall shape of the deformable "upper wall" of the Hele Shaw
  /// cell from Foeppl von Karman.
  /// Shape function and its deriv. w.r.t. to global coordinate
  /// is passed in for efficiency to allow computation
  /// of mesh velocity without having to recompute it.
  inline virtual void get_upper_wall_data(const Vector<double>& s,
                                          const Vector<double>& x,
                                          const Shape& psi,
                                          const DShape& dpsidx,
                                          double& b,
                                          double& dbdt,
                                          Vector<double>& dbdx) const
   {
    //Find number of nodes
    const unsigned n_node = this->nnode();
    
    //Get the index at which the vertical displacement is stored
    const unsigned w_nodal_index = this->nodal_index_fvk();
    
    // Mesh velocity
    Vector<double> mesh_velocity(2,0.0);

    // Gradient of wall displacement
    Vector<double> dwdx(2,0.0);
    
    //this->get_gradient_of_deflection(s, dwdx);

    //Initialise value of w and dw/dt
    double w = 0.0;
    double w_prev=0.0;
    double dwdt = 0.0;
    
    //Loop over the local nodes and sum
    for(unsigned l=0;l<n_node;l++)
     {
      // Build up fvk wall displacement
      w += this->nodal_value(l,w_nodal_index)*psi[l];
      w_prev += this->nodal_value(1,l,w_nodal_index)*psi[l];

      // Get the data's timestepper
      TimeStepper* time_stepper_pt = this->node_pt(l)->time_stepper_pt();
      
      // Number of timsteps (past & present)
      const unsigned n_time = time_stepper_pt->ntstorage();
      
      //Add the contributions to the time derivative
      for(unsigned t=0;t<n_time;t++)
       {
        dwdt += time_stepper_pt->weight(1,t)*
         this->nodal_value(t,l,w_nodal_index)*psi[l];
       }
     }
    
    // Mesh velocity and gradient
    double ale_bit=0.0;
    for(unsigned l=0;l<n_node;l++)
     {
      for(unsigned j=0;j<2;j++)
       {
        mesh_velocity[j] += this->raw_dnodal_position_dt(l,j)*psi(l);
        dwdx[j] += this->nodal_value(l,w_nodal_index)*dpsidx(l,j);
       }
     }
    ale_bit=mesh_velocity[0]*dwdx[0]+mesh_velocity[1]*dwdx[1];
    

    if (Junk::Do_junk_string)
     {
      std::stringstream junk;
      junk << dwdt << " " << ale_bit << " " << dwdt-ale_bit  << " " 
           << mesh_velocity[0] << " " << dwdx[0] << " " 
           << mesh_velocity[1] << " " << dwdx[1] << " "
           << w_prev << " ";
      Junk::Junk_string=junk.str();
     }
      
    //oomph_info << "ALE bit: " << ale_bit << std::endl;

    // Subtract ALE bits...
    dwdt-=ale_bit;
    
    // Rescale on Hele - Shaw transverse non-dim
    b=1.0+w/aspect_ratio(); 
    dbdt=dwdt/aspect_ratio();
    dbdx[0] = dwdx[0]/aspect_ratio();
    dbdx[1] = dwdx[1]/aspect_ratio();
   }
  
  // Overload the pressure gradient function to make sure we don't
  // get spurious velocities inside the bubble
  void get_pressure_gradient(const Vector<double>& s,
                              Vector<double>& gradient) const
   {
    // We're in the bubble, the pressure is constant and the gradient zero    
    if(this->Hele_shaw_disabled)
     {
      gradient[0] = 0.0;
      gradient[1] = 0.0;
     }
    else
     {
      // Find out how many nodes there are in the element
      const unsigned n_node = this->nnode();
      
      // Get the index at which the unknown is stored
      const unsigned p_nodal_index = this->p_index_hele_shaw();;
      
      // Set up memory for the shape and test functions
      Shape psi(n_node);
      DShape dpsidx(n_node,2);
      
      // Call the derivatives of the shape and test functions
      this->dshape_eulerian(s,psi,dpsidx);
      
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
   }
  
  /// The current nondimensionalisation has velocity[i] = -h^2 *dp/dx_i
  void get_velocity(const Vector<double>& s, Vector<double>& velocity) const
  {
   
   /// To find the velocity, we multiply the pressure gradient by b^2. We need
   /// to interpolate to find x(s) in order to call b(x,t) via 
   /// get_upper_wall_data.
   
   // Find out how many nodes there are in the element
   const unsigned n_node = this->nnode();
   
   // Set up memory for the shape and test functions
   Shape psi(n_node);
   DShape dpsidx(n_node,2);
   
   // Call the derivatives of the shape and test functions
   this->dshape_eulerian(s,psi,dpsidx);
   
   // Calculate local values of unknown
   // Allocate and initialise to zero
   Vector<double> interpolated_x(2,0.0);
   Vector<double> pressure_gradient(2,0.0);
   
   // Initialise
   double dbdt = 0.0;
   double b = 0.0;
   Vector<double> dbdx(2,0.0);
   
   // Loop over nodes to assemble the coordinate
   for(unsigned l=0; l<n_node; l++)
    {
     // Loop over coordinate directions
     for(unsigned j=0; j<2; j++)
      {
       interpolated_x[j] += this->raw_nodal_position(l,j)*psi(l);
      }
    }
   
   // Now get the gap width; local coordinate of integration point
   // is only used in multiphysics
   get_upper_wall_data(s,interpolated_x,psi,dpsidx,b,dbdt,dbdx);
   this->get_pressure_gradient(s,pressure_gradient);
   
   /// Now assemble the velocity components.
   for(unsigned j=0; j<2; j++)
    {
     velocity[j] = -b*b*pressure_gradient[j];
    }
   
  }

  /// Aspect ratio: Reference gap width / in-plane lengthscale
  double aspect_ratio() const
  {
#ifdef PARANOID
   if (Aspect_ratio_pt==0)
    {
     throw OomphLibError(
      "Aspect_ratio_pt not set yet for elastic-walled HeleShaw elements",
      "THeleShawFoepplvonKarmanElement::aspect_ratio()",
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
    
  /// Thin film homotopy parameter
  double thin_film_homotopy() const
  {
   if (Thin_film_homotopy_pt==0)
    {
     return 1.0;
    }
   else
    {
     return *Thin_film_homotopy_pt;
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
  
  /// \short Number of nodal values required at node n: sum of Hele Shaw
  /// and Foeppl von Karman ones.
  unsigned required_nvalue(const unsigned &n) const
  {
   return THeleShawElement<NNODE_1D>::required_nvalue(n)
    + TFoepplvonKarmanDisplacementElement<NNODE_1D>::required_nvalue(n);
  }

  /// \short Overload function that returns index of nodal value that
  /// represents FvK wall displacement (if no argument or argument 
  /// is 0; otherwise get the auxiliary fields required in 
  /// FvK; check that element for details.
 virtual inline unsigned nodal_index_fvk(const unsigned& i=0) const {return i;}

 /// \short Overload function that returns index of nodal value that
 /// stores Hele Shaw pressure
 virtual inline unsigned p_index_hele_shaw() const {return 4;}



////////// Joaos alteration
//
//void fill_in_contribution_to_mass_matrix(Vector<double> &residuals,
//                                      DenseMatrix<double> &mass_matrix)
//     {
//      TFoepplvonKarmanDisplacementElement<NNODE_1D>::
//        fill_in_contribution_to_mass_matrix(residuals,mass_matrix);
//      THeleShawElement<NNODE_1D>::
//        fill_in_contribution_to_mass_matrix(residuals,mass_matrix);       
//     }                                 
////////// Joaos alteration



////////// Joaos alteration
  void fill_in_contribution_to_jacobian_and_mass_matrix(
     Vector<double> &residuals,
     DenseMatrix<double> &jacobian, DenseMatrix<double> &mass_matrix)
     {
     TFoepplvonKarmanDisplacementElement<NNODE_1D>::fill_in_contribution_to_jacobian(residuals,jacobian);
     THeleShawElement<NNODE_1D>::fill_in_contribution_to_jacobian(residuals,jacobian); 
     }                                 
////////// Joaos alteration


 /// Assemble residuals from contributions of underlying elements
 void fill_in_contribution_to_residuals(Vector<double> &residuals)
 {
  TFoepplvonKarmanDisplacementElement<NNODE_1D>::
   fill_in_contribution_to_residuals(residuals);
  THeleShawElement<NNODE_1D>::fill_in_contribution_to_residuals(residuals);
 }
 
 ///\short Compute the element's residual vector and the Jacobian matrix.
 /// Jacobian is computed by finite-differencing.
 void fill_in_contribution_to_jacobian(Vector<double> &residuals,
                                      DenseMatrix<double> &jacobian)
 {
 
  oomph_info << "hierher shouldn't really call this -- overloaded by\n"
             << "brutal global FD routine in wrapper class in driver code\n"
             << "tidy this at some point...\n";
  exit(0);
  
  // these are only the diagonal entries... 
  TFoepplvonKarmanDisplacementElement<NNODE_1D>::
   fill_in_contribution_to_jacobian(residuals,jacobian);
  THeleShawElement<NNODE_1D>::
   fill_in_contribution_to_jacobian(residuals,jacobian);
 }
 
 ///  Overload the standard output function with the broken default
 void output(std::ostream &outfile)
 {
  this->output(outfile,3);
 }
 

 /// Output
 void output(std::ostream &outfile,
             const unsigned &nplot)
 {
  // Vector of local coordinates and velocity
  Vector<double> s(2);
  Vector<double> velocity(2);
  
  // Tecplot header info
  outfile << this->tecplot_zone_string(nplot);
  
  // Loop over plot points
  unsigned num_plot_points = this->nplot_points(nplot);
  for (unsigned iplot=0; iplot<num_plot_points; iplot++)
   {    
    // Get local coordinates and velocity at plot point
    this->get_s_plot(iplot,nplot,s);
    this->get_velocity(s,velocity);
    
    // We're in the bubble, so get the pressure from the external data    
    double pressure_hs=0.0;
    if(this->Hele_shaw_disabled)
     {
      if (this->nexternal_data()!=0)
       {
        pressure_hs = this->external_data_pt(0)->value(0);
       }
      else
       {
        pressure_hs=0.0;
       }
     }
    else
     {
      pressure_hs=this->interpolated_p_hele_shaw(s);
     }
    
    // output
    Vector<double> xx(2);
    for(unsigned i=0; i<2; i++)
     {
      xx[i]=this->interpolated_x(s,i);
      outfile << xx[i] << " ";
     }
    outfile << velocity[0] << " "
            << velocity[1] << " "
            << pressure_hs << " "
            << this->interpolated_w_fvk(s) << " "; //"\n";
    
    
    // Get gap width and (Eulerian!) wall velocity
    // (Pass in psi and dpsidx to avoid their recomputation
    // in fsi context)
    unsigned n_dim=this->dim();
    unsigned n_node=this->nnode();
    Shape psi(n_node);
    DShape dpsidx(n_node,n_dim);
    this->dshape_eulerian(s,psi,dpsidx);
    double h = 1.0;
    double dhdt = 0.0;
    Vector<double> dhdx(2,0.0);
    get_upper_wall_data(s,xx,psi,dpsidx,h,dhdt,dhdx);
    outfile << h << " " << dhdt << " ";
    
    // Get in-plane stress
    DenseMatrix<double> sigma(2,2,0.0);
    DenseMatrix<double> strain(2,2,0.0);
    this->get_stress_and_strain_for_output(s,sigma,strain);
    
    outfile << sigma(0,0) << " " 
            << sigma(1,1) << " " 
            << sigma(0,1) << " "
            << strain(0,0) << " "
            << strain(1,1) << " "
            << strain(0,1) << " ";
    
    /* // Convert to polars */
    /* if ((xx[0]==0.0)&&(xx[1]==0.0)) */
    /*  { */
    /*   outfile << "0.0 0.0 0.0 "; */
    /*  } */
    /* else */
    /*  { */
    /*   double theta=atan2(xx[1],xx[0]); */
    /*   double sigma_rr= */
    /*    sigma_xx*cos(theta)*cos(theta)+ */
    /*    sigma_yy*sin(theta)*sin(theta)+ */
    /*    sigma_xy*sin(2.0*theta); */
      
    /*   double sigma_tt= */
    /*    sigma_xx*sin(theta)*sin(theta)+ */
    /*    sigma_yy*cos(theta)*cos(theta)- */
    /*    sigma_xy*sin(2.0*theta); */
      
    /*   double sigma_rt= */
    /*    (sigma_yy-sigma_xx)*sin(theta)*cos(theta)+ */
    /*    sigma_xy*cos(2.0*theta); */
      
    /*   outfile << sigma_rr << " "  */
    /*           << sigma_tt << " "  */
    /*           << sigma_rt << " ";  */
    /*  } */


    outfile << std::endl;
    
   }
  
  // Write tecplot footer (e.g. FE connectivity lists)
  this->write_tecplot_zone_footer(outfile,nplot);  
 }

 
 /// \short C-style output function: Broken default
 void output(FILE* file_pt)
 {FiniteElement::output(file_pt);}
 
 ///  \short C-style output function: Broken default
 void output(FILE* file_pt, const unsigned &n_plot)
 {FiniteElement::output(file_pt,n_plot);}
 
 /// \short Output function for an exact solution: Broken default
 void output_fct(std::ostream &outfile, const unsigned &Nplot,
                 FiniteElement::SteadyExactSolutionFctPt 
                 exact_soln_pt)
 {FiniteElement::output_fct(outfile,Nplot,exact_soln_pt);}
 
 
 /// \short Output function for a time-dependent exact solution:
 /// Broken default.
 void output_fct(std::ostream &outfile, const unsigned &Nplot,
                 const double& time,
                 FiniteElement::UnsteadyExactSolutionFctPt 
                 exact_soln_pt)
 {
  FiniteElement::
   output_fct(outfile,Nplot,time,exact_soln_pt);
 }
 
 /// \short Validate against exact solution at given time
 /// Solution is provided via function pointer.
 /// Plot at a given number of plot points and compute L2 error
 /// and L2 norm of velocity solution over element
 /// Call the broken default
 void compute_error(std::ostream &outfile,
                    FiniteElement::UnsteadyExactSolutionFctPt exact_soln_pt,
                    const double& time,
                    double& error, double& norm)
  {FiniteElement::compute_error(outfile,exact_soln_pt,
                                time,error,norm);}
 
 /// \short Validate against exact solution.
 /// Solution is provided via function pointer.
 /// Plot at a given number of plot points and compute L2 error
 /// and L2 norm of velocity solution over element
 /// Call the broken default
 void compute_error(std::ostream &outfile,
                    FiniteElement::SteadyExactSolutionFctPt exact_soln_pt,
                    double& error, double& norm)
 {FiniteElement::compute_error(outfile,exact_soln_pt,error,norm);}
 

 /// Call self tests of underlying elements
 unsigned self_test()
 {
   unsigned hs_test_result = THeleShawElement<NNODE_1D>::self_test();
   unsigned fvk_test_result = 
    TFoepplvonKarmanDisplacementElement<NNODE_1D>::self_test();

   if (hs_test_result == 0 && fvk_test_result == 0)
    {
     return 0;
    }
   else
    {
     return 1;
    }
  }


 /// \short Pin all Foeppl von Karmann dofs -- turns problem into
 /// single physics Hele Shaw problem
 void pin_fvk()
 {
  unsigned n_node = this->nnode();
  for(unsigned i_node = 0; i_node < n_node; i_node++)
   {
    Node *nod_pt = this->node_pt(i_node);
    for(unsigned index = 0;index < 4; index++)
     {
      nod_pt->pin(index);
     }
   }
 }


 /// \short Unpin all Foeppl von Karmann dofs -- turns problem into
 /// single physics Hele Shaw problem
 void unpin_fvk()
 {
  unsigned n_node = this->nnode();
  for(unsigned i_node = 0; i_node < n_node; i_node++)
   {
    Node *nod_pt = this->node_pt(i_node);
    for(unsigned index = 0;index < 4; index++)
     {
      nod_pt->unpin(index);
     }
   }
 }
 
 
 /// \short Pin all Hele Shaw dofs -- turns problem into single
 /// physics Foeppl von Karmann problem
 void pin_hs()
 {
  unsigned n_node = this->nnode();
  for(unsigned i_node = 0; i_node < n_node; i_node++)
   {
    Node *nod_pt = this->node_pt(i_node);
    unsigned n_value = nod_pt->nvalue();
    for(unsigned index = 4; index < n_value; index++)
     {
      nod_pt->pin(index);
     }
   }
 }
 
 
 /// \short Unpin all Hele Shaw dofs -- turns problem into single
 /// physics Foeppl von Karmann problem
 void unpin_hs()
 {
  unsigned n_node = this->nnode();
  for(unsigned i_node = 0; i_node < n_node; i_node++)
   {
    Node *nod_pt = this->node_pt(i_node);
    unsigned n_value = nod_pt->nvalue();
    for(unsigned index = 4; index < n_value; index++)
     {
      nod_pt->unpin(index);
     }
   }
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
    
 /// The value of the Capillary number
 const double& ca() const 
 {
  return *Ca_pt;
 }
 
 /// Pointer to the inverse Capillary number
 double*& ca_pt() {return Ca_pt;}
 
  protected:
 
 /// Pointer to FSI parameter Q
 double *Q_pt;
 
 /// \short Order of recovery shape functions for Z2 error estimation:
 /// Same order as shape functions.
 unsigned nrecovery_order() {return (NNODE_1D-1);}
 
 /// \short Number of 'flux' terms for Z2 error estimation -- currently
 /// based on Hele Shaw
 unsigned num_Z2_flux_terms()
 {
  return THeleShawElement<NNODE_1D>::num_Z2_flux_terms();
 }
 
 /// Get 'flux' for Z2 error recovery -- currently based on Hele Shaw
 void get_Z2_flux(const Vector<double>& s, Vector<double>& flux)
 {
  Vector<double> fluid_contribution(2,0.0);
  THeleShawElement<NNODE_1D>::get_Z2_flux(s, fluid_contribution);
  // obacht
  Vector<double> sheet_contribution(2,0.0);
  TFoepplvonKarmanDisplacementElement<NNODE_1D>::get_Z2_flux(
   s, sheet_contribution);
  for(unsigned i =0;i<2;i++)
   {
    // obacht
    flux[i] = fluid_contribution[i] + sheet_contribution[i];
   }
 }
 
 /// \short Number of vertex nodes in the element
 unsigned nvertex_node() const
 {return TElement<2,NNODE_1D>::nvertex_node();}
 
 /// \short Pointer to the j-th vertex node in the element
 Node* vertex_node_pt(const unsigned& j) const
  {return TElement<2,NNODE_1D>::vertex_node_pt(j);}
 
private:


 /// \short Pointer to the aspect_ratio ratio: reference gap width / 
 /// in-plane lengthscale
 double* Aspect_ratio_pt;
 
 /// Default value for physical constants
 static double Default_Physical_Constant_Value;

 /// Flag to indicate whether we apply thin-film effects or not
 bool Apply_thin_film_effects;
   
 /// Pointer to the Capillary number
 double *Ca_pt;
 
 /// Thin film effect contribution when computing the effective film thickness
 double thin_film_effect(const double& Ca) const
 {
  if(Apply_thin_film_effects)
   {
    return thin_film_homotopy()*pow(Ca,2.0/3.0)/(0.76+2.16*pow(Ca,2.0/3.0));
   }
  else
   {
    return 0.0;
   }
 }

 // Pointer to the bubble tip node
 Node* Bubble_tip_node_pt;

 // Pointer to the bubble tip velocity
 double* Bubble_tip_velocity_pt;

 /// Pointer to a homotopy parameter to control thin film effects
 double* Thin_film_homotopy_pt;

};



///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////


//==========================================================
/// Instantiation of default value for FSI parameter
//==========================================================
template<unsigned NNODE_1D>
double THeleShawFoepplvonKarmanDisplacementElement<NNODE_1D>::
Default_Physical_Constant_Value = 0.0;

}

#endif

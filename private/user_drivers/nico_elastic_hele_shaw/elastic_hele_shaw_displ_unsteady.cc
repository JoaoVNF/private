//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC//    
//LIC//           Version 0.85. June 9, 2008.
//LIC// 
//LIC// Copyright (C) 2006-2008 Matthias Heil and Andrew Hazel
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
#include <fenv.h> 

//Generic routines 
#include "generic.h"
#include "solid.h" 

#ifdef USE_USR_SRC_LIBRARY

// The equations
#include "hele_shaw_and_fvk.h"

#else

#include "hs_displ_fvk_elements.h" 
#include "hele_shaw_interface_elements.h"
#include "hele_shaw_flux_elements.h"

#endif

// The mesh
#include "meshes/triangle_mesh.h"

namespace oomph
{

//=================================================================
/// Custom element inheriting from wrapped pseudo-solid version 
/// of Hele-Shaw and Foeppl von Karman multiphysics elements
/// so that the jacobian and residual calculations can be fine-tuned
//=================================================================
template <unsigned NNODE_1D>
class MyElement : public PseudoSolidNodeUpdateElement<
                  THeleShawFoepplvonKarmanDisplacementElement<NNODE_1D>,
                  TPVDElement<2, NNODE_1D> >
{

public:
 
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
   this->fill_in_jacobian_from_external_by_fd(full_residuals,jacobian,true);

   //There could also be nodal data
   this->fill_in_jacobian_from_nodal_by_fd(full_residuals,jacobian);
  }

};

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=================================================================
/// MyElement<3> upgraded to be projectable for mesh refinement
//=================================================================
 class ProjectableMyElement :
  public virtual ProjectableElement<MyElement<3> >
 {
  
 public:
  
  /// \short Specify the values associated with field fld.
  /// The information is returned in a vector of pairs which comprise
  /// the Data object and the value within it, that correspond to field fld.
  Vector<std::pair<Data*,unsigned> > data_values_of_field(const unsigned& fld)
   {
    
#ifdef PARANOID
    if (fld > 4)
     {
      std::stringstream error_stream;
      error_stream
       << "MyElement elements only store nine fields so fld must be"
       << "0 to 4 rather than " << fld << std::endl;
      throw OomphLibError(
       error_stream.str(),
       "ProjectableMyElement::data_values_of_field()",
       OOMPH_EXCEPTION_LOCATION);
     }
#endif
    
    // Create the vector
    unsigned nnod = this->nnode();
    Vector<std::pair<Data*,unsigned> > data_values(nnod);
    
    // Loop over all nodes
    for (unsigned j = 0; j < nnod; j++)
     {
      // Add the data value associated field: The node itself
      data_values[j] = std::make_pair(this->node_pt(j), fld);
     }
    
    // Return the vector
    return data_values;
   }
  
  /// \short Number of fields to be projected
  unsigned nfields_for_projection()
   {
    return 5;
   }
  
  /// \short Number of history values to be stored for fld-th field.
  unsigned nhistory_values_for_projection(const unsigned &fld)
   {
#ifdef PARANOID
    if (fld > 4)
     {
      std::stringstream error_stream;
      error_stream
       << "MyElement elements only store nine fields so fld must be"
       << "0 to 4 rather than " << fld << std::endl;
      throw OomphLibError(
       error_stream.str(),
       "ProjectableMyElement::nhistory_values_for_projection()",
       OOMPH_EXCEPTION_LOCATION);
     }
#endif
    return this->node_pt(0)->ntstorage();
   }
  
  ///\short Number of positional history values
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
    if (fld > 4)
     {
      std::stringstream error_stream;
      error_stream
       << "MyElement elements only store nine fields so fld must be"
       << "0 to 4 rather than " << fld << std::endl;
      throw OomphLibError(
       error_stream.str(),
       "ProjectableMyElement::jacobian_and_shape_of_field()",
       OOMPH_EXCEPTION_LOCATION);
     }
#endif
    unsigned n_dim = this->dim();
    unsigned n_node = this->nnode();
    Shape test(n_node);
    DShape dpsidx(n_node, n_dim), dtestdx(n_node, n_dim);
    double J = this->dshape_and_dtest_eulerian_fvk(s, psi, dpsidx,
                                                   test, dtestdx);
    return J;
   }
  
  
  
  /// \short Return interpolated field fld at local coordinate s, at time level
  /// t (t=0: present; t>0: history values)
  double get_field(const unsigned &t,
                   const unsigned &fld,
                   const Vector<double>& s)
   {
#ifdef PARANOID
    if (fld > 4)
     {
      std::stringstream error_stream;
      error_stream
       << "MyElement elements only store nine fields so fld must be"
       << "0 to 4 rather than " << fld << std::endl;
      throw OomphLibError(
       error_stream.str(),
       "ProjectableMyElement::get_field()",
       OOMPH_EXCEPTION_LOCATION);
     }
#endif
    
    //Find the index at which the variable is stored
    unsigned w_nodal_index = this->nodal_index_fvk(fld);
    
   //Local shape function
    unsigned n_node = this->nnode();
    Shape psi(n_node);
    
    //Find values of shape function
    this->shape(s,psi);
    
    //Initialise value of u
    double interpolated_w = 0.0;
    
    //Sum over the local nodes
    for(unsigned l = 0; l < n_node; l++)
     {
      interpolated_w += this->nodal_value(t,l,w_nodal_index)*psi[l];
     }
    return interpolated_w;
   }
  
  
  
  
  ///Return number of values in field fld: One per node
  unsigned nvalue_of_field(const unsigned &fld)
   {
#ifdef PARANOID
    if (fld > 4)
     {
      std::stringstream error_stream;
      error_stream
       << "MyElement elements only store nine fields so fld must be"
       << "0 to 4 rather than " << fld << std::endl;
      throw OomphLibError(
       error_stream.str(),
       "ProjectableMyElement::nvalue_of_field()",
       OOMPH_EXCEPTION_LOCATION);
     }
#endif
    return this->nnode();
   }
  
  
  ///Return local equation number of field fld of node j.
  int local_equation(const unsigned &fld,
                     const unsigned &j)
   {
#ifdef PARANOID
    if (fld > 4)
     {
      std::stringstream error_stream;
      error_stream
       << "MyElement elements only store nine fields so fld must be"
       << "0 to 4 rather than " << fld << std::endl;
      throw OomphLibError(
       error_stream.str(),
       "ProjectableMyElement::local_equation()",
       OOMPH_EXCEPTION_LOCATION);
     }
#endif
    const unsigned w_nodal_index = this->nodal_index_fvk(fld);
    return this->nodal_local_eqn(j, w_nodal_index);
   }
  
 };


///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////


//=======================================================================
/// Face geometry for element is the same as that for the underlying
/// wrapped element
//=======================================================================
template <>
class FaceGeometry<MyElement<3> >
  : public virtual FaceGeometry<SolidTElement<2,3> >
 {
   public:
  FaceGeometry() : FaceGeometry<SolidTElement<2,3> >() {}
 };


//=======================================================================
/// Face geometry for element is the same as that for the underlying
/// wrapped element
//=======================================================================
template <>
class FaceGeometry<ProjectableMyElement >
  : public virtual FaceGeometry<SolidTElement<2,3> >
 {
   public:
  FaceGeometry() : FaceGeometry<SolidTElement<2,3> >() {}
 };
} // end of oomph namespace



//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

/// Start of the driver code itself

using namespace std;
using namespace oomph;
using MathematicalConstants::Pi;


//=======================================================================
/// Namespace for problem parameters
//=======================================================================
namespace Problem_Parameter
{

 /// Adaptation frequency (negative: never; 1: always)
 int Adaptation_frequency=1;

 /// Initial/max element area
 double Max_el_area=0.05; // 0.005;

 /// Min element area for refinement
 double Min_el_area=0.005; //0.01; // 0.0001;

 /// Max permitted error (spatial)
 double Max_permitted_error=0.0005;

 /// Min permitted error (spatial)
 double Min_permitted_error=0.0001;

 /// Number of polyline segments on bubble. Negative: infer from
 /// max element area
 int N_segment_bubble=50;

 /// Max length of polyline segments on bubble -- negative: don't limit.
 double Max_polyline_length_bubble=0.03;

 /// Output directory
 std::string Directory = "RESLT";

 /// Name of restart file
 std::string Restart_file="";

 /// Hacky: Pointer to problem to get acccess to global time
 Problem* Problem_pt = 0;

 ///Physical parameters in the problem

 ///Sheet thickness (in m)
 double h_sheet = 0.34*1.0e-3; //3*1e-5;

 ///Youngs modulus of the sheet (in Pa)
 double E_sheet = 1.44*1.0e6; //3.71*1e9;

 ///Poissons ratio of the sheet 
 double nu_sheet = 0.5; //0.44;

 ///Bending stiffness of the sheet
 double K_sheet = E_sheet*h_sheet*h_sheet*h_sheet/
  (12.0*(1.0-nu_sheet*nu_sheet));

 ///Flow rate (in m^3/s)
 double V_rate = 189.0*1.0e-6/60.0; //45 //550*1e-6/60;

 ///Initial cell depth (in m)
 double b0_cell = 1.05e-3; 

 ///The cell width (in m)
 double W_cell = 30.0*1.0e-3; 

 ///Fluid viscosity (Pa.s)
 double mu_fluid = 0.099;

 ///Surface tension (N/m)
 double gamma_fluid = 21.0*1.0e-3;

 // Hele Shaw parameters
 //---------------------
 /// \short Storage for pointer to Data item whose single value stores
 /// pressure in the bubble
 Data *P_bubble_pt = 0;

 /// \short The volume of the bubble region under the membrane (dependent
 /// variable; needs to be updated before each timestep.
 double Prescribed_volume = 0.0;

 /// Non-dimensional gap width 
 double Aspect_ratio = b0_cell/W_cell;

 /// Time-scale used in non-dimnesionalisation
 double Time_scale = W_cell*W_cell*b0_cell/V_rate;

 /// Capillary number 
 double Ca = mu_fluid*W_cell/(gamma_fluid*Time_scale);

 ///Inverse Capillary number
 double Ca_inv = 1.0/Ca;

 /// FSI parameter
 double Target_q = 12.0*mu_fluid*V_rate/
  (Aspect_ratio*Aspect_ratio*Aspect_ratio*K_sheet);
 double Q=Target_q;
 
 /// \short Backup for bubble pressure (to retain value during spatial
 /// adaptation during which volume control element (which stores
 /// this value in its internal Data) is temporarily deleted.
 double P_bubble_backup=0.0;

 // Foeppl von Karman parameters
 //-----------------------------

 /// Foeppl von Karman parameter 
 double Eta = 12.0*(1.0-nu_sheet*nu_sheet)*(W_cell/h_sheet)*(W_cell/h_sheet); 

 /// Initial bubble width (in m)
 double R_bubble_initial_dim = 7.5*1.0e-3;       // joao==> used before 7.5*1.0e-3

 /// Initial bubble width
 double R_bubble_initial = R_bubble_initial_dim/W_cell; 
 
 /// Amplitude of radius perturbation of the bubble interface 
 double Radius_perturb = R_bubble_initial/10.0;           // joao==> used before  R_bubble_initial/10.0
 
 /// Wavenumber of pressure perturbation at interface
 unsigned N_perturb = 12; 

 /// Initial bubble length (in m)
 double L_bubble_initial_dim = 90.0*1.0e-3;

 /// Initial bubble length
 double L_bubble_initial = L_bubble_initial_dim/W_cell; 
 
 /// Pseudo-solid Poisson ratio
 double Nu = 0.3;
 
 /// Constitutive law used to determine the mesh deformation
 ConstitutiveLaw *Constitutive_law_pt = 0;

 /// Length of channel 
 double Length=20.0; //2.0;

 /// \short Function describing the prescribed motion of the upper wall. 
 /// Remember that h modulates the non-dimensional thickness specified 
 /// as the aspect ratio
 // void get_upper_wall_data(const Vector<double>& x,
 //                          double& h,
 //                          double& dhdt)
 // {
 //  // Get time and radius
 //  double t = Problem_pt->time_stepper_pt()->time();
 //  double r = sqrt(x[0]*x[0] + x[1]*x[1]);

 //  // Gap width changes linearly from conical into parabolic profile
 //  h = 1.0 + Prescribed_parabolic_gap_width_growth_factor*(1.0 - r*r)*t + 
 //   Prescribed_conical_gap_width_slope*r;
 //  dhdt = Prescribed_parabolic_gap_width_growth_factor*(1.0 - r*r);
 // }
  

 /// Amplitude of pressure perturbation at interface 
 double Pressure_perturb = 0.0;

 /// External pressure (in Pa)
 double Target_p_ext_dim = -81.0; //-96.0; //-153.2; //-324.5; //-67.0
 double Target_p_ext=W_cell*W_cell*W_cell/K_sheet*
  Target_p_ext_dim;
  //b0_cell*b0_cell*b0_cell/(12.0*V_rate*mu_fluid)*Target_p_ext_dim; 

 /// External pressure that is actually applied
 double P_ext=0.0;

 /// Assigns the value of pressure depending on the position
 void external_pressure(const Vector<double>& x, const double& w, 
                        double& pressure)
  {
   pressure=P_ext;
  }

 // Pre stress in the sheet (Pa)
 double Pre_stress = 30.0e3; //40.0e3;     //joao==> used before 160.0e3

 void get_pre_stress(DenseMatrix<double>& sigma_0)
 {
  sigma_0(0,0) = 0.0;
  sigma_0(1,1) = Pre_stress/E_sheet;
  sigma_0(0,1) = 0.0;
  sigma_0(1,0) = sigma_0(0,1);
 }
 
 double Thin_film_homotopy_parameter=1.0;

 // /// Thin film effect contribution when computing the effective film thickness
 double thin_film_effect()
 {
  if(CommandLineArgs::command_line_flag_has_been_set("--inc_thin_films"))
   {
    return pow(Ca,2.0/3.0)/(0.76+2.16*pow(Ca,2.0/3.0));
   }
  else
   {
    return 0.0;
   }
 }

 /// Update the dependent parameters
 void update_dependent_parameters()
 {
  Pressure_perturb=0.0;
   //Radius_perturb*Aspect_ratio*Aspect_ratio*Ca_inv/12.0*
   //float(N_perturb*N_perturb-1)/(R_bubble_initial*R_bubble_initial);
 }

 Node* Tip_node_pt = 0;
 double Bubble_tip_velocity = 0.0;
 Node* Top_node_pt = 0;
 Node* Bottom_node_pt = 0;

 /// Doc parameters
 void doc_parameters()
 {
  oomph_info 
   << "Initial radius: " 
   << Problem_Parameter::R_bubble_initial << std::endl
   << "Initial radius perturbation: " 
   << Problem_Parameter::Radius_perturb << std::endl
   << "Pressure perturbation: " 
   << Problem_Parameter::Pressure_perturb << std::endl
   << "Initial number of fingers: " 
   << Problem_Parameter::N_perturb << std::endl
   << "Aspect ratio: "
   << Problem_Parameter::Aspect_ratio << std::endl
   << "Inverse capillary number: "
   << Problem_Parameter::Ca_inv << std::endl
   << "Foeppl-von Karman coefficient, Eta: "
   << Problem_Parameter::Eta << std::endl
   << "FSI coefficient, Q: "
   << Problem_Parameter::Q << std::endl
   << "External pressure, P_ext: "
   << Problem_Parameter::P_ext << std::endl
   << "Initial/max element area: "
   << Problem_Parameter::Max_el_area << "\n"
   << "Min element area for refinement: "
   << Problem_Parameter::Min_el_area << "\n"
   << "Max length of polyline segments on bubble (neg: not constrained): "
   << Problem_Parameter::Max_polyline_length_bubble << "\n"
   << "Adaptation frequency (neg: never): "
   << Problem_Parameter::Adaptation_frequency << "\n"
   << std::endl;
 }


 /// \short Pressure acting on air-liquid interface -- basically the
 /// bubble pressure but here modified by sinusoidal perturbation
 /// to force fingering with a particular wavenumber.
 void bubble_pressure_function(const Vector<double>& x, double& pressure)
 {
  if (0 == P_bubble_pt)
   {
    pressure = 0.0;
   }
  else
   {
    pressure = P_bubble_pt->value(0);
    double phi = 2.0*MathematicalConstants::Pi*x[1]; //atan2(x[1], x[0]);
    pressure += Pressure_perturb*sin(double(N_perturb)*phi);
   }
 }


}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////


//==start_of_problem_class============================================
/// Class definition
//====================================================================
template<class ELEMENT>
class UnstructuredHSFvKProblem : public virtual Problem
{

public:

 /// Constructor
 UnstructuredHSFvKProblem();
    
 /// Destructor
 ~UnstructuredHSFvKProblem()
  {
   delete Bulk_mesh_pt;
  };

 /// \short hierher kill
 void actions_before_newton_step()
  {
  }


 /// \short Actions before conv check hierher
 void actions_before_newton_convergence_check()
  {
   if (CommandLineArgs::command_line_flag_has_been_set("--pin_fvk"))
    {
     // Update pinned fvk displacements (mesh moves!)
     update_pinned_fvk_displacements();
    }

   // doc_solution();
   // pause("done doc solution before newton conv check");
  }

 /// \short Actions before implicit timestep: Update the target volume
 void actions_before_implicit_timestep()
  {
   // Current time
   double t = time_stepper_pt()->time();

   // Intial volume for bubble
   // obacht need to add thin film effect here
   Problem_Parameter::Prescribed_volume = Pi*0.5*
    Problem_Parameter::R_bubble_initial*
    Problem_Parameter::R_bubble_initial+
    Problem_Parameter::L_bubble_initial*2.0*Problem_Parameter::R_bubble_initial;
   
   if (CommandLineArgs::command_line_flag_has_been_set("--inc_thin_films"))
    {
     Problem_Parameter::Prescribed_volume *= 
      (1.0-Problem_Parameter::thin_film_effect());
    }

   Problem_Parameter::Bubble_tip_velocity = 
    Problem_Parameter::Tip_node_pt->dposition_dt(0);

   // Increase
   Problem_Parameter::Prescribed_volume += t;

   if(t > 1.0e-4 && t < 0.025)
    {
     Problem_Parameter::update_dependent_parameters();
    }
   else
    {
     //Problem_Parameter::Radius_perturb=0.0;
     Problem_Parameter::Pressure_perturb=0.0;
    }

   if (CommandLineArgs::command_line_flag_has_been_set("--pin_fvk"))
    {
     // Set all Foeppl von Karmann displacements
     update_pinned_fvk_displacements();
    }
  }

 /// \short Actions before adapt. Delete old volume constraint element and
 /// update global mesh
 void actions_before_adapt()
  {


   // // hierher
   // {
   //  ofstream some_file;
   //  char filename[100];
   //  sprintf(filename,"%s/gauss_point_data_before_adapt%i.dat",
   //          Doc_info.directory().c_str(),
   //          Doc_info.number());
   //  some_file.open(filename);
   //  unsigned nel = Bulk_mesh_pt->nelement();
   //  for (unsigned e = 0; e < nel; e++)
   //   {
   //    dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e))->
   //     output_at_gauss_points(some_file);
   //   }
   //  some_file.close();

   //  sprintf(filename,"%s/soln_before_adapt%i.dat",Doc_info.directory().c_str(),
   //          Doc_info.number());
   //  some_file.open(filename);
   //  this->Bulk_mesh_pt->output(some_file,5); 
   //  some_file.close();
   // }


   // Make backup of surface mesh
   Backed_up_surface_mesh_pt=new BackupMeshForProjection<TElement<1,3> >
    (Free_surface_bc_mesh_pt,Inner_boundary);


   if (!CommandLineArgs::command_line_flag_has_been_set("--pin_hs"))
    {
     // Keep a temporary copy of the pressure from the old volume contraint
     // element
     Problem_Parameter::P_bubble_backup =
      Volume_constraint_element_pt->pressure_data_pt()->value(0);
     
     // Empty volume constraint mesh
     Volume_constraint_mesh_pt->flush_element_and_node_storage();
     
     // Delete the volume constraint element
     delete Volume_constraint_element_pt;
     
     // Kill the Hele Shaw free surface elements
     delete_surface_elements();

     //delete_cross_section_elements();
    }

   // The constraint element is gone and therefore also the pointer to
   // the bubble pressure (it's internal data). Therefore we need to
   // manually set the nodal pressures in the bubble region to the
   // backed up bubble pressure, otherwise we'll be in trouble when
   // adapting
   {
    unsigned n_bubble_regions = Bubble_regions.size();
    for(unsigned r = 0; r < n_bubble_regions; r++)
     {
      unsigned n_element = Bulk_mesh_pt->nregion_element(Bubble_regions[r]);
      for(unsigned e = 0; e < n_element; e++)
       {
        ELEMENT* el_pt = dynamic_cast<ELEMENT*>
         (Bulk_mesh_pt->region_element_pt(Bubble_regions[r], e));
        
        unsigned nnod = el_pt->nnode();
        for(unsigned j = 0; j < nnod; j++)
         {
          Node* nod_pt = el_pt->node_pt(j);
          if( (!(nod_pt->is_on_boundary(Inner_boundary)  )) )
           {
            // Pin HS pressure
            //nod_pt->pin(4);
            // obacht
            nod_pt->set_value(4,Problem_Parameter::P_bubble_backup);
           }
         }  
       }
     }
   }

   // Kill mesh as geom object
   delete Mesh_as_geom_object_pt;
   Mesh_as_geom_object_pt=0;

   // ...and rebuild the global mesh
   rebuild_global_mesh();
  }
 
 /// \short Actions after adapt: 
 /// Setup the problem again -- remember that the mesh has been
 /// completely rebuilt and its elements don't have any
 /// pointers to source fcts etc. yet
 /// Also create new volume constraint element with previous pressure value
 void actions_after_adapt()
  {   

   set_tip_node_pt();

   if (!CommandLineArgs::command_line_flag_has_been_set("--pin_hs"))
    {
     // Create new volume constraint element -- recycle any previous
     // assignment of bubble pressure
     Volume_constraint_element_pt =
      new FoepplvonKarmanVolumeConstraintElement
      <ELEMENT, RefineableSolidTriangleMesh>(
       Bulk_mesh_pt,
       Bubble_regions,
       Problem_Parameter::P_bubble_backup);
     
     // Set pointer to prescribed volume
     Volume_constraint_element_pt->set_prescribed_volume(
      &Problem_Parameter::Prescribed_volume);
     
     // Add it to the mesh
     Volume_constraint_mesh_pt->add_element_pt(Volume_constraint_element_pt);
     
     // Set the pressure pointer to allow external access to the pressure
     Problem_Parameter::P_bubble_pt = 
      Volume_constraint_element_pt->pressure_data_pt();
     
     // Create the interface elements along the bubble boundary
     create_surface_elements();
     //create_cross_section_elements();
    }

   // Now project from backup of original contact mesh to new one
   Backed_up_surface_mesh_pt->project_onto_new_mesh(
    Free_surface_bc_mesh_pt);

   // Set the Lagrange multiplier value to zero
   unsigned n_element = Free_surface_bc_mesh_pt->nelement();
   for(unsigned e=0;e<n_element;e++)
    {
     HeleShawInterfaceElement<ELEMENT>* el_pt =
      dynamic_cast<HeleShawInterfaceElement<ELEMENT>*>(
       Free_surface_bc_mesh_pt->element_pt(e));
     el_pt->set_lagrange_multipliers_to_zero();
    }

   // Rebuild the mesh and re-apply pinning and boundary conditions
   rebuild_global_mesh();
   complete_problem_setup();


   double bubble_volume=Volume_constraint_element_pt->get_bounded_volume();

   oomph_info<<"Bubble volume: Actual: "<<bubble_volume<<", Prescribed: "
             <<Problem_Parameter::Prescribed_volume<<std::endl;

   // oomph_info << "obacht resetting prescribed to actual bubble volume\n";
   // Problem_Parameter::Prescribed_volume=bubble_volume;
   
   // // hierher
   // {
   //  ofstream some_file;
   //  char filename[100];
   //  sprintf(filename,"%s/gauss_point_data_after_adapt%i.dat",
   //          Doc_info.directory().c_str(),
   //          Doc_info.number());
   //  some_file.open(filename);
   //  unsigned nel = Bulk_mesh_pt->nelement();
   //  for (unsigned e = 0; e < nel; e++)
   //   {
   //    dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e))->
   //     output_at_gauss_points(some_file);
   //   }
   //  some_file.close();

   //  sprintf(filename,"%s/soln_after_adapt%i.dat",Doc_info.directory().c_str(),
   //          Doc_info.number());
   //  some_file.open(filename);
   //  this->Bulk_mesh_pt->output(some_file,5); 
   //  some_file.close();
   // }

   if (CommandLineArgs::command_line_flag_has_been_set("--pin_fvk"))
    {
     // Set all Foeppl von Karmann displacements
     update_pinned_fvk_displacements();
    }

   // Now reset Lagrangian coordinates of the pseudo-solid
   if (!CommandLineArgs::command_line_flag_has_been_set(
        "--suppress_reset_lagr_coords"))
    {
     Bulk_mesh_pt->set_lagrangian_nodal_coordinates();
    }

   // (Re-)create mesh as geom object representation of the mesh
   Mesh_as_geom_object_pt=new MeshAsGeomObject(Bulk_mesh_pt);
   
   // Kill backed up mesh
   delete Backed_up_surface_mesh_pt;
   Backed_up_surface_mesh_pt=0;
  }
 
 /// Update after solve (empty)
 void actions_after_newton_solve() {}

 /// Update the problem specs before solve (empty)
 void actions_before_newton_solve()
  {
  }

 /// Restart
 void restart()
  {
   // Pointer to restart file
   ifstream* restart_file_pt=0;
   
   // Open restart file from stem
   restart_file_pt=new ifstream(Problem_Parameter::Restart_file.c_str(),
                                ios_base::in);
   if (restart_file_pt!=0)
    {
     oomph_info << "Have opened "
                << Problem_Parameter::Restart_file.c_str() 
                << " for restart. " << std::endl;
    }
   else
    {
     std::ostringstream error_stream;
     error_stream
      << "ERROR while trying to open " 
      << Problem_Parameter::Restart_file.c_str()
      << " for restart." << std::endl;
     
     throw OomphLibError(
      error_stream.str(),
      "restart()",
      OOMPH_EXCEPTION_LOCATION);
    }
   
   
   // Read restart data:
   //-------------------
   if (restart_file_pt!=0)
    {  
     // Read line up to termination sign
     string input_string;
     getline(*restart_file_pt,input_string,'#');
     
     // Ignore rest of line
     restart_file_pt->ignore(80,'\n');
     
     // Doc number
     Doc_info.number()=unsigned(atoi(input_string.c_str()));
     
     // Refine the mesh and read in the generic problem data
     Problem::read(*restart_file_pt);
    }
   
   // Clean up
   delete restart_file_pt;
   
   // Bump step number back 
   Doc_info.number()++; 
  }

 /// Dump problem data to allow for later restart
 void dump_it(ofstream& dump_file)
  {
   // Write doc number
   dump_file << Doc_info.number() << " # current doc number" << std::endl;

   // Dump the refinement pattern and the generic problem data
   Problem::dump(dump_file);
  }

  
 /// Doc the solution
 void doc_solution();

 /// \short Global error norm for adaptive timestepping
 double global_temporal_error_norm()
  {
   double global_error = 0.0;
   
   // Base error estimate on motion of air-liquid interface
   //------------------------------------------------------
   if (Use_interface_motion_for_temporal_error)
    {
     unsigned count = 0;
     unsigned ibound = Inner_boundary;
     unsigned num_nod = Bulk_mesh_pt->nboundary_node(ibound);
     for (unsigned inod = 0; inod < num_nod; inod++)
      {
       // Get node
       Node* nod_pt = Bulk_mesh_pt->boundary_node_pt(ibound, inod);
       double errx = nod_pt->position_time_stepper_pt()->
        temporal_error_in_position(nod_pt,0);
       double erry = nod_pt->position_time_stepper_pt()->
        temporal_error_in_position(nod_pt,1);
       
       // Add the square of the individual error to the global error
       count++;
       global_error += errx*errx + erry*erry;
      }
     
     
     // Divide by the number of nodes
     global_error /= double(count);
     
     // Return square root...
     global_error=sqrt(global_error);
    }
   // Base error estimate on norm of fvk deflection
   //----------------------------------------------
   else
    {
     unsigned count=0;
     unsigned num_nod = Bulk_mesh_pt->nnode();
     for (unsigned inod = 0; inod < num_nod; inod++)
      {
       // Get node
       Node* nod_pt=Bulk_mesh_pt->node_pt(inod);
       
       // Error based on FvK displacement
       double err=nod_pt->time_stepper_pt()->temporal_error_in_value(nod_pt,0);
       
       // Add the square of the individual error to the global error
       count++;
       global_error += err*err;
      }

     // Divide by the number of nodes
     global_error /= double(count);
     
     // Return square root...
     global_error=sqrt(global_error);
    }
   return global_error;
  }

 /// Helper function to pin all Foeppl von Karman dofs (for Hele Shaw only)
 void pin_fvk()
  {
   // Loop over all elements
   unsigned nel = Bulk_mesh_pt->nelement();
   for (unsigned e = 0; e < nel; e++)
    {
     ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));
     el_pt->pin_fvk();
    }
  }

 /// Helper function to pin all Foeppl von Karman dofs (for Hele Shaw only)
 void pin_hs()
  {
   // Loop over all elements
   unsigned nel = Bulk_mesh_pt->nelement();
   for (unsigned e = 0; e < nel; e++)
    {
     ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));
     el_pt->pin_hs();
    }
  }

 /// Helper function to pin all Foeppl von Karman dofs (for Hele Shaw only)
 void unpin_hs()
  {
   // Loop over all elements
   unsigned nel = Bulk_mesh_pt->nelement();
   for (unsigned e = 0; e < nel; e++)
    {
     ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));
     el_pt->unpin_hs();
    }
  }

 /// Update pinned fvk displacements (mesh moves!)
 void update_pinned_fvk_displacements()
  {
   // hierher uncomment if we want to do validation again

   // oomph_info << "This is wrong: Need to set history value of displacement\n"
   //            << " to the displacement at the point where the node was\n"
   //            << " at the previous timestep\n";
   // //abort();

   // // Loop over current and previous timesteps
   // unsigned ntime=0; //time_stepper_pt()->ndt();
   // for (int it=ntime;it>=0;it--)
   //  {
   //   // Get time
   //   double t = time_pt()->time(unsigned(it));

   //   oomph_info << "Assigning pinned fvk displacement for t=" << t << std::endl;
     
   //   // Set all Foeppl von Karmann displacements
   //   unsigned n_node = Bulk_mesh_pt->nnode();
   //   for(unsigned i = 0; i < n_node; i++)
   //    {
   //     Node *nod_pt = Bulk_mesh_pt->node_pt(i);
   //     double x=nod_pt->x(0);
   //     double y=nod_pt->x(1);
   //     double r=sqrt(x*x+y*y);
   //     double w=Problem_Parameter::Pinned_fvk_constant+
   //      (1.0-r*r)*(1.0-r*r)*Problem_Parameter::Pinned_fvk_amplitude;
   //     if (Problem_Parameter::Pinned_fvk_ref_time>0.0)
   //      {
   //       w*=(t*t)/(Problem_Parameter::Pinned_fvk_ref_time*
   //                 Problem_Parameter::Pinned_fvk_ref_time);
   //      }
   //     nod_pt->set_value(it,0,w);
   //    }
   //  }
  }

 /// Apply boundary conditions
 void apply_boundary_conditions(const bool& pin_bubble_ends)
  {
   Problem_Parameter::Top_node_pt->pin(6);
   Problem_Parameter::Top_node_pt->set_value(6,0.0);
   Problem_Parameter::Bottom_node_pt->pin(6);
   Problem_Parameter::Bottom_node_pt->set_value(6,0.0);
   
   // Set the boundary conditions for problem: All nodes are
   // free by default -- just pin the ones that have Dirichlet conditions
   // here. 
   unsigned nbound = Bulk_mesh_pt->nboundary();
   for(unsigned ibound = 0; ibound < nbound; ibound++)
    {
     unsigned num_nod = Bulk_mesh_pt->nboundary_node(ibound);
     for (unsigned inod = 0; inod < num_nod; inod++)
      {
       // Get node
       Node* nod_pt = Bulk_mesh_pt->boundary_node_pt(ibound, inod);
       SolidNode* solid_node_pt = dynamic_cast<SolidNode*>(nod_pt);
       if(ibound == Bottom_boundary || 
          ibound == Outflow_boundary ||
          ibound == Top_boundary ||
          ibound == Top_inflow_boundary ||
          ibound == Actual_inflow_boundary ||
          ibound == Bottom_inflow_boundary )
        {
         // Loop over current and previous timesteps
         unsigned ntime=time_stepper_pt()->ndt();
         for (unsigned t=0;t<ntime;t++)
          {         
           // Pin FvK unknown values
           // Transverse Displacement
           if(ibound == Bottom_boundary || 
              ibound == Top_boundary)
            {
             nod_pt->pin(0);
             nod_pt->set_value(t,0,0.0);
            }
           else if(!CommandLineArgs::command_line_flag_has_been_set("--moving_frame_bc"))
            {
             nod_pt->pin(0);
             nod_pt->set_value(t,0,0.0);
            }
              
           // In-plane displacements
           nod_pt->pin(2);
           nod_pt->set_value(t,2,0.0);
           nod_pt->pin(3);
           nod_pt->set_value(t,3,0.0);
           
           // Pin HS pressure
           // if (ibound == Outflow_boundary && inod == 0)
           //  {
           //   nod_pt->pin(4);
           //   nod_pt->set_value(t,4,Problem_Parameter::Outlet_pressure);     
           //   oomph_info << "Pinning pressure at " 
           //              << nod_pt->x(0) << " " 
           //              << nod_pt->x(1) << " " 
           //              << std::endl;
           //  }
          }
         
         // Pin positions of nodes on the outer boundary
         solid_node_pt->pin_position(0);
         if ( ibound == Top_inflow_boundary ||
              ibound == Actual_inflow_boundary ||
              ibound == Bottom_inflow_boundary )
          {
           if  ( solid_node_pt->is_on_boundary(Top_inflow_boundary) &&
                 solid_node_pt->is_on_boundary(Actual_inflow_boundary) )
            {
             if (pin_bubble_ends)
              {
               solid_node_pt->pin_position(1);
              }
             else            
              {
               solid_node_pt->unpin_position(1);
              }
            }
           if  ( solid_node_pt->is_on_boundary(Bottom_inflow_boundary) &&
                 solid_node_pt->is_on_boundary(Actual_inflow_boundary) )
            {
             if (pin_bubble_ends)
              {
               solid_node_pt->pin_position(1);
              }
             else            
              {
               solid_node_pt->unpin_position(1);
              }
            }
          }
         else
          {
           solid_node_pt->pin_position(1);
          }
        }
       
      }   
    } // end loop over boundaries
  }

 /// \short Helper function to (re-)set boundary condition
 /// and complete the build of  all elements
 void complete_problem_setup();

private:

 /// Doc info object for labeling output
 DocInfo Doc_info;

 /// Use interface motion to estimate temporal error
 bool Use_interface_motion_for_temporal_error;

 /// Mesh as geom object representation (updated after every adaptation)
 MeshAsGeomObject* Mesh_as_geom_object_pt;

 /// Pointers to bulk mesh
 RefineableSolidTriangleMesh<ELEMENT>* Bulk_mesh_pt;

 /// Mesh for the volume constraint element
 Mesh* Volume_constraint_mesh_pt;

 /// Pointers to mesh of flux BC elements
 Mesh* Free_surface_bc_mesh_pt;

// Mesh* Cross_section_mesh_pt;

 /// Single volume constraint element instance
 FoepplvonKarmanVolumeConstraintElement<ELEMENT, RefineableSolidTriangleMesh>
  *Volume_constraint_element_pt;

 /// Vector of the regions comprising the bubble
 Vector<unsigned> Bubble_regions;

 /// \short Backup of Free_surface_bc_mesh_pt so the Lagrange multipliers
 /// and tangents can be projected across
 BackupMeshForProjection<TElement<1,3> >*  Backed_up_surface_mesh_pt;

 /// Trace file to document norm of solution
 ofstream Trace_file;

 /// Enum for boundary ids
 enum
 {
  Bottom_boundary=0,
  Outflow_boundary=1,
  Top_boundary=2,
  Top_inflow_boundary=3,
  Actual_inflow_boundary=4,
  Bottom_inflow_boundary=5,
  Inner_boundary = 6//,
  //Cross_section_boundary = 7
 };

 /// Snap nodes onto curvilinear boundary
 void snap_onto_boundary()
  {
   unsigned ibound = Inner_boundary;
   unsigned num_nod = Bulk_mesh_pt->nboundary_node(ibound);
   for (unsigned inod = 0; inod < num_nod; inod++)
    {
     // Get node
     Node* nod_pt = Bulk_mesh_pt->boundary_node_pt(ibound,inod);
     
     if(nod_pt->x(0) > Problem_Parameter::L_bubble_initial)
      {
       double phi = atan2(nod_pt->x(1), 
                          nod_pt->x(0)-Problem_Parameter::L_bubble_initial);
       nod_pt->x(0) = Problem_Parameter::R_bubble_initial*cos(phi)+
        Problem_Parameter::L_bubble_initial;
       nod_pt->x(1) = Problem_Parameter::R_bubble_initial*sin(phi);
      }
    }
   // Now reset Lagrangian coordinates
   Bulk_mesh_pt->set_lagrangian_nodal_coordinates();
  }

 /// \short Create surface elements
 void create_surface_elements()
  {
   // Loop over the appropriate boundaries
   unsigned b = Inner_boundary;
   // How many bulk fluid elements are adjacent to boundary b in region 0?
   unsigned n_element = Bulk_mesh_pt->nboundary_element_in_region(b, 0);
   for(unsigned e = 0; e < n_element; e++)
    {
     // Get pointer to the bulk fluid element that is 
     // adjacent to boundary b in region 0
     ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>(
      Bulk_mesh_pt->boundary_element_in_region_pt(b, 0, e));
     
     //Find the index of the face of element e along boundary b in region 0
     int face_index = Bulk_mesh_pt->
      face_index_at_boundary_in_region(b, 0, e);
     
     // Build it
     HeleShawInterfaceElement<ELEMENT>* surface_element_pt = new
      HeleShawInterfaceElement<ELEMENT>(bulk_elem_pt, face_index);
     
     // Set the inverse capillary number
     surface_element_pt->ca_inv_pt() = 
      &Problem_Parameter::Ca_inv;
     
     // Set the non-dim reference gap width
     surface_element_pt->aspect_ratio_pt() = 
      &Problem_Parameter::Aspect_ratio;
     
     surface_element_pt->bubble_pressure_fct_pt() = 
      &Problem_Parameter::bubble_pressure_function;
     
     // Add bubble pressure (computed from volume constraint)
     // as external data for these elements
     surface_element_pt->add_external_data(Problem_Parameter::P_bubble_pt);
     
     if (CommandLineArgs::command_line_flag_has_been_set("--inc_thin_films"))
      {
       surface_element_pt->enable_thin_film_effects();
       surface_element_pt->thin_film_homotopy_pt() = 
        &Problem_Parameter::Thin_film_homotopy_parameter;
      }
     
     // Add the prescribed-flux element to the surface mesh
     Free_surface_bc_mesh_pt->add_element_pt(surface_element_pt);
     
    } // end of loop over bulk elements adjacent to boundary b
  }
 
 /// \short Delete surface elements
 void delete_surface_elements()
  {
   // How many surface elements are in the surface mesh
   unsigned n_element = Free_surface_bc_mesh_pt->nelement();
   
   // Loop over the surface elements
   for(unsigned e = 0; e < n_element; e++)
    {
     // Kill surface element
     delete Free_surface_bc_mesh_pt->element_pt(e);
    }
   
   // Wipe the mesh
   Free_surface_bc_mesh_pt->flush_element_and_node_storage();
  }

 void set_tip_node_pt()
  {
   Problem_Parameter::Tip_node_pt = 
    Bulk_mesh_pt->boundary_node_pt(Inner_boundary,0);
   
   /// Now check all nodes on bubble boundaries.. 
   /// Is this the one with the largest x-coordinate?
   unsigned num_nod=Bulk_mesh_pt->nboundary_node(Inner_boundary);
   for(unsigned inod=0; inod<num_nod; inod++)
    {
     // Get node
     Node* nod_pt=Bulk_mesh_pt->boundary_node_pt(Inner_boundary,inod);
     if (nod_pt->x(0) > Problem_Parameter::Tip_node_pt->x(0))
      {
       Problem_Parameter::Tip_node_pt = nod_pt;
      }
     if(nod_pt->is_on_boundary(Top_inflow_boundary))
      {
       Problem_Parameter::Top_node_pt = nod_pt;
      }
     if(nod_pt->is_on_boundary(Bottom_inflow_boundary))
      {
       Problem_Parameter::Bottom_node_pt = nod_pt;
      }
    }
   
  }

}; // end_of_problem_class


//========================================================================
/// Constructor of the problem class - build mesh and set up elements with
/// the appropriate parameters
//========================================================================
template<class ELEMENT>
UnstructuredHSFvKProblem<ELEMENT>::UnstructuredHSFvKProblem()
 : Mesh_as_geom_object_pt(0), Volume_constraint_mesh_pt(0), 
   Volume_constraint_element_pt(0)   
{

 Use_interface_motion_for_temporal_error=true;

 if (CommandLineArgs::command_line_flag_has_been_set("--error_estimate_use_displacement"))
  {
      
   oomph_info 
    << "///////////////////////////////////////////////////////////////////////"
    << std::endl;
   oomph_info << "Using displacement for error estimate"  << std::endl;
   oomph_info 
    << "///////////////////////////////////////////////////////////////////////"
    << std::endl << std::endl;
   
   Use_interface_motion_for_temporal_error=false;
  }

 if (CommandLineArgs::command_line_flag_has_been_set("--use_mumps"))
  {
   
#ifdef OOMPH_HAS_MUMPS
   
   
   oomph_info 
    << "///////////////////////////////////////////////////////////////////////"
    << std::endl;
   oomph_info << "Using MUMPS solver"  << std::endl;
   oomph_info 
    << "///////////////////////////////////////////////////////////////////////"
    << std::endl << std::endl;
   
   linear_solver_pt() = new MumpsSolver;

#else

   oomph_info 
    << "///////////////////////////////////////////////////////////////////////"
    << std::endl;
   oomph_info 
    << "Ignoring request to use MUMPS solver because we don't have it ..."  
    << std::endl;
   oomph_info 
    << "///////////////////////////////////////////////////////////////////////"
    << std::endl << std::endl;
   
   
#endif
  }



 // Allow for rough startup
 Max_newton_iterations=15;
 this->Problem::Max_residuals = 1000.0; // 1.0e20;
 
 // Output directory
 Doc_info.set_directory(Problem_Parameter::Directory);
 
 Vector<double> zeta(1);
 Vector<double> posn(2);
 
 //Outer boundary
 //--------------

 // Outer boundary: 6 distinct boundaries
 TriangleMeshClosedCurve* outer_boundary_pt = 0;
 Vector<TriangleMeshCurveSection*> outer_polyline_boundary_pt(6);

 // Segment vertex boundaries
 Vector<Vector<double> > bound_vertex(2);
 bound_vertex[0].resize(2);
 bound_vertex[1].resize(2);


 // Bottom boundary
 bound_vertex[0][0] = 0.0;
 bound_vertex[0][1] = -0.5;


 bound_vertex[1][0] = Problem_Parameter::Length;
 bound_vertex[1][1] = -0.5;

 outer_polyline_boundary_pt[0] = new TriangleMeshPolyLine(bound_vertex,
   Bottom_boundary);

 // Outflow boundary
 bound_vertex[0][0] = Problem_Parameter::Length;
 bound_vertex[0][1] = -0.5;

 bound_vertex[1][0] = Problem_Parameter::Length;
 bound_vertex[1][1] = 0.5;

 outer_polyline_boundary_pt[1] = new TriangleMeshPolyLine(bound_vertex,
   Outflow_boundary);

 // Top boundary
 bound_vertex[0][0] = Problem_Parameter::Length;
 bound_vertex[0][1] = 0.5;

 bound_vertex[1][0] = 0.0;
 bound_vertex[1][1] = 0.5;

 outer_polyline_boundary_pt[2] = new TriangleMeshPolyLine(bound_vertex,
   Top_boundary);

 bound_vertex.resize(2);
 bound_vertex[0].resize(2);
 bound_vertex[1].resize(2);

 // Top inflow boundary
 bound_vertex[0][0] = 0.0;
 bound_vertex[0][1] = 0.5;

 bound_vertex[1][0] = 0.0;
 bound_vertex[1][1] = Problem_Parameter::R_bubble_initial; 

 outer_polyline_boundary_pt[3] = new TriangleMeshPolyLine(bound_vertex,
   Top_inflow_boundary);

 // Actual inflow boundary
 bound_vertex[0][0] = 0.0;
 bound_vertex[0][1] = Problem_Parameter::R_bubble_initial;

 bound_vertex[1][0] = 0.0;
 bound_vertex[1][1] = -Problem_Parameter::R_bubble_initial;

 outer_polyline_boundary_pt[4] = new TriangleMeshPolyLine(bound_vertex,
   Actual_inflow_boundary);

 // Bottom inflow boundary
 bound_vertex[0][0] = 0.0;
 bound_vertex[0][1] = -Problem_Parameter::R_bubble_initial;

 bound_vertex[1][0] = 0.0;
 bound_vertex[1][1] = -0.5;

 
 outer_polyline_boundary_pt[5] = new TriangleMeshPolyLine(bound_vertex,
    Bottom_inflow_boundary);

 outer_boundary_pt =
  new TriangleMeshClosedCurve(outer_polyline_boundary_pt);

 //Inner boundary
 //--------------
 
 double A = Problem_Parameter::R_bubble_initial;
 double B = Problem_Parameter::R_bubble_initial;
 Ellipse* inner_boundary_ellipse_pt = new Ellipse(A, B);

 Vector<TriangleMeshOpenCurve*> inner_boundary_pt(1);

 Vector<TriangleMeshCurveSection*> inner_polyline_boundary_pt(1);

 // Number of segments on bubble
 unsigned nsegment=10;
 if (Problem_Parameter::N_segment_bubble>0)
  {
   nsegment = Problem_Parameter::N_segment_bubble; 
  }
 // Infer from target element area
 else
  {
   nsegment=unsigned(double(nsegment)*
                     sqrt(0.005/Problem_Parameter::Max_el_area));
  }

 double unit_zeta = -MathematicalConstants::Pi/double(nsegment);
 
 unsigned storage=nsegment+1;
 unsigned start_point=0;
 unsigned end_point=storage;

 if(Problem_Parameter::L_bubble_initial > 0.0)
  {
   storage+=2;
   start_point=1;
   end_point=storage-1;
  }

 Vector<Vector<double> > bound_hole(storage);

 //First part
 for(unsigned ipoint = start_point; ipoint < end_point; ipoint++)
  {
   bound_hole[ipoint].resize(2);

   zeta[0] = 0.5*MathematicalConstants::Pi+unit_zeta*double(ipoint-start_point);
   inner_boundary_ellipse_pt->position(zeta,posn);
   bound_hole[ipoint][0] = posn[0]+Problem_Parameter::L_bubble_initial;
   bound_hole[ipoint][1] = posn[1];
  }
 
 if(Problem_Parameter::L_bubble_initial > 0.0)
  {
   bound_hole[0].resize(2);
   bound_hole[0][0] = 0.0;
   bound_hole[0][1] = Problem_Parameter::R_bubble_initial;
   bound_hole[storage-1].resize(2);
   bound_hole[storage-1][0] = 0.0;
   bound_hole[storage-1][1] = -Problem_Parameter::R_bubble_initial;
  }

 inner_polyline_boundary_pt[0] = new TriangleMeshPolyLine(bound_hole,
   Inner_boundary);

 // Limit the size of the boundary segments 
 double max_length = Problem_Parameter::Max_polyline_length_bubble; 
 if (max_length>0.0)
  {
   inner_polyline_boundary_pt[0]->set_maximum_length(max_length);
  }

 inner_boundary_pt[0] = new TriangleMeshOpenCurve(inner_polyline_boundary_pt);

 // Connect initial vertex to final vertex on upper inlet boundary
 inner_polyline_boundary_pt[0]->connect_initial_vertex_to_polyline(
  dynamic_cast<TriangleMeshPolyLine*>
  (outer_polyline_boundary_pt[3]),1);
 
 // Connect final vertex to first vertex on lower outlet boundary
 inner_polyline_boundary_pt[0]->connect_final_vertex_to_polyline(
  dynamic_cast<TriangleMeshPolyLine*>
  (outer_polyline_boundary_pt[5]),0);

 Vector<double> bubble_region_coords1(2);
 bubble_region_coords1[0] = Problem_Parameter::R_bubble_initial/2.0;
 bubble_region_coords1[1] = Problem_Parameter::R_bubble_initial/2.0;


 double tol = 0.001;
 dynamic_cast<TriangleMeshCurve*>(inner_boundary_pt[0])->
 set_polyline_refinement_tolerance(tol);

 // Create timestepper
 bool adapt = true;
 this->add_time_stepper_pt(new BDF<2>(adapt));


 //Create the mesh
 //---------------

 //Create mesh parameters object
 TriangleMeshParameters mesh_parameters(outer_boundary_pt);

 mesh_parameters.internal_open_curves_pt() = inner_boundary_pt;

 mesh_parameters.add_region_coordinates(1, bubble_region_coords1);
 
 // Add the bubble region to the list of regions
 Bubble_regions.push_back(1);

 // Use attributes to enable identification of regions
 mesh_parameters.enable_use_attributes();

 // Set target area for elements in the initial mesh obacht
 mesh_parameters.element_area() = Problem_Parameter::Max_el_area; //0.005 

 // Build mesh
 Bulk_mesh_pt = new RefineableSolidTriangleMesh<ELEMENT>(mesh_parameters,
                                                  this->time_stepper_pt(0));

 // Snap nodes onto curvilinear boundary
 snap_onto_boundary();

 Z2ErrorEstimator* error_estimator_pt = new Z2ErrorEstimator;
 Bulk_mesh_pt->spatial_error_estimator_pt() = error_estimator_pt;
 Bulk_mesh_pt->max_permitted_error() = Problem_Parameter::Max_permitted_error; 
 Bulk_mesh_pt->min_permitted_error() = Problem_Parameter::Min_permitted_error;
 Bulk_mesh_pt->max_element_size() = Problem_Parameter::Max_el_area; 
 Bulk_mesh_pt->min_element_size() = Problem_Parameter::Min_el_area;


 set_tip_node_pt();

 if (!CommandLineArgs::command_line_flag_has_been_set("--pin_hs"))
  {
   /// Create the initial volume constraint element
   Volume_constraint_element_pt
    = new FoepplvonKarmanVolumeConstraintElement
    <ELEMENT, RefineableSolidTriangleMesh> (Bulk_mesh_pt,Bubble_regions);
   
   // Set volume constraint
   Volume_constraint_element_pt->set_prescribed_volume(
    &Problem_Parameter::Prescribed_volume);
   
   /// Set the problem pressure to be that of the volume constraint element
   Problem_Parameter::P_bubble_pt = 
    Volume_constraint_element_pt->pressure_data_pt();
   
   /// Create the volume constraint mesh and add the element to the mesh
   Volume_constraint_mesh_pt = new Mesh;
   Volume_constraint_mesh_pt->add_element_pt(Volume_constraint_element_pt);
   
   /// New mesh for the interface elements, and create them
   Free_surface_bc_mesh_pt = new Mesh;
   create_surface_elements();

   // Cross_section_mesh_pt = new Mesh;
   // create_cross_section_elements();
  }

 /// Pin the appropriate unknowns and set boundary conditions
 complete_problem_setup();

 // Bulk_mesh_pt->output("mesh.dat");
 // Bulk_mesh_pt->output_boundaries("boundaries.dat");
 // exit(0);


 // Create mesh as geom object representation of the mesh
 Mesh_as_geom_object_pt=new MeshAsGeomObject(Bulk_mesh_pt);

 // Set all Foeppl von Karmann displacements to zero
 unsigned n_node = Bulk_mesh_pt->nnode();
 for(unsigned i = 0; i < n_node; i++)
  {
   Node *nod_pt = Bulk_mesh_pt->node_pt(i);
   double w=0.0;
   nod_pt->set_value(0,w);
  }
if (CommandLineArgs::command_line_flag_has_been_set("--pin_fvk"))
 {
  update_pinned_fvk_displacements();
 }

 /// Add the sub meshes and build the global mesH
 add_sub_mesh(Bulk_mesh_pt);
 if (!CommandLineArgs::command_line_flag_has_been_set("--pin_hs"))
  {
   add_sub_mesh(Free_surface_bc_mesh_pt);
   add_sub_mesh(Volume_constraint_mesh_pt);
  }
 build_global_mesh();
 

 bool open_trace=true;
#ifdef OOMPH_HAS_MPI
 if (communicator_pt()->my_rank()!=0)
  {
   open_trace=false;
  }
#endif

 // Open trace file
 if (open_trace)
  {
   char filename[100];
   sprintf(filename, "%s/trace.dat", Doc_info.directory().c_str());
   Trace_file.open(filename);
  }

 // Assign boundary conditions
 oomph_info << "Number of equations: "
            << this->assign_eqn_numbers() << '\n';
 
}



//==start_of_complete======================================================
 /// Set boundary condition exactly, and complete the build of 
 /// all elements
//========================================================================
template<class ELEMENT>
void UnstructuredHSFvKProblem<ELEMENT>::complete_problem_setup()
{   
 
 // Apply boundary conditions obacht
 bool pin_bubble_ends=false; //true;
 apply_boundary_conditions(pin_bubble_ends);
 

 // Fix the pressure inside the bubble and disable HS equations
 {
  unsigned n_bubble_regions = Bubble_regions.size();
  for(unsigned r = 0; r < n_bubble_regions; r++)
   {
    unsigned n_element = Bulk_mesh_pt->nregion_element(Bubble_regions[r]);
    for(unsigned e = 0; e < n_element; e++)
     {
      ELEMENT* el_pt = dynamic_cast<ELEMENT*>
       (Bulk_mesh_pt->region_element_pt(Bubble_regions[r], e));
      
      // Switch off the element's contribution to the Hele Shaw eqns
      el_pt->disable_hele_shaw();

      // Pin pressure in the proper interior of the bubble
      // (pressure on the bubble surface is computed!)
      unsigned nnod = el_pt->nnode();
      for(unsigned j = 0; j < nnod; j++)
       {
        Node* nod_pt = el_pt->node_pt(j);
        if( (!(nod_pt->is_on_boundary(Inner_boundary)  )) )
         {
          // Pin HS pressure
          nod_pt->pin(4);
         }
       }  
      el_pt->add_external_data(Problem_Parameter::P_bubble_pt);
     }
   }
 }

 
 // Complete the build of all elements so they are fully functional
 unsigned n_element = Bulk_mesh_pt->nelement();
 for(unsigned e = 0; e < n_element; e++)
  {
   // Upcast from GeneralisedElement to the present element
   ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));
   
   //Set the pressure function pointers and the physical constants
   el_pt->eta_pt() = &Problem_Parameter::Eta;

   // Validation pressure (fluid everywhere -- overwritten below for bubble)
   //el_pt->pressure_fct_pt()=&Problem_Parameter::get_validation_pressure_fluid;

   // External pressure
   el_pt->pressure_fct_pt()=&Problem_Parameter::external_pressure;

   // Put the sheet under tension
   el_pt->pre_stress_fct_pt() = &Problem_Parameter::get_pre_stress;

   //Set the FSI parameter
   el_pt->q_pt() = &Problem_Parameter::Q;

   // Set the non-dim reference gap width
   el_pt->aspect_ratio_pt() =  &Problem_Parameter::Aspect_ratio;

   // Set the constitutive law for pseudo-elastic mesh deformation
   el_pt->constitutive_law_pt() = Problem_Parameter::Constitutive_law_pt;

   // Linear wall?
   if (CommandLineArgs::command_line_flag_has_been_set("--linear_wall"))
    {
     el_pt->use_linear_bending_model();
    }


   if (CommandLineArgs::command_line_flag_has_been_set("--pin_hs"))
    {
     // Switch off the element's contribution to the Hele Shaw eqns
     el_pt->disable_hele_shaw();
     
     // Pin pressure 
     unsigned nnod = el_pt->nnode();
     for(unsigned j = 0; j < nnod; j++)
      {
       Node* nod_pt = el_pt->node_pt(j);
       nod_pt->pin(4);
      }  
    }
   
  }


 // Pin Foeppl von Karman?
 if (CommandLineArgs::command_line_flag_has_been_set("--pin_fvk"))
  {
   pin_fvk();
  }

 // Add the bubble pressure as external data to the nodes in the bubble
 // region
 unsigned n_bubble_regions = Bubble_regions.size();
 for(unsigned r = 0; r < n_bubble_regions; r++)
  {
   n_element = Bulk_mesh_pt->nregion_element(Bubble_regions[r]);
   for(unsigned e = 0; e < n_element; e++)
    {
     
     // Upcast from GeneralisedElement to the present element
     ELEMENT* el_pt = dynamic_cast<ELEMENT*>(
      Bulk_mesh_pt->region_element_pt(Bubble_regions[r],e));
     
     if (!CommandLineArgs::command_line_flag_has_been_set("--pin_hs"))
      {
       el_pt->set_volume_constraint_pressure_data_as_external_data(
         Volume_constraint_element_pt->pressure_data_pt());
      }

     // // Validation pressure overwrite for bubble)
     // el_pt->pressure_fct_pt()=
     //  &Problem_Parameter::get_validation_pressure_bubble;

     // External pressure
     el_pt->pressure_fct_pt()=&Problem_Parameter::external_pressure;


     //el_pt->add_external_data(Problem_Parameter::P_bubble_pt);
      
     if(CommandLineArgs::command_line_flag_has_been_set("--inc_thin_films"))
      {
       el_pt->enable_thin_film_effects();
       el_pt->ca_pt() = &Problem_Parameter::Ca;
       el_pt->thin_film_homotopy_pt() = 
        &Problem_Parameter::Thin_film_homotopy_parameter;
      }
    }
  }

 // Add the bubble pressure as external data to the interface elements
 n_element = Free_surface_bc_mesh_pt->nelement();
 for(unsigned e = 0; e < n_element; e++)
  {
   // Upcast from GeneralisedElement to the present element
   HeleShawInterfaceElement<ELEMENT>* el_pt = 
    dynamic_cast<HeleShawInterfaceElement<ELEMENT>* >(
     Free_surface_bc_mesh_pt->element_pt(e));
   
   el_pt->add_external_data(Problem_Parameter::P_bubble_pt,true);

   // if we add in thin film effects we need to compute pressure gradients
   // and sheet deflections and hence need to add in external data
   if (CommandLineArgs::command_line_flag_has_been_set("--inc_thin_films"))
    {
     ELEMENT* bulk_el_pt = dynamic_cast<ELEMENT*>(el_pt->bulk_element_pt());
     unsigned bulk_nodes = bulk_el_pt->nnode();
     for(unsigned inod=0; inod<bulk_nodes; inod++)
      {
       Node* nod_pt = dynamic_cast<Node*>(bulk_el_pt->node_pt(inod));
       SolidNode* solid_nod_pt = dynamic_cast<SolidNode*>(nod_pt);
       if( !(nod_pt->is_on_boundary(Inner_boundary)) )
        {
         el_pt->add_external_data(nod_pt, true);
         el_pt->add_external_data(solid_nod_pt->variable_position_pt(),true);
        }
      }
    }
  }

}



//==start_of_doc_solution=================================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void UnstructuredHSFvKProblem<ELEMENT>::doc_solution()
{ 
 
 // Doc it?
 bool do_it=true;
 
 // No need for full output from all processors
#ifdef OOMPH_HAS_MPI
 if (communicator_pt()->my_rank()!=0)
  {
   do_it=false;
   oomph_info << "Omitting (most) output from proc: " 
              << communicator_pt()->my_rank() << std::endl;
  }
#endif
 
 // Full doc?
 if (do_it)
  {
   Problem_Parameter::doc_parameters();
   
   oomph_info << "Outputting for step: " << Doc_info.number() 
              << " ; time= " 
              << Problem_Parameter::Problem_pt->time_stepper_pt()->time()
              << std::endl;
   
   ofstream some_file;
   char filename[100];
 
   // Number of plot points
   unsigned npts = 5;
 
   // Full soln
   //----------
   sprintf(filename,"%s/soln%i.dat",Doc_info.directory().c_str(),
           Doc_info.number());
   some_file.open(filename);
   this->Bulk_mesh_pt->output(some_file,npts); 
   some_file << "TEXT X = 22, Y = 92, CS=FRAME T = \"" << "\"\n";
   some_file.close();
 
   // Output boundaries
   //------------------
   sprintf(filename,"%s/boundaries%i.dat",Doc_info.directory().c_str(),
           Doc_info.number());
   some_file.open(filename);
   Bulk_mesh_pt->output_boundaries(some_file);
   some_file.close();
 
 
   // // Output boundary coordinates
   // //----------------------------
   // oomph_info << "about to do nboundary\n";
   // oomph_info.stream_pt()->flush();
   // unsigned nb=Bulk_mesh_pt->nboundary();
   // oomph_info << "done nboundary: " << nb << std::endl;
   // oomph_info.stream_pt()->flush();
   // for (unsigned b=0;b<nb;b++)
   //  {

   //   oomph_info << "outputting bound coords for bound: " << b << std::endl;
   //   oomph_info.stream_pt()->flush();
   //   sprintf(filename,"%s/boundary_coords_for_boundary%i_%i.dat",
   //           Doc_info.directory().c_str(),
   //           b,
   //           Doc_info.number());
   //   some_file.open(filename);
   //   Bulk_mesh_pt->template doc_boundary_coordinates<ELEMENT>(b,some_file);
   //   some_file.close();


   //   oomph_info << "done outputting bound coords for bound: " << b << std::endl;
   //   oomph_info.stream_pt()->flush();
   //  } 


   // Quantities at Gauss points
   //---------------------------
   // sprintf(filename,"%s/gauss_point_data%i.dat",Doc_info.directory().c_str(),
   //         Doc_info.number());
   // some_file.open(filename);
   // unsigned nel = Bulk_mesh_pt->nelement();
   // for (unsigned e = 0; e < nel; e++)
   //  {
   //   dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e))->
   //    output_at_gauss_points(some_file);
   //  }
   // some_file.close();

   // Fluid region
   //-------------
   // sprintf(filename,"%s/fluid%i.dat",Doc_info.directory().c_str(),
   //         Doc_info.number());
   // some_file.open(filename);
   // unsigned r = 0;
   // unsigned nel = Bulk_mesh_pt->nregion_element(r);
   // for (unsigned e = 0; e < nel; e++)
   //  {
   //   Bulk_mesh_pt->region_element_pt(r, e)->output(some_file, npts);
   //  }
   // some_file.close();


   // Fluid region
   //------------- 
   sprintf(filename,"%s/fluid_coarse%i.dat",Doc_info.directory().c_str(),
           Doc_info.number());
   some_file.open(filename);
   unsigned r = 0;
   unsigned nel = Bulk_mesh_pt->nregion_element(r);
   for (unsigned e = 0; e < nel; e++)
    {
     Bulk_mesh_pt->region_element_pt(r, e)->output(some_file, 2);
    }
   some_file.close();



   // Bubble region
   //--------------
   // sprintf(filename,"%s/bubble%i.dat",Doc_info.directory().c_str(),
   //         Doc_info.number());
   // some_file.open(filename);
   // unsigned n_bubble_regions = Bubble_regions.size();
   // for(unsigned r = 0; r < n_bubble_regions; r++)
   //  {
   //   nel = Bulk_mesh_pt->nregion_element(Bubble_regions[r]);
   //   for (unsigned e = 0; e < nel; e++)
   //    {
   //     Bulk_mesh_pt->
   //      region_element_pt(Bubble_regions[r],e)->output(some_file, npts);
   //    }
   //  }
   // some_file.close();

   // Bubble region
   //--------------
   sprintf(filename,"%s/bubble_coarse%i.dat",Doc_info.directory().c_str(),
           Doc_info.number());
   some_file.open(filename);
   unsigned n_bubble_regions = Bubble_regions.size();
   for(unsigned r = 0; r < n_bubble_regions; r++)
    {
     nel = Bulk_mesh_pt->nregion_element(Bubble_regions[r]);
     for (unsigned e = 0; e < nel; e++)
      {
       Bulk_mesh_pt->
        region_element_pt(Bubble_regions[r],e)->output(some_file,2);
      }
    }
   some_file.close();


   // Interface
   //----------
   {
    sprintf(filename,"%s/interface%i.dat",Doc_info.directory().c_str(),
            Doc_info.number());
    some_file.open(filename);
    
    // How many surface elements are in the surface mesh
    unsigned n_element = Free_surface_bc_mesh_pt->nelement();
    
    // Loop over the surface elements
    for(unsigned e = 0; e < n_element; e++)
     {
      HeleShawInterfaceElement<ELEMENT>* surface_element_pt = 
       dynamic_cast<HeleShawInterfaceElement<ELEMENT>*>(
        Free_surface_bc_mesh_pt->element_pt(e));
      surface_element_pt->output(some_file,npts);
     }
    some_file.close();
   }


   // Extract solution along horizontal line across mesh
   sprintf(filename,"%s/horizontal_line%i.dat",Doc_info.directory().c_str(),
           Doc_info.number());
   some_file.open(filename);
   double w_0=0.0;
   {
    // Cooordinate we're trying to find
    Vector<double> x(2);
    Vector<double> s(2);
    unsigned nplot=1001;
    for (unsigned j=0;j<nplot;j++)
     {
      // Point we're trying to find
      //x[0]=double(j)/double(nplot-1);
      x[0]=double(j)*(Problem_Parameter::Length)/double(nplot-1);
      x[1]=0.0;
    
      // Sub-geomobject (FE) that contains that point
      GeomObject* sub_geom_object_pt=0;
    
      // Find it
      Mesh_as_geom_object_pt->locate_zeta(x, sub_geom_object_pt, s);

      if (sub_geom_object_pt!=0)
       {
        // Cast to our element
        ELEMENT* el_pt=dynamic_cast<ELEMENT*>(sub_geom_object_pt);
        if (el_pt!=0)
         {
          double w  = el_pt->interpolated_w_fvk(s,0);

          if(j==0)
           {
            w_0 = w;
           }
    
          // Get in-plane stress
          DenseMatrix<double> sigma(2,2,0.0);
          DenseMatrix<double> strain(2,2,0.0);
          el_pt->get_stress_and_strain_for_output(s,sigma,strain);


          some_file << el_pt->interpolated_x(s,0) << " " 
                    << el_pt->interpolated_x(s,1) << " " 
                    << el_pt->interpolated_w_fvk(s,0) << " "
                    << el_pt->interpolated_p_hele_shaw(s) << " "
                    << el_pt->interpolated_w_fvk(s,1) << " "
                    << el_pt->interpolated_w_fvk(s,2) << " "
                    << el_pt->interpolated_w_fvk(s,3) << " "
                    << sigma(0,0) << " " 
                    << sigma(1,1) << " " 
                    << sigma(0,1) << " "
                    << strain(0,0) << " "
                    << strain(1,1) << " "
                    << strain(0,1) << " "
                    << std::endl;
         }
        else
         {
          oomph_info << "Cast failed for point at " << x[0] << " " 
                     << x[1] << "\n"; 
         }
       }
      else
       {
        oomph_info << "Point at " << x[0] << " " 
                   << x[1] << " not found.\n"; 
       }
     }
   }
   some_file.close();


   // Extract solution along vertical line across mesh
   sprintf(filename,"%s/vertical_line%i.dat",Doc_info.directory().c_str(),
           Doc_info.number());
   some_file.open(filename);
   {
    // Cooordinate we're trying to find
    Vector<double> x(2);
    Vector<double> s(2);
    unsigned nplot=101;
    for (unsigned j=0;j<nplot;j++)
     {
      // Point we're trying to find
      x[0]=Problem_Parameter::Length/2.0;
      x[1]=-0.5+double(j)/double(nplot-1);
    
      // Sub-geomobject (FE) that contains that point
      GeomObject* sub_geom_object_pt=0;
    
      // Find it
      Mesh_as_geom_object_pt->locate_zeta(x, sub_geom_object_pt, s);

      if (sub_geom_object_pt!=0)
       {
        // Cast to our element
        ELEMENT* el_pt=dynamic_cast<ELEMENT*>(sub_geom_object_pt);
        if (el_pt!=0)
         {
          double w  = el_pt->interpolated_w_fvk(s,0);

          if(j==0)
           {
            w_0 = w;
           }
    
          // Get in-plane stress
          DenseMatrix<double> sigma(2,2,0.0);
          DenseMatrix<double> strain(2,2,0.0);
          el_pt->get_stress_and_strain_for_output(s,sigma,strain);


          some_file << el_pt->interpolated_x(s,0) << " " 
                    << el_pt->interpolated_x(s,1) << " " 
                    << el_pt->interpolated_w_fvk(s,0) << " "
                    << el_pt->interpolated_p_hele_shaw(s) << " "
                    << el_pt->interpolated_w_fvk(s,1) << " "
                    << el_pt->interpolated_w_fvk(s,2) << " "
                    << el_pt->interpolated_w_fvk(s,3) << " "
                    << sigma(0,0) << " " 
                    << sigma(1,1) << " " 
                    << sigma(0,1) << " "
                    << strain(0,0) << " "
                    << strain(1,1) << " "
                    << strain(0,1) << " "
                    << std::endl;
         }
        else
         {
          oomph_info << "Cast failed for point at " << x[0] << " " 
                     << x[1] << "\n"; 
         }
       }
      else
       {
        oomph_info << "Point at " << x[0] << " " 
                   << x[1] << " not found.\n"; 
       }
     }
   }
   some_file.close();



   // Get the volume of the bubble region
   double bubble_volume = 0.0;
   if (!CommandLineArgs::command_line_flag_has_been_set("--pin_hs"))
    {
     bubble_volume=Volume_constraint_element_pt->get_bounded_volume();
    }

   // Get mean radius
   double mean_radius = 0.0;
   double arc_length = 0.0;
   double arc_length_contribution;
   double mean_radius_nod_contribution; 
   double deviation = 0.0;

   if (!CommandLineArgs::command_line_flag_has_been_set("--pin_hs"))
    {
     // How many surface elements are in the surface mesh
     unsigned n_element = Free_surface_bc_mesh_pt->nelement();
   
     // Loop over the surface elements
     for(unsigned e = 0; e < n_element; e++)
      {
       HeleShawInterfaceElement<ELEMENT>* surface_element_pt = 
        dynamic_cast<HeleShawInterfaceElement<ELEMENT>*>(
         Free_surface_bc_mesh_pt->element_pt(e));
     
       surface_element_pt->compute_mean_radius(
        arc_length_contribution, mean_radius_nod_contribution);
     
       mean_radius += mean_radius_nod_contribution;
       arc_length += arc_length_contribution;
      }
     mean_radius /= arc_length;
   
   
     // Get deviation of the radius
     double deviation_contribution=0.0;
   
     // Loop over the surface elements
     for(unsigned e = 0; e < n_element; e++)
      {
       HeleShawInterfaceElement<ELEMENT>* surface_element_pt = 
        dynamic_cast<HeleShawInterfaceElement<ELEMENT>*>(
         Free_surface_bc_mesh_pt->element_pt(e));
     
       surface_element_pt->compute_deviation_from_mean_radius(
        mean_radius, deviation_contribution);
     
       deviation += deviation_contribution;
     
      }
     deviation = sqrt(2.0*deviation/arc_length);
    }

   double p_bubble=0.0;
   if (Problem_Parameter::P_bubble_pt!=0)
    {
     p_bubble=Problem_Parameter::P_bubble_pt->value(0);
    }
   // Output time, bubble pressure and bubble volume etc to the trace file
   Trace_file << time_stepper_pt()->time() << " "
              << mean_radius << " "
              << deviation << " "
              << p_bubble << " "
              << Problem_Parameter::Tip_node_pt->x(0) << " "
              << Problem_Parameter::Tip_node_pt->dposition_dt(0) << " "
              << bubble_volume << " "
              << global_temporal_error_norm() << " "
              << time_stepper_pt()->time_pt()->dt() << " "
              << Bulk_mesh_pt->nnode() << " "
              << Bulk_mesh_pt->nelement() << " "
              << ndof() << " "
              << Doc_info.number() << '\n';

   Trace_file.flush();
   
  } // end of full output
 
 
 // Write restart file on all procs -- awkward but this call is
 // needed because dump re-orders the nodes in the mesh!
 int rank=communicator_pt()->my_rank();
 ofstream some_file;
 char filename[100];
 sprintf(filename,"%s/restart_proc%i_%i.dat",
         Doc_info.directory().c_str(),
         rank,
         Doc_info.number());
 some_file.open(filename);
 some_file.precision(20); 
 dump_it(some_file);
 some_file.close();
 
 // Increment the doc_info number
 Doc_info.number()++;

} // end of doc


//=======start_of_main========================================
///Driver code for demo of inline triangle mesh generation
//============================================================
int main(int argc, char **argv)
{

#ifdef OOMPH_HAS_MPI
 MPI_Helpers::init(argc,argv,false);
#endif


 //feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW);

 // Store command line arguments
 CommandLineArgs::setup(argc,argv);
 
 // Sometimes the Newton method produces inverted elements while
 // it's iterating towards a solution with non-inverted elements:
 // accept these!
 FiniteElement::Accept_negative_jacobian = true;
 
 // Define possible command line arguments and parse the ones that
 // were actually specified

 // Name of restart file
 CommandLineArgs::specify_command_line_flag("--restart_file",
                                            &Problem_Parameter::Restart_file);
  
 // Pin FvK, i.e. Hele Shaw only?
 CommandLineArgs::specify_command_line_flag("--pin_fvk");
 
 // Pin Hele Shaw i.e. FvK only?
 CommandLineArgs::specify_command_line_flag("--pin_hs");

 // Perturbation wavenumber
 CommandLineArgs::specify_command_line_flag("--n_perturb",
                                            &Problem_Parameter::N_perturb);
 
 // Initial bubble radius
 CommandLineArgs::specify_command_line_flag(
  "--r_bubble",
  &Problem_Parameter::R_bubble_initial);
 
 // Initial bubble length
 CommandLineArgs::specify_command_line_flag(
  "--l_bubble",
  &Problem_Parameter::L_bubble_initial);

 // Initial perturbation of the bubble radius
 CommandLineArgs::specify_command_line_flag(
  "--r_perturb",
  &Problem_Parameter::Radius_perturb);

 // Initial number of fingers
 CommandLineArgs::specify_command_line_flag(
  "--n_perturb",
  &Problem_Parameter::N_perturb);

 // Perturbation of the bubble pressure
 CommandLineArgs::specify_command_line_flag(
  "--p_perturb",
  &Problem_Parameter::Pressure_perturb);

 // Aspect ratio
 CommandLineArgs::specify_command_line_flag(
  "--aspect_ratio",
  &Problem_Parameter::Aspect_ratio);

 // Inverse capillary number
 CommandLineArgs::specify_command_line_flag(
  "--ca_inv",
  &Problem_Parameter::Ca_inv);

 // FvK coefficient : Eta
 CommandLineArgs::specify_command_line_flag(
  "--eta",
  &Problem_Parameter::Eta);

 // fsi parameter
 CommandLineArgs::specify_command_line_flag(
  "--q",
  &Problem_Parameter::Q);

 // What sort of pre-stress are we using?
 CommandLineArgs::specify_command_line_flag("--pre_stress", 
                                            &Problem_Parameter::Pre_stress);

 // Linear wall?
 CommandLineArgs::specify_command_line_flag("--linear_wall");

 // Use additional fvk pressure?
 CommandLineArgs::specify_command_line_flag("--use_additional_fvk_pressure");

 // Adaptation frequency (negative: never; 1: always)
 CommandLineArgs::specify_command_line_flag(
  "--adapt_frequency",
  &Problem_Parameter::Adaptation_frequency);

 // Suppress update of Lagrangian coords after adapt?
 CommandLineArgs::specify_command_line_flag(
  "--suppress_reset_lagr_coords");

 // Timestep
 double dt = 1.0e-4;
 CommandLineArgs::specify_command_line_flag("--dt",&dt);

 // Suppress temporal adaptativity?
 CommandLineArgs::specify_command_line_flag("--suppress_temporal_adaptivity");
 
 // Number of steps
 unsigned nstep=1000;
 CommandLineArgs::specify_command_line_flag("--nstep",&nstep);

 // Target error for adaptive timestepping
 double epsilon_t = 1.0e-1; //1.0e-5; // obacht
 CommandLineArgs::specify_command_line_flag("--epsilon_t",&epsilon_t);

 // Output directory
 CommandLineArgs::specify_command_line_flag("--dir",
                                            &Problem_Parameter::Directory);
 

 // Initial/max element area
 CommandLineArgs::specify_command_line_flag("--max_el_area",
                                            &Problem_Parameter::Max_el_area);


 // Min element area during refinement
 CommandLineArgs::specify_command_line_flag("--min_el_area",
                                            &Problem_Parameter::Min_el_area);


 // Max spatial error
 CommandLineArgs::specify_command_line_flag(
  "--max_permitted_spatial_error",
  &Problem_Parameter::Max_permitted_error);

 // Min spatial error
 CommandLineArgs::specify_command_line_flag(
  "--min_permitted_spatial_error",
  &Problem_Parameter::Min_permitted_error);
 
 // Max length of polylines on bubble
 CommandLineArgs::specify_command_line_flag(
  "--max_polyline_length_bubble",
  &Problem_Parameter::Max_polyline_length_bubble);

 // Estimating the time step based on the displacement of the sheet
 CommandLineArgs::specify_command_line_flag(
  "--error_estimate_use_displacement");

 // Use mumps
 CommandLineArgs::specify_command_line_flag("--use_mumps");

 // Are we applying thin film corrections?
 CommandLineArgs::specify_command_line_flag("--inc_thin_films");

 // Are we using the same boundary conditions as in the moving frame?
 CommandLineArgs::specify_command_line_flag("--moving_frame_bc");

 // Parse command line
 CommandLineArgs::parse_and_assign(); 
 
#ifdef OOMPH_HAS_MPI

  // Switch off output modifier
 oomph_info.output_modifier_pt() = &default_output_modifier;

 // Define processor-labeled output file for all on-screen stuff
 std::ofstream output_stream;
 char filename[1000];
 sprintf(filename,"%s/OUTPUT.%i",Problem_Parameter::Directory.c_str(),
         MPI_Helpers::communicator_pt()->my_rank());
 output_stream.open(filename);
 oomph_info.stream_pt() = &output_stream;
 OomphLibWarning::set_stream_pt(&output_stream);
 OomphLibError::set_stream_pt(&output_stream);   

#endif

 // Doc what has actually been specified on the command line
 CommandLineArgs::doc_specified_flags();
 
 // Build constititive eqn for pseudo-solid
 Problem_Parameter::Constitutive_law_pt =
  new GeneralisedHookean(&Problem_Parameter::Nu);
 
 
 // Create Problem instance
 UnstructuredHSFvKProblem<ProjectableMyElement> problem;
 
 // Set pointer to problem (so we have access to global time)
 Problem_Parameter::Problem_pt = &problem;
 
 // Restart file specified via command line 
 if (CommandLineArgs::command_line_flag_has_been_set("--restart_file"))
  {
   problem.restart();
   oomph_info << "Done restart\n";
   //problem.doc_solution();
   //pause("done");
  }
 // No restart
 else
  {
   problem.initialise_dt(dt);
   problem.actions_before_implicit_timestep();
   problem.assign_initial_values_impulsive(dt);

   // Set initial volume for steady solve
   Problem_Parameter::Prescribed_volume = Pi*0.5*
    Problem_Parameter::R_bubble_initial*
    Problem_Parameter::R_bubble_initial+
    Problem_Parameter::L_bubble_initial*2.0*Problem_Parameter::R_bubble_initial;

   if (CommandLineArgs::command_line_flag_has_been_set("--inc_thin_films"))
    {
     Problem_Parameter::Prescribed_volume *= 
      (1.0-Problem_Parameter::thin_film_effect());
    }
   
   // Output initial conditions
   problem.doc_solution();
   
   bool do_steady_presolve=true; //true; // obacht
   if (do_steady_presolve)
    {
     // Steady solve without FSI
     //double q_backup = Problem_Parameter::Q;
     //double pre_stress_backup = Problem_Parameter::Pre_stress;

     //Problem_Parameter::Q = 1.0e-7;
     //Problem_Parameter::Pre_stress = 0.0;
     //Problem_Parameter::update_dependent_parameters();
     
     oomph_info << "Pre-solving for prescribed volume: " 
                << Problem_Parameter::Prescribed_volume  
                << " and pressure perturbation: " 
                << Problem_Parameter::Pressure_perturb
                << " and P_ext: " 
                << Problem_Parameter::P_ext
                << " and Q : "
                << Problem_Parameter::Q
                << std::endl;
     
     oomph_info << "Pinning Hele Shaw dofs." << std::endl;
     problem.pin_hs();
     //Problem_Parameter::P_bubble_pt->pin(0);
     //Problem_Parameter::P_bubble_pt->set_value(0,0.0);
     oomph_info << "Number of equations: " 
                <<problem.assign_eqn_numbers() << std::endl;

     double increment = -50.0;
     while(Problem_Parameter::P_ext > 
           Problem_Parameter::Target_p_ext-increment)
     //while(Problem_Parameter::Outlet_pressure > 
     //      Problem_Parameter::Target_outlet_pressure-increment)
      {
       // Do steady Newton solve first to get steady equilibrium shape
       // with perturbed interface shape
       problem.steady_newton_solve();
       
       // Output steady solution
       problem.doc_solution();

       Problem_Parameter::P_ext += increment;
       //Problem_Parameter::Outlet_pressure += increment;
       //problem.apply_boundary_conditions(false);
      }
     
     Problem_Parameter::P_ext=Problem_Parameter::Target_p_ext;
     //Problem_Parameter::Outlet_pressure = 
     // Problem_Parameter::Target_outlet_pressure;
     //problem.apply_boundary_conditions(false);
     problem.steady_newton_solve();
     problem.doc_solution();
     // Re-enable FSI
     //Problem_Parameter::Q = q_backup;  
     //Problem_Parameter::Pre_stress = pre_stress_backup;   
     
     oomph_info << "Unpinning Hele Shaw dofs." << std::endl;
     problem.unpin_hs();
     //Problem_Parameter::P_bubble_pt->unpin(0);
     problem.complete_problem_setup();
     oomph_info << "Number of equations: " 
                <<problem.assign_eqn_numbers() << std::endl;
    }
   else
    {
     oomph_info << "Bypassing steady pre-solve\n";
    }

  }

 //exit(1);

 // hierher 
 //Problem_Parameter::P_ext=1.0e-5;

 // Reset the pressure perturbation
 //Problem_Parameter::Radius_perturb = 0.0; 
 Problem_Parameter::update_dependent_parameters();

 unsigned max_adapt = 0;
 bool first = true;
 //problem.unsteady_newton_solve(1.0e-12,max_adapt,first);
 for(unsigned i = 0; i < nstep; i++)
  {   
   // Allow adaptation every few steps
   if( (Problem_Parameter::Adaptation_frequency>0) &&
       (i%Problem_Parameter::Adaptation_frequency == 0) && i>0 )
    {
     problem.adapt();//max_adapt = 1;
     problem.unsteady_newton_solve(1.0e-8,max_adapt,first);
    }
   if (Problem_Parameter::Adaptation_frequency<0) max_adapt = 0;
   
   // After the first step reset the pressure perturbation
   // if(i==20)
   //  Problem_Parameter::Radius_perturb = 0.0;

   // if(i==50)
   //  {
   //   Problem_Parameter::Thin_film_homotopy_parameter=1.0;
   //  }

   // Upgrade the parameters
   Problem_Parameter::update_dependent_parameters();
   

   //problem.unsteady_newton_solve(dt);


   double next_dt = dt;
   if (CommandLineArgs::command_line_flag_has_been_set
       ("--suppress_temporal_adaptivity"))
    {
     problem.unsteady_newton_solve(dt,max_adapt,first);
    }
   else
    {
     next_dt=
      //problem.doubly_adaptive_unsteady_newton_solve(dt,epsilon_t,
      //                                              max_adapt,first);
      problem.adaptive_unsteady_newton_solve(dt,epsilon_t,first);
     oomph_info << "Suggested next dt: " << next_dt << std::endl;
    }
   first=false;
   
   dt = next_dt; 
   
   problem.doc_solution();

   // Revert to no adaptation
   max_adapt = 0;

  }

 oomph_info << "done\n";

#ifdef OOMPH_HAS_MPI
 MPI_Helpers::finalize();
#endif

} //End of main


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
 
 // Joao defining the fill in jacobian and mass, in order to use in the eigensolver
 void fill_in_contribution_to_jacobian_and_mass_matrix(Vector<double> &residuals,
                                       DenseMatrix<double> &jacobian,DenseMatrix<double> &mass_matrix)
  {
  this->fill_in_contribution_to_jacobian(residuals,jacobian);   
  }

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
   
   // obacht new
   // Now fill in the off-diagonal entries (the shape derivatives),
   this->fill_in_shape_derivatives(jacobian); 
  }

 inline void update_in_external_fd(const unsigned &i)
  {
   //oomph_info << "Aloha!"<<std::endl;
   // if we are in the bubble we need to update the nodal pressures to
   // make sure they match the (updated) bubble pressure
   if(this->Hele_shaw_disabled)
    {
     // Get the current bubble pressure
     double bubble_pressure = this->external_data_pt(0)->value(0);
     // Loop over all nodes
     unsigned n_node = this->nnode();
     for(unsigned inod=0; inod<n_node; inod++)
      {
       Node* nod_pt = this->node_pt(inod);
       // Only update when the node is in the interior of the bubble where the
       // pressure dof is pinned
       if(nod_pt->is_pinned(4))
        {
         //oomph_info << "Updating pressure!"<<std::endl;
         nod_pt->set_value(4,bubble_pressure);
        }
      }
    }
  }

 inline void update_before_external_fd()
  {
   const unsigned i=0;
   update_in_external_fd(i);
  }

 inline void update_before_nodal_fd()
  {
   const unsigned i=0;
   update_in_external_fd(i);
  }

 inline void update_before_solid_position_fd()
  {
   const unsigned i=0;
   update_in_external_fd(i);
  }

 inline void update_in_solid_position_fd(const unsigned &i)
  {
   update_in_external_fd(i);
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
       << "MyElement elements only store five fields so fld must be"
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
       << "MyElement elements only store five fields so fld must be"
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
       << "MyElement elements only store five fields so fld must be"
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
       << "MyElement elements only store five fields so fld must be"
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
       << "MyElement elements only store five fields so fld must be"
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
       << "MyElement elements only store five fields so fld must be"
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


//=====================================================================
/// A constraint element that constrains the inflation upstream
/// and based on that computes the flux downstream.
//=====================================================================
 class DisplacementConstraintElement : public GeneralisedElement
 {
 public:

  DisplacementConstraintElement(double& pressure_gradient,
                                double& bubble_pressure,
                                double& bubble_tip_speed,
                                double& downstream_area,
                                Mesh* downstream_flux_mesh_pt,
                                Mesh* flow_rate_flux_mesh_pt,
                                SolidNode* bubble_tip_node_pt,
                                double* prescribed_x_position_pt,
                                Node* displacement_control_node_pt,
                                double* prescribed_displacement_pt)
   {
    // Create internal data object, which is the unknown pressure
    // gradient ahead of the bubble
    this->add_internal_data(new Data(1),true);
    this->internal_data_pt(0)->set_value(0,pressure_gradient);
    
    // The unknown bubble pressure is also internal data of this element
    this->add_internal_data(new Data(1),true);
    this->internal_data_pt(1)->set_value(0,bubble_pressure);
    
    // The unknown bubble tip velocity is also internal data of this element
    this->add_internal_data(new Data(1),true);
    this->internal_data_pt(2)->set_value(0,bubble_tip_speed);
    
    // The unknown flux downstream is also internal data of this element
    this->add_internal_data(new Data(1),true);
    this->internal_data_pt(3)->set_value(0,downstream_area);

    // Set pointer to the flux mesh ahead of the bubble
    Downstream_flux_mesh_pt = downstream_flux_mesh_pt;
    // Set pointer to the flux mesh that constitutes the bubble inflow
    Flow_rate_flux_mesh_pt = flow_rate_flux_mesh_pt;
    // Set pointer to solid node of which we track x position
    Bubble_tip_node_pt = bubble_tip_node_pt;
    // Set pointer to the displacement control node
    Displacement_control_node_pt = displacement_control_node_pt;
    // Set the prescribed x position of the bubble tip
    Prescribed_x_position_pt = prescribed_x_position_pt;
    // Set the prescribed displacement
    Prescribed_displacement_pt = prescribed_displacement_pt;

    // Now we need to add external data

    // Looks like we need to add all the nodes in the bulk elements
    // adjacent to the flux mesh
    // obacht WHY?? because we compute gradients!
    unsigned n_el=Downstream_flux_mesh_pt->nelement();
    for(unsigned e=0;e<n_el;e++)
     {     
      HeleShawFluxElement<ProjectableMyElement>* el_pt = 
       dynamic_cast<HeleShawFluxElement<ProjectableMyElement>*>(
        Downstream_flux_mesh_pt->element_pt(e));
      ProjectableMyElement* bulk_el_pt = dynamic_cast<ProjectableMyElement*>(
       el_pt->bulk_element_pt());
      unsigned n_node = bulk_el_pt->nnode();
      for(unsigned inod=0;inod<n_node;inod++)
       {
        this->add_external_data(bulk_el_pt->node_pt(inod),true);
       }
     }

    n_el=Flow_rate_flux_mesh_pt->nelement();
    for(unsigned e=0;e<n_el;e++)
     {     
      HeleShawFluxElement<ProjectableMyElement>* el_pt = 
       dynamic_cast<HeleShawFluxElement<ProjectableMyElement>*>(
        Flow_rate_flux_mesh_pt->element_pt(e));
      unsigned n_node = el_pt->nnode();
      for(unsigned inod=0;inod<n_node;inod++)
       {
        this->add_external_data(el_pt->node_pt(inod),true);
        // Need to add the position pt because changes in bubble width
        // affect flow rate
        SolidNode* solid_nod_pt = dynamic_cast<SolidNode*>(
         el_pt->node_pt(inod));
        this->add_external_data(solid_nod_pt->variable_position_pt(),true);
       }
     }

    // add the node we track as external data (but not the node itself,
    // instead the position data of the node)
    this->add_external_data(Bubble_tip_node_pt->variable_position_pt(),true);
    // add the displacement control node
    this->add_external_data(Displacement_control_node_pt,true);
   }
        
  ~DisplacementConstraintElement(){}

  double calculate_volume_flux_ahead_of_bubble()
   {
    // get the volume flux from the flux elements
    double volume_flux = 0.0;
    unsigned n_element=Downstream_flux_mesh_pt->nelement();
    for(unsigned el=0; el<n_element; el++)
     {
      // Get pointer to flux element
      HeleShawFluxElement<ProjectableMyElement>* el_pt = 
       dynamic_cast<HeleShawFluxElement<ProjectableMyElement>*>(
        Downstream_flux_mesh_pt->element_pt(el));
    
      // Add its contribution
      volume_flux += el_pt->compute_volume_flux();
     }
    return volume_flux;
   }

  double calculate_bubble_inflow_cross_section()
   {
    // get the cross section from the flux elements
    double cross_section = 0.0;
    unsigned n_element=Flow_rate_flux_mesh_pt->nelement();
    for(unsigned el=0; el<n_element; el++)
     {
      // Get pointer to flux element
      HeleShawFluxElement<ProjectableMyElement>* el_pt = 
       dynamic_cast<HeleShawFluxElement<ProjectableMyElement>*>(
        Flow_rate_flux_mesh_pt->element_pt(el));
    
      // Add its contribution
      cross_section += el_pt->compute_cross_section();
     }
    return cross_section;
   }

  void fill_in_contribution_to_residuals(Vector<double>& residuals)
   {
    // Calculate the total volume flux across all sub flux meshes
    double volume_flux_ahead_of_bubble =
     calculate_volume_flux_ahead_of_bubble();

    // Calculate the cross section at the bubble inflow
    double bubble_inflow_cross_section =
     calculate_bubble_inflow_cross_section();

    // Get the current x position of the bubble tip
    double current_x_position = Bubble_tip_node_pt->position(0);

    // Get the current transverse displacement of the control node
    double current_w = Displacement_control_node_pt->value(0);

    // Equation for pressure gradient
    int pressure_gradient_eqn=internal_local_eqn(0,0);
    // Equation for bubble pressure
    int bubble_pressure_eqn=internal_local_eqn(1,0);
    // Equation for bubble tip speed
    int bubble_tip_speed_eqn=internal_local_eqn(2,0);
    // Equation for downstream flux
    int downstream_area_eqn=internal_local_eqn(3,0);

    
    if(pressure_gradient_eqn>=0)
     {
      // Get the residuals
      residuals[pressure_gradient_eqn] += 
       volume_flux_ahead_of_bubble
       + this->internal_data_pt(2)->value(0) * 
       this->internal_data_pt(3)->value(0);
     }
    if(bubble_pressure_eqn>=0)
     {
      // Get the residuals
      residuals[bubble_pressure_eqn] += 
       current_x_position - *Prescribed_x_position_pt;
     }
    if(bubble_tip_speed_eqn>=0)
     {
      // Get the residuals
      residuals[bubble_tip_speed_eqn] +=
       this->internal_data_pt(2)->value(0)*bubble_inflow_cross_section - 1.0;
     }
    if(downstream_area_eqn>=0)
     {
      //oomph_info << "Flux dof not pinned!" << std::endl;
      // Get the residuals
      residuals[downstream_area_eqn] +=
       current_w - *Prescribed_displacement_pt;
     }

   }


  // Joao defining the fill in jacobian and mass, in order to use in the eigensolver
 void fill_in_contribution_to_jacobian_and_mass_matrix(Vector<double> &residuals,
                                       DenseMatrix<double> &jacobian,DenseMatrix<double> &mass_matrix)
  { 
  this->fill_in_contribution_to_jacobian(residuals,jacobian);   
  } 
  // Joao defining the fill in jacobian and mass, in order to use in the eigensolver


 private:

  Mesh* Downstream_flux_mesh_pt;
  Mesh* Flow_rate_flux_mesh_pt;
  SolidNode* Bubble_tip_node_pt;
  double* Prescribed_x_position_pt;
  Node* Displacement_control_node_pt;
  double* Prescribed_displacement_pt;
 };    


//=====================================================================
/// A constraint element that constrains the flux downstream (and
/// optionally the inflation upstream).
//=====================================================================
 class FluxConstraintElement : public GeneralisedElement
 {
 public:

  FluxConstraintElement(double& pressure_gradient,
                        double& bubble_pressure,
                        double& downstream_area,                               //Joao: added to implement pressure control
                        double& bubble_tip_speed,
                        double& flow_rate_var,                                 ///Joao: added to keep the buble speed at the desired constant value
                        Mesh* downstream_flux_mesh_pt,
                        Mesh* flow_rate_flux_mesh_pt,
                        SolidNode* bubble_tip_node_pt,
                        double* prescribed_x_position_pt,
                        double* prescribed_area_pt,
                        double* prescribed_bubble_pressure_pt,                 //Joao: added to implement pressure control
                        bool* impose_area_pt,                                  //Joao: added to implement pressure control
                        bool* impose_bubble_speed_pt,                          ///Joao: flag to switch imposed parameter: bubble speed or flow rate
                        double* dimensional_desired_bubble_speed_pt,           ///Joao: added to keep the buble speed at the desired constant value
                        double* dimensional_desired_flow_rate_pt,              ///Joao: added to keep the flow rate at the desired constant value
                        double* undeformed_rectangular_dimensional_area_pt)    ///Joao: added to keep the buble speed at the desired constant value
                                                                        
   {
    // Create internal data object, which is the unknown pressure
    // gradient ahead of the bubble
    this->add_internal_data(new Data(1),true);
    this->internal_data_pt(0)->set_value(0,pressure_gradient);
    
    // The unknown bubble pressure is also internal data of this element
    this->add_internal_data(new Data(1),true);
    this->internal_data_pt(1)->set_value(0,bubble_pressure);
    
    // The unknown bubble tip velocity is also internal data of this element
    this->add_internal_data(new Data(1),true);
    this->internal_data_pt(2)->set_value(0,bubble_tip_speed);

    // The unknown flow rate is also internal data of this element           ///Joao: added to keep the buble speed at the desired constant value
    this->add_internal_data(new Data(1),true);                               ///Joao: added to keep the buble speed at the desired constant value
    this->internal_data_pt(3)->set_value(0,flow_rate_var);                   ///Joao: added to keep the buble speed at the desired constant value

    // The unknown flux downstream is also internal data of this element     
    this->add_internal_data(new Data(1),true);                               //Joao: added to implement pressure control                              
    this->internal_data_pt(4)->set_value(0,downstream_area);                 //Joao: added to implement pressure control

    // Set pointer to the flux mesh ahead of the bubble
    Downstream_flux_mesh_pt = downstream_flux_mesh_pt;
    // Set pointer to the flux mesh that constitutes the bubble inflow
    Flow_rate_flux_mesh_pt = flow_rate_flux_mesh_pt;
    // Set pointer to solid node of which we track x position
    Bubble_tip_node_pt = bubble_tip_node_pt;
    // Set the prescribed x position of the bubble tip
    Prescribed_x_position_pt = prescribed_x_position_pt;
    // Set the prescribed flux ahead of the bubble
    Prescribed_area_pt = prescribed_area_pt;

    Prescribed_bubble_pressure_pt = prescribed_bubble_pressure_pt;                           //Joao: added to implement pressure control          
    Impose_area_pt = impose_area_pt;                                                         //Joao: added to implement pressure control

    Impose_bubble_speed_pt = impose_bubble_speed_pt;                                         ///Joao: flag to switch imposed parameter: bubble speed or flow rate
    Dimensional_desired_bubble_speed_pt = dimensional_desired_bubble_speed_pt;               ///Joao: added to keep the buble speed at the desired constant value
    Dimensional_desired_flow_rate_pt = dimensional_desired_flow_rate_pt;                     ///Joao: added to keep the flow rate at the desired constant value
    Undeformed_rectangular_dimensional_area_pt = undeformed_rectangular_dimensional_area_pt; ///Joao: added to keep the buble speed at the desired constant value

    // Now we need to add external data

    // Looks like we need to add all the nodes in the bulk elements
    // adjacent to the flux mesh
    // obacht WHY?? because we compute gradients!
    unsigned n_el=Downstream_flux_mesh_pt->nelement();
    for(unsigned e=0;e<n_el;e++)
     {     
      HeleShawFluxElement<ProjectableMyElement>* el_pt = 
       dynamic_cast<HeleShawFluxElement<ProjectableMyElement>*>(
        Downstream_flux_mesh_pt->element_pt(e));
      ProjectableMyElement* bulk_el_pt = dynamic_cast<ProjectableMyElement*>(
       el_pt->bulk_element_pt());
      unsigned n_node = bulk_el_pt->nnode();
      for(unsigned inod=0;inod<n_node;inod++)
       {
        this->add_external_data(bulk_el_pt->node_pt(inod),true);
       }
     }

    n_el=Flow_rate_flux_mesh_pt->nelement();
    for(unsigned e=0;e<n_el;e++)
     {     
      HeleShawFluxElement<ProjectableMyElement>* el_pt = 
       dynamic_cast<HeleShawFluxElement<ProjectableMyElement>*>(
        Flow_rate_flux_mesh_pt->element_pt(e));
      unsigned n_node = el_pt->nnode();
      for(unsigned inod=0;inod<n_node;inod++)
       {
        this->add_external_data(el_pt->node_pt(inod),true);
        // Need to add the position pt because changes in bubble width
        // affect flow rate
        SolidNode* solid_nod_pt = dynamic_cast<SolidNode*>(
         el_pt->node_pt(inod));
        this->add_external_data(solid_nod_pt->variable_position_pt(),true);
       }
     }

    // add the node we track as external data (but not the node itself,
    // instead the position data of the node)
    this->add_external_data(Bubble_tip_node_pt->variable_position_pt(),true);
   }
        
  ~FluxConstraintElement(){}

  double calculate_volume_flux_ahead_of_bubble()
   {
    // get the volume flux from the flux elements
    double volume_flux = 0.0;
    unsigned n_element=Downstream_flux_mesh_pt->nelement();
    for(unsigned el=0; el<n_element; el++)
     {
      // Get pointer to flux element
      HeleShawFluxElement<ProjectableMyElement>* el_pt = 
       dynamic_cast<HeleShawFluxElement<ProjectableMyElement>*>(
        Downstream_flux_mesh_pt->element_pt(el));
    
      // Add its contribution
      volume_flux += el_pt->compute_volume_flux();
     }
    return volume_flux;
   }

  double calculate_bubble_inflow_cross_section()
   {
    // get the cross section from the flux elements
    double cross_section = 0.0;
    unsigned n_element=Flow_rate_flux_mesh_pt->nelement();
    for(unsigned el=0; el<n_element; el++)
     {
      // Get pointer to flux element
      HeleShawFluxElement<ProjectableMyElement>* el_pt = 
       dynamic_cast<HeleShawFluxElement<ProjectableMyElement>*>(
        Flow_rate_flux_mesh_pt->element_pt(el));
    
      // Add its contribution
      cross_section += el_pt->compute_cross_section();
     }
    return cross_section;
   }

  void fill_in_contribution_to_residuals(Vector<double>& residuals)
   {
    // Calculate the total volume flux across all sub flux meshes
    double volume_flux_ahead_of_bubble =
     calculate_volume_flux_ahead_of_bubble();

    // Calculate the cross section at the bubble inflow
    double bubble_inflow_cross_section =
     calculate_bubble_inflow_cross_section();

    // Get the current x position of the bubble tip
    double current_x_position = Bubble_tip_node_pt->position(0);

    // Equation for pressure gradient
    int pressure_gradient_eqn=internal_local_eqn(0,0);
    // Equation for bubble pressure
    int bubble_pressure_eqn=internal_local_eqn(1,0);
    // Equation for bubble tip speed
    int bubble_tip_speed_eqn=internal_local_eqn(2,0);
    // Equation for flow rate                                       ///Joao: added to keep the buble speed at the desired constant value
    int flow_rate_var_eqn=internal_local_eqn(3,0);                  ///Joao: added to keep the buble speed at the desired constant value
    // Equation for area downstream                                 
    int downstream_area_eqn=internal_local_eqn(4,0);                //Joao: added to implement pressure control

    
    if(pressure_gradient_eqn>=0)
     {
      // Get the residuals
      residuals[pressure_gradient_eqn] += 
       volume_flux_ahead_of_bubble + this->internal_data_pt(4)->value(0)* 
       this->internal_data_pt(2)->value(0);
     }
    if(bubble_pressure_eqn>=0)
     {
      // Get the residuals
      if((*Impose_area_pt))
       {
        residuals[bubble_pressure_eqn] += 
         current_x_position - *Prescribed_x_position_pt;
       }
      else                                                       //Joao: added to implement pressure control
       {
        residuals[bubble_pressure_eqn] +=
         this->internal_data_pt(1)->value(0) - (*Prescribed_bubble_pressure_pt);
       }
     }
    if(bubble_tip_speed_eqn>=0)
     {
      //oomph_info << "Flux dof not pinned!" << std::endl;
      // Get the residuals
      residuals[bubble_tip_speed_eqn] +=
       this->internal_data_pt(2)->value(0)*bubble_inflow_cross_section - 1.0;
     }
    if(flow_rate_var_eqn>=0)                                    ///Joao: added to keep the buble speed or the flow rate at the desired constant value
     {                                                                               
      if((*Impose_bubble_speed_pt))
       {
        // Get the residuals             ///Imposing bubble speed                                                         
        residuals[flow_rate_var_eqn] +=                                              
         this->internal_data_pt(2)->value(0)-(*Dimensional_desired_bubble_speed_pt)*(*Undeformed_rectangular_dimensional_area_pt)/(this->internal_data_pt(3)->value(0));
        }
      else
       {
        // Get the residuals             ///Imposing flow rate
        residuals[flow_rate_var_eqn] += 
         this->internal_data_pt(3)->value(0) - (*Dimensional_desired_flow_rate_pt);
        }  
      }  
    if(downstream_area_eqn>=0)               //Joao: added to implement pressure control
     {
      // Get the residuals
      if((*Impose_area_pt))
       {
        residuals[downstream_area_eqn] +=
         this->internal_data_pt(4)->value(0) - (*Prescribed_area_pt);
       }
      else
       {
        residuals[downstream_area_eqn] +=
         current_x_position - *Prescribed_x_position_pt;
       }
     }
                                                                         
 }


  // Joao defining the fill in jacobian and jacobian_and_mass, in order to use in the eigensolver
 void fill_in_contribution_to_jacobian_and_mass_matrix(Vector<double> &residuals,
                                       DenseMatrix<double> &jacobian,DenseMatrix<double> &mass_matrix)
  {
   this->fill_in_contribution_to_jacobian(residuals,jacobian);   
//      oomph_info<<"Done fill_in_contribution_to_jacobian_and_mass_matrix at FluxConstraintElement."<<std::endl;
//    exit(1); 
   } 
  // Joao defining the fill in jacobian and jacobian_and_mass, in order to use in the eigensolver


 private:

  Mesh* Downstream_flux_mesh_pt;
  Mesh* Flow_rate_flux_mesh_pt;
  SolidNode* Bubble_tip_node_pt;
  double* Prescribed_x_position_pt;
  double* Prescribed_area_pt;
  double* Prescribed_bubble_pressure_pt;                         //Joao: added to implement pressure control
  bool* Impose_area_pt;                                          //Joao: added to implement pressure control
  bool* Impose_bubble_speed_pt;                                                    ///Joao: flag to switch imposed parameter: bubble speed or flow rate
  double* Dimensional_desired_bubble_speed_pt;                                     ///Joao: added to keep the buble speed at the desired constant value
  double* Dimensional_desired_flow_rate_pt;                                        ///Joao: added to keep the flow rate at the desired constant value
  double* Undeformed_rectangular_dimensional_area_pt;                              ///Joao: added to keep the buble speed at the desired constant value
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
 int Adaptation_frequency=-1;

 /// Initial/max element area
 double Max_el_area=0.04; //0.012 // 0.005;

 /// Min element area for refinement
 double Min_el_area=0.0012; //0.0012 //0.01; // 0.0001;

 /// Max permitted error (spatial)
 double Max_permitted_error=0.0005; //0.003 //0.0005

 /// Min permitted error (spatial)
 double Min_permitted_error=0.00005; //0.0003 //0.0001

 /// Number of polyline segments on bubble. Negative: infer from
 /// max element area
 int N_segment_bubble=600; // 200

 /// Max length of polyline segments on bubble -- negative: don't limit.
 double Max_polyline_length_bubble=0.01;

 /// Output directory
 std::string Directory = "RESLT";

 /// Name of restart file
 std::string Restart_file="";

 ///Physical parameters in the problem

 ///Sheet thickness (in m)
 double h_sheet = 0.34*1.0e-3; ///joao===>>> Lucie's parameter: 0.34*1.0e-3 --- Callum's parameter: 0.46*1.0e-3

 ///Youngs modulus of the sheet (in Pa)
 double E_sheet = 1.44*1.0e6; //3.24*1.0e5; //1.44*1.0e6;

 ///Poissons ratio of the sheet 
 double nu_sheet = 0.5; //0.44;

 ///Bending stiffness of the sheet
 double K_sheet = E_sheet*h_sheet*h_sheet*h_sheet/
  (12.0*(1.0-nu_sheet*nu_sheet));

 ///Initial cell depth (in m)
 double b0_cell = 1.05e-3; //1.05e-3; //90.0e-6; //3.0e-3;

 ///The cell width (in m)
 double W_cell = 30.0*1.0e-3; //0.17;

 ///Fluid viscosity (Pa.s)
 double mu_fluid = 0.099;                 

 ///Surface tension (N/m)
 double gamma_fluid = 21.0*1.0e-3;

 // Hele Shaw parameters
 //---------------------
 /// \short Storage for pointer to Data item whose single value stores
 /// pressure in the bubble
 Data* P_bubble_pt = 0;

 /// Non-dimensional gap width 
 double Aspect_ratio = b0_cell/W_cell;

 /// Storage for pointer to Data item whose single value stores
 /// the bubble tip speed
 Data* U_bubble_tip_pt = 0;

 /// Storage for pointer to Data item whose single value stores
 /// the flow rate variable
 Data* flow_rate_var_pt = 0;                             ///Joao: added to keep the buble speed at the desired constant value      

 // Decide if we will use Flux_constraint_element or Displacement_constraint_element.
 // If we decide to use Flux_constraint_element we still have to decide if we will prescribe 
 // the bubble pressure or the downstream area.
 //Joao: added to implement pressure control
 double Flux_constraint = true;

 //Joao: added to implement pressure control
 double Prescribed_bubble_pressure = 0.0818928;

 /// \short The prescribed volumetric flow rate (in m^3/s) 
 double Prescribed_flow_rate =2.40012*1.0e-6;//177.0156*1.0e-6/60.0; // flow_rate_var_pt->value(0);  ///Joao: added to keep the buble speed at the desired constant value  //////Joao ====>>> Before 187.0*1.0e-6/60.0 

 double Flow_rate_auxiliar_temporary = Prescribed_flow_rate;

 /// Time-scale used in non-dimnesionalisation
 double Time_scale = W_cell*W_cell*b0_cell/Prescribed_flow_rate;

 /// Capillary number 
 double Ca = mu_fluid*W_cell/(gamma_fluid*Time_scale);          /////Joao ==>> using flow rate as a measure of velocity

 ///Inverse Capillary number
 double Ca_inv = 1.0/Ca;

 /// FSI parameter
 double Q = 12.0*mu_fluid*Prescribed_flow_rate/
  (Aspect_ratio*Aspect_ratio*Aspect_ratio*K_sheet);                   /////Joao ==>> using flow rate as a measure of velocity

 /// \short Backup for bubble pressure (to retain value during spatial
 /// adaptation during which volume control element (which stores
 /// this value in its internal Data) is temporarily deleted.
 double P_bubble_backup=0.0;

 double U_bubble_tip_backup=1.0;

 //boolean flag:
 //True  -> Dimensional bubble speed is kept constant
 //False -> Dimensional flow rate is kept constant
 bool Impose_U_bubble_tip = true;

 /// Imposed bubble tip velocity (in m/s)
 /// when it's kept constant
 double U_bubble_tip_dimensional_desired = 0.1;  ///Joao: added to keep the buble speed at the desired constant value

 /// Imposed flow rate (in m^3/s)
 /// when it's kept constant
 double Flow_rate_dimensional_desired = 2.40012*1.0e-6;  ///Joao: added to keep the flow rate at the desired constant value

 /// Area of the undeformed rectangular cross_section (in m^2)
 double undeformed_dimensional_rectangular_cross_section = W_cell*b0_cell;

 double flow_rate_var_backup=2.40012*1.0e-6;        ///Joao: added to keep the buble speed at the desired constant value

 // Foeppl von Karman parameters
 //-----------------------------

 /// Foeppl von Karman parameter 
 double Eta = 12.0*(1.0-nu_sheet*nu_sheet)*(W_cell/h_sheet)*(W_cell/h_sheet);

 /// Initial bubble width (in m)
 double W_bubble_initial_dim = 15.0*1.0e-3;

 /// Initial bubble width
 double W_bubble_initial = W_bubble_initial_dim/W_cell;
 
 /// Amplitude of radius perturbation of the bubble interface 
 double Radius_perturb = W_bubble_initial/100.0;
 
 /// Wavenumber of pressure perturbation at interface
 unsigned N_perturb = 12; 
 
 /// Pseudo-solid Poisson ratio
 double Nu = 0.3; // obacht was 0.3
 
 /// Constitutive law used to determine the mesh deformation
 ConstitutiveLaw *Constitutive_law_pt = 0;

 /// Length of channel downstream
 double L_downstream=15.0;

 /// Length of channel upstream
 double L_upstream=10.0;

 /// Amplitude of pressure perturbation at interface 
 double Pressure_perturb = 0.0;

 /// External pressure
 double P_ext=0.0;

 double Prescribed_bubble_tip_x_position = 0.0;

 Node* Tip_node_pt = 0;
 Node* Top_node_pt = 0;
 Node* Bottom_node_pt = 0;

 Data* Downstream_area_pt = 0;
 double Downstream_area_backup = 1.0;

 double Prescribed_area = 1.0;               /// Joao ===>> before 1.0                          
 bool Prescribe_area = false;                /// Joao ===>> before false                      

 Node* Displacement_control_node_pt = 0;
 double Prescribed_displacement = 0.0;

 Node* Pressure_control_node_pt = 0;
 double Prescribed_pressure = 0.0;

 bool Include_transmural_pressure = false;

 /// Assigns the value of pressure depending on the position
 void external_pressure(const Vector<double>& x, const double& w, 
                        double& pressure)
  {
   double AiA0=-Prescribed_area;
   if(!Prescribe_area)
    {
     AiA0=Downstream_area_pt->value(0);
    }
   double transmural_pressure_dim=
    83.5*pow(AiA0,3.0)-292.4*AiA0*AiA0+469.4*AiA0-270.9;
   if( (Q > 0.0) && Include_transmural_pressure)
    {
     P_ext=W_cell*W_cell*W_cell/K_sheet*transmural_pressure_dim;
    }
   else
    {
     P_ext=0.0;
    }
   pressure=P_ext;
   //oomph_info << "obacht "<<pressure<<" "<<AiA0<<std::endl;
  }

 // Pre stress in the sheet (Pa)
 double Pre_stress = 30.0e3;  ///joao===>>> Lucie's parameter: 30.0e3 --- Callum's parameter: 35.0e3/110.0e3

 void get_pre_stress(DenseMatrix<double>& sigma_0)
 {
  sigma_0(0,0) = 0.0;
  sigma_0(1,1) = Pre_stress/E_sheet;
  sigma_0(0,1) = 0.0;
  sigma_0(1,0) = sigma_0(0,1);
 }

 bool Update_ca_only=true;

 /// Update the dependent parameters
 void update_dependent_parameters()
 {
  double flow_rate_aux;                                            ///Joao: added to keep the buble speed at the desired constant value
  if((flow_rate_var_pt == 0)||(!Flux_constraint))                                        ///Joao: added to keep the buble speed at the desired constant value
  {                                                                ///Joao: added to keep the buble speed at the desired constant value
  flow_rate_aux=flow_rate_var_backup;                              ///Joao: added to keep the buble speed at the desired constant value
   }                                                               ///Joao: added to keep the buble speed at the desired constant value
  else                                                             ///Joao: added to keep the buble speed at the desired constant value
  {                                                                ///Joao: added to keep the buble speed at the desired constant value
  flow_rate_aux=flow_rate_var_pt->value(0);                        ///Joao: added to keep the buble speed at the desired constant value
   }                                                               ///Joao: added to keep the buble speed at the desired constant value
  oomph_info<<"Updating parameters."<<std::endl;                   ///Joao: added to keep the buble speed at the desired constant value
  Time_scale = W_cell*W_cell*b0_cell/flow_rate_aux;                ///Joao: added to keep the buble speed at the desired constant value
  Ca = mu_fluid*W_cell/(gamma_fluid*Time_scale);
  Ca_inv = 1.0/Ca;
  if(!Update_ca_only)
   {
    Q = 12.0*mu_fluid*flow_rate_aux/                               ///Joao: added to keep the buble speed at the desired constant value
     (Aspect_ratio*Aspect_ratio*Aspect_ratio*K_sheet);
   }
  undeformed_dimensional_rectangular_cross_section = W_cell*b0_cell;
 }
 /// \short Pressure acting on air-liquid interface -- basically the
 /// bubble pressure but here modified by sinusoidal perturbation
 /// to force fingering with a particular wavenumber.
 void bubble_pressure_function(const Vector<double>& x, double& pressure)
 {
  if (P_bubble_pt == 0)
   {
    pressure = 0.0;
   }
  else
   {
    pressure = P_bubble_pt->value(0);
    double phi = atan2(x[1], x[0]);
    pressure += Pressure_perturb*sin(double(N_perturb)*phi);
    //oomph_info << "Bubble pressure: "<<pressure<<std::endl;
   }
 }

 void moving_frame_velocity_fct(const Vector<double>&x, Vector<double>& U_frame)
 {
  U_frame[0] = U_bubble_tip_pt->value(0);
  U_frame[1] = 0.0;
 }

 // Pressure gradient ahead of bubble
 Data* Pressure_gradient_downstream_pt = 0;
 double Pressure_gradient_downstream_backup=0.0;

 void normal_flux_ahead_of_bubble(const double& b, const Vector<double>& dbdx,
                                  double& flux)
 {
  /// Ahead of the bubble we have the boundary condition dp/dx = -G.
  /// We need to supply the function flux = b^3 n.grad p = -G b^3.
  double G = 0.0;
  if (Pressure_gradient_downstream_pt!=0)
   {
    G = Pressure_gradient_downstream_pt->value(0);
   }
  flux = -G * (b*b*b);
 }

 double Thin_film_homotopy_parameter = 1.0;

 // /// Thin film effect contribution to the kinematic BC
 // double thin_film_effect_kinematic_bc()
 // {
 //  //if(Apply_thin_film_effects)
 //  // {
 //    return pow(Ca,2.0/3.0)/(0.76+2.16*pow(Ca,2.0/3.0));
 //    // }
 //    //else
 //    // {
 //    //return 0.0;
 //    //}
 // }

 // void thin_film_flux_behind_bubble(const double& b, const Vector<double>& dbdx,
 //                                   double& flux)
 // {
 //  /// Ahead of the bubble we have the boundary condition dp/dx = -G.
 //  /// We need to supply the function flux = b^3 n.grad p = -G b^3.
 //  double combined_gradient = dbdx[0]*dbdx[1];
 //  double sum_square_gradients = dbdx[0]*dbdx[0]+dbdx[1]*dbdx[1];
 //  double G = 0.0;
 //  Vector<double> dummy_x(2,0.0);
 //  Vector<double> U_wall(2,0.0);
 //  moving_frame_velocity_fct(dummy_x,U_wall);
 //  flux = 0.0; //-U_wall[0]*b*thin_film_effect_kinematic_bc();
 // }

 // void normal_flux_behind_bubble(const double& b, const Vector<double>& dbdx,
 //                                  double& flux)
 // {
 //  /// Ahead of the bubble we have the boundary condition dp/dx = -G.
 //  /// We need to supply the function flux = b^3 n.grad p = -G b^3.
 //  Vector<double> dummy_x(2,0.0);
 //  Vector<double> U_wall(2,0.0);
 //  flux = U_wall[0]*b*(1.0-(dbdx[0]*dbdx[0])/
 //                      (2.0*(1.0+(dbdx[0]*dbdx[0]))));
 // }

 // Shape parameter in Saffman-Taylor solution
 double Lambda=W_bubble_initial;

 // Shape of Saffman-Taylor finger
 double x_as_function_of_y(double& y)
 {
  return Prescribed_bubble_tip_x_position + 
   (1.0-Lambda)/MathematicalConstants::Pi * 
   log(1.0/2.0*(1.0+cos(MathematicalConstants::Pi*2.0*y/Lambda)));
 }

 double y_as_function_of_x(double& x)
 {
  return Lambda/(2.0*MathematicalConstants::Pi) * 
   acos(2.0*exp(MathematicalConstants::Pi*
                (x-Prescribed_bubble_tip_x_position)/(1.0-Lambda)) - 1.0);
 }

 bool Area_continuation = false;                               
 double Target_prescribed_area=0.0;
 double Prescribed_area_prev=0.0;

 bool Displacement_continuation = false;
 double Target_prescribed_displacement=0.0;
 double Prescribed_displacement_prev=0.0;

 bool FSI_continuation = false;
 double Target_q = Q;
 double Q_prev=0.0;

 bool Aspect_ratio_continuation = false;
 double Target_aspect_ratio = Aspect_ratio;
 double Aspect_ratio_prev =0.0;

 bool Flow_rate_continuation = false;
 double Target_prescribed_flow_rate=0.0;
 double Prescribed_flow_rate_prev=0.0;

 double Max_increment = 0.0;

 double Amplitude=0.0;

 Vector<double> Vector_of_eigenvalues_rp(20,2.0);
 Vector<double> Vector_of_eigenvalues_ip(20,2.0);

 /// Doc parameters
 void doc_parameters()
 {
  double downstream_area = Prescribed_area;
  if(!Prescribe_area)
   {
    downstream_area = Downstream_area_pt->value(0);
   }
  oomph_info 
   << "Initial width: " 
   << Problem_Parameter::W_bubble_initial << std::endl
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
   << "Number of polyline segments on half of bubble (neg: inferred from element area): "
   << Problem_Parameter::N_segment_bubble << "\n"
   << "Max length of polyline segments on bubble (neg: not constrained): "
   << Problem_Parameter::Max_polyline_length_bubble << "\n"
   << "Adaptation frequency (neg: never): "
   << Problem_Parameter::Adaptation_frequency << "\n"
   << "Pressure gradient downstream: "
   << Problem_Parameter::Pressure_gradient_downstream_pt->value(0) << "\n"
   << "Area downstream: "
   << downstream_area << "\n"
   << "U bubble tip: "<<U_bubble_tip_pt->value(0) << "\n"
   << std::endl;
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

 Problem* make_copy() 
  {
   return (new UnstructuredHSFvKProblem<ELEMENT>());
  }

 /// \short Actions before conv check hierher
 void actions_before_newton_convergence_check()
  {

   if (CommandLineArgs::command_line_flag_has_been_set("--pin_fvk"))
    {
     // Update pinned fvk displacements (mesh moves!)
     update_pinned_fvk_displacements();
    }

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
        
        // Pin pressure in the proper interior of the bubble
        // (pressure on the bubble surface is computed!)
        unsigned nnod = el_pt->nnode();
        for(unsigned j = 0; j < nnod; j++)
         {
          Node* nod_pt = el_pt->node_pt(j);
          if( (!(nod_pt->is_on_boundary(Inner_boundary)  )) )
           {
            // obacht
            nod_pt->set_value(4,Problem_Parameter::P_bubble_pt->value(0));
           }
         }  
       }
     }
   }
   // DenseDoubleMatrix m(ndof());
   // DoubleVector residuals;
   // get_jacobian(residuals, m);
   // char filename[100];
   // ofstream some_file;
   // sprintf(filename,"%s/jacobian%i.dat",
   //        Doc_info.directory().c_str(),
   //        Counter);
   // some_file.open(filename);
   
   // //m.sparse_indexed_output("junk_hsfvk.dat");
   // m.sparse_indexed_output(filename);
   // //m.output("junk_fvk.dat");
   
   // //ofstream dump_residuals("residuals_hsfvk.dat");
   // ofstream dump_residuals;
   // sprintf(filename,"%s/residuals%i.dat",
   //         Doc_info.directory().c_str(),
   //         Counter);
   // dump_residuals.open(filename);
   
   // unsigned n_residuals = ndof();
   // std::cout  << "in main, ndof: " << ndof() << "Number "<< Counter << '\n';
   // for(unsigned i=0;i<n_residuals;i++)
   //  {
   //   dump_residuals << i << ' ' << residuals[i] << '\n';
   //  }
   // dump_residuals.flush();
   // dump_residuals.close();

   // ofstream descr_file;
   // //descr_file.open("most_recent_description.dat");
   // sprintf(filename,"%s/description%i.dat",
   //         Doc_info.directory().c_str(),
   //         Counter);
   // descr_file.open(filename);
   // describe_dofs(descr_file);
   // descr_file.close();
   // Counter++;

   // if(Problem_Parameter::Prescribed_flux > 0.0)
   //  Pressure_gradient_downstream_constraint_element_pt->
   //  set_prescribed_flux(Problem_Parameter::Prescribed_flux);

   // DenseDoubleMatrix m(ndof());
   // DoubleVector residuals;
   // get_jacobian(residuals, m);
   
   // m.sparse_indexed_output("junk_hsfvk.dat");
   // //m.output("junk_fvk.dat");
   
   // ofstream dump_residuals("residuals_hsfvk.dat");
   
   // unsigned n_residuals = ndof();
   // std::cout  << "in main, ndof: " << ndof() << '\n';
   // for(unsigned i=0;i<n_residuals;i++)
   //  {
   //   dump_residuals << i << ' ' << residuals[i] << '\n';
   //  }
   // dump_residuals.flush();
   
   // ofstream descr_file;
   // descr_file.open("most_recent_description.dat");
   // describe_dofs(descr_file);
   // descr_file.close();
   // doc_solution();

   //doc_solution();
   // pause("done doc solution before newton conv check");
  }

 /// \short Actions before adapt. Delete old volume constraint element and
 /// update global mesh
 void actions_before_adapt()
  {
   oomph_info << " Beginning actions before adapt "<<std::endl;
   oomph_info << "Top tail node: "
              <<Problem_Parameter::Top_node_pt->x(0)<<", "
              <<Problem_Parameter::Top_node_pt->x(1)<<", "
              <<Problem_Parameter::Top_node_pt->value(4)<<", "
              <<Problem_Parameter::Top_node_pt->value(0)<<std::endl;
   oomph_info << "Bottom tail node: "
              <<Problem_Parameter::Bottom_node_pt->x(0)<<", "
              <<Problem_Parameter::Bottom_node_pt->x(1)<<", "
              <<Problem_Parameter::Bottom_node_pt->value(4)<<", "
              <<Problem_Parameter::Bottom_node_pt->value(0)<<std::endl;


   // Make backup of surface mesh
   Backed_up_surface_mesh_pt=new BackupMeshForProjection<TElement<1,3> >
    (Free_surface_bc_mesh_pt,Inner_boundary);

   // // hierher
   // {
   // ofstream some_file;
   // char filename[100];
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

   // sprintf(filename,"%s/soln_before_adapt%i.dat",Doc_info.directory().c_str(),
   //         Doc_info.number());
   // some_file.open(filename);
   // some_file.precision(20);
   // this->Bulk_mesh_pt->output(some_file,5); 
   // some_file.close();
   // }


   if (!CommandLineArgs::command_line_flag_has_been_set("--pin_hs"))
    {
     // Kill the Hele Shaw free surface elements
     delete_surface_elements();

     // Kill the Hele Shaw flux elements
     delete_outflow_elements();
     delete_top_inflow_elements();
     delete_bottom_inflow_elements();
     delete_actual_inflow_elements();
    }

   // Kill the element that enforces constraints
   delete_constraint_element();

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
   // reset the pointer to the tip node
   set_tip_node_pt();
   set_displacement_control_node_pt();
   set_pressure_control_node_pt();

   if (!CommandLineArgs::command_line_flag_has_been_set("--pin_hs"))
    {
     // Create the interface elements along the bubble boundary
     create_surface_elements();

     // Create the flux elements
     create_actual_inflow_elements();
     create_top_inflow_elements();
     create_bottom_inflow_elements();
     create_outflow_elements();
    }
   // create the element that enforces constraints
   create_constraint_element();

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
   complete_problem_setup();
   rebuild_global_mesh();

   // // hierher
   // {
   //  ofstream some_file;
   //  char filename[100];
   //  // sprintf(filename,"%s/gauss_point_data_after_adapt%i.dat",
   //  //         Doc_info.directory().c_str(),
   //  //         Doc_info.number());
   //  // some_file.open(filename);
   //  // unsigned nel = Bulk_mesh_pt->nelement();
   //  // for (unsigned e = 0; e < nel; e++)
   //  //  {
   //  //   dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e))->
   //  //    output_at_gauss_points(some_file);
   //  //  }
   //  // some_file.close();
    
   //  sprintf(filename,"%s/soln_after_adapt%i.dat",Doc_info.directory().c_str(),
   //          Doc_info.number());
   //  some_file.open(filename);
   //  some_file.precision(20);
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

   if(Is_steady)
    {
     // Compute the actual flux ahead of the bubble
     double area_ahead_of_bubble = 0.0;
     {
      unsigned n_element = Outflow_mesh_pt->nelement();
      for(unsigned e = 0; e < n_element; e++)
       {
        // Upcast from GeneralisedElement to the present element
        HeleShawFluxElement<ELEMENT>* el_pt = 
         dynamic_cast<HeleShawFluxElement<ELEMENT>* >(
          Outflow_mesh_pt->element_pt(e));
        
        area_ahead_of_bubble += el_pt->compute_cross_section();
       }
     }
     
     // Compute the actual flux behind the bubble top
     double flux_behind_bubble_top = 0.0;
     {
      unsigned n_element = Top_inflow_mesh_pt->nelement();
      for(unsigned e = 0; e < n_element; e++)
       {
        // Upcast from GeneralisedElement to the present element
        HeleShawFluxElement<ELEMENT>* el_pt = 
         dynamic_cast<HeleShawFluxElement<ELEMENT>* >(
          Top_inflow_mesh_pt->element_pt(e));
        
        flux_behind_bubble_top += el_pt->compute_volume_flux();
       }
     }
     
     // Compute the actual flux behind the bubble bottom
     double flux_behind_bubble_bottom = 0.0;
     {
      unsigned n_element = Bottom_inflow_mesh_pt->nelement();
      for(unsigned e = 0; e < n_element; e++)
       {
        // Upcast from GeneralisedElement to the present element
        HeleShawFluxElement<ELEMENT>* el_pt = 
         dynamic_cast<HeleShawFluxElement<ELEMENT>* >(
          Bottom_inflow_mesh_pt->element_pt(e));
        
        flux_behind_bubble_bottom += el_pt->compute_volume_flux();
       }
     }
     
     // Compute the actual flow rate
     double flow_rate = 0.0;
     {
      unsigned n_element = Actual_inflow_mesh_pt->nelement();
      for(unsigned e = 0; e < n_element; e++)
       {
        // Upcast from GeneralisedElement to the present element
        HeleShawFluxElement<ELEMENT>* el_pt = 
         dynamic_cast<HeleShawFluxElement<ELEMENT>* >(
          Actual_inflow_mesh_pt->element_pt(e));
        
        flow_rate += el_pt->compute_cross_section();
       }
     }
     flow_rate *= Problem_Parameter::U_bubble_tip_pt->value(0); 
     
     double prescribed_area = Problem_Parameter::Prescribed_area;
     if(!Problem_Parameter::Prescribe_area)
      {
       oomph_info << "Before read downstream area pointer"<<std::endl;
       prescribed_area = Problem_Parameter::Downstream_area_pt->value(0);
       oomph_info << "After read downstream area pointer"<<std::endl;
      }
     
     oomph_info << "Displacement of control node: Prescribed: "
                << Problem_Parameter::Prescribed_displacement 
                << ", actual: "
                <<Problem_Parameter::Displacement_control_node_pt->value(0)
                <<std::endl;
     oomph_info << "Area ahead of bubble: Prescribed: "
                << prescribed_area
                << ", actual: "<<area_ahead_of_bubble<<std::endl;
     oomph_info << "Flow rate into bubble: Prescribed: "
                << Problem_Parameter::Prescribed_flow_rate
                << ", actual: "<<flow_rate<<std::endl;
     oomph_info << "Flux behind bubble: Top: "
                << flux_behind_bubble_top
                << ", bottom: "<<flux_behind_bubble_bottom<<std::endl;
     oomph_info << "Tip node position: "
                << Problem_Parameter::Tip_node_pt->x(0)<<", "
                << Problem_Parameter::Tip_node_pt->x(1)<<std::endl;
     oomph_info << "Prescribed x position of tip node: "
                << Problem_Parameter::Prescribed_bubble_tip_x_position<<std::endl;
     oomph_info << "Pressure gradient downstream: "
                << Problem_Parameter::Pressure_gradient_downstream_pt->value(0)
                <<std::endl;
     oomph_info << "Bubble pressure: "
                << Problem_Parameter::P_bubble_pt->value(0) <<std::endl;
     oomph_info << "Top tail node: "
                <<Problem_Parameter::Top_node_pt->x(0)<<", "
                <<Problem_Parameter::Top_node_pt->x(1)<<", "
                <<Problem_Parameter::Top_node_pt->value(4)<<", "
                <<Problem_Parameter::Top_node_pt->value(0)<<std::endl;
     oomph_info << "Bottom tail node: "
                <<Problem_Parameter::Bottom_node_pt->x(0)<<", "
                <<Problem_Parameter::Bottom_node_pt->x(1)<<", "
                <<Problem_Parameter::Bottom_node_pt->value(4)<<", "
                <<Problem_Parameter::Bottom_node_pt->value(0)<<std::endl;
     
     // obacht
//     oomph_info << "obacht setting prescribed area ahead of the bubble to the "
//                << "actual value!" << std::endl;    

//     if(Problem_Parameter::Prescribe_area)
//      {
//       Problem_Parameter::Prescribed_area = area_ahead_of_bubble;
//      }
//     else
//      {
//       Problem_Parameter::Downstream_area_pt->set_value(0,area_ahead_of_bubble);
//       Problem_Parameter::Prescribed_displacement = 
//        Problem_Parameter::Displacement_control_node_pt->value(0);
//      }
     
      Problem_Parameter::Prescribed_displacement = 
        Problem_Parameter::Displacement_control_node_pt->value(0);

//     oomph_info << "obacht setting prescribed flow rate into the bubble to the "
//                << "actual value!" << std::endl;
//     Problem_Parameter::Prescribed_flow_rate = flow_rate;
    }
   //doc_solution();
  }
 
 /// Update after solve (empty)
 void actions_after_newton_solve() 
  {}

 /// Update the problem specs before solve (empty)
 void actions_before_newton_solve()                         
  {                                                         
   Newton_steps_taken = 0;                                  
   if(Update_dependent_parameters)                          
    {                                                       //
     Problem_Parameter::update_dependent_parameters();      ///Joao added line
    }                                                       // 
  }

 void actions_after_newton_step()
  {
   // Newton_steps_taken+=1;

   // DenseDoubleMatrix m(ndof());
   // DoubleVector residuals;
   // get_jacobian(residuals, m);
   // char filename[100];
   // ofstream some_file;
   // sprintf(filename,"%s/jacobian%i.dat",
   //        Doc_info.directory().c_str(),
   //        Counter);
   // some_file.open(filename);
   
   // //m.sparse_indexed_output("junk_hsfvk.dat");
   // m.sparse_indexed_output(filename);
   // //m.output("junk_fvk.dat");
   // some_file.close();
   
   // //ofstream dump_residuals("residuals_hsfvk.dat");
   // ofstream dump_residuals;
   // sprintf(filename,"%s/residuals%i.dat",
   //         Doc_info.directory().c_str(),
   //         Counter);
   // dump_residuals.open(filename);
   
   // unsigned n_residuals = ndof();
   // std::cout  << "in main, ndof: " << ndof() << "Number "<< Counter << '\n';
   // for(unsigned i=0;i<n_residuals;i++)
   //  {
   //   dump_residuals << i << ' ' << residuals[i] << '\n';
   //  }
   // dump_residuals.flush();
   // dump_residuals.close();

   // ofstream descr_file;
   // //descr_file.open("most_recent_description.dat");
   // sprintf(filename,"%s/description%i.dat",
   //         Doc_info.directory().c_str(),
   //         Counter);
   // descr_file.open(filename);
   // describe_dofs(descr_file);
   // descr_file.close();
   // Counter++;

   // doc_solution();
  }

 void actions_after_parameter_increase(double* const &parameter_pt)
  {
   if(Problem_Parameter::Displacement_continuation)
    {
     if(std::fabs(*parameter_pt - 
                  Problem_Parameter::Prescribed_displacement_prev) > 
        std::fabs(Problem_Parameter::Max_increment))
      {
       *parameter_pt = 
        Problem_Parameter::Prescribed_displacement_prev+
        Problem_Parameter::Max_increment;
      }
     if( (Problem_Parameter::Max_increment > 0.0) &&
         (*parameter_pt > 
          Problem_Parameter::Target_prescribed_displacement) )
      {
       *parameter_pt = 
        Problem_Parameter::Target_prescribed_displacement;
      }
     if( (Problem_Parameter::Max_increment < 0.0) &&
         (*parameter_pt < 
          Problem_Parameter::Target_prescribed_displacement) )
      {
       *parameter_pt = 
        Problem_Parameter::Target_prescribed_displacement;
      }
    }

   if(Problem_Parameter::Area_continuation)
    {
     if(std::fabs(*parameter_pt - 
                  Problem_Parameter::Prescribed_area_prev) > 
        std::fabs(Problem_Parameter::Max_increment))
      {
       *parameter_pt = 
        Problem_Parameter::Prescribed_area_prev+
        Problem_Parameter::Max_increment;
      }
     if( (Problem_Parameter::Max_increment > 0.0) &&
         (*parameter_pt > 
          Problem_Parameter::Target_prescribed_area) )
      {
       *parameter_pt = 
        Problem_Parameter::Target_prescribed_area;
      }
     if( (Problem_Parameter::Max_increment < 0.0) &&
         (*parameter_pt < 
          Problem_Parameter::Target_prescribed_area) )
      {
       *parameter_pt = 
        Problem_Parameter::Target_prescribed_area;
      }
    }

   if(Problem_Parameter::Aspect_ratio_continuation)
    {
     if(std::fabs(*parameter_pt - 
                  Problem_Parameter::Aspect_ratio_prev) > 
        std::fabs(Problem_Parameter::Max_increment))
      {
       *parameter_pt = Problem_Parameter::Aspect_ratio_prev+
        Problem_Parameter::Max_increment;
      }
     if( (Problem_Parameter::Max_increment > 0.0) &&
         (*parameter_pt > 
          Problem_Parameter::Target_aspect_ratio) )
      {
       *parameter_pt = Problem_Parameter::Target_aspect_ratio;
      }
     if( (Problem_Parameter::Max_increment < 0.0) &&
         (*parameter_pt < 
          Problem_Parameter::Target_aspect_ratio) )
      {
       *parameter_pt = Problem_Parameter::Target_aspect_ratio;
      }
    }

   if(Problem_Parameter::FSI_continuation)
    {
     if(std::fabs(*parameter_pt - 
                  Problem_Parameter::Q_prev) > 
        std::fabs(Problem_Parameter::Max_increment))
      {
       *parameter_pt = Problem_Parameter::Q_prev+
        Problem_Parameter::Max_increment;
      }
     if( (Problem_Parameter::Max_increment > 0.0) &&
         (*parameter_pt > Problem_Parameter::Target_q) )
      {
       *parameter_pt = Problem_Parameter::Target_q;
      }
     if( (Problem_Parameter::Max_increment < 0.0) &&
         (*parameter_pt < Problem_Parameter::Target_q) )
      {
       *parameter_pt = Problem_Parameter::Target_q;
      }
    }

   if(Problem_Parameter::Flow_rate_continuation)
    {
     if(std::fabs(*parameter_pt - 
                  Problem_Parameter::Prescribed_flow_rate_prev) > 
        std::fabs(Problem_Parameter::Max_increment))
      {
       *parameter_pt = Problem_Parameter::Prescribed_flow_rate_prev+
        Problem_Parameter::Max_increment;
      }
     if( (Problem_Parameter::Max_increment > 0.0) &&
         (*parameter_pt > 
          Problem_Parameter::Target_prescribed_flow_rate) )
      {
       *parameter_pt = 
        Problem_Parameter::Target_prescribed_flow_rate;
      }
     if( (Problem_Parameter::Max_increment < 0.0) &&
         (*parameter_pt < 
          Problem_Parameter::Target_prescribed_flow_rate) )
      {
       *parameter_pt = 
        Problem_Parameter::Target_prescribed_flow_rate;
      }
    }

   if(Update_dependent_parameters)
    {
     Problem_Parameter::update_dependent_parameters();
    }

   oomph_info << "obacht "<<*parameter_pt<<std::endl;
  }

  void actions_before_newton_step()                          ///Joao: added to keep the buble speed at the desired constant value
   {                                                         ///Joao: added to keep the buble speed at the desired constant value
    if(Update_dependent_parameters)                          ///Joao: added to keep the buble speed at the desired constant value
     {                                                       ///Joao: added to keep the buble speed at the desired constant value  
      Problem_Parameter::update_dependent_parameters();      ///Joao: added to keep the buble speed at the desired constant value
     }                                                       ///Joao: added to keep the buble speed at the desired constant value
    }                                                        ///Joao: added to keep the buble speed at the desired constant value
 // void actions_before_newton_step()
 //  {
 //   if(Problem_Parameter::Displacement_continuation)
 //    {
 //     if(std::fabs(Problem_Parameter::Prescribed_displacement - 
 //                  Problem_Parameter::Prescribed_displacement_prev) > 
 //        std::fabs(Problem_Parameter::Max_increment))
 //      {
 //       Problem_Parameter::Prescribed_displacement = 
 //        Problem_Parameter::Prescribed_displacement_prev+
 //        Problem_Parameter::Max_increment;
 //      }
 //     if( (Problem_Parameter::Max_increment > 0.0) &&
 //         (Problem_Parameter::Prescribed_displacement > 
 //          Problem_Parameter::Target_prescribed_displacement) )
 //      {
 //       Problem_Parameter::Prescribed_displacement = 
 //        Problem_Parameter::Target_prescribed_displacement;
 //      }
 //     if( (Problem_Parameter::Max_increment < 0.0) &&
 //         (Problem_Parameter::Prescribed_displacement < 
 //          Problem_Parameter::Target_prescribed_displacement) )
 //      {
 //       Problem_Parameter::Prescribed_displacement = 
 //        Problem_Parameter::Target_prescribed_displacement;
 //      }
 //    }

 //   if(Problem_Parameter::Aspect_ratio_continuation)
 //    {
 //     if(std::fabs(Problem_Parameter::Aspect_ratio - 
 //                  Problem_Parameter::Aspect_ratio_prev) > 
 //        std::fabs(Problem_Parameter::Max_increment))
 //      {
 //       Problem_Parameter::Aspect_ratio = Problem_Parameter::Aspect_ratio_prev+
 //        Problem_Parameter::Max_increment;
 //      }
 //     if( (Problem_Parameter::Max_increment > 0.0) &&
 //         (Problem_Parameter::Aspect_ratio > 
 //          Problem_Parameter::Target_aspect_ratio) )
 //      {
 //       Problem_Parameter::Aspect_ratio = Problem_Parameter::Target_aspect_ratio;
 //      }
 //     if( (Problem_Parameter::Max_increment < 0.0) &&
 //         (Problem_Parameter::Aspect_ratio < 
 //          Problem_Parameter::Target_aspect_ratio) )
 //      {
 //       Problem_Parameter::Aspect_ratio = Problem_Parameter::Target_aspect_ratio;
 //      }
 //    }

 //   if(Problem_Parameter::FSI_continuation)
 //    {
 //     if(std::fabs(Problem_Parameter::Q - 
 //                  Problem_Parameter::Q_prev) > 
 //        std::fabs(Problem_Parameter::Max_increment))
 //      {
 //       Problem_Parameter::Q = Problem_Parameter::Q_prev+
 //        Problem_Parameter::Max_increment;
 //      }
 //     if( (Problem_Parameter::Max_increment > 0.0) &&
 //         (Problem_Parameter::Q > Problem_Parameter::Target_q) )
 //      {
 //       Problem_Parameter::Q = Problem_Parameter::Target_q;
 //      }
 //     if( (Problem_Parameter::Max_increment < 0.0) &&
 //         (Problem_Parameter::Q < Problem_Parameter::Target_q) )
 //      {
 //       Problem_Parameter::Q = Problem_Parameter::Target_q;
 //      }
 //    }

 //   if(Problem_Parameter::Bubble_speed_continuation)
 //    {
 //     if(std::fabs(Problem_Parameter::U_bubble_tip - 
 //                  Problem_Parameter::U_bubble_tip_prev) > 
 //        std::fabs(Problem_Parameter::Max_increment))
 //      {
 //       Problem_Parameter::U_bubble_tip = Problem_Parameter::U_bubble_tip_prev+
 //        Problem_Parameter::Max_increment;
 //      }
 //     if( (Problem_Parameter::Max_increment > 0.0) &&
 //         (Problem_Parameter::U_bubble_tip > 
 //          Problem_Parameter::Target_U_bubble_tip) )
 //      {
 //       Problem_Parameter::U_bubble_tip = 
 //        Problem_Parameter::Target_U_bubble_tip;
 //      }
 //     if( (Problem_Parameter::Max_increment < 0.0) &&
 //         (Problem_Parameter::U_bubble_tip < 
 //          Problem_Parameter::Target_U_bubble_tip) )
 //      {
 //       Problem_Parameter::U_bubble_tip = 
 //        Problem_Parameter::Target_U_bubble_tip;
 //      }
 //    }

 //   if(Update_dependent_parameters)
 //    {
 //     Problem_Parameter::update_dependent_parameters();
 //    }
 //  }

 /// Restart
 void restart(bool unsteady_initial_config_steady)
  {
  
   // obacht required for restart

   if(unsteady_initial_config_steady)
   {
//   this->add_time_stepper_pt(&Problem::Continuation_time_stepper);                //Joao added line
   this->set_timestepper_for_all_data(&Problem::Continuation_time_stepper);       //Joao added line
//   this->set_timestepper_for_all_data(this->time_stepper_pt(0));
    }
  
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
 void unpin_fvk()
  {
   // Loop over all elements
   unsigned nel = Bulk_mesh_pt->nelement();
   for (unsigned e = 0; e < nel; e++)
    {
     ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));
     el_pt->unpin_fvk();
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
   
   // // Set all Foeppl von Karmann displacements
   // double fraction=0.1;
   // double shift=-5.0;
   // unsigned n_node = Bulk_mesh_pt->nnode();
   // for(unsigned i = 0; i < n_node; i++)
   //  {
   //   Node* nod_pt = Bulk_mesh_pt->node_pt(i);
   //   double x=nod_pt->x(0);
   //   double y=nod_pt->x(1);
   //   //double w=Problem_Parameter::Slope*(x-Problem_Parameter::L_downstream);
   //   double w=-Problem_Parameter::Amplitude*
   //    (cos(MathematicalConstants::Pi*
   //         (-x-shift-fraction*Problem_Parameter::L_downstream)/
   //         (fraction*Problem_Parameter::L_upstream+
   //          fraction*Problem_Parameter::L_downstream))-1.0);
   //   if(-x > fraction*Problem_Parameter::L_downstream+shift)
   //    {
   //     w=0.0;
   //    }
   //   if(-x < -fraction*Problem_Parameter::L_upstream+shift)
   //    {
   //     w=2.0*Problem_Parameter::Amplitude;
   //    }
   //   w*=(1.0-4.0*y*y);
   //   nod_pt->set_value(0,w);

   //   if(nod_pt == Problem_Parameter::Displacement_control_node_pt)
   //    {
   //     oomph_info<<"obacht setting displacement"<<std::endl;
   //     Problem_Parameter::Prescribed_displacement=w;
   //    }
   //   else
   //    {
   //     nod_pt->pin(0);
   //    }
   //  }

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

 /// \short Helper function to (re-)set boundary condition
 /// and complete the build of  all elements
 void complete_problem_setup();

 /// Apply boundary conditions
 void apply_boundary_conditions()
  {
   // dynamic_cast<SolidNode*>(Problem_Parameter::Displacement_control_node_pt)
   //  ->pin_position(0);
   // dynamic_cast<SolidNode*>(Problem_Parameter::Displacement_control_node_pt)
   //  ->pin_position(1);

   // dynamic_cast<SolidNode*>(Problem_Parameter::Pressure_control_node_pt)
   //  ->pin_position(0);
   // dynamic_cast<SolidNode*>(Problem_Parameter::Pressure_control_node_pt)
   //  ->pin_position(1);

   // Problem_Parameter::Top_node_pt->pin(6);
   // Problem_Parameter::Top_node_pt->set_value(6,0.0);
   // Problem_Parameter::Bottom_node_pt->pin(6);
   // Problem_Parameter::Bottom_node_pt->set_value(6,0.0);

   // Problem_Parameter::Pressure_control_node_pt->pin(4);
   // Problem_Parameter::Pressure_control_node_pt->set_value(4,0.0);

   // Set the boundary conditions for problem: All nodes are
   // free by default -- just pin the ones that have Dirichlet conditions
   // here. 
   //bool pressure_dof_pinned = false;
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
         // Pin all FvK displacement values on top and bottom boundaries
         if(ibound == Bottom_boundary ||
            ibound == Top_boundary)
          { 
           // Loop over current and previous timesteps
           unsigned ntime=time_stepper_pt()->ndt();
           for (unsigned t=0;t<ntime;t++)
            {   
             // Displacement
             nod_pt->pin(0);
             nod_pt->set_value(t,0,0.0);
             nod_pt->pin(2);
             nod_pt->set_value(t,2,0.0);
             nod_pt->pin(3);
             nod_pt->set_value(t,3,0.0);
            }
           
           solid_node_pt->pin_position(0);
           solid_node_pt->pin_position(1);
          }

         // Pin normal and tangential FvK displacements on in- and outflow
         // boundaries
         if ( ibound == Top_inflow_boundary ||
              ibound == Actual_inflow_boundary ||
              ibound == Bottom_inflow_boundary  || 
              ibound == Outflow_boundary )
          {
           // Loop over current and previous timesteps
           unsigned ntime=time_stepper_pt()->ndt();
           for (unsigned t=0;t<ntime;t++)
            {   
             // Displacement
             nod_pt->pin(2);
             nod_pt->set_value(t,2,0.0);
             nod_pt->pin(3);
             nod_pt->set_value(t,3,0.0);
            }

           solid_node_pt->pin_position(0);

           if(ibound == Outflow_boundary)
            {
             solid_node_pt->pin_position(1);
            }

           // obacht
           // if(ibound == Top_inflow_boundary ||
           //    ibound == Bottom_inflow_boundary)
           //  {
           //   nod_pt->pin(4);
           //  }
          }   
        }
       if(ibound == Dummy_boundary)
        {
         solid_node_pt->pin_position(0);
         solid_node_pt->pin_position(1);
        }
      }   
    } // end loop over boundaries
  }

 void enable_parameter_update()
  {
   Update_dependent_parameters=true;
  }

 void disable_parameter_update()
  {
   Update_dependent_parameters=false;
  }

 void make_steady()
  {
   Is_steady = true;
   this->time_stepper_pt()->make_steady();
   Use_continuation_timestepper=true;
   //Use_finite_differences_for_continuation_derivatives = true;
   Bifurcation_detection = true;
  }

 void make_unsteady()
  {
   Is_steady = false;
   Use_continuation_timestepper=false;
   //Use_finite_differences_for_continuation_derivatives = true;
   Bifurcation_detection = false;
  }


// joao  
 void doc_eigen_values(unsigned n_eval, Vector<complex<double> >  eigenvalues, int sensible_eigenvalues, int Number_of_negative_eigenvalues, unsigned perturbation_evec, bool retain_perturbed_solution_at_function_exit)
 {
    ofstream some_file;
    char filename[100];
    
    sprintf(filename,"%s/eigenvalues%i.dat",Doc_info.directory().c_str(),
            Doc_info.number());
    some_file.open(filename);
    
    some_file << "There are " << sensible_eigenvalues << " reasonable eigenvalues" << '\n';
    some_file << "There are " << Number_of_negative_eigenvalues << " negative eigenvalues" << '\n';
    
    // Loop over the eigenvalues
    for(unsigned e = 0; e < n_eval; e++)
     {
        some_file << eigenvalues[e] << " " << '\n';
     }
    if(retain_perturbed_solution_at_function_exit)
     {
        some_file << " " << '\n'; 
        some_file << eigenvalues[perturbation_evec] << " Kept as a perturbation " << '\n';
     }
    some_file.close();
 }
//joao



 void solve_for_eigenproblem()
  {
   /// Will doc perturbed solutions for 10 eigenvectors
   bool doc_eigenvector = false;                                        
   /// System state at function exit will be perturbed by eigenvector 0
   /// Good for finite amplitude calculations
   bool retain_perturbed_solution_at_function_exit = false; 
   
   /// If we want to output solutions, we should backup the initial state.
   unsigned n_dof = ndof();
   Vector<double> backup(n_dof);
   // Store original state.
   for(unsigned n=0; n<n_dof; n++)
    {
     backup[n] = dof(n);
    }
   
   if (doc_eigenvector == true)
    {
     std::cout << "Doc initial solution" << std::endl;
     doc_solution();
    }
   /// Reset eigenvalue guesses.
   for(unsigned n=0; n<20; n++)
    {
     Problem_Parameter::Vector_of_eigenvalues_rp[n]=20.0;
     Problem_Parameter::Vector_of_eigenvalues_ip[n]=20.0;
    }
   
   //Set external storage for the eigenvalues
   Vector<complex<double> > eigenvalues;
   Vector<DoubleVector> eigenvectors;
   
   //Desired number of eigenvalues
   unsigned n_eval = 20;
   int sensible_eigenvalues =0;
   Vector<unsigned> position_of_unstable_eigenvalues(n_eval,n_eval+1);
   //Solve the eigenproblem
   std::cout << "Now attempt to solve eigenproblem, n_eval = " 
             << n_eval << std::endl;


//   exit(1);     //Joao testing --->>> fine       
          
   this->eigen_solver_pt()->solve_eigenproblem(this,n_eval,eigenvalues,eigenvectors);
//Joao   ===>>> before: this->eigen_solver_pt().solve_eigenproblem(this,n_eval,eigenvalues,eigenvectors);

//   exit(1);     //Joao testing --->>> fine

   std::cout << "N_eval is " << n_eval << std::endl;
   std::cout << "Eigenvalues_size is " << eigenvalues.size() << std::endl;
   
   /// Describe eigenvalues
   std::cout << "And describe eigenvalues" << std::endl;
   n_eval = eigenvalues.size();
   int Number_of_negative_eigenvalues = 0;
   for (unsigned n=0; n<n_eval; n++)
    {
     if (isinf(real(eigenvalues[n]))!=true && 
         isnan(real(eigenvalues[n]))!=true)
      {
       std::cout << "Eigenvalue " << eigenvalues[n]  << std::endl;
       sensible_eigenvalues++;
       if (real(eigenvalues[n])<0)
        {
         Number_of_negative_eigenvalues++;
         position_of_unstable_eigenvalues[n]=n;
        }
      }
     if (doc_eigenvector == true)
      {
       std::cout << "Doc eigenvector with eigenvalue: " 
                 << real(eigenvalues[n]) << " " 
                 << imag(eigenvalues[n]) << " " << std::endl;
       /// The function assign_eigenvector_to_dofs(eigenvectors[n]) 
       /// is not helpful for finding the perturbed interface position.
       
       /// Perturb original solution. Choice of factor 10 is arbitrary.
       for(unsigned i=0; i<n_dof; i++)
        {
         dof(i) = backup[i] + 1.0*eigenvectors[n][i];
        }
       actions_after_change_in_bifurcation_parameter();
       doc_solution();
      }
    }
   
   for(unsigned n=0; n<n_dof; n++)
    {
     dof(n) = backup[n];
    }
   actions_after_change_in_bifurcation_parameter();
   
   for (unsigned n=0; n<20; n++)
    {
     if (n_eval>=n+1)
      {
       Problem_Parameter::Vector_of_eigenvalues_rp[n]=real(eigenvalues[n]);
       Problem_Parameter::Vector_of_eigenvalues_ip[n]=imag(eigenvalues[n]);
      }
    }
   
   std::cout << "There are " << sensible_eigenvalues 
             << " reasonable eigenvalues" << std::endl;

   unsigned perturbation_evec = 0;
   if (retain_perturbed_solution_at_function_exit==true)
    {
     std::cout << "Keep perturbation with eigenvectors" << std::endl; 
     
     Vector<double> Sum_of_eigenvectors(n_dof,0.0);

     bool keep_real_stable_eigenvector=false;
     bool keep_real_unstable_eigenvector=true;
     bool keep_unstable_eigenvector=false;

     if(keep_real_stable_eigenvector==true)
       {
        for(unsigned n=0; n<=n_eval; n++)
           {
            if( imag(eigenvalues[n])==0.0 )
              {
              perturbation_evec=n;
              }
           }
        for(unsigned n=0; n<=n_eval; n++)
           {
            if( (imag(eigenvalues[n])==0.0)&&(real(eigenvalues[n])>real(eigenvalues[perturbation_evec])) )
              {
              perturbation_evec=n;
              }
           }
       }

     if(keep_real_unstable_eigenvector==true)
       {
        for(unsigned n=0; n<=n_eval; n++)
           {
            if( imag(eigenvalues[n])==0.0 )
              {
              perturbation_evec=n;
              }
           }
        for(unsigned n=0; n<=n_eval; n++)
           {
            if( (imag(eigenvalues[n])==0.0)&&(real(eigenvalues[n])<real(eigenvalues[perturbation_evec])) )
              {
              perturbation_evec=n;
              }
           }
       }

     bool first_position = true;
     if(keep_unstable_eigenvector==true)
       {
        for(unsigned n=0; n<=n_eval; n++)
         {
          if( (position_of_unstable_eigenvalues[n]<n_eval)&&(first_position) )
            {
             perturbation_evec = n;
             first_position = false;
            }
         }
        }



     for(unsigned i=0; i<n_dof; i++)
      {
//If we want to add all the eigenvalues to the solution
//       for(unsigned n_evec=0; n_evec<n_eval; n_evec++)
//        {
//          Sum_of_eigenvectors[i] += eigenvectors[n_evec][i];      
//         }
       dof(i) = backup[i] + 2.0*eigenvectors[perturbation_evec][i];
      }
     actions_after_change_in_bifurcation_parameter();
     doc_solution();
     
    }
   doc_eigen_values(n_eval,eigenvalues,sensible_eigenvalues, Number_of_negative_eigenvalues,perturbation_evec,retain_perturbed_solution_at_function_exit);
   doc_solution();
  }



private:

 unsigned Counter;

 double Next_dt;

 /// Doc info object for labeling output
 DocInfo Doc_info;

 /// Use interface motion to estimate temporal error
 bool Use_interface_motion_for_temporal_error;

 /// Mesh as geom object representation (updated after every adaptation)
 MeshAsGeomObject* Mesh_as_geom_object_pt;

 /// Pointers to bulk mesh
 RefineableSolidTriangleMesh<ELEMENT>* Bulk_mesh_pt;

 /// Mesh for the position constraint element
 Mesh *Constraint_mesh_pt;

 /// Pointers to mesh of flux BC elements
 Mesh* Free_surface_bc_mesh_pt;

 /// Single position constraint element instance
 DisplacementConstraintElement* Displacement_constraint_element_pt;
 FluxConstraintElement* Flux_constraint_element_pt;

 /// Vector of the regions comprising the bubble
 Vector<unsigned> Bubble_regions;

 /// Trace file to document norm of solution
 ofstream Trace_file;

 Mesh* Outflow_mesh_pt;

 Mesh* Actual_inflow_mesh_pt;

 Mesh* Top_inflow_mesh_pt;

 Mesh* Bottom_inflow_mesh_pt;

 /// \short Backup of Free_surface_bc_mesh_pt so the Lagrange multipliers
 /// and tangents can be projected across
 BackupMeshForProjection<TElement<1,3> >*  Backed_up_surface_mesh_pt;

 // The actual Newton steps
 unsigned Newton_steps_taken;

 // Flag that indicates whether we want to update parameters or not
 bool Update_dependent_parameters;

 // Flag indicating whether we are after steady solutions or not
 bool Is_steady;

 /// Enum for boundary ids
 enum
 {
  Bottom_boundary=0,
  Outflow_boundary=1,
  Top_boundary=2,
  Top_inflow_boundary=3,
  Actual_inflow_boundary=4,
  Bottom_inflow_boundary=5,
  Inner_boundary = 6,
  Dummy_boundary = 7
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
     
     double x = nod_pt->x(0);
     double y = nod_pt->x(1);

     if(std::fabs(y)<0.05 && std::fabs(y)>1.0e-5)
      {
       double new_x = Problem_Parameter::x_as_function_of_y(y);
       nod_pt->x(0) = new_x;
      }
     else if(std::fabs(y)<=1.0e-5)
      {
       continue;
      }
     else if(y<0)
      {
       double new_y = -Problem_Parameter::y_as_function_of_x(x);
       nod_pt->x(1) = new_y;
      }
     else if(y>0)
      {
       double new_y = Problem_Parameter::y_as_function_of_x(x);
       nod_pt->x(1) = new_y;
      }
     else
      {
       oomph_info << "obacht Trouble!"<<std::endl;
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
     
     // Set the wall speed function
     surface_element_pt->wall_speed_fct_pt() = 
      &Problem_Parameter::moving_frame_velocity_fct;
     
     if (CommandLineArgs::command_line_flag_has_been_set("--inc_thin_films"))  
      {
       surface_element_pt->enable_thin_film_effects();
       surface_element_pt->thin_film_homotopy_pt() =
        &Problem_Parameter::Thin_film_homotopy_parameter;
      }

     // Add bubble pressure (computed from position constraint)
     // as external data for these elements
     //surface_element_pt->add_external_data(Problem_Parameter::P_bubble_pt);

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
   
   Problem_Parameter::Prescribed_bubble_tip_x_position = 
    Problem_Parameter::Tip_node_pt->x(0);
   
  }

 void set_displacement_control_node_pt()
  {
   Problem_Parameter::Displacement_control_node_pt =
    Bulk_mesh_pt->boundary_node_pt(Actual_inflow_boundary,0);

   unsigned num_nod=Bulk_mesh_pt->nboundary_node(Actual_inflow_boundary);
   for(unsigned inod=0; inod<num_nod; inod++)
    {
     // Get node
     Node* nod_pt=Bulk_mesh_pt->boundary_node_pt(Actual_inflow_boundary,inod);
     if(std::fabs(nod_pt->x(1)) < 1.0e-8)
      {
       Problem_Parameter::Displacement_control_node_pt = nod_pt;
       break;
      }
    }
  }

 void set_pressure_control_node_pt()
  {
   Problem_Parameter::Pressure_control_node_pt =
    Bulk_mesh_pt->boundary_node_pt(Outflow_boundary,0);

   unsigned num_nod=Bulk_mesh_pt->nboundary_node(Outflow_boundary);
   for(unsigned inod=0; inod<num_nod; inod++)
    {
     // Get node
     Node* nod_pt=Bulk_mesh_pt->boundary_node_pt(Outflow_boundary,inod);
     if(std::fabs(nod_pt->x(1)) < 1.0e-8)
      {
       Problem_Parameter::Pressure_control_node_pt = nod_pt;
       break;
      }
    }
  }

 void create_constraint_element()
  { 
   oomph_info << "Construct constraints"<<std::endl;

   // Create new position constraint element -- recycle any previous
   // assignment of bubble pressure
   oomph_info << "Tip node pt = "<<Problem_Parameter::Tip_node_pt<<std::endl;

   if(Problem_Parameter::Flux_constraint)
    {
     Flux_constraint_element_pt =
      new FluxConstraintElement(
       Problem_Parameter::Pressure_gradient_downstream_backup,
       Problem_Parameter::P_bubble_backup,
       Problem_Parameter::Downstream_area_backup,                                     //Joao: added to implement pressure control
       Problem_Parameter::U_bubble_tip_backup,
       Problem_Parameter::flow_rate_var_backup,                                       ///Joao: added to keep the buble speed at the desired constant value
       Outflow_mesh_pt,
       Actual_inflow_mesh_pt,
       dynamic_cast<SolidNode*>(Problem_Parameter::Tip_node_pt),
       &Problem_Parameter::Prescribed_bubble_tip_x_position,
       &Problem_Parameter::Prescribed_area,
       &Problem_Parameter::Prescribed_bubble_pressure,                                //Joao: added to implement pressure control
       &Problem_Parameter::Prescribe_area,                                            //Joao: added to implement pressure control
       &Problem_Parameter::Impose_U_bubble_tip,                                       ///Joao: flag to switch imposed parameter: bubble speed or flow rate
       &Problem_Parameter::U_bubble_tip_dimensional_desired,                          ///Joao: added to keep the buble speed at the desired constant value
       &Problem_Parameter::Flow_rate_dimensional_desired,                             ///Joao: added to keep the flow rate at the desired constant value
       &Problem_Parameter::undeformed_dimensional_rectangular_cross_section);         ///Joao: added to keep the buble speed at the desired constant value

     
     // Set the pressure pointer to allow external access to the pressure
     Problem_Parameter::Pressure_gradient_downstream_pt = 
      Flux_constraint_element_pt->internal_data_pt(0);
     Problem_Parameter::P_bubble_pt = 
      Flux_constraint_element_pt->internal_data_pt(1);
     Problem_Parameter::U_bubble_tip_pt = 
      Flux_constraint_element_pt->internal_data_pt(2);
     Problem_Parameter::flow_rate_var_pt =                                ///Joao: added to keep the buble speed at the desired constant value
      Flux_constraint_element_pt->internal_data_pt(3);                    ///Joao: added to keep the buble speed at the desired constant value
     Problem_Parameter::Downstream_area_pt =                              //Joao: added to implement pressure control                                 
      Flux_constraint_element_pt->internal_data_pt(4);                    //Joao: added to implement pressure control

     // Add it to the mesh
     Constraint_mesh_pt->add_element_pt(Flux_constraint_element_pt);
    }
   else
    {
     Displacement_constraint_element_pt =
      new DisplacementConstraintElement(
       Problem_Parameter::Pressure_gradient_downstream_backup,
       Problem_Parameter::P_bubble_backup,
       Problem_Parameter::U_bubble_tip_backup,
       Problem_Parameter::Downstream_area_backup,
       Outflow_mesh_pt,
       Actual_inflow_mesh_pt,
       dynamic_cast<SolidNode*>(Problem_Parameter::Tip_node_pt),
       &Problem_Parameter::Prescribed_bubble_tip_x_position,
       Problem_Parameter::Displacement_control_node_pt,
       &Problem_Parameter::Prescribed_displacement);
     
     // Set the pressure pointer to allow external access to the pressure
     Problem_Parameter::Pressure_gradient_downstream_pt = 
      Displacement_constraint_element_pt->internal_data_pt(0);
     Problem_Parameter::P_bubble_pt = 
      Displacement_constraint_element_pt->internal_data_pt(1);
     Problem_Parameter::U_bubble_tip_pt = 
      Displacement_constraint_element_pt->internal_data_pt(2);
     Problem_Parameter::Downstream_area_pt =
      Displacement_constraint_element_pt->internal_data_pt(3);
     // Add it to the mesh
     Constraint_mesh_pt->add_element_pt(Displacement_constraint_element_pt);
     Problem_Parameter::flow_rate_var_pt = 0; //Joao: So we reset the pointer to the next adaptation
    }
   oomph_info << "Constraint elements created"<<std::endl;
  }

 void delete_constraint_element()
  {
   oomph_info << "Remove constraints"<<std::endl;

   Problem_Parameter::Pressure_gradient_downstream_backup = 
    Problem_Parameter::Pressure_gradient_downstream_pt->value(0);

   oomph_info << "Flag 1"<<std::endl;
   
   // Keep a temporary copy of the pressure from the old position contraint
   // element
   Problem_Parameter::P_bubble_backup =
    Problem_Parameter::P_bubble_pt->value(0);

   oomph_info << "Flag 2"<<std::endl;

   Problem_Parameter::U_bubble_tip_backup =
    Problem_Parameter::U_bubble_tip_pt->value(0);

   oomph_info << "Flag 3"<<std::endl;


   if((Problem_Parameter::flow_rate_var_pt == 0)||(!Problem_Parameter::Flux_constraint))
   {

   oomph_info << "Flag 4"<<std::endl;

   Problem_Parameter::flow_rate_var_backup =                                          ///Joao: added to keep the buble speed at the desired constant value
    Problem_Parameter::Prescribed_flow_rate;                                          ///Joao: added to keep the buble speed at the desired constant value    
    }
   else
   {

   oomph_info << "Flag 5"<<std::endl;

   Problem_Parameter::flow_rate_var_backup =                                          ///Joao: added to keep the buble speed at the desired constant value
    Problem_Parameter::flow_rate_var_pt->value(0);                                    
    }

   oomph_info << "Flag 6"<<std::endl;



     // Compute the actual flux ahead of the bubble
     double area_ahead_of_bubble = 0.0;
     {
      unsigned n_element = Outflow_mesh_pt->nelement();
      for(unsigned e = 0; e < n_element; e++)
       {
        // Upcast from GeneralisedElement to the present element
        HeleShawFluxElement<ELEMENT>* el_pt = 
         dynamic_cast<HeleShawFluxElement<ELEMENT>* >(
          Outflow_mesh_pt->element_pt(e));
        
        area_ahead_of_bubble += el_pt->compute_cross_section();
       }
     }

   oomph_info << "Flag 7"<<std::endl;

   Problem_Parameter::Downstream_area_backup = area_ahead_of_bubble;

   oomph_info << "Flag 8"<<std::endl;

//   Problem_Parameter::Downstream_area_pt->set_value(0,area_ahead_of_bubble);



//   if(!Problem_Parameter::Prescribe_area)
//    {
//     Problem_Parameter::Downstream_area_backup =
//      Problem_Parameter::Downstream_area_pt->value(0);
//    }    

   oomph_info << "Flag 9"<<std::endl;

   unsigned n_el = Constraint_mesh_pt->nelement();
   for(unsigned el=0;el<n_el;el++)
    {
     delete Constraint_mesh_pt->element_pt(el);
    }

   oomph_info << "Flag 10"<<std::endl;

   Constraint_mesh_pt->flush_element_and_node_storage();
   
   oomph_info << "Done removing constraints"<<std::endl;
  }

 void create_outflow_elements()
  {
   // Loop over the appropriate boundaries
   unsigned b = Outflow_boundary;
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
     HeleShawFluxElement<ELEMENT>* flux_element_pt = new
      HeleShawFluxElement<ELEMENT>(bulk_elem_pt, face_index);

     flux_element_pt->flux_fct_pt() = 
      &Problem_Parameter::normal_flux_ahead_of_bubble;
     
     // Set the wall speed function
     flux_element_pt->wall_speed_fct_pt() = 
      &Problem_Parameter::moving_frame_velocity_fct;
     
     // Set the non-dim reference gap width
     flux_element_pt->aspect_ratio_pt() = 
      &Problem_Parameter::Aspect_ratio;
     
     // Add pressure gradient (computed from flux constraint)
     // as external data for these elements
     //flux_element_pt->add_external_data(
     // Problem_Parameter::Pressure_gradient_downstream_pt);

     Outflow_mesh_pt->add_element_pt(flux_element_pt);
    }
  }

 void delete_outflow_elements()
  {
   // How many flux elements are in the flux mesh
   unsigned n_element = Outflow_mesh_pt->nelement();
   
   // Loop over the flux elements
   for(unsigned e = 0; e < n_element; e++)
    {
     // Kill flux element
     delete Outflow_mesh_pt->element_pt(e);
    }
   
   // Wipe the mesh
   Outflow_mesh_pt->flush_element_and_node_storage();
  }

 void create_top_inflow_elements()
  {
   // Loop over the appropriate boundaries
   unsigned b = Top_inflow_boundary;
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
     HeleShawFluxElement<ELEMENT>* flux_element_pt = new
      HeleShawFluxElement<ELEMENT>(bulk_elem_pt, face_index);
     
     // Set the non-dim reference gap width
     flux_element_pt->aspect_ratio_pt() = 
      &Problem_Parameter::Aspect_ratio;
     
     // Set the wall speed function
     flux_element_pt->wall_speed_fct_pt() = 
      &Problem_Parameter::moving_frame_velocity_fct;

     // flux_element_pt->flux_fct_pt() = 
     //  &Problem_Parameter::normal_flux_behind_bubble;
     
     // if (CommandLineArgs::command_line_flag_has_been_set("--inc_thin_films"))
     //  {
     //   flux_element_pt->enable_thin_film_effects();
     //   flux_element_pt->ca_inv_pt() = &Problem_Parameter::Ca_inv;
     //   flux_element_pt->set_bubble_tail_nodes(
     //    Problem_Parameter::Top_node_pt,Problem_Parameter::Bottom_node_pt);
     //  }

     Top_inflow_mesh_pt->add_element_pt(flux_element_pt);
    }
  }

 void delete_top_inflow_elements()
  {
   // How many flux elements are in the flux mesh
   unsigned n_element = Top_inflow_mesh_pt->nelement();
   
   // Loop over the flux elements
   for(unsigned e = 0; e < n_element; e++)
    {
     // Kill flux element
     delete Top_inflow_mesh_pt->element_pt(e);
    }
   
   // Wipe the mesh
   Top_inflow_mesh_pt->flush_element_and_node_storage();
  }

 void create_bottom_inflow_elements()
  {
   // Loop over the appropriate boundaries
   unsigned b = Bottom_inflow_boundary;
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
     HeleShawFluxElement<ELEMENT>* flux_element_pt = new
      HeleShawFluxElement<ELEMENT>(bulk_elem_pt, face_index);
     
     // Set the non-dim reference gap width
     flux_element_pt->aspect_ratio_pt() = 
      &Problem_Parameter::Aspect_ratio;
     
     // Set the wall speed function
     flux_element_pt->wall_speed_fct_pt() = 
      &Problem_Parameter::moving_frame_velocity_fct;
     
     // if (CommandLineArgs::command_line_flag_has_been_set("--inc_thin_films"))
     //  {
     //   flux_element_pt->enable_thin_film_effects();
     //   flux_element_pt->ca_inv_pt() = &Problem_Parameter::Ca_inv;
     //   flux_element_pt->set_bubble_tail_nodes(
     //    Problem_Parameter::Top_node_pt,Problem_Parameter::Bottom_node_pt);
     //  }

     // flux_element_pt->flux_fct_pt() = 
     //  &Problem_Parameter::normal_flux_behind_bubble;

     Bottom_inflow_mesh_pt->add_element_pt(flux_element_pt);
    }
  }

 void delete_bottom_inflow_elements()
  {
   // How many flux elements are in the flux mesh
   unsigned n_element = Bottom_inflow_mesh_pt->nelement();
   
   // Loop over the flux elements
   for(unsigned e = 0; e < n_element; e++)
    {
     // Kill flux element
     delete Bottom_inflow_mesh_pt->element_pt(e);
    }
   
   // Wipe the mesh
   Bottom_inflow_mesh_pt->flush_element_and_node_storage();
  }

 void create_actual_inflow_elements()
  {
   // Loop over the appropriate boundaries
   unsigned b = Actual_inflow_boundary;
   // How many bulk fluid elements are adjacent to boundary b in region 1?
   unsigned n_element = Bulk_mesh_pt->nboundary_element_in_region(b, 1);
   for(unsigned e = 0; e < n_element; e++)
    {
     // Get pointer to the bulk fluid element that is 
     // adjacent to boundary b in region 0
     ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>(
      Bulk_mesh_pt->boundary_element_in_region_pt(b, 1, e));
     
     //Find the index of the face of element e along boundary b in region 0
     int face_index = Bulk_mesh_pt->
      face_index_at_boundary_in_region(b, 1, e);
     
     // Build it
     HeleShawFluxElement<ELEMENT>* flux_element_pt = new
      HeleShawFluxElement<ELEMENT>(bulk_elem_pt, face_index);

     // flux_element_pt->flux_fct_pt() = 
     //  &Problem_Parameter::thin_film_flux_behind_bubble;
     
     // Set the non-dim reference gap width
     flux_element_pt->aspect_ratio_pt() = 
      &Problem_Parameter::Aspect_ratio;
     
     // Set the wall speed function
     flux_element_pt->wall_speed_fct_pt() = 
      &Problem_Parameter::moving_frame_velocity_fct;
     
     if (CommandLineArgs::command_line_flag_has_been_set("--inc_thin_films")) 
      {
       flux_element_pt->enable_thin_film_effects();
       flux_element_pt->ca_pt() = &Problem_Parameter::Ca;
       flux_element_pt->thin_film_homotopy_pt() =
        &Problem_Parameter::Thin_film_homotopy_parameter;
      }

     Actual_inflow_mesh_pt->add_element_pt(flux_element_pt);
    }
  }

 void delete_actual_inflow_elements()
  {
   // How many flux elements are in the flux mesh
   unsigned n_element = Actual_inflow_mesh_pt->nelement();
   
   // Loop over the flux elements
   for(unsigned e = 0; e < n_element; e++)
    {
     // Kill flux element
     delete Actual_inflow_mesh_pt->element_pt(e);
    }
   
   // Wipe the mesh
   Actual_inflow_mesh_pt->flush_element_and_node_storage();
  }

}; // end_of_problem_class


//========================================================================
/// Constructor of the problem class - build mesh and set up elements with
/// the appropriate parameters
//========================================================================
template<class ELEMENT>
UnstructuredHSFvKProblem<ELEMENT>::UnstructuredHSFvKProblem()
 : Counter(0), Mesh_as_geom_object_pt(0), Constraint_mesh_pt(0), 
   Displacement_constraint_element_pt(0), Flux_constraint_element_pt(0), 
   Backed_up_surface_mesh_pt(0),
   Update_dependent_parameters(true)               //// Joao ===>>> Before false
{

 // obacht
 //enable_globally_convergent_newton_method();
 Max_newton_iterations = 40;
 Always_take_one_newton_step = true;
 //Use_continuation_timestepper=true;
 //Newton_solver_tolerance = 1.0e-10;

// Create the eigen solver 
 this->eigen_solver_pt() = new ANASAZI;

 //Use_finite_differences_for_continuation_derivatives = true;
 //Bifurcation_detection = true;

 if (CommandLineArgs::command_line_flag_has_been_set(
      "--inc_transmural_pressure"))
  {
   Problem_Parameter::Include_transmural_pressure = true;
   oomph_info << "Including transmural pressure as external load on the sheet"
              << std::endl;
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

 //linear_solver_pt() = new FD_LU;

 // Allow for rough startup
 this->Problem::Max_residuals = 2000.0; // 1.0e20;
 
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
 bound_vertex[0][0] = -Problem_Parameter::L_upstream;
 bound_vertex[0][1] = -0.5;

 bound_vertex[1][0] = Problem_Parameter::L_downstream;
 bound_vertex[1][1] = -0.5;

 outer_polyline_boundary_pt[0] = new TriangleMeshPolyLine(bound_vertex,
   Bottom_boundary);

 bound_vertex.resize(3);
 bound_vertex[0].resize(2);
 bound_vertex[1].resize(2);
 bound_vertex[2].resize(2);
 

 // Outflow boundary
 bound_vertex[0][0] = Problem_Parameter::L_downstream;
 bound_vertex[0][1] = -0.5;

 bound_vertex[1][0] = Problem_Parameter::L_downstream;
 bound_vertex[1][1] = 0.0;

 bound_vertex[2][0] = Problem_Parameter::L_downstream;
 bound_vertex[2][1] = 0.5;

 outer_polyline_boundary_pt[1] = new TriangleMeshPolyLine(bound_vertex,
   Outflow_boundary);

 bound_vertex.resize(2);
 bound_vertex[0].resize(2);
 bound_vertex[1].resize(2);

 // Top boundary
 bound_vertex[0][0] = Problem_Parameter::L_downstream;
 bound_vertex[0][1] = 0.5;

 bound_vertex[1][0] = -Problem_Parameter::L_upstream;
 bound_vertex[1][1] = 0.5;

 outer_polyline_boundary_pt[2] = new TriangleMeshPolyLine(bound_vertex,
   Top_boundary);

 // Top inflow boundary
 bound_vertex[0][0] = -Problem_Parameter::L_upstream;
 bound_vertex[0][1] = 0.5;

 bound_vertex[1][0] = -Problem_Parameter::L_upstream;
 bound_vertex[1][1] = Problem_Parameter::W_bubble_initial/2.0; 

 outer_polyline_boundary_pt[3] = new TriangleMeshPolyLine(bound_vertex,
   Top_inflow_boundary);

 bound_vertex.resize(3);
 bound_vertex[0].resize(2);
 bound_vertex[1].resize(2);
 bound_vertex[2].resize(2);

 // Actual inflow boundary
 bound_vertex[0][0] = -Problem_Parameter::L_upstream;
 bound_vertex[0][1] = Problem_Parameter::W_bubble_initial/2.0;

 bound_vertex[1][0] = -Problem_Parameter::L_upstream;
 bound_vertex[1][1] = 0.0;

 bound_vertex[2][0] = -Problem_Parameter::L_upstream;
 bound_vertex[2][1] = -Problem_Parameter::W_bubble_initial/2.0;

 outer_polyline_boundary_pt[4] = new TriangleMeshPolyLine(bound_vertex,
   Actual_inflow_boundary);

 bound_vertex.resize(2);
 bound_vertex[0].resize(2);
 bound_vertex[1].resize(2);

 // Bottom inflow boundary
 bound_vertex[0][0] = -Problem_Parameter::L_upstream;
 bound_vertex[0][1] = -Problem_Parameter::W_bubble_initial/2.0;

 bound_vertex[1][0] = -Problem_Parameter::L_upstream;
 bound_vertex[1][1] = -0.5;

 
 outer_polyline_boundary_pt[5] = new TriangleMeshPolyLine(bound_vertex,
    Bottom_inflow_boundary);

 outer_boundary_pt =
  new TriangleMeshClosedCurve(outer_polyline_boundary_pt);
 


 //Inner boundary
 //--------------
 
 double A = Problem_Parameter::W_bubble_initial/2.0;
 double B = Problem_Parameter::W_bubble_initial/2.0;
 Ellipse* inner_boundary_ellipse_pt = new Ellipse(A, B);

 Vector<TriangleMeshOpenCurve*> inner_boundary_pt(2);

 Vector<TriangleMeshCurveSection*> inner_polyline_boundary_pt(1);

 Vector<TriangleMeshCurveSection*> dummy_polyline_boundary_pt(1);

 // Number of segments on bubble
 unsigned nsegment=30;
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
 Vector<Vector<double> > bound_hole(nsegment + 3);

 bound_hole[0].resize(2);
 bound_hole[0][0] = -Problem_Parameter::L_upstream;
 bound_hole[0][1] = Problem_Parameter::W_bubble_initial/2.0;

 //First part
 for(unsigned ipoint = 0; ipoint < nsegment+2; ipoint++)
  {
   bound_hole[ipoint+1].resize(2);

   zeta[0] = 0.5*MathematicalConstants::Pi+unit_zeta*double(ipoint);
   inner_boundary_ellipse_pt->position(zeta,posn);
   bound_hole[ipoint+1][0] = posn[0]-Problem_Parameter::W_bubble_initial/2.0;
   bound_hole[ipoint+1][1] = posn[1];
  }

 bound_hole[nsegment+2].resize(2);
 bound_hole[nsegment+2][0] = -Problem_Parameter::L_upstream;
 bound_hole[nsegment+2][1] = -Problem_Parameter::W_bubble_initial/2.0;

 if (CommandLineArgs::command_line_flag_has_been_set("--saffman_taylor"))
  { 
   // Snap to Saffman-Taylor finger shape
   for(unsigned i=1;i<bound_hole.size()-1;i++)
    {
     double x = bound_hole[i][0];
     double y = bound_hole[i][1];
     
     //std::cout<<x<<" "<<y<<std::endl;
     
     if(std::fabs(y)<0.05 && std::fabs(y)>1.0e-5)
      {
       double new_x = Problem_Parameter::x_as_function_of_y(y);
       bound_hole[i][0] = new_x;
      }
     else if(std::fabs(y)<=1.0e-5)
      {
       continue;
      }
     else if(y<0)
      {
       double new_y = -Problem_Parameter::y_as_function_of_x(x);
       bound_hole[i][1] = new_y;
      }
     else if(y>0)
      {
       double new_y = Problem_Parameter::y_as_function_of_x(x);
       bound_hole[i][1] = new_y;
      }
     else
      {
       oomph_info << "obacht Trouble!"<<std::endl;
      }
    }
  }

 ofstream file;
 file.open("bubble_boundary.dat");
 inner_polyline_boundary_pt[0] = new TriangleMeshPolyLine(bound_hole,
   Inner_boundary);
 inner_polyline_boundary_pt[0]->output(file);
 file.close();

 // Limit the size of the boundary segments 
 double max_length = Problem_Parameter::Max_polyline_length_bubble; 
 if (max_length>0.0)
  {
   inner_polyline_boundary_pt[0]->set_maximum_length(max_length);
  }

 // for(unsigned i=0; i<bound_hole.size();i++)
 //  {
 //   std::cout<<bound_hole[i][0] << " " <<bound_hole[i][1]<<std::endl;
 //  }

 bound_vertex[0][0] = -Problem_Parameter::L_upstream;
 bound_vertex[0][1] = 0.0;

 bound_vertex[1][0] = -Problem_Parameter::L_upstream+1.0;
 bound_vertex[1][1] = 0.0;

 dummy_polyline_boundary_pt[0] = new TriangleMeshPolyLine(bound_vertex,
   Dummy_boundary);

 inner_boundary_pt[0] = new TriangleMeshOpenCurve(inner_polyline_boundary_pt);

 inner_boundary_pt[1] = new TriangleMeshOpenCurve(dummy_polyline_boundary_pt);

 // Connect initial vertex to final vertex on upper inlet boundary
 inner_polyline_boundary_pt[0]->connect_initial_vertex_to_polyline(
  dynamic_cast<TriangleMeshPolyLine*>
  (outer_polyline_boundary_pt[3]),1);
 
 // Connect final vertex to first vertex on lower outlet boundary
 inner_polyline_boundary_pt[0]->connect_final_vertex_to_polyline(
  dynamic_cast<TriangleMeshPolyLine*>
  (outer_polyline_boundary_pt[5]),0);

 // Connect initial vertex to middle vertex on actual inflow boundary
 dummy_polyline_boundary_pt[0]->connect_initial_vertex_to_polyline(
  dynamic_cast<TriangleMeshPolyLine*>
  (outer_polyline_boundary_pt[4]),1);

 Vector<double> bubble_region_coords1(2);
 bubble_region_coords1[0] = -Problem_Parameter::W_bubble_initial;
 bubble_region_coords1[1] = 0.0;


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

 // Set target area for elements in the initial mesh 
 mesh_parameters.element_area() = Problem_Parameter::Max_el_area; 

 // Build mesh
 Bulk_mesh_pt = new RefineableSolidTriangleMesh<ELEMENT>(mesh_parameters,
                                                  this->time_stepper_pt(0));

 Bulk_mesh_pt->set_lagrangian_nodal_coordinates();

 //Bulk_mesh_pt->output("mesh_before_snapping.dat");
 Bulk_mesh_pt->output_boundaries("boundaries_before_snapping.dat");

 if (CommandLineArgs::command_line_flag_has_been_set("--saffman_taylor"))
  {
   // Snap nodes onto curvilinear boundary
   snap_onto_boundary();
  }
 set_tip_node_pt();
 set_displacement_control_node_pt();
 set_pressure_control_node_pt();


 Z2ErrorEstimator* error_estimator_pt = new Z2ErrorEstimator;
 Bulk_mesh_pt->spatial_error_estimator_pt() = error_estimator_pt;
 Bulk_mesh_pt->max_permitted_error() = Problem_Parameter::Max_permitted_error; 
 Bulk_mesh_pt->min_permitted_error() = Problem_Parameter::Min_permitted_error;
 Bulk_mesh_pt->max_element_size() = Problem_Parameter::Max_el_area; 
 Bulk_mesh_pt->min_element_size() = Problem_Parameter::Min_el_area;


 if (!CommandLineArgs::command_line_flag_has_been_set("--pin_hs"))
  { 
   /// New mesh for the interface elements, and create them
   Free_surface_bc_mesh_pt = new Mesh;
   create_surface_elements();

   Actual_inflow_mesh_pt = new Mesh;
   create_actual_inflow_elements();

   /// New mesh for the flux elements, and create them
   Top_inflow_mesh_pt = new Mesh;
   create_top_inflow_elements();

   Bottom_inflow_mesh_pt = new Mesh;
   create_bottom_inflow_elements();
   
   Outflow_mesh_pt = new Mesh;
   create_outflow_elements();
  }

 oomph_info << "Iip_node_pt = "<<Problem_Parameter::Tip_node_pt<<std::endl;
 oomph_info << "P_bubble_pt = "<<Problem_Parameter::P_bubble_pt<<std::endl;
 oomph_info << "Pressure_gradient_downstream_pt = "
            <<Problem_Parameter::Pressure_gradient_downstream_pt<<std::endl;

 Constraint_mesh_pt = new Mesh;
 create_constraint_element();

 oomph_info << "P_bubble_pt = "<<Problem_Parameter::P_bubble_pt<<std::endl;
 oomph_info << "Pressure_gradient_downstream_pt = "
            <<Problem_Parameter::Pressure_gradient_downstream_pt<<std::endl;

 /// Pin the appropriate unknowns and set boundary conditions
 complete_problem_setup();

 // oomph_info << "Pointer "<<Free_surface_bc_mesh_pt->element_pt(0)->external_data_pt(0)<<std::endl;

 Bulk_mesh_pt->output("mesh_after_snapping.dat");
 Bulk_mesh_pt->output_boundaries("boundaries_after_snapping.dat");
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
   add_sub_mesh(Outflow_mesh_pt);
   // if (CommandLineArgs::command_line_flag_has_been_set("--inc_thin_films"))
   //  {
   //   //add_sub_mesh(Top_inflow_mesh_pt);
   //   //add_sub_mesh(Bottom_inflow_mesh_pt);
   //   add_sub_mesh(Actual_inflow_mesh_pt);
   //  }
  }
 add_sub_mesh(Constraint_mesh_pt);
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
 
 // Apply boundary conditions
 apply_boundary_conditions();
 
 // Complete the build of all elements so they are fully functional
 unsigned n_element = Bulk_mesh_pt->nelement();
 for(unsigned e = 0; e < n_element; e++)
  {
   // Upcast from GeneralisedElement to the present element
   ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));
   
   //Set the pressure function pointers and the physical constants
   el_pt->eta_pt() = &Problem_Parameter::Eta;

   // External pressure
   el_pt->pressure_fct_pt()=&Problem_Parameter::external_pressure;

   // Put the sheet under tension
   el_pt->pre_stress_fct_pt() = &Problem_Parameter::get_pre_stress;

   //Set the FSI parameter
   el_pt->q_pt() = &Problem_Parameter::Q;

   // Set the non-dim reference gap width
   el_pt->aspect_ratio_pt() =  &Problem_Parameter::Aspect_ratio;
     
   // Set the wall speed function
   el_pt->wall_speed_fct_pt() = 
    &Problem_Parameter::moving_frame_velocity_fct;

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
       // el_pt->set_volume_constraint_pressure_data_as_external_data(
       //   Volume_constraint_element_pt->pressure_data_pt());
       el_pt->add_external_data(Problem_Parameter::P_bubble_pt,true);
      }

     // External pressure
     el_pt->pressure_fct_pt()=&Problem_Parameter::external_pressure;
      
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
         // obacht
         nod_pt->set_value(4,Problem_Parameter::P_bubble_pt->value(0));
         //*nod_pt->value_pt(4) = *Problem_Parameter::P_bubble_pt->value_pt(0);
        }
      }  

    }
  }

 n_element = Bulk_mesh_pt->nelement();
 for(unsigned e = 0; e < n_element; e++)
  {
   // Upcast from GeneralisedElement to the present element
   ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));

   el_pt->add_external_data(Problem_Parameter::U_bubble_tip_pt,true);
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

   el_pt->add_external_data(Problem_Parameter::U_bubble_tip_pt,true);

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

 // Add the pressure gradients as external data to the flux elements
 n_element = Outflow_mesh_pt->nelement();
 for(unsigned e = 0; e < n_element; e++)
  {
   // Upcast from GeneralisedElement to the present element
   HeleShawFluxElement<ELEMENT>* el_pt = 
    dynamic_cast<HeleShawFluxElement<ELEMENT>* >(
     Outflow_mesh_pt->element_pt(e));
   
   el_pt->add_external_data(
    Problem_Parameter::Pressure_gradient_downstream_pt,true);

   el_pt->add_external_data(Problem_Parameter::U_bubble_tip_pt,true);

   // Because we not only need the pressure gradient G to compute the flux,
   // but also the deflection of the sheet dbdx, we need to add
   // external data
   ELEMENT* bulk_el_pt = dynamic_cast<ELEMENT*>(el_pt->bulk_element_pt());
   unsigned bulk_nodes = bulk_el_pt->nnode();
   for(unsigned inod=0; inod<bulk_nodes; inod++)
    {
     Node* nod_pt = dynamic_cast<Node*>(bulk_el_pt->node_pt(inod));
     if( !(nod_pt->is_on_boundary(Outflow_boundary)) )
      {
       el_pt->add_external_data(nod_pt, true);
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
   
   oomph_info << "Outputting for step: " << Doc_info.number();
   if(!Is_steady)   
    {
     oomph_info << " ; time= " << this->time_stepper_pt()->time();
    }        
   oomph_info << std::endl;
   
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
   // nel = Bulk_mesh_pt->nregion_element(r);
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
      x[0]=-Problem_Parameter::L_upstream+
       double(j)*(Problem_Parameter::L_downstream + 
                  Problem_Parameter::L_upstream)/double(nplot-1);
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
        
          // // Convert to polars
          // double sigma_rr=0.0;
          // double sigma_tt=0.0;
          // double sigma_rt=0.0;
          // if ((x[0]==0.0)&&(x[1]==0.0))
          //  {
          //  }
          // else
          //  {
          //   double theta=atan2(x[1],x[0]);
          //   sigma_rr=
          //    sigma_xx*cos(theta)*cos(theta)+
          //    sigma_yy*sin(theta)*sin(theta)+
          //    sigma_xy*sin(2.0*theta);
          
          //   sigma_tt=
          //    sigma_xx*sin(theta)*sin(theta)+
          //    sigma_yy*cos(theta)*cos(theta)-
          //    sigma_xy*sin(2.0*theta);
          
          //   sigma_rt=
          //    (sigma_yy-sigma_xx)*sin(theta)*cos(theta)+
          //    sigma_xy*cos(2.0*theta);
          
          //  }


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
                    // << sigma_rr << " " 
                    // << sigma_tt << " " 
                    // << sigma_rt << " "
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

/////////////////////////////////////////////////////////////////////////////////                     /////////Joao
/////////////////////////////////////////////////////////////////////////////////
//   double asymmetry_level = 0.0;                                            
//   double mean_y_nod;
//
//   
//     // How many surface elements are in the surface mesh
//     unsigned n_element = Free_surface_bc_mesh_pt->nelement();
//   
//   Vector<double> mean_y_vector(n_element,0.0);  
//
//     // Loop over the surface elements
//     for(unsigned e = 0; e < n_element; e++)
//      {
//
//       HeleShawInterfaceElement<ELEMENT>* surface_element_pt = 
//        dynamic_cast<HeleShawInterfaceElement<ELEMENT>*>(
//         Free_surface_bc_mesh_pt->element_pt(e));
//     
//       surface_element_pt->compute_mean_y(arc_length_contribution,mean_y_nod);
//         
//         mean_y_vector[e] = mean_y_nod;  
//
//      }
//     for(unsigned e = 0; e < n_element; e++)
//      {
//       if(e==0)
//        {
//       asymmetry_level = 0.0;
//        }
//      asymmetry_level += (mean_y_vector[e] + mean_y_vector[n_element-1-e])*(mean_y_vector[e] + mean_y_vector[n_element-1-e])/double(n_element);
//      }

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////



   double p_bubble=0.0;
   if (Problem_Parameter::P_bubble_pt!=0)
    {
     p_bubble = Problem_Parameter::P_bubble_pt->value(0);
    }
   // Output time, bubble pressure and bubble volume etc to the trace file
   // obacht overwrite
   w_0 = Problem_Parameter::Displacement_control_node_pt->value(0);

   // Compute the actual flux ahead of the bubble
   double flux_ahead_of_bubble = 0.0;
   {
    unsigned n_element = Outflow_mesh_pt->nelement();
    for(unsigned e = 0; e < n_element; e++)
     {
      // Upcast from GeneralisedElement to the present element
      HeleShawFluxElement<ELEMENT>* el_pt = 
       dynamic_cast<HeleShawFluxElement<ELEMENT>* >(
        Outflow_mesh_pt->element_pt(e));
      
      flux_ahead_of_bubble += el_pt->compute_volume_flux();
     }
   }

   // Compute the actual flux behind the bubble top
   double flux_behind_bubble_top = 0.0;
   {
    unsigned n_element = Top_inflow_mesh_pt->nelement();
    for(unsigned e = 0; e < n_element; e++)
     {
      // Upcast from GeneralisedElement to the present element
      HeleShawFluxElement<ELEMENT>* el_pt = 
       dynamic_cast<HeleShawFluxElement<ELEMENT>* >(
        Top_inflow_mesh_pt->element_pt(e));
      
      flux_behind_bubble_top += el_pt->compute_volume_flux();
     }
   }

   // Compute the actual flux behind the bubble bottom
   double flux_behind_bubble_bottom = 0.0;
   {
    unsigned n_element = Bottom_inflow_mesh_pt->nelement();
    for(unsigned e = 0; e < n_element; e++)
     {
      // Upcast from GeneralisedElement to the present element
      HeleShawFluxElement<ELEMENT>* el_pt = 
       dynamic_cast<HeleShawFluxElement<ELEMENT>* >(
        Bottom_inflow_mesh_pt->element_pt(e));
      
      flux_behind_bubble_bottom += el_pt->compute_volume_flux();
     }
   }

   // Compute mass flow due to thin films
   double thin_film_mass_flow = 0.0;
   // {
   //  unsigned n_element = Actual_inflow_mesh_pt->nelement();
   //  for(unsigned e = 0; e < n_element; e++)
   //   {
   //    // Upcast from GeneralisedElement to the present element
   //    HeleShawFluxElement<ELEMENT>* el_pt = 
   //     dynamic_cast<HeleShawFluxElement<ELEMENT>* >(
   //      Actual_inflow_mesh_pt->element_pt(e));
      
   //    thin_film_mass_flow += el_pt->compute_thin_film_mass_flow();
   //   }
   // }
   {
    unsigned n_element = Free_surface_bc_mesh_pt->nelement();
    for(unsigned e = 0; e < n_element; e++)
     {
      // Upcast from GeneralisedElement to the present element
      HeleShawInterfaceElement<ELEMENT>* el_pt = 
       dynamic_cast<HeleShawInterfaceElement<ELEMENT>* >(
        Free_surface_bc_mesh_pt->element_pt(e));
      
      thin_film_mass_flow += el_pt->compute_flow_into_thin_films();
     }
   }

   double cross_section_air=0.0;
   {
    unsigned n_element = Actual_inflow_mesh_pt->nelement();
    for(unsigned e = 0; e < n_element; e++)
     {
      // Upcast from GeneralisedElement to the present element
      HeleShawFluxElement<ELEMENT>* el_pt = 
       dynamic_cast<HeleShawFluxElement<ELEMENT>* >(
        Actual_inflow_mesh_pt->element_pt(e));
      
      cross_section_air += el_pt->compute_cross_section();
     }
   }

   double prescribed_downstream_area = Problem_Parameter::Prescribed_area;
   if((Problem_Parameter::Downstream_area_pt == 0)||((Problem_Parameter::Prescribe_area)&&(Problem_Parameter::Flux_constraint)))
    {

     oomph_info << "Flag 1b"<<std::endl;
     prescribed_downstream_area = Problem_Parameter::Prescribed_area;
    }
   else
    {
     prescribed_downstream_area = 
      Problem_Parameter::Downstream_area_pt->value(0);
     oomph_info << "Flag 2b"<<std::endl;

    } 


   double inflow_rate;
   if((Problem_Parameter::flow_rate_var_pt == 0)||(!Problem_Parameter::Flux_constraint))
    {
     oomph_info << "Flag 3b"<<std::endl;     

     inflow_rate=Problem_Parameter::flow_rate_var_backup;

     oomph_info << "Flag 4b"<<std::endl;
     }
   else
    {
     oomph_info << "Flag 5b"<<std::endl;

     inflow_rate=Problem_Parameter::flow_rate_var_pt->value(0);    

     oomph_info << "Flag 6b"<<std::endl;
     }

   Trace_file << time_stepper_pt()->time() << " "
              << mean_radius << " "
              << deviation << " "
              << Problem_Parameter::Top_node_pt->x(1)-
    Problem_Parameter::Bottom_node_pt->x(1) << " "
              << p_bubble << " "
              << w_0 << " "
              << Problem_Parameter::U_bubble_tip_pt->value(0) << " "
              << Problem_Parameter::Q << " "
              << Problem_Parameter::Ca << " "
              << Problem_Parameter::Aspect_ratio << " "
              << prescribed_downstream_area << " "
              << Problem_Parameter::Pressure_gradient_downstream_pt->value(0) 
              << " "
              << flux_ahead_of_bubble << " "
              << flux_behind_bubble_top << " "
              << flux_behind_bubble_bottom << " "
              << thin_film_mass_flow << " "
              << cross_section_air << " "
              << Bulk_mesh_pt->nnode() << " "
              << Bulk_mesh_pt->nelement() << " "
              << Problem_Parameter::Pre_stress << " "
              << ndof() << " "
              << Doc_info.number() << " " 
              << inflow_rate << '\n';

   Trace_file.flush();
   
  } // end of full output
 
 
 // obacht
 //Write restart file on all procs -- awkward but this call is
 //needed because dump re-orders the nodes in the mesh!
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
 
 // Initial bubble width
 CommandLineArgs::specify_command_line_flag(
  "--w_bubble",
  &Problem_Parameter::W_bubble_initial);

 // Initial perturbation of the bubble radius
 CommandLineArgs::specify_command_line_flag(
  "--r_perturb",
  &Problem_Parameter::Radius_perturb);

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

 // Start with Saffman Taylor finger as initial shape?
 CommandLineArgs::specify_command_line_flag("--saffman_taylor");

 // Are we applying thin film corrections?
 CommandLineArgs::specify_command_line_flag("--inc_thin_films");

 // Are we performing a steady solve?
 CommandLineArgs::specify_command_line_flag("--steady");

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

 // Number of segments on half of bubble surface
 CommandLineArgs::specify_command_line_flag(
  "--n_segment_bubble",
  &Problem_Parameter::N_segment_bubble);
 
 // Max length of polylines on bubble
 CommandLineArgs::specify_command_line_flag(
  "--max_polyline_length_bubble",
  &Problem_Parameter::Max_polyline_length_bubble);

 // Estimating the time step based on the displacement of the sheet
 CommandLineArgs::specify_command_line_flag(
  "--error_estimate_use_displacement");
 
 // Max length of polylines on bubble
 // double flux=-1.0;
 // CommandLineArgs::specify_command_line_flag(
 //  "--flux",&flux);

 // Use mumps
 CommandLineArgs::specify_command_line_flag("--use_mumps");

 // Include transmural pressure as external load
 CommandLineArgs::specify_command_line_flag("--inc_transmural_pressure");

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
 
 // Restart file specified via command line 
 if (CommandLineArgs::command_line_flag_has_been_set("--restart_file"))
  {
//   problem.restart();
//   oomph_info << "Done restart\n";
//   problem.doc_solution();
   
   bool unsteady_initial_config_steady = true;
   Problem_Parameter::Displacement_continuation=false;
   Problem_Parameter::Aspect_ratio_continuation=false;
   Problem_Parameter::FSI_continuation=false;
   Problem_Parameter::Flow_rate_continuation=false;

/////////////////////////////////////////////////// choose which continuation to make
   bool stability_solver=false;
   bool pre_stress_continuation=false;
   bool continuation_over_Aspect_ratio=false;
   bool continuation_over_FSI=false;
   bool continuation_over_Displacement=false;
   bool continuation_over_area=false;
   bool continuation_over_bubble_pressure=true;
   bool continuation_over_U_bubble_tip=false;
   bool continuation_over_Flow_rate=false;
   bool time_evolution=false;
/////////////////////////////////////////////////// choose which continuation to make   


    if (!CommandLineArgs::command_line_flag_has_been_set("--inc_thin_films"))
      {
       Problem_Parameter::Thin_film_homotopy_parameter = 0.0;
       }

  //Start time evolution
  if(time_evolution)
    {

     Problem_Parameter::Impose_U_bubble_tip=false;
     double flow_rate_desired_backup = 158.4822*1.0e-6/60.0;
     double U_bubble_backup = 0.1;
     if(Problem_Parameter::Impose_U_bubble_tip)
      {
       Problem_Parameter::U_bubble_tip_dimensional_desired=U_bubble_backup;
       }
     else
      {
       Problem_Parameter::Flow_rate_dimensional_desired = flow_rate_desired_backup;
       }

     unsteady_initial_config_steady = false;

     problem.enable_parameter_update();
     Problem_Parameter::Update_ca_only=false;  
     double area_backup = 0.92; //Problem_Parameter::Downstream_area_pt->value(0);
     epsilon_t = 1.0e-3;
     if(unsteady_initial_config_steady)
     {
     Problem_Parameter::Prescribe_area = true;
     problem.restart(unsteady_initial_config_steady);
     oomph_info << "Done restart\n";
     Problem_Parameter::Prescribed_area = area_backup;

     oomph_info << "Before make_unsteady\n";
     problem.make_unsteady();
     oomph_info << "After make_unsteady\n";

     problem.set_timestepper_for_all_data(problem.time_stepper_pt(0));
     problem.initialise_dt(1.0e-4);
     problem.assign_initial_values_impulsive();
//     problem.set_initial_condition();

     oomph_info << "Before adapt\n";
     problem.adapt();
     oomph_info << "After adapt\n";
     Problem_Parameter::Prescribed_area = area_backup;
//     problem.initialise_dt(1.0e-4); 
     problem.doc_solution();
     problem.adaptive_unsteady_newton_solve(1.0e-4,epsilon_t,false);
     problem.doc_solution();
//     problem.adapt();
//     problem.unsteady_newton_solve(1.0e-4,epsilon_t,false);//joao change the boolean flag
//     problem.doc_solution();
      }
     else
     {
      Problem_Parameter::Prescribe_area = true;
      problem.make_unsteady();
      problem.restart(unsteady_initial_config_steady);
      oomph_info << "Done restart\n";   
      Problem_Parameter::Prescribed_area = area_backup;
      problem.doc_solution();

     //problem.set_timestepper_for_all_data(problem.time_stepper_pt(0));
     //problem.initialise_dt(1.0e-4);
     //problem.assign_initial_values_impulsive();

      problem.adapt();
//      problem.unsteady_newton_solve(1.0e-4,0,false);
      problem.adaptive_unsteady_newton_solve(1.0e-4,epsilon_t,false);
      problem.doc_solution();
      }

     bool first = false; //true; in the first step
     bool perturb_flow_rate = false;
     double perturb_delay = 0.2;
     double flow_rate_perturb_parameter = 0.0;
     double flow_rate_perturb_attenuation = 2.0;
     unsigned max_adapt = 2;
     for(unsigned i = 0; i < nstep; i++)
      {   
       // Allow adaptation every few steps
       if( (i>0) && (Problem_Parameter::Adaptation_frequency>0) &&
           (i%Problem_Parameter::Adaptation_frequency == 0) )
        {
         problem.adapt();//max_adapt = 1;
         Problem_Parameter::Prescribed_area = area_backup;
//         problem.adaptive_unsteady_newton_solve(1.0e-4,epsilon_t,first);    
         problem.unsteady_newton_solve(1.0e-4,0,first);    
         first = false;
        }

       double next_dt = dt;

       if (CommandLineArgs::command_line_flag_has_been_set
           ("--suppress_temporal_adaptivity"))
        {
         problem.unsteady_newton_solve(dt,0,first);
         }
       else
        {
            if(perturb_flow_rate)
             {              
              if(double(problem.time_stepper_pt()->time())>perturb_delay)
               {
                flow_rate_perturb_parameter = 8.0*(double(problem.time_stepper_pt()->time())-perturb_delay);
                Problem_Parameter::Flow_rate_dimensional_desired = 
                       flow_rate_desired_backup*(1.0+(flow_rate_perturb_parameter/flow_rate_perturb_attenuation)*pow(2.718281828,-flow_rate_perturb_parameter));
                Problem_Parameter::update_dependent_parameters();
                }
              }
         next_dt = problem.adaptive_unsteady_newton_solve(dt,epsilon_t,first);
         //next_dt = problem.doubly_adaptive_unsteady_newton_solve(dt,epsilon_t,max_adapt,first);
         oomph_info << "Suggested next dt: " << next_dt << std::endl;
         }
       
       dt = next_dt;
       if(dt>0.1)
       {
        dt=0.1;
       }  
       
       problem.doc_solution();

      }
     exit(1);
     }  
  //End time evolution



   if(stability_solver)
   {// stability

   unsteady_initial_config_steady = true;
   Problem_Parameter::Prescribe_area=true;
   problem.enable_parameter_update();
   Problem_Parameter::Update_ca_only=false;
   Problem_Parameter::Impose_U_bubble_tip=false;
   
     Problem_Parameter::U_bubble_tip_dimensional_desired=0.1;
     Problem_Parameter::Flow_rate_dimensional_desired = 194.277*1.0e-6/60.0;
     
   problem.restart(unsteady_initial_config_steady);
   oomph_info << "Done restart\n";   
   double area_backup = 0.54; //Problem_Parameter::Downstream_area_pt->value(0);
   Problem_Parameter::Prescribed_area = area_backup;
//   problem.doc_solution();
   
   bool solve_for_flow_rate_control = true;
   bool time_evolution_for_perturbed_solution = true; //Only for flow rate control
   bool solve_for_bubble_speed_control = false;

//Imposing the bubble speed
   if(solve_for_bubble_speed_control)
     {
      Problem_Parameter::Prescribe_area=true;
      Problem_Parameter::Impose_U_bubble_tip=true;
      problem.make_steady();
      problem.adapt();
      Problem_Parameter::Prescribed_area = area_backup;
      problem.steady_newton_solve(0);
//      problem.doc_solution();
      oomph_info << "Start eigenproblem\n";   
      problem.solve_for_eigenproblem();
      oomph_info << "Done eigenproblem\n";
     }

//Imposing the flow rate      
//Optionally time evolution
   if(solve_for_flow_rate_control)
     {
      Problem_Parameter::Prescribe_area=true;
      Problem_Parameter::Impose_U_bubble_tip=false;
      problem.make_steady();
      problem.adapt();
      Problem_Parameter::Prescribed_area = area_backup;
      problem.steady_newton_solve(0);
//      problem.doc_solution();
      oomph_info << "Start eigenproblem\n";   
      problem.solve_for_eigenproblem();
      oomph_info << "Done eigenproblem\n";

      //Time evolution
      if(time_evolution_for_perturbed_solution)
        {
         epsilon_t = 1.0e-3;
         problem.make_unsteady();
         problem.set_timestepper_for_all_data(problem.time_stepper_pt(0));
         problem.initialise_dt(1.0e-4);
         problem.assign_initial_values_impulsive();
         problem.adapt();
         Problem_Parameter::Prescribed_area = area_backup;

         problem.doc_solution();
         problem.adaptive_unsteady_newton_solve(1.0e-4,epsilon_t,true);
         problem.doc_solution();

     bool first = false; //true; in the first step
     for(unsigned i = 0; i < nstep; i++)
      {   
       // Allow adaptation every few steps
       if( (i>0) && (Problem_Parameter::Adaptation_frequency>0) &&
           (i%Problem_Parameter::Adaptation_frequency == 0) )
        {
         problem.adapt();
         Problem_Parameter::Prescribed_area = area_backup;
         problem.adaptive_unsteady_newton_solve(1.0e-4,epsilon_t,false);
//         problem.unsteady_newton_solve(1.0e-4,0,first);          
         first = false;
        }

       double next_dt = dt;
       if (CommandLineArgs::command_line_flag_has_been_set
           ("--suppress_temporal_adaptivity"))
        {
         problem.unsteady_newton_solve(dt,0,first);
         }
       else
        {
         next_dt=
          problem.adaptive_unsteady_newton_solve(dt,epsilon_t,first);
         oomph_info << "Suggested next dt: " << next_dt << std::endl;
         }
       
       dt = next_dt;
       if(dt>0.1)
       {
        dt=0.1;
       }    
       problem.doc_solution();
      }//End of time evolution
        }
     }

   exit(1);
   }//end of stability



  //Start of continuation over Flow rate
   if(continuation_over_Flow_rate)
     {

   Problem_Parameter::Displacement_continuation = false;
   Problem_Parameter::Prescribe_area=true; 
      
   problem.restart(unsteady_initial_config_steady);
   oomph_info << "Done restart\n";   
   problem.doc_solution();

   problem.enable_parameter_update();
   Problem_Parameter::Update_ca_only=false;
   
//   Problem_Parameter::Prescribed_displacement=Problem_Parameter::Displacement_control_node_pt->value(0);
   Problem_Parameter::Prescribed_area=0.6;//Problem_Parameter::Downstream_area_pt->value(0);
   double area_backup = Problem_Parameter::Prescribed_area;
//   double displacement_backup = Problem_Parameter::Prescribed_displacement;
   double flow_rate_backup = 250.0*1.0e-6/60.0;//Problem_Parameter::Prescribed_flow_rate;
   problem.adapt();
   problem.steady_newton_solve(0);
   problem.doc_solution();
//   Problem_Parameter::Downstream_area_pt->set_value(0,area_backup);
//   Problem_Parameter::Prescribed_displacement = displacement_backup;
   Problem_Parameter::Prescribed_flow_rate = flow_rate_backup;
   Problem_Parameter::update_dependent_parameters();
   problem.steady_newton_solve(0);
   problem.doc_solution();

     {
      Problem_Parameter::Target_prescribed_flow_rate = 500.0*1.0e-6/60.0;
      Problem_Parameter::Flow_rate_continuation = true;
      Problem_Parameter::Prescribed_flow_rate_prev = 
       Problem_Parameter::Prescribed_flow_rate;
      Problem_Parameter::Max_increment = 20.0*1.0e-6/60.0;
      double desired_increment=5.0*1.0e-6/60.0; // minus: decrease Ca, plus: increase
//      if (CommandLineArgs::command_line_flag_has_been_set("--inc_thin_films"))            
//       {
//        desired_increment *= -1.0;
//       }
      double actual_increment=desired_increment;
      int count=0;
      while(Problem_Parameter::Prescribed_flow_rate < 
            Problem_Parameter::Target_prescribed_flow_rate)
       {
        if( (count!=Problem_Parameter::Adaptation_frequency) ||
            (Problem_Parameter::Adaptation_frequency < 0) )
         {
//          actual_increment = problem.arc_length_step_solve(
//           &Problem_Parameter::Prescribed_flow_rate,desired_increment,0);
          Problem_Parameter::Prescribed_flow_rate += desired_increment;
          Problem_Parameter::update_dependent_parameters();
          problem.steady_newton_solve(0);
          problem.doc_solution();
//          desired_increment=actual_increment;
          count++;
         }
        else
         {
          double area_backup = Problem_Parameter::Prescribed_area;
          double flow_rate_backup = Problem_Parameter::Prescribed_flow_rate;
          problem.adapt();
          problem.steady_newton_solve(0);
          Problem_Parameter::Prescribed_area = area_backup;
          Problem_Parameter::Prescribed_flow_rate = flow_rate_backup;
          problem.steady_newton_solve(0);
          count=0;
         }
        Problem_Parameter::Prescribed_flow_rate_prev = 
         Problem_Parameter::Prescribed_flow_rate;
       }
     }
     Problem_Parameter::Prescribed_flow_rate = 
      Problem_Parameter::Target_prescribed_flow_rate;
     Problem_Parameter::Flow_rate_continuation = false;
     problem.steady_newton_solve(0);
     problem.doc_solution();
     problem.disable_parameter_update();
     exit(1);      
     }

  ////continuation over bubble pressure
  if(continuation_over_bubble_pressure)
  {
   unsteady_initial_config_steady = true;
   problem.enable_parameter_update();
   Problem_Parameter::Update_ca_only=false;

   double U_bubble_backup = 0.1;
   double flow_rate_desired_backup = 3.16541*1.0e-6;
   double area_backup = 0.71385;
   double pressure_backup = 0.0602722;

   double Prescribed_bubble_pressure_target = 0.1;
   double Prescribed_area_target = 1.0;
   double first_increment = 0.0005; //0.01;//-0.0004;
   double desired_increment = first_increment; 
   double actual_increment = first_increment;

   Problem_Parameter::Prescribed_displacement = 
     Problem_Parameter::Displacement_control_node_pt->value(0);
   Problem_Parameter::Prescribed_flow_rate = flow_rate_desired_backup;
   Problem_Parameter::flow_rate_var_backup = flow_rate_desired_backup;
   Problem_Parameter::Prescribed_bubble_pressure = pressure_backup;
   Problem_Parameter::Prescribed_area = area_backup;
   Problem_Parameter::Downstream_area_backup = area_backup;

   Problem_Parameter::U_bubble_tip_dimensional_desired = U_bubble_backup;
   Problem_Parameter::Flow_rate_dimensional_desired = flow_rate_desired_backup;


   Problem_Parameter::Flux_constraint = false;
   Problem_Parameter::Prescribe_area = false;
   Problem_Parameter::Impose_U_bubble_tip=true;

   if(Problem_Parameter::Directory == "Pb_1")
    {
     U_bubble_backup = 0.1;
     flow_rate_desired_backup = 2.58783*1.0e-6;
     area_backup = 0.94;
     pressure_backup = 0.0759193;

     Problem_Parameter::Restart_file="Ainf_0.94_SF_pressure_control/restart_proc0_366.dat";
     problem.restart(unsteady_initial_config_steady);
     oomph_info << "Done restart\n";
     Problem_Parameter::Prescribed_area = area_backup;   
     problem.doc_solution();

     Prescribed_bubble_pressure_target = 0.08;
     Prescribed_area_target = 1.0;
     first_increment = 0.0001; //0.01;//-0.0004;
     desired_increment = first_increment; 
     actual_increment = first_increment;
     }

   if(Problem_Parameter::Directory == "Pb_2")
    {
     U_bubble_backup = 0.1;
     flow_rate_desired_backup = 3.16541*1.0e-6;
     area_backup = 0.71385;
     pressure_backup = 0.0602722;

     Problem_Parameter::Restart_file="Ainf_0.713_3_tip_displacement_control/restart_proc0_459.dat";
     problem.restart(unsteady_initial_config_steady);
     oomph_info << "Done restart\n";
     Problem_Parameter::Prescribed_area = area_backup;   
     problem.doc_solution();

     Prescribed_bubble_pressure_target = 0.09;
     Prescribed_area_target = 1.0;
     first_increment = 0.0005; //0.01;//-0.0004;
     desired_increment = first_increment; 
     actual_increment = first_increment;
     }
      
   Problem_Parameter::Flux_constraint = true;
   Problem_Parameter::Prescribe_area=true;
   Problem_Parameter::Impose_U_bubble_tip=true;

   Problem_Parameter::flow_rate_var_backup = flow_rate_desired_backup;

   problem.adapt();

   Problem_Parameter::Prescribed_displacement = 
     Problem_Parameter::Displacement_control_node_pt->value(0);
   Problem_Parameter::U_bubble_tip_dimensional_desired = U_bubble_backup;
   Problem_Parameter::Flow_rate_dimensional_desired = flow_rate_desired_backup;
   Problem_Parameter::Prescribed_flow_rate = flow_rate_desired_backup;
   Problem_Parameter::Prescribed_bubble_pressure = pressure_backup;
   Problem_Parameter::Prescribed_area = area_backup;
   Problem_Parameter::Downstream_area_backup = area_backup;

   problem.steady_newton_solve(0);
   problem.doc_solution();

   Problem_Parameter::Flux_constraint = true;
   Problem_Parameter::Prescribe_area=false;
   Problem_Parameter::Impose_U_bubble_tip=true;
   problem.adapt();

   Problem_Parameter::Prescribed_displacement = 
     Problem_Parameter::Displacement_control_node_pt->value(0);
   Problem_Parameter::U_bubble_tip_dimensional_desired = U_bubble_backup;
   Problem_Parameter::Flow_rate_dimensional_desired = flow_rate_desired_backup;
   Problem_Parameter::Prescribed_flow_rate = flow_rate_desired_backup;
   Problem_Parameter::Prescribed_bubble_pressure = pressure_backup;
   Problem_Parameter::Prescribed_area = area_backup;
   Problem_Parameter::Downstream_area_backup = area_backup;

   problem.steady_newton_solve(0);
   problem.doc_solution();


   double direction = 1.0;
   if(pressure_backup>Prescribed_bubble_pressure_target){direction = -1.0;}

   int count=0;
   while(direction*(Problem_Parameter::Prescribed_bubble_pressure - Prescribed_bubble_pressure_target) < 0.0 
         )
//   while(Problem_Parameter::Prescribed_area < //
//         Prescribed_area_target)              //
        {
        oomph_info<<"obacht "<<Problem_Parameter::Prescribed_bubble_pressure<<" out of "
                  <<Prescribed_bubble_pressure_target<<std::endl;       
        if( (count!=Problem_Parameter::Adaptation_frequency) ||
            (Problem_Parameter::Adaptation_frequency < 0) )
         {
          actual_increment = problem.arc_length_step_solve(
           &Problem_Parameter::Prescribed_bubble_pressure,desired_increment,0); 

//          Problem_Parameter::Prescribed_bubble_pressure += desired_increment; 
//          Problem_Parameter::Prescribed_area+= desired_increment;//
//          problem.steady_newton_solve(0);
          problem.doc_solution();

          if((actual_increment)*(actual_increment)>(first_increment)*(first_increment))
           {
            desired_increment = 
              sqrt((first_increment)*(first_increment))*
                sqrt((actual_increment)*(actual_increment))/(actual_increment);            
           }
          else
           {
            desired_increment = actual_increment;
           }

          count++;
         }
        else
         {
//          double area_backup = Problem_Parameter::Prescribed_area;//Problem_Parameter::Downstream_area_pt->value(0);
//          double displacement_backup = 
//           Problem_Parameter::Prescribed_displacement;
//          double flow_rate_backup = 
//           Problem_Parameter::Prescribed_flow_rate;
//          Problem_Parameter::Prescribed_area = area_backup;
          pressure_backup = Problem_Parameter::Prescribed_bubble_pressure;
          problem.adapt(); //preblem during the adaptation while in Area_continuation
          Problem_Parameter::Prescribed_bubble_pressure = pressure_backup;
//          Problem_Parameter::Prescribed_area = area_backup;
          problem.steady_newton_solve(0);
//          Problem_Parameter::Downstream_area_pt->set_value(0,area_backup);
//          Problem_Parameter::Prescribed_displacement = displacement_backup;
//          Problem_Parameter::Prescribed_flow_rate = flow_rate_backup;
//          problem.steady_newton_solve(0);
          count=0;
         }

       }
      Problem_Parameter::Prescribed_bubble_pressure = Prescribed_bubble_pressure_target;
//      Problem_Parameter::Prescribed_area = Prescribed_area_target;//
      problem.steady_newton_solve(0);
      problem.doc_solution();
      exit(1);





  }///End of continuation over bubble pressure

  //Start of continuation over area
  if(continuation_over_area)
  {
   unsteady_initial_config_steady = true;
   Problem_Parameter::Prescribe_area=true;
   problem.enable_parameter_update();
   Problem_Parameter::Update_ca_only=false;
   Problem_Parameter::Impose_U_bubble_tip=true;
   double U_bubble_backup = 0.1;
   double flow_rate_desired_backup = 2.40012*1.0e-6;
   double initial_area = 1.0;
   double area_backup = initial_area;
   Problem_Parameter::Prescribed_area = area_backup;

   Problem_Parameter::Target_prescribed_area=1.36;
   Problem_Parameter::Area_continuation=true;
   Problem_Parameter::Prescribed_area_prev = 
       Problem_Parameter::Prescribed_area;
   Problem_Parameter::Max_increment = 0.02;
   double first_increment = 0.002;
   double desired_increment = first_increment;
   double actual_increment = first_increment;
   double solve_eigen_problem = false;

   if(Problem_Parameter::Directory == "JOAO_1")
    {
     Problem_Parameter::Restart_file="RESLT_2/restart_proc0_414.dat";
     initial_area = 0.59;
     area_backup = initial_area;
     U_bubble_backup = 0.1;
     flow_rate_desired_backup = 3.29571*1.0e-6;
     Problem_Parameter::Prescribed_area = area_backup;
     Problem_Parameter::Target_prescribed_area=1.36;
     Problem_Parameter::Area_continuation=true;
     Problem_Parameter::Prescribed_area_prev = 
     Problem_Parameter::Prescribed_area;
     Problem_Parameter::Max_increment = 0.02;
     first_increment = 0.005;
     solve_eigen_problem = true;
    }

   if(Problem_Parameter::Directory == "JOAO_2")
    {
     Problem_Parameter::Restart_file="RESLT_2/restart_proc0_413.dat";
     initial_area = 0.59;
     area_backup = initial_area;
     U_bubble_backup = 0.1;
     flow_rate_desired_backup = 3.11803*1.0e-6;
     Problem_Parameter::Prescribed_area = area_backup;
     Problem_Parameter::Target_prescribed_area=1.36;
     Problem_Parameter::Area_continuation=true;
     Problem_Parameter::Prescribed_area_prev = 
     Problem_Parameter::Prescribed_area;
     Problem_Parameter::Max_increment = 0.02;
     first_increment = 0.005;
     solve_eigen_problem = true;
    }

   if(Problem_Parameter::Directory == "JOAO_3")
    {
//     Problem_Parameter::Restart_file="Lucie_values/Flow_rate_and_bubble_speed_control/Ub=0.1/2_tipped/RESLT_3/restart_proc0_453.dat";
     Problem_Parameter::Restart_file="JOAO_3_old/restart_proc0_473.dat";
     initial_area = 0.945082;
     area_backup = initial_area;
     U_bubble_backup = 0.1;
     flow_rate_desired_backup = 2.7516*1.0e-6;
     Problem_Parameter::Prescribed_area = area_backup;
     Problem_Parameter::Target_prescribed_area=0.36;
     Problem_Parameter::Area_continuation=true;
     Problem_Parameter::Prescribed_area_prev = 
     Problem_Parameter::Prescribed_area;
     Problem_Parameter::Max_increment = -0.02;
     first_increment = 0.002;
     solve_eigen_problem = false;
    }

   if(Problem_Parameter::Directory == "JOAO_4")
    {
     Problem_Parameter::Restart_file="Lucie_values/Flow_rate_and_bubble_speed_control/Ub=0.1/Max_el_area_0.03/Ainf_from_0.58_to_1.885_upper_branch/restart_proc0_364.dat";
     initial_area = 0.94;
     area_backup = initial_area;
     U_bubble_backup = 0.1;
     flow_rate_desired_backup = 2.58783*1.0e-6;
     Problem_Parameter::Prescribed_area = area_backup;
     Problem_Parameter::Target_prescribed_area=0.36;
     Problem_Parameter::Area_continuation=true;
     Problem_Parameter::Prescribed_area_prev = 
     Problem_Parameter::Prescribed_area;
     Problem_Parameter::Max_increment = -0.02;
     first_increment = 0.002;
     solve_eigen_problem = false;
    }

   if(Problem_Parameter::Directory == "JOAO_5")
    {
     Problem_Parameter::Restart_file="Lucie_values/Flow_rate_and_bubble_speed_control/Ub=0.1/Max_el_area_0.03/Ainf_from_0.43_to_0.725_lower_branch/restart_proc0_278.dat";
     initial_area = 0.72;
     area_backup = initial_area;
     U_bubble_backup = 0.1;
     flow_rate_desired_backup = 2.98017*1.0e-6/60.0;
     Problem_Parameter::Prescribed_area = area_backup;
     Problem_Parameter::Target_prescribed_area=1.36;
     Problem_Parameter::Area_continuation=true;
     Problem_Parameter::Prescribed_area_prev = 
     Problem_Parameter::Prescribed_area;
     Problem_Parameter::Max_increment = 0.02;
     first_increment = -0.002;
     solve_eigen_problem = false;
    }


   desired_increment = first_increment;
   actual_increment = first_increment;
   double direction = 1.0;
   if(initial_area>Problem_Parameter::Target_prescribed_area){direction = -1.0;}

   Problem_Parameter::U_bubble_tip_dimensional_desired = U_bubble_backup;
   Problem_Parameter::Flow_rate_dimensional_desired = flow_rate_desired_backup;

   Problem_Parameter::Prescribe_area=false;
   problem.restart(unsteady_initial_config_steady);
   oomph_info << "Done restart\n";
   Problem_Parameter::Prescribed_area = area_backup;   
   problem.doc_solution();
      
   Problem_Parameter::Prescribe_area=false;
   problem.make_steady();
   problem.adapt();

   Problem_Parameter::Prescribed_area = area_backup;
   problem.steady_newton_solve(0);
   problem.doc_solution();

   Problem_Parameter::Prescribe_area=true;
   problem.adapt();
   problem.doc_solution();
   oomph_info << "Done doc\n";
   problem.steady_newton_solve(0);
   problem.doc_solution();
   exit(1);

   if(Problem_Parameter::Prescribe_area)
   {
//   Problem_Parameter::Prescribed_displacement=Problem_Parameter::Displacement_control_node_pt->value(0);
   Problem_Parameter::Prescribed_area = area_backup;
   problem.steady_newton_solve(0);
   problem.doc_solution();

   Problem_Parameter::Prescribe_area=true;
   problem.adapt();
   problem.doc_solution();
   exit(1);

//   problem.adapt();
//   Problem_Parameter::Prescribed_area = area_backup;
//   problem.steady_newton_solve(0);
//   problem.doc_solution();

//   Problem_Parameter::Downstream_area_pt->set_value(0,area_backup);
//   Problem_Parameter::Prescribed_displacement = displacement_backup;
//   Problem_Parameter::Prescribed_flow_rate = flow_rate_backup;
//   problem.steady_newton_solve(0);
//   problem.doc_solution();   
   }
   else
   {
   Problem_Parameter::Displacement_continuation = false;
   Problem_Parameter::Prescribe_area = true;
//   Problem_Parameter::Prescribed_displacement=Problem_Parameter::Displacement_control_node_pt->value(0);
   Problem_Parameter::Prescribed_area=Problem_Parameter::Downstream_area_pt->value(0);
   problem.doc_solution();   
   area_backup = Problem_Parameter::Prescribed_area;
//   double displacement_backup = Problem_Parameter::Prescribed_displacement;
//   double flow_rate_backup = Problem_Parameter::Prescribed_flow_rate;
   problem.adapt();

   problem.steady_newton_solve(0);
   problem.doc_solution();
//   Problem_Parameter::Downstream_area_pt->set_value(0,area_backup);
//   Problem_Parameter::Prescribed_displacement = displacement_backup;
//   Problem_Parameter::Prescribed_flow_rate = flow_rate_backup;
//   problem.steady_newton_solve(0);
//   problem.doc_solution();
   
//   Problem_Parameter::Prescribe_area=true;
//   Problem_Parameter::Displacement_continuation = false;     
//   Problem_Parameter::Prescribed_area=Problem_Parameter::Downstream_area_pt->value(0);
//   problem.doc_solution();
     
//   problem.adapt();
//   problem.steady_newton_solve(0);
//   Problem_Parameter::Prescribed_area = area_backup;
//   Problem_Parameter::Prescribed_flow_rate = flow_rate_backup;     
   //problem.complete_problem_setup();
   //oomph_info<<"Number of equations: "<<problem.assign_eqn_numbers()<<std::endl;
//   problem.steady_newton_solve(0); //1
//   problem.doc_solution();      
   }     

   int count=0;
   while( direction*(Problem_Parameter::Target_prescribed_area-Problem_Parameter::Prescribed_area) > 0.0 )
        {
        oomph_info<<"obacht "<<Problem_Parameter::Prescribed_area<<" out of "
                  <<Problem_Parameter::Target_prescribed_area<<std::endl;       
        if( (count!=Problem_Parameter::Adaptation_frequency) ||
            (Problem_Parameter::Adaptation_frequency < 0) )
         {

//          actual_increment = problem.arc_length_step_solve(
//           &Problem_Parameter::Prescribed_area,desired_increment,0);
//
//          if(actual_increment > sqrt(first_increment*first_increment) )
//           {
//            desired_increment = sqrt(first_increment*first_increment);
//           }
//          else
//           {
//             if(actual_increment < -sqrt(first_increment*first_increment) )
//              {
//               desired_increment = -sqrt(first_increment*first_increment);
//              }
//             else
//              {
//               desired_increment=actual_increment;
//              } 
//           }


//          desired_increment=actual_increment;

          Problem_Parameter::Prescribed_area += desired_increment; 
          problem.steady_newton_solve(0);
          problem.doc_solution();
          count++;

          if(solve_eigen_problem)
           {
            Problem_Parameter::Impose_U_bubble_tip=false;            
            if( !(Problem_Parameter::flow_rate_var_pt == 0) )
             {
              Problem_Parameter::Flow_rate_dimensional_desired = 
                Problem_Parameter::flow_rate_var_pt->value(0);
             }
            area_backup = Problem_Parameter::Prescribed_area;
            problem.adapt();
            Problem_Parameter::Prescribed_area = area_backup;
            problem.steady_newton_solve(0);
            
            oomph_info << "Start eigenproblem\n";   
            problem.solve_for_eigenproblem();
            oomph_info << "Done eigenproblem\n";
            
            Problem_Parameter::Impose_U_bubble_tip=true;
            Problem_Parameter::U_bubble_tip_dimensional_desired = U_bubble_backup;
            problem.adapt();
            Problem_Parameter::Prescribed_area = area_backup;
           }

         }
        else
         {
          if(!solve_eigen_problem)
           {
            area_backup = Problem_Parameter::Prescribed_area;
//            Problem_Parameter::Downstream_area_pt->value(0);
//            double displacement_backup = Problem_Parameter::Prescribed_displacement;
//            Problem_Parameter::Downstream_area_pt->set_value(0,area_backup);
            problem.adapt(); 
            Problem_Parameter::Prescribed_area = area_backup;
            problem.steady_newton_solve(0);
           }
          count=0;
         }
        Problem_Parameter::Prescribed_area_prev = 
         Problem_Parameter::Prescribed_area;
       }
      Problem_Parameter::Prescribed_area = 
       Problem_Parameter::Target_prescribed_area;
      Problem_Parameter::Area_continuation = false;
      problem.steady_newton_solve(0);
      problem.doc_solution();
      exit(1);  
  
  }//End of continuation over area


  //Start of continuation over bubble tip velocity
  if(continuation_over_U_bubble_tip)
   {
   unsteady_initial_config_steady = true;
   Problem_Parameter::Prescribe_area=true;
   problem.enable_parameter_update();
   Problem_Parameter::Update_ca_only=false;
   Problem_Parameter::Impose_U_bubble_tip=true;
   Problem_Parameter::U_bubble_tip_dimensional_desired=0.1;

   problem.restart(unsteady_initial_config_steady);
   oomph_info << "Done restart\n";   
   problem.doc_solution();
      
   Problem_Parameter::Prescribe_area=true;
   problem.make_steady();
   problem.adapt();

   double area_backup = 0.55; //Problem_Parameter::Downstream_area_pt->value(0);

   if(Problem_Parameter::Prescribe_area)
    {
     Problem_Parameter::Prescribed_area = area_backup;
     problem.steady_newton_solve(0);
     problem.doc_solution();

     problem.adapt();
     Problem_Parameter::Prescribed_area = area_backup;
     problem.steady_newton_solve(0);
     problem.doc_solution();
     }


   Problem_Parameter::Max_increment = -0.01;
   double Target_U_bubble_tip_dimensional_desired=0;
   double U_bubble_prev=Problem_Parameter::U_bubble_tip_dimensional_desired;
   double desired_increment=-0.001;
   double actual_increment=desired_increment;

   int count=0;
   while(Problem_Parameter::U_bubble_tip_dimensional_desired > 
         Target_U_bubble_tip_dimensional_desired)
        {
        oomph_info<<"obacht "<<Problem_Parameter::U_bubble_tip_dimensional_desired<<" out of "
                  <<Target_U_bubble_tip_dimensional_desired<<std::endl;       
        if( (count!=Problem_Parameter::Adaptation_frequency) ||
            (Problem_Parameter::Adaptation_frequency < 0) )
         {
//          actual_increment = problem.arc_length_step_solve(
//           &Problem_Parameter::U_bubble_tip_dimensional_desired,desired_increment,0); 
          Problem_Parameter::U_bubble_tip_dimensional_desired += desired_increment; 
          problem.steady_newton_solve(0);
          problem.doc_solution();
//          desired_increment=actual_increment;
          count++;
         }
        else
         {
//          double area_backup = Problem_Parameter::Prescribed_area;//Problem_Parameter::Downstream_area_pt->value(0);
//          double displacement_backup = 
//           Problem_Parameter::Prescribed_displacement;
//          double flow_rate_backup = 
//           Problem_Parameter::Prescribed_flow_rate;
          Problem_Parameter::Prescribed_area = area_backup;
          problem.adapt(); //preblem during the adaptation while in Area_continuation
          Problem_Parameter::Prescribed_area = area_backup;
          problem.steady_newton_solve(0);
//          Problem_Parameter::Downstream_area_pt->set_value(0,area_backup);
//          Problem_Parameter::Prescribed_displacement = displacement_backup;
//          Problem_Parameter::Prescribed_flow_rate = flow_rate_backup;
//          problem.steady_newton_solve(0);
          count=0;
         }
        U_bubble_prev = 
         Problem_Parameter::U_bubble_tip_dimensional_desired;
       }
      Problem_Parameter::U_bubble_tip_dimensional_desired = 
        Target_U_bubble_tip_dimensional_desired;
      problem.steady_newton_solve(0);
      problem.doc_solution();
      exit(1);    

    }//End of continuation over bubble tip velocity
 
  //Start of continuation over Displacement
  if(continuation_over_Displacement)
  {
//   Problem_Parameter::Q=200.0;
   problem.restart(unsteady_initial_config_steady);
   oomph_info << "Done restart\n";   
   problem.doc_solution();

   Problem_Parameter::Prescribe_area=false;
   Problem_Parameter::Prescribed_displacement=Problem_Parameter::Displacement_control_node_pt->value(0);
   Problem_Parameter::Prescribed_area=Problem_Parameter::Downstream_area_pt->value(0);
   double area_backup = Problem_Parameter::Prescribed_area;
   double displacement_backup = 
       Problem_Parameter::Prescribed_displacement;
   double flow_rate_backup = 
       Problem_Parameter::Prescribed_flow_rate;
   problem.adapt();
   problem.steady_newton_solve(0);
   Problem_Parameter::Downstream_area_pt->set_value(0,area_backup);
   Problem_Parameter::Prescribed_displacement = displacement_backup;
   Problem_Parameter::Prescribed_flow_rate = flow_rate_backup;
   problem.steady_newton_solve(0);
   problem.doc_solution();
   
   Problem_Parameter::Target_prescribed_displacement=0.1;
   Problem_Parameter::Target_prescribed_area=0.36;   
   Problem_Parameter::Displacement_continuation = true;
   Problem_Parameter::Prescribed_displacement_prev = 
       Problem_Parameter::Prescribed_displacement;
   Problem_Parameter::Max_increment = -1.5e-2;
   double desired_increment=-1.5e-3; // obacht usually negative
   double actual_increment=desired_increment;
   int count=0;
//   while(Problem_Parameter::Prescribed_displacement < 
//         Problem_Parameter::Target_prescribed_displacement)
   while(Problem_Parameter::Downstream_area_pt->value(0) > 
         Problem_Parameter::Target_prescribed_area)

        {
        oomph_info<<"obacht "<<Problem_Parameter::Prescribed_displacement<<" out of "
                  <<Problem_Parameter::Target_prescribed_displacement<<std::endl;       
        if( (count!=Problem_Parameter::Adaptation_frequency) ||
            (Problem_Parameter::Adaptation_frequency < 0) )
         {
//          actual_increment = problem.arc_length_step_solve(
//           &Problem_Parameter::Prescribed_displacement,desired_increment,0);
          Problem_Parameter::Prescribed_displacement += desired_increment; 
          problem.steady_newton_solve(0);          
          problem.doc_solution();
          desired_increment=actual_increment;
          count++;
         }
        else
         {
          double area_backup = Problem_Parameter::Downstream_area_pt->value(0);
          double displacement_backup = 
           Problem_Parameter::Prescribed_displacement;
          double flow_rate_backup = 
           Problem_Parameter::Prescribed_flow_rate;
          problem.adapt();
          problem.steady_newton_solve(0);
          Problem_Parameter::Downstream_area_pt->set_value(0,area_backup);
          Problem_Parameter::Prescribed_displacement = displacement_backup;
          Problem_Parameter::Prescribed_flow_rate = flow_rate_backup;
          problem.steady_newton_solve(0);
          count=0;
         }
        Problem_Parameter::Prescribed_displacement_prev = 
         Problem_Parameter::Prescribed_displacement;
       }
      Problem_Parameter::Prescribed_displacement = 
       Problem_Parameter::Target_prescribed_displacement;
      Problem_Parameter::Displacement_continuation = false;
      problem.steady_newton_solve(0);
      problem.doc_solution();
      exit(1);  
  }//End of continuation over Displacement
  
  //Start of continuation over Aspect ratio
  if(continuation_over_Aspect_ratio)
  {
   Problem_Parameter::Q=200.0;
   Problem_Parameter::Aspect_ratio = 0.055;
   problem.restart(unsteady_initial_config_steady);
   oomph_info << "Done restart\n";   
   problem.doc_solution();

   Problem_Parameter::Prescribe_area=false;
   Problem_Parameter::Prescribed_displacement=Problem_Parameter::Displacement_control_node_pt->value(0);
   Problem_Parameter::Prescribed_area=Problem_Parameter::Downstream_area_pt->value(0);
   double area_backup = Problem_Parameter::Prescribed_area;
   double displacement_backup = 
       Problem_Parameter::Prescribed_displacement;
   double flow_rate_backup = 
       Problem_Parameter::Prescribed_flow_rate;
   problem.adapt();
   problem.steady_newton_solve(0);
   problem.doc_solution();
   Problem_Parameter::Downstream_area_pt->set_value(0,area_backup);
   Problem_Parameter::Prescribed_displacement = displacement_backup;
   Problem_Parameter::Prescribed_flow_rate = flow_rate_backup;
   problem.steady_newton_solve(0);
   problem.doc_solution();
      
   Problem_Parameter::Aspect_ratio_continuation = true;
   Problem_Parameter::Aspect_ratio_prev = Problem_Parameter::Aspect_ratio;
   Problem_Parameter::Max_increment = -1.0;
   double desired_increment=-1.0e-2; // obacht usually negative
   double actual_increment=desired_increment;
   int count=0;
   while(Problem_Parameter::Aspect_ratio > 
            Problem_Parameter::Target_aspect_ratio)
       //while(Problem_Parameter::Aspect_ratio_continuation)
       {
        if( (count!=Problem_Parameter::Adaptation_frequency) ||
            (Problem_Parameter::Adaptation_frequency < 0) )
         {
//          actual_increment = problem.arc_length_step_solve(
//           &Problem_Parameter::Aspect_ratio,desired_increment,0);
          Problem_Parameter::Aspect_ratio += desired_increment; 
          problem.steady_newton_solve(0);
          problem.doc_solution();
          desired_increment=actual_increment;
          count++;
         }
        else
         {
          double area_backup = Problem_Parameter::Downstream_area_pt->value(0);
          double displacement_backup = 
           Problem_Parameter::Prescribed_displacement;
          double flow_rate_backup = 
           Problem_Parameter::Prescribed_flow_rate;
          problem.adapt();
          problem.steady_newton_solve(0);
          Problem_Parameter::Downstream_area_pt->set_value(0,area_backup);
          Problem_Parameter::Prescribed_displacement = displacement_backup;
          Problem_Parameter::Prescribed_flow_rate = flow_rate_backup;
          problem.steady_newton_solve(0);
          count=0;
         }
        Problem_Parameter::Aspect_ratio_prev = Problem_Parameter::Aspect_ratio;
       }
      
      Problem_Parameter::Aspect_ratio = Problem_Parameter::Target_aspect_ratio;
      Problem_Parameter::Aspect_ratio_continuation = false;
      problem.steady_newton_solve(0);
      problem.doc_solution();
      exit(1);
  }//End of continuation over Aspect ratio

  //Start of continuation over FSI
  if(continuation_over_FSI)
  {
   Problem_Parameter::Q=200.0;
   problem.restart(unsteady_initial_config_steady);
   oomph_info << "Done restart\n";   
   problem.doc_solution();

   Problem_Parameter::Prescribe_area=false;
   Problem_Parameter::Prescribed_displacement=Problem_Parameter::Displacement_control_node_pt->value(0);
   Problem_Parameter::Prescribed_area=Problem_Parameter::Downstream_area_pt->value(0);
   double area_backup = Problem_Parameter::Prescribed_area;
   double displacement_backup = 
       Problem_Parameter::Prescribed_displacement;
   double flow_rate_backup = 
       Problem_Parameter::Prescribed_flow_rate;
   problem.adapt();
   problem.steady_newton_solve(0);
   problem.doc_solution();
   Problem_Parameter::Downstream_area_pt->set_value(0,area_backup);
   Problem_Parameter::Prescribed_displacement = displacement_backup;
   Problem_Parameter::Prescribed_flow_rate = flow_rate_backup;
   problem.steady_newton_solve(0);
   problem.doc_solution();
   
   Problem_Parameter::FSI_continuation = true;
   Problem_Parameter::Q_prev = Problem_Parameter::Q;
   Problem_Parameter::Max_increment = 1.0e6;
   double desired_increment=3.0e1;
   double actual_increment=desired_increment;
   int count=0;
   while(Problem_Parameter::Q < Problem_Parameter::Target_q)
     //while(Problem_Parameter::FSI_continuation)
       {
        oomph_info<<"obacht "<<Problem_Parameter::Q<<" out of "
                  <<Problem_Parameter::Target_q<<std::endl;
        if( (count!=Problem_Parameter::Adaptation_frequency) ||
            (Problem_Parameter::Adaptation_frequency < 0) )
         {
          //actual_increment = problem.arc_length_step_solve(
           //&Problem_Parameter::Q,desired_increment,0);
          Problem_Parameter::Q += desired_increment;
          problem.steady_newton_solve(0);
          problem.doc_solution();
          desired_increment *= 1.5;
          count++;
         }
        else
         {
          double area_backup = Problem_Parameter::Downstream_area_pt->value(0);
          double displacement_backup = 
           Problem_Parameter::Prescribed_displacement;
          double flow_rate_backup = 
           Problem_Parameter::Prescribed_flow_rate;
          problem.adapt();
          problem.steady_newton_solve(0);
          Problem_Parameter::Downstream_area_pt->set_value(0,area_backup);
          Problem_Parameter::Prescribed_displacement = displacement_backup;
          Problem_Parameter::Prescribed_flow_rate = flow_rate_backup;
          problem.steady_newton_solve(0);
          count=0;
         }
        Problem_Parameter::Q_prev = Problem_Parameter::Q;
       }
     Problem_Parameter::Q = Problem_Parameter::Target_q;
     Problem_Parameter::FSI_continuation = false;
     problem.steady_newton_solve(0);
     problem.doc_solution();
     exit(1);   
  }//End of continuation over FSI

   //pause("done");
  }
 // No restart
 else
  {
   if(CommandLineArgs::command_line_flag_has_been_set("--steady"))
    {
     problem.make_steady();
     unsigned max_adapt=0;
     double tol=1.0e-4;
     
     oomph_info << "Target FSI parameter: "
                <<Problem_Parameter::Target_q<<std::endl;
     
     double pre_stress_backup = Problem_Parameter::Pre_stress;
     
     // Steady solve without FSI
     Problem_Parameter::Q = 1.0e-2;
     Problem_Parameter::Pre_stress = 0.0;
     Problem_Parameter::Aspect_ratio = 0.1;
     Problem_Parameter::Prescribed_flow_rate = 81.0*1.0e-6/60.0;
     Problem_Parameter::update_dependent_parameters();
     
     //if (CommandLineArgs::command_line_flag_has_been_set("--saffman_taylor"))
     // { 
     // problem.adapt();
     // }
     
     //Problem_Parameter::P_bubble_pt->pin(0);
     //Problem_Parameter::P_bubble_pt->set_value(0,0.0);
     //dynamic_cast<SolidNode*>(Problem_Parameter::Tip_node_pt)->pin_position(0);
     // Problem_Parameter::Downstream_flux_pt->pin(0);
     // Problem_Parameter::Downstream_flux_pt->set_value(0,flux);
     
     // Output initial conditions
     problem.doc_solution();
     
     // Do steady Newton solve first to get steady equilibrium shape
     // with perturbed interface shape
     // problem.unsteady_newton_solve(0.1);
     //problem.complete_problem_setup();
     //problem.pin_pressure_dof();
     // problem.pin_fvk();
     // oomph_info<<"Number of equations: "<<problem.assign_eqn_numbers()<<std::endl;
     
     // Problem_Parameter::Amplitude += 1.0e-2;
     // problem.update_pinned_fvk_displacements();
     // oomph_info<<"Number of equations: "<<problem.assign_eqn_numbers()<<std::endl;
     // problem.doc_solution();
     
     problem.steady_newton_solve(max_adapt); //1
     
     // Output steady solution
     problem.doc_solution();
     //problem.adapt();
     
     // for(unsigned i=0; i<10;i++)
     //  {
     //   //Problem_Parameter::Slope -= 0.01;
     //   Problem_Parameter::Amplitude += 1.0e-2;
     //   problem.update_pinned_fvk_displacements();
     //   problem.steady_newton_solve(0);
     //   problem.doc_solution();
     //  }
     
     //Problem_Parameter::P_bubble_pt->pin(0);
     //Problem_Parameter::P_bubble_pt->set_value(0,0.0);
     //oomph_info<<"Number of equations: "<<problem.assign_eqn_numbers()<<std::endl;
     
     //problem.steady_newton_solve(0);
     
     // Output steady solution
     //problem.doc_solution();
     //problem.set_initial_guess_for_flux();
     //problem.pin_flux();
     //problem.complete_problem_setup();
     //oomph_info<<"Number of equations: "<<problem.assign_eqn_numbers()<<std::endl;
     
     //problem.steady_newton_solve(0); //1
     
     // Output steady solution
     //problem.doc_solution();
     
     // problem.adapt();
     // problem.steady_newton_solve(0);
     // problem.doc_solution();
     
     //exit(1);
     
     if (CommandLineArgs::command_line_flag_has_been_set("--inc_thin_films")) 
      {
       double desired_increment=1.0e-1;
       double actual_increment=desired_increment;
       int count=0;
       while(Problem_Parameter::Thin_film_homotopy_parameter < 1.0)
        {
         if( (count!=Problem_Parameter::Adaptation_frequency) ||
             (Problem_Parameter::Adaptation_frequency < 0) )
          {
           actual_increment = problem.arc_length_step_solve(
            &Problem_Parameter::Thin_film_homotopy_parameter,desired_increment,0);
           problem.doc_solution();
           desired_increment=actual_increment;
           count++;
          }
         else
          {
           double area_backup = Problem_Parameter::Downstream_area_pt->value(0);
           double displacement_backup = 
            Problem_Parameter::Prescribed_displacement;
           double flow_rate_backup = 
            Problem_Parameter::Prescribed_flow_rate;
           problem.adapt();
           problem.steady_newton_solve(0);
           Problem_Parameter::Downstream_area_pt->set_value(0,area_backup);
           Problem_Parameter::Prescribed_displacement = displacement_backup;
           Problem_Parameter::Prescribed_flow_rate = flow_rate_backup;
           problem.steady_newton_solve(0);
           count=0;
          }
        }
       Problem_Parameter::Thin_film_homotopy_parameter = 1.0;
       problem.steady_newton_solve(max_adapt);
       problem.doc_solution();
      }
     //exit(1);
     //========================================================================
     problem.enable_parameter_update();
     {
      Problem_Parameter::Target_prescribed_flow_rate = 187.0*1.0e-6/60.0;                 //////Joao ====>>> Before 187.0*1.0e-6/60.0
      Problem_Parameter::Flow_rate_continuation = true;
      Problem_Parameter::Prescribed_flow_rate_prev = 
       Problem_Parameter::Prescribed_flow_rate;
      Problem_Parameter::Max_increment = 10.0;
      double desired_increment=1.0e-1; // minus: decrease Ca, plus: increase
      // if (CommandLineArgs::command_line_flag_has_been_set("--inc_thin_films"))
      //  {
      //   desired_increment *= -1.0;
      //  }
      double actual_increment=desired_increment;
      int count=0;
      while(Problem_Parameter::Prescribed_flow_rate < 
            Problem_Parameter::Target_prescribed_flow_rate)
       //while(Problem_Parameter::Bubble_speed_continuation)
       {
        if( (count!=Problem_Parameter::Adaptation_frequency) ||
            (Problem_Parameter::Adaptation_frequency < 0) )
         {
          actual_increment = problem.arc_length_step_solve(
           &Problem_Parameter::Prescribed_flow_rate,desired_increment,0);
          problem.doc_solution();
          desired_increment=actual_increment;
          count++;
         }
        else
         {
          double area_backup = Problem_Parameter::Downstream_area_pt->value(0);
          double displacement_backup = 
           Problem_Parameter::Prescribed_displacement;
          double flow_rate_backup = 
           Problem_Parameter::Prescribed_flow_rate;
          problem.adapt();
          problem.steady_newton_solve(0);
          Problem_Parameter::Downstream_area_pt->set_value(0,area_backup);
          Problem_Parameter::Prescribed_displacement = displacement_backup;
          Problem_Parameter::Prescribed_flow_rate = flow_rate_backup;
          problem.steady_newton_solve(0);
          count=0;
         }
        Problem_Parameter::Prescribed_flow_rate_prev = 
         Problem_Parameter::Prescribed_flow_rate;
       }
     }
     Problem_Parameter::Prescribed_flow_rate = 
      Problem_Parameter::Target_prescribed_flow_rate;
     Problem_Parameter::Flow_rate_continuation = false;
     problem.steady_newton_solve(max_adapt);
     problem.doc_solution();
     problem.disable_parameter_update();
     
     //========================================================================
     
     // Re-enable FSI
     Problem_Parameter::Q = 1.0e1; //1.0e-3;
     Problem_Parameter::Pre_stress = pre_stress_backup; 
     
     // Need to reset the BC
     //problem.complete_problem_setup();
     //oomph_info<<"Number of equations: "<<problem.assign_eqn_numbers()<<std::endl;
     
     //problem.adapt();
     problem.steady_newton_solve(max_adapt);
     
     // Output steady solution
     problem.doc_solution();
     //exit(1);
     
     //========================================================================
     
     double q_backup = Problem_Parameter::Target_q;
     {
      Problem_Parameter::Target_q = 2.0e2;
      Problem_Parameter::FSI_continuation = true;
      Problem_Parameter::Q_prev = Problem_Parameter::Q;
      Problem_Parameter::Max_increment = 1.0e6;
      double desired_increment=1.0e1;
      double actual_increment=desired_increment;
      int count=0;
      while(Problem_Parameter::Q < Problem_Parameter::Target_q-tol)
       //while(Problem_Parameter::FSI_continuation)
       {
        oomph_info<<"obacht "<<Problem_Parameter::Q<<" out of "
                  <<Problem_Parameter::Target_q<<std::endl;
        if( (count!=Problem_Parameter::Adaptation_frequency) ||
            (Problem_Parameter::Adaptation_frequency < 0) )
         {
          actual_increment = problem.arc_length_step_solve(
           &Problem_Parameter::Q,desired_increment,0);
          problem.doc_solution();
          desired_increment=actual_increment;
          count++;
         }
        else
         {
         //problem.steady_newton_solve(max_adapt);
          double area_backup = Problem_Parameter::Downstream_area_pt->value(0);
          double displacement_backup = 
           Problem_Parameter::Prescribed_displacement;
          double flow_rate_backup = 
           Problem_Parameter::Prescribed_flow_rate;
          problem.adapt();
          problem.steady_newton_solve(0);
          Problem_Parameter::Downstream_area_pt->set_value(0,area_backup);
          Problem_Parameter::Prescribed_displacement = displacement_backup;
          Problem_Parameter::Prescribed_flow_rate = flow_rate_backup;
          problem.steady_newton_solve(0);
          //problem.doc_solution();
          //exit(1);
          count=0;
         }
        Problem_Parameter::Q_prev = Problem_Parameter::Q;
       }
     }
     Problem_Parameter::Q = Problem_Parameter::Target_q;
     Problem_Parameter::FSI_continuation = false;
     problem.steady_newton_solve(max_adapt);
     problem.doc_solution();
     //exit(1);
     
     //========================================================================
     { 
      Problem_Parameter::Target_prescribed_displacement=0.02;
      Problem_Parameter::Displacement_continuation = true;
      Problem_Parameter::Prescribed_displacement_prev = 
       Problem_Parameter::Prescribed_displacement;
      Problem_Parameter::Max_increment = 1.0;
      double desired_increment=2.0e-3; // obacht usually negative
      double actual_increment=desired_increment;
      int count=0;
      while(Problem_Parameter::Prescribed_displacement < 
            Problem_Parameter::Target_prescribed_displacement-tol)
       //while(Problem_Parameter::Displacement_continuation)
       {
        oomph_info<<"obacht "<<Problem_Parameter::Prescribed_displacement<<" out of "
                  <<Problem_Parameter::Target_prescribed_displacement<<std::endl;       
        if( (count!=Problem_Parameter::Adaptation_frequency) ||
            (Problem_Parameter::Adaptation_frequency < 0) )
         {
//          actual_increment = problem.arc_length_step_solve(
//           &Problem_Parameter::Prescribed_displacement,desired_increment,0);
          Problem_Parameter::Prescribed_displacement += desired_increment; 
          problem.steady_newton_solve(0);          
          problem.doc_solution();
          desired_increment=actual_increment;
          count++;
         }
        else
         {
          double area_backup = Problem_Parameter::Downstream_area_pt->value(0);
          double displacement_backup = 
           Problem_Parameter::Prescribed_displacement;
          double flow_rate_backup = 
           Problem_Parameter::Prescribed_flow_rate;
          problem.adapt();
          problem.steady_newton_solve(0);
          Problem_Parameter::Downstream_area_pt->set_value(0,area_backup);
          Problem_Parameter::Prescribed_displacement = displacement_backup;
          Problem_Parameter::Prescribed_flow_rate = flow_rate_backup;
          problem.steady_newton_solve(0);
          count=0;
         }
        Problem_Parameter::Prescribed_displacement_prev = 
         Problem_Parameter::Prescribed_displacement;
       }
      Problem_Parameter::Prescribed_displacement = 
       Problem_Parameter::Target_prescribed_displacement;
      Problem_Parameter::Displacement_continuation = false;
      problem.steady_newton_solve(max_adapt);
      problem.doc_solution();
      //exit(1);
     }
     //========================================================================
     {
      Problem_Parameter::Aspect_ratio_continuation = true;
      Problem_Parameter::Aspect_ratio_prev = Problem_Parameter::Aspect_ratio;
      Problem_Parameter::Max_increment = -1.0;
      double desired_increment=-1.0e-2; // obacht usually negative
      double actual_increment=desired_increment;
      int count=0;
      while(Problem_Parameter::Aspect_ratio > 
            Problem_Parameter::Target_aspect_ratio+tol)
       //while(Problem_Parameter::Aspect_ratio_continuation)
       {
        if( (count!=Problem_Parameter::Adaptation_frequency) ||
            (Problem_Parameter::Adaptation_frequency < 0) )
         {
          actual_increment = problem.arc_length_step_solve(
           &Problem_Parameter::Aspect_ratio,desired_increment,0);
          problem.doc_solution();
          desired_increment=actual_increment;
          count++;
         }
        else
         {
          double area_backup = Problem_Parameter::Downstream_area_pt->value(0);
          double displacement_backup = 
           Problem_Parameter::Prescribed_displacement;
          double flow_rate_backup = 
           Problem_Parameter::Prescribed_flow_rate;
          problem.adapt();
          problem.steady_newton_solve(0);
          Problem_Parameter::Downstream_area_pt->set_value(0,area_backup);
          Problem_Parameter::Prescribed_displacement = displacement_backup;
          Problem_Parameter::Prescribed_flow_rate = flow_rate_backup;
          problem.steady_newton_solve(0);
          count=0;
         }
        Problem_Parameter::Aspect_ratio_prev = Problem_Parameter::Aspect_ratio;
       }
      
      Problem_Parameter::Aspect_ratio = Problem_Parameter::Target_aspect_ratio;
      Problem_Parameter::Aspect_ratio_continuation = false;
      problem.steady_newton_solve(max_adapt);
      problem.doc_solution();
     }
     //========================================================================
     {
      Problem_Parameter::Target_prescribed_displacement=0.1;
      Problem_Parameter::Displacement_continuation = true;
      Problem_Parameter::Prescribed_displacement_prev = 
       Problem_Parameter::Prescribed_displacement;
      Problem_Parameter::Max_increment = 1.0;
      double desired_increment=1.0e-2; // obacht usually negative
      double actual_increment=desired_increment;
      int count=0;
      while(Problem_Parameter::Prescribed_displacement < 
          Problem_Parameter::Target_prescribed_displacement-tol)
       //while(Problem_Parameter::Displacement_continuation)
       {
        oomph_info<<"obacht "<<Problem_Parameter::Prescribed_displacement<<" out of "
                   <<Problem_Parameter::Target_prescribed_displacement<<std::endl;
        if( (count!=Problem_Parameter::Adaptation_frequency) ||
            (Problem_Parameter::Adaptation_frequency < 0) )
         {
          actual_increment = problem.arc_length_step_solve(
           &Problem_Parameter::Prescribed_displacement,desired_increment,0);
          problem.doc_solution();
          desired_increment=actual_increment;
          count++;
         }
        else
         {
          double area_backup = Problem_Parameter::Downstream_area_pt->value(0);
          double displacement_backup = 
           Problem_Parameter::Prescribed_displacement;
          double flow_rate_backup = 
           Problem_Parameter::Prescribed_flow_rate;
          problem.adapt();
          problem.steady_newton_solve(0);
          Problem_Parameter::Downstream_area_pt->set_value(0,area_backup);
          Problem_Parameter::Prescribed_displacement = displacement_backup;
          Problem_Parameter::Prescribed_flow_rate = flow_rate_backup;
          problem.steady_newton_solve(0);
          count=0;
         }
        Problem_Parameter::Prescribed_displacement_prev = 
         Problem_Parameter::Prescribed_displacement;
       }
      Problem_Parameter::Prescribed_displacement = 
       Problem_Parameter::Target_prescribed_displacement;
      Problem_Parameter::Displacement_continuation = false;
      problem.steady_newton_solve(max_adapt);
      problem.doc_solution();
      //exit(1);
     }
  
     //========================================================================
     {
      Problem_Parameter::Target_q = q_backup; //2.0e4;
      Problem_Parameter::FSI_continuation = true;
      Problem_Parameter::Q_prev = Problem_Parameter::Q;
      Problem_Parameter::Max_increment = 1.0e6;
      double desired_increment=1.0e-2;
      double actual_increment=desired_increment;
      int count=0;
      while(Problem_Parameter::Q < Problem_Parameter::Target_q-tol)
       //while(Problem_Parameter::FSI_continuation)
       {
        oomph_info<<"obacht "<<Problem_Parameter::Q<<" out of "
                  <<Problem_Parameter::Target_q<<std::endl;
        if( (count!=Problem_Parameter::Adaptation_frequency) ||
            (Problem_Parameter::Adaptation_frequency < 0) )
         {
          actual_increment = problem.arc_length_step_solve(
           &Problem_Parameter::Q,desired_increment,max_adapt);
          problem.doc_solution();
          desired_increment=actual_increment;
          count++;
         }
        else
         {
          double area_backup = Problem_Parameter::Downstream_area_pt->value(0);
          double displacement_backup = 
           Problem_Parameter::Prescribed_displacement;
          double flow_rate_backup = 
           Problem_Parameter::Prescribed_flow_rate;
          problem.adapt();
          problem.steady_newton_solve(0);
          Problem_Parameter::Downstream_area_pt->set_value(0,area_backup);
          Problem_Parameter::Prescribed_displacement = displacement_backup;
          Problem_Parameter::Prescribed_flow_rate = flow_rate_backup;
          problem.steady_newton_solve(0);
          count=0;
         }
        Problem_Parameter::Q_prev = Problem_Parameter::Q;
       }
      
      Problem_Parameter::Q = Problem_Parameter::Target_q;
      Problem_Parameter::FSI_continuation = false;
      problem.steady_newton_solve(max_adapt);
      problem.doc_solution();
      //exit(1);
     }
     
     //========================================================================
     Problem_Parameter::Target_prescribed_area = 1.5; //-0.475; //  -0.36     ////// Joao ==>> before 1.0                                  
     Problem_Parameter::Displacement_continuation = true;
     Problem_Parameter::Prescribed_displacement_prev = 
      Problem_Parameter::Prescribed_displacement;
     Problem_Parameter::Target_prescribed_displacement = 0.0;
     Problem_Parameter::Max_increment = -2.5e-3;
     {
      double desired_increment=-5.0e-3;
      // if (CommandLineArgs::command_line_flag_has_been_set("--inc_thin_films"))
      //   {
      //    desired_increment *= -1.0;
      //   }
      double actual_increment=desired_increment;
      int count=0;
      while(Problem_Parameter::Downstream_area_pt->value(0) > 
            Problem_Parameter::Target_prescribed_area-tol)
       {
        if( (count!=Problem_Parameter::Adaptation_frequency) ||
            (Problem_Parameter::Adaptation_frequency < 0) )
         {
          actual_increment = problem.arc_length_step_solve(
           &Problem_Parameter::Prescribed_displacement,desired_increment,0);
          problem.doc_solution();
          desired_increment=actual_increment;
          count++;
         }
        else
         {
          double area_backup = Problem_Parameter::Downstream_area_pt->value(0);
          double displacement_backup = 
           Problem_Parameter::Prescribed_displacement;
          double flow_rate_backup = 
           Problem_Parameter::Prescribed_flow_rate;
          problem.adapt();
          problem.steady_newton_solve(0);
          Problem_Parameter::Downstream_area_pt->set_value(0,area_backup);
          Problem_Parameter::Prescribed_displacement = displacement_backup;
          Problem_Parameter::Prescribed_flow_rate = flow_rate_backup;
          problem.steady_newton_solve(0);
          count=0;
         }
        Problem_Parameter::Prescribed_displacement_prev = 
         Problem_Parameter::Prescribed_displacement;
       }
      //exit(1);
     }
     
     Problem_Parameter::Displacement_continuation = false;
     
     Problem_Parameter::Prescribed_area = 
      Problem_Parameter::Target_prescribed_area;
     //Problem_Parameter::Downstream_area_pt->value(0);
     
     double area_backup = Problem_Parameter::Prescribed_area;
     double flow_rate_backup = Problem_Parameter::Prescribed_flow_rate;
     Problem_Parameter::Prescribe_area = true;
     problem.adapt();
     problem.steady_newton_solve(0);
     Problem_Parameter::Prescribed_area = area_backup;
     Problem_Parameter::Prescribed_flow_rate = flow_rate_backup;
     
     //problem.complete_problem_setup();
     //oomph_info<<"Number of equations: "<<problem.assign_eqn_numbers()<<std::endl;
     problem.steady_newton_solve(max_adapt); //1
     problem.doc_solution();
     
     //========================================================================
     // {
     //  Problem_Parameter::Target_prescribed_displacement=1.0;
     //  Problem_Parameter::Max_increment = 2.5e-3;
     //  int count=0;
     //  while(Problem_Parameter::Prescribed_displacement < 
     //        Problem_Parameter::Target_prescribed_displacement)
     //  //while(Problem_Parameter::Displacement_continuation)
     //   {
     //    if( (count!=Problem_Parameter::Adaptation_frequency) ||
     //        (Problem_Parameter::Adaptation_frequency < 0) )
     //     {
     //      Problem_Parameter::Prescribed_displacement+=
     //       Problem_Parameter::Max_increment;
     //      Problem_Parameter::Displacement_control_node_pt->set_value(
     //       0, Problem_Parameter::Prescribed_displacement);
     //      problem.steady_newton_solve(0);
     //      problem.doc_solution();
     //      count++;
     //     }
     //    else
     //     {
     //      double flux_backup = Problem_Parameter::Prescribed_flux;
     //      problem.adapt();
     //      problem.steady_newton_solve(0);
     //      Problem_Parameter::Prescribed_flux = flux_backup;
     //      problem.steady_newton_solve(0);
     //      count=0;
     //     }
     //   }
     //  Problem_Parameter::Prescribed_displacement = 
     //   Problem_Parameter::Target_prescribed_displacement;
     //  Problem_Parameter::Displacement_control_node_pt->set_value(
     //   0,Problem_Parameter::Prescribed_displacement);
     //  problem.steady_newton_solve(max_adapt);
     //  problem.doc_solution();
     //  exit(1);
     // }
     
     //========================================================================
     problem.enable_parameter_update();
     Problem_Parameter::Update_ca_only=false;
     {
      Problem_Parameter::Target_prescribed_flow_rate = 1.0*1.0e-6/60.0;
      Problem_Parameter::Flow_rate_continuation = true;
      Problem_Parameter::Prescribed_flow_rate_prev = 
       Problem_Parameter::Prescribed_flow_rate;
      Problem_Parameter::Max_increment = -10.0*1.0e-6/60.0;
      double desired_increment=-1.0*1.0e-6/60.0; // minus: decrease Ca, plus: increase
//      if (CommandLineArgs::command_line_flag_has_been_set("--inc_thin_films"))            
//       {
//        desired_increment *= -1.0;
//       }
      double actual_increment=desired_increment;
      int count=0;
      while(Problem_Parameter::Prescribed_flow_rate > 
            Problem_Parameter::Target_prescribed_flow_rate)
       {
        if( (count!=Problem_Parameter::Adaptation_frequency) ||
            (Problem_Parameter::Adaptation_frequency < 0) )
         {
//          actual_increment = problem.arc_length_step_solve(
//           &Problem_Parameter::Prescribed_flow_rate,desired_increment,0);
          Problem_Parameter::Prescribed_flow_rate += desired_increment;
          problem.steady_newton_solve(0);
          problem.doc_solution();
//          desired_increment=actual_increment;
          count++;
         }
        else
         {
          double area_backup = Problem_Parameter::Prescribed_area;
          double flow_rate_backup = Problem_Parameter::Prescribed_flow_rate;
          problem.adapt();
          problem.steady_newton_solve(0);
          Problem_Parameter::Prescribed_area = area_backup;
          Problem_Parameter::Prescribed_flow_rate = flow_rate_backup;
          problem.steady_newton_solve(0);
          count=0;
         }
        Problem_Parameter::Prescribed_flow_rate_prev = 
         Problem_Parameter::Prescribed_flow_rate;
       }
     }
     Problem_Parameter::Prescribed_flow_rate = 
      Problem_Parameter::Target_prescribed_flow_rate;
     Problem_Parameter::Flow_rate_continuation = false;
     problem.steady_newton_solve(max_adapt);
     problem.doc_solution();
     problem.disable_parameter_update();
     
    }
   else
    {
     if (!CommandLineArgs::command_line_flag_has_been_set("--inc_thin_films"))
       {
        Problem_Parameter::Thin_film_homotopy_parameter = 0.0;
        }
     problem.make_unsteady();
     problem.enable_parameter_update();
     Problem_Parameter::Update_ca_only=false;

     Problem_Parameter::Impose_U_bubble_tip=false;
     double flow_rate_desired_backup = 144.0072*1.0e-6/60.0;
     double U_bubble_backup = 0.1;
     double area_backup = 1.0;
     if(Problem_Parameter::Impose_U_bubble_tip)
      {
       Problem_Parameter::U_bubble_tip_dimensional_desired=U_bubble_backup;
       }
     else
      {
       Problem_Parameter::Flow_rate_dimensional_desired = flow_rate_desired_backup;
       }

     Problem_Parameter::Prescribed_area = area_backup; 
     Problem_Parameter::Prescribe_area = true;
     problem.adapt();
     bool first = true;
     for(unsigned i = 0; i < nstep; i++)
      {   
       // Allow adaptation every few steps
       if( (i>0) && (Problem_Parameter::Adaptation_frequency>0) &&
           (i%Problem_Parameter::Adaptation_frequency == 0) )
        {
         problem.adapt();//max_adapt = 1;
         problem.unsteady_newton_solve(1.0e-8,0,first);
         first = false;
        }

       double next_dt = dt;
       bool perturb_flow_rate=false;
       if (CommandLineArgs::command_line_flag_has_been_set
           ("--suppress_temporal_adaptivity"))
        {
         problem.unsteady_newton_solve(dt,0,first);
        }
       else
        {

//            if(perturb_flow_rate)
//             {
//              double flow_rate_perturb_parameter = 8.0*double(problem.time_stepper_pt()->time());
//              Problem_Parameter::Prescribed_flow_rate = 
//                     flow_rate_backup*(1.0+flow_rate_perturb_parameter*pow(2.718281828,-flow_rate_perturb_parameter));
//              }

         next_dt=
          //problem.doubly_adaptive_unsteady_newton_solve(dt,epsilon_t,
          //                                              max_adapt,first);
          problem.adaptive_unsteady_newton_solve(dt,epsilon_t,first);
         oomph_info << "Suggested next dt: " << next_dt << std::endl;
        }
       
       dt = next_dt; 
       if(dt>0.1)
       {
        dt=0.1;
       }
       
       problem.doc_solution();

      }
    }
   oomph_info << "done\n";
  }
#ifdef OOMPH_HAS_MPI
 MPI_Helpers::finalize();
#endif

} //End of main


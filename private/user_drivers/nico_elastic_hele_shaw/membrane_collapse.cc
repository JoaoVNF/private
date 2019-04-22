//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC//    Version 1.0; svn revision $LastChangedRevision: 1097 $
//LIC//
//LIC// $LastChangedDate: 2015-12-17 11:53:17 +0000 (Thu, 17 Dec 2015) $
//LIC// 
//LIC// Copyright (C) 2006-2016 Matthias Heil and Andrew Hazel
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

// The equations
//#include "foeppl_von_karman.h"
#include "nico_elastic_hele_shaw.h"

// The mesh
#include "meshes/triangle_mesh.h"

// obacht are we using stress or displacement formulation
#define STRESS_FORMULATION

using namespace std;
using namespace oomph;
using MathematicalConstants::Pi;

namespace Problem_Parameter
{
 /// Width of the channel (m)
 double W_channel=30.0*1.0e-3;

 /// Length of the channel (m)
 double L_channel=210.0*1.0e-3; //180.0*1.0e-3; //90.0*1.0e-3;

 /// Depth of the channel (m)
 double H_channel=1.05e-3;

 /// Gap width when collapsed (m)
 double Gap_size=50.0*1.0e-6; //1.0e-4;

 /// Thickness of the sheet (m)
 double d=0.34*1.0e-3;

 /// Young's modulus of the sheet (Pa)
 double E=1.44*1.0e6;

 /// Poisson ratio of the sheet
 double nu=0.5;

 /// Pre-factor in FvK equations
 double eta=12.0*(1.0-nu*nu)*(W_channel/d)*(W_channel/d);

 /// Material constant
 double K=E*d*d*d/(12.0*(1.0-nu*nu));

 /// Pressure on the membrane (Pa), relative to atmospheric pressure 
 /// negative because we want to push the membrane down
 double P= 0.0; //-50.0;

 // Storage for the dimensional pressure on the membrane
 Data *P_dim_pt = 0;
  
 /// Heavy-side function activating the impact function
 double Heavy(const double& z)
 {
  if(z>=-H_channel/W_channel+Gap_size/W_channel) return 0.0;
  else return 1.0;
 }
  
 /// Spring stiffness
 double C=100000.0;

 // Assigns the value of pressure (spatially constant)
 void get_pressure(const Vector<double>& x, const double& w, double& pressure)
  {
   pressure=P_dim_pt->value(0)*W_channel*W_channel*W_channel/K
    +Heavy(w)*pow(std::fabs(w+H_channel/W_channel-Gap_size/W_channel),1.0/2.0)
    *C;
  }

 void get_airy_forcing(const Vector<double>& x, double& airy_forcing)
  {
   airy_forcing = 0.0;
  }

 // Pre stress in the sheet (Pa)
 double Pre_stress = 30.0e3; //40.0e3;

 void get_pre_stress(DenseMatrix<double>& sigma_0)
 {
  sigma_0(0,0) = 0.0;
  sigma_0(1,1) = Pre_stress/E;
  sigma_0(0,1) = 0.0;
  sigma_0(1,0) = sigma_0(0,1);
 }

 // Maximum area in the mesh
 double Uniform_element_area=0.005;

}


///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////


//==start_of_problem_class============================================
/// Class definition
//====================================================================
template<class ELEMENT>
class MembraneCollapseFvKProblem : public virtual Problem
{

public:

 /// Constructor
 MembraneCollapseFvKProblem();
    
 /// Destructor
 ~MembraneCollapseFvKProblem()
  {
   delete Sheet_mesh_pt;
  };

 /// Actions before adapt. Delete old volume constraint element and update
 /// global mesh
 void actions_before_adapt()
  {
   //rebuild_global_mesh();
  }
 
 /// \short Actions after adapt: 
 /// Setup the problem again -- remember that the mesh has been
 /// completely rebuilt and its element's don't have any
 /// pointers to source fcts etc. yet
 /// Also create new volume constraint element with previous pressure value
 void actions_after_adapt()
  {
   //rebuild_global_mesh();
   complete_problem_setup();
  }
 
 /// Update after solve (empty)
 void actions_after_newton_solve()
  {
  }

 /// Update the problem specs before solve: Re-apply boundary conditons
 void actions_before_newton_solve()
  {
   //apply_boundary_conditions();
  }
  
 /// Doc the solution
 void doc_solution();
 

 // Get the displacement in the centre of the domain
 double displacement_at_domain_centre()
  {
   return Centre_nod_pt->value(0);
  }

private:

 /// Doc info object for labeling output
 DocInfo Doc_info;

 /// \short Helper function to (re-)set boundary condition
 /// and complete the build of  all elements
 void complete_problem_setup();

 /// Pointers to specific mesh
 RefineableTriangleMesh<ELEMENT>* Sheet_mesh_pt;

 /// Trace file to document norm of solution
 ofstream Trace_file;

 /// Pointer to the node in the centre of the domain
 Node* Centre_nod_pt;

 // Keep track of boundary ids
 enum
  {
   Outer_boundary0 = 0,
   Outer_boundary1 = 1,
   Outer_boundary2 = 2,
   Outer_boundary3 = 3,
   Inner_boundary0 = 4,
   Inner_boundary1 = 5
  };

}; // end_of_problem_class


template<class ELEMENT>
MembraneCollapseFvKProblem<ELEMENT>::MembraneCollapseFvKProblem()
{
 // We'll do pseudo-timestepping to adjust the spring stiffness
 // Allocate the timestepper -- this constructs the Problem's 
 // time object with a sufficient amount of storage to store the
 // previous timsteps. 
 //this->add_time_stepper_pt(new BDF<2>(true));

 Max_newton_iterations = 20;
 Always_take_one_newton_step = true;

 Vector<Vector<double> > vertex_coord(3);
 for(unsigned i=0;i<3;i++)
  {
   vertex_coord[i].resize(2);
  }

 //Outer boundary
 //--------------

 TriangleMeshClosedCurve* outer_boundary_pt = 0;

 Vector<TriangleMeshCurveSection*> outer_boundary_polyline_pt(4);
 
 vertex_coord[0][0]=0.0;
 vertex_coord[0][1]=0.5;
 vertex_coord[1][0]=0.0;
 vertex_coord[1][1]=0.0;
 vertex_coord[2][0]=0.0;
 vertex_coord[2][1]=-0.5;

 outer_boundary_polyline_pt[0] = new TriangleMeshPolyLine(vertex_coord,
                                                          Outer_boundary0);

 vertex_coord[0][0]=0.0;
 vertex_coord[0][1]=-0.5;
 vertex_coord[1][0]=Problem_Parameter::L_channel/
  (2.0*Problem_Parameter::W_channel);
 vertex_coord[1][1]=-0.5;
 vertex_coord[2][0]=Problem_Parameter::L_channel/
  Problem_Parameter::W_channel;
 vertex_coord[2][1]=-0.5;

 outer_boundary_polyline_pt[1] = new TriangleMeshPolyLine(vertex_coord,
                                                          Outer_boundary1);

 vertex_coord[0][0]=Problem_Parameter::L_channel/
  Problem_Parameter::W_channel;
 vertex_coord[0][1]=-0.5;
 vertex_coord[1][0]=Problem_Parameter::L_channel/
  Problem_Parameter::W_channel;
 vertex_coord[1][1]=0.0;
 vertex_coord[2][0]=Problem_Parameter::L_channel/
  Problem_Parameter::W_channel;
 vertex_coord[2][1]=0.5;

 outer_boundary_polyline_pt[2] = new TriangleMeshPolyLine(vertex_coord,
                                                          Outer_boundary2);

 vertex_coord[0][0]=Problem_Parameter::L_channel/
  Problem_Parameter::W_channel;
 vertex_coord[0][1]=0.5;
 vertex_coord[1][0]=Problem_Parameter::L_channel/
  (2.0*Problem_Parameter::W_channel);
 vertex_coord[1][1]=0.5;
 vertex_coord[2][0]=0.0;
 vertex_coord[2][1]=0.5;

 outer_boundary_polyline_pt[3] = new TriangleMeshPolyLine(vertex_coord,
                                                          Outer_boundary3);

 outer_boundary_pt =
  new TriangleMeshClosedCurve(outer_boundary_polyline_pt);


 //Inner boundaries
 //----------------

 //Cross-sectional boundary lines
 Vector<TriangleMeshOpenCurve*> inner_open_boundary_pt(2);

 vertex_coord.resize(2);
 for(unsigned i=0;i<2;i++)
  {
   vertex_coord[i].resize(2);
  }

 //Line along the channel
 vertex_coord[0][0] = 0.0;
 vertex_coord[0][1] = 0.0;
 vertex_coord[1][0] = Problem_Parameter::L_channel/
  Problem_Parameter::W_channel;
 vertex_coord[1][1] = 0.0;

 TriangleMeshPolyLine *inner_open_polyline1_pt =
  new TriangleMeshPolyLine(vertex_coord, Inner_boundary0);

 inner_open_polyline1_pt->connect_initial_vertex_to_polyline(
   dynamic_cast<TriangleMeshPolyLine *>(outer_boundary_polyline_pt[0]),1);

 inner_open_polyline1_pt->connect_final_vertex_to_polyline(
   dynamic_cast<TriangleMeshPolyLine *>(outer_boundary_polyline_pt[2]),1);

 Vector<TriangleMeshCurveSection *> inner_boundary_line1_pt(1);
 inner_boundary_line1_pt[0] = inner_open_polyline1_pt;

 inner_open_boundary_pt[0] =
  new TriangleMeshOpenCurve(inner_boundary_line1_pt);

 //Line across the channel
 vertex_coord[0][0] = Problem_Parameter::L_channel/
  (2.0*Problem_Parameter::W_channel);
 vertex_coord[0][1] = -0.5;
 vertex_coord[1][0] = Problem_Parameter::L_channel/
  (2.0*Problem_Parameter::W_channel);
 vertex_coord[1][1] = 0.5;

 TriangleMeshPolyLine *inner_open_polyline2_pt =
  new TriangleMeshPolyLine(vertex_coord, Inner_boundary1);

 inner_open_polyline2_pt->connect_initial_vertex_to_polyline(
   dynamic_cast<TriangleMeshPolyLine *>(outer_boundary_polyline_pt[1]),1);

 inner_open_polyline2_pt->connect_final_vertex_to_polyline(
   dynamic_cast<TriangleMeshPolyLine *>(outer_boundary_polyline_pt[3]),1);

 Vector<TriangleMeshCurveSection *> inner_boundary_line2_pt(1);
 inner_boundary_line2_pt[0] = inner_open_polyline2_pt;

 inner_open_boundary_pt[1] =
  new TriangleMeshOpenCurve(inner_boundary_line2_pt);

 //Create the mesh
 //---------------
 //Create mesh parameters object
 TriangleMeshParameters mesh_parameters(outer_boundary_pt);

 mesh_parameters.internal_open_curves_pt() = inner_open_boundary_pt;
 mesh_parameters.element_area() = Problem_Parameter::Uniform_element_area;

 Sheet_mesh_pt = new RefineableTriangleMesh<ELEMENT>(mesh_parameters);

 Z2ErrorEstimator* error_estimator_pt = new Z2ErrorEstimator;
 Sheet_mesh_pt->spatial_error_estimator_pt() = error_estimator_pt;

 Sheet_mesh_pt->max_element_size() = Problem_Parameter::Uniform_element_area;
 Sheet_mesh_pt->min_element_size() = 0.002;

 //Sheet_mesh_pt->max_permitted_error()=0.00005;
 //Sheet_mesh_pt->min_permitted_error()=0.00001;

 Problem_Parameter::P_dim_pt = new Data(1);
 Problem_Parameter::P_dim_pt->set_value(0,Problem_Parameter::P);
// Problem_Parameter::P_dim_pt->pin(0); ////Joao: before --->>> Problem_Parameter::P_dim_pt->pin(0);

 complete_problem_setup();

 add_sub_mesh(Sheet_mesh_pt);

 build_global_mesh();

 char filename[100];
 sprintf(filename, "RESLT/trace.dat");
 Trace_file.open(filename);

 oomph_info << "Number of equations: "
            << this->assign_eqn_numbers() << '\n';
}



//==start_of_complete======================================================
 /// Set boundary condition exactly, and complete the build of 
 /// all elements
//========================================================================
template<class ELEMENT>
void MembraneCollapseFvKProblem<ELEMENT>::complete_problem_setup()
{   

 // Set the boundary conditions for problem: All nodes are
 // free by default -- just pin the ones that have Dirichlet conditions
 // here. 
 //unsigned nbound=Sheet_mesh_pt->nboundary();

 for(unsigned ibound=0;ibound<Inner_boundary0;ibound++)
  {
   unsigned num_nod=Sheet_mesh_pt->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
    {
     // Get node
     Node* nod_pt=Sheet_mesh_pt->boundary_node_pt(ibound,inod);
     
#ifndef STRESS_FORMULATION

     // pin displacement
     nod_pt->pin(0);

     // pin tangential displacement on remaining boundaries
     if( (ibound == Outer_boundary0) || (ibound == Outer_boundary2) )
      {
       nod_pt->pin(3);
      }

     if( (ibound == Outer_boundary1) || (ibound == Outer_boundary3) )
      {
       nod_pt->pin(2);
      }
#else
     // Pin unknown values (displacement w and Airy stress function \phi)
     // This BC corresponds to stress-free clamping
     nod_pt->pin(0);
     nod_pt->pin(2);
     nod_pt->pin(3); // need this one as well
#endif

    }   
  } // end loop over boundaries
 

 /// Create higher order Gauss schemes
 TGauss<2,5>* int1_pt = new TGauss<2,5>;
 TGauss<2,13>* int2_pt = new TGauss<2,13>;
 
 // Complete the build of all elements so they are fully functional
 unsigned n_element = Sheet_mesh_pt->nelement();
 for(unsigned e=0;e<n_element;e++)
  {
   // Upcast from GeneralisedElement to the present element
   ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Sheet_mesh_pt->element_pt(e));
   
   // obacht set integration scheme
   //el_pt->set_integration_scheme(int1_pt);
   el_pt->set_integration_scheme(int2_pt);

   //Set the pressure function pointers and the physical constants
   el_pt->eta_pt() = &Problem_Parameter::eta;
   el_pt->pressure_fct_pt() = &Problem_Parameter::get_pressure;  
   el_pt->pre_stress_fct_pt() = &Problem_Parameter::get_pre_stress;
#ifdef STRESS_FORMULATION
   el_pt->airy_forcing_fct_pt() = &Problem_Parameter::get_airy_forcing;
#else 
   el_pt->nu_pt() = &Problem_Parameter::nu;
#endif

   // Linear wall?
   if (CommandLineArgs::command_line_flag_has_been_set("--linear_wall"))
    {
     el_pt->use_linear_bending_model();
    }
  }

 // Updating the centre node pointer
 // Loop over the nodes on the inner boundaries
 Centre_nod_pt = 0;
 unsigned n_node=this->Sheet_mesh_pt->nboundary_node(Inner_boundary0);
 for(unsigned i=0;i<n_node;i++)
  {
   // Get node
   Node* nod_pt=this->Sheet_mesh_pt->boundary_node_pt(Inner_boundary0,i);

   if(nod_pt->is_on_boundary(Inner_boundary1))
    {
     if(Centre_nod_pt!=0)
      {
       oomph_info<<"Strange.. More than one node in the centre\n";
       abort();
      }
     else
      {
       Centre_nod_pt=nod_pt;
       oomph_info<<"Centre node now at "<<nod_pt->x(0)<<", "
                 <<nod_pt->x(1)<<"\n";
      }
    }
  }

}

//==start_of_doc_solution=================================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void MembraneCollapseFvKProblem<ELEMENT>::doc_solution()
{ 
 oomph_info << "Docing step: "<<Doc_info.number()<<"\n";

 ofstream some_file;
 char filename[100];
 
 // Number of plot points
 unsigned npts = 5;
 
 sprintf(filename,"RESLT/soln%i.dat",Doc_info.number());
 some_file.open(filename);
 this->Sheet_mesh_pt->output(some_file,npts);
 some_file.close();
 
 // Output boundaries
 //------------------
 // sprintf(filename,"RESLT/boundaries%i.dat",Doc_info.number());
 // some_file.open(filename);
 // Sheet_mesh_pt->output_boundaries(some_file);
 // some_file.close();


 // Output along channel
 // --------------------
 sprintf(filename,"RESLT/soln_along%i.dat",Doc_info.number());
 some_file.open(filename);
 unsigned n_boundary_nodes = Sheet_mesh_pt->nboundary_node(Inner_boundary0);
 for(unsigned inod = 0; inod < n_boundary_nodes; inod++)
  {
   some_file << Sheet_mesh_pt->boundary_node_pt(Inner_boundary0,inod)->x(0)
             << " " 
             << Sheet_mesh_pt->boundary_node_pt(Inner_boundary0,inod)->value(0)
             << std::endl;
  }
 some_file.close();


 // Output across channel
 // --------------------
 sprintf(filename,"RESLT/soln_across%i.dat",Doc_info.number());
 some_file.open(filename);
 n_boundary_nodes = Sheet_mesh_pt->nboundary_node(Inner_boundary1);
 for(unsigned inod = 0; inod < n_boundary_nodes; inod++)
  {
   some_file << Sheet_mesh_pt->boundary_node_pt(Inner_boundary1,inod)->x(1)
             << " " 
             << Sheet_mesh_pt->boundary_node_pt(Inner_boundary1,inod)->value(0)
             << std::endl;
  }
 some_file.close();


 // Doc error and return of the square of the L2 error
 //---------------------------------------------------
 //double error,norm,dummy_error,zero_norm;
 // double dummy_error,zero_norm;
 // sprintf(filename,"RESLT/error%i.dat",Doc_info.number());
 // some_file.open(filename);
 
 // Sheet_mesh_pt->compute_error(some_file,TestSoln::zero,
 //                              dummy_error,zero_norm);
 // some_file.close();

 // // Doc L2 error and norm of solution
 // oomph_info << "Norm of computed solution: " << sqrt(dummy_error) << std::endl;

 Trace_file << Problem_Parameter::P_dim_pt->value(0) << " "
            << Problem_Parameter::C << " "
            << displacement_at_domain_centre() << std::endl;

 //Trace_file.flush();

 // Increment the doc_info number
 Doc_info.number()++;

} // end of doc


//=======start_of_main========================================
///Driver code for demo of inline triangle mesh generation
//============================================================
int main(int argc, char **argv)
{

 feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW);


 // Store command line arguments
 CommandLineArgs::setup(argc,argv);

 // Define possible command line arguments and parse the ones that
 // were actually specified

 // Validation?
 CommandLineArgs::specify_command_line_flag("--validation");

 // Linear wall?
 CommandLineArgs::specify_command_line_flag("--linear_wall");

 // What sort of pre-stress are we using?
 CommandLineArgs::specify_command_line_flag("--pre_stress", 
                                            &Problem_Parameter::Pre_stress);

 // The initial pressure value
 CommandLineArgs::specify_command_line_flag("--pressure", 
                                            &Problem_Parameter::P);

 // The element area
 CommandLineArgs::specify_command_line_flag("--uniform_element_area", 
                                            &Problem_Parameter::Uniform_element_area);

 // Tolerance for opposite wall contact
 double impact_tol = 0.01;
 CommandLineArgs::specify_command_line_flag("--impact_tol", &impact_tol);

 // Parse command line
 CommandLineArgs::parse_and_assign(); 

 // Doc what has actually been specified on the command line
 CommandLineArgs::doc_specified_flags();

 // Problem instance
#ifdef STRESS_FORMULATION
 MembraneCollapseFvKProblem<
 ProjectableFoepplvonKarmanElement<TFoepplvonKarmanElement<3> > >
  problem;
#else
 MembraneCollapseFvKProblem<
 ProjectableFoepplvonKarmanDisplacementElement<
 TFoepplvonKarmanDisplacementElement<3> > >
  problem;
#endif



/////////////
bool vary_spring_stifness = false;
bool vary_pressure = true;
/////////////

 if(vary_pressure)
 {
   for(unsigned i=0;i<200;i++)
    {
     Problem_Parameter::C=1000.0;
     Problem_Parameter::P = 0.0;
     Problem_Parameter::P_dim_pt->set_value(0,Problem_Parameter::P);
     double increment=-0.5;
     //problem.newton_solve();
     increment = problem.arc_length_step_solve(&Problem_Parameter::P,increment,0);
     problem.doc_solution();
     Problem_Parameter::P+=increment;
    }
   exit(1);
  }


 //problem.doc_solution();
 if(vary_spring_stifness)
 {
   Problem_Parameter::P = -2.0;
   Problem_Parameter::C=0.0;
   double increment=1000.0;
   for(unsigned i=0;i<200;i++)
    {
     //problem.newton_solve();
     increment = problem.arc_length_step_solve(&Problem_Parameter::C,increment,0);
     problem.doc_solution();
     Problem_Parameter::C+=increment;

     // Calculate the deviation from the minimum gap width
     double deviation = problem.displacement_at_domain_centre()+
      Problem_Parameter::H_channel/Problem_Parameter::W_channel-
      Problem_Parameter::Gap_size/Problem_Parameter::W_channel;

     double relative_deviation = std::fabs(deviation)/
      (Problem_Parameter::Gap_size/Problem_Parameter::W_channel);

     if( (deviation > 0.0) || (relative_deviation <= impact_tol) ) break;

    }
   exit(1);
  }
 

} //End of main


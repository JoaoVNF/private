
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
//Non-inline functions for HeleShaw elements
#include "hele_shaw_elements.h"


namespace oomph
{
 
 
//======================================================================
/// Set the data for the number of Variables at each node, always one
/// in every case
//======================================================================
 template<unsigned NNODE_1D>
  const unsigned QHeleShawElement<NNODE_1D>::Initial_Nvalue = 1;
 
 
//======================================================================
/// Compute element residual Vector and/or element Jacobian matrix
///
/// flag=1: compute both
/// flag=0: compute only residual Vector
///
/// Pure version without hanging nodes
//======================================================================
 void  HeleShawEquations::
  fill_in_generic_residual_contribution_hele_shaw(Vector<double> &residuals,
                                                  DenseMatrix<double> &jacobian,
                                                  const unsigned& flag)
 {

  // Do nothing (I.e. make no contribution to residuals or
  // Jacobian) if we're not solving the Hele Shaw equations
  if (Hele_shaw_disabled) 
   {
    return; 
   }

  // Find out how many nodes there are
  const unsigned n_node = nnode();
  
  // Set up memory for the shape and test functions
  Shape psi(n_node), test(n_node);
  DShape dpsidx(n_node,2), dtestdx(n_node,2);
  
  // Index at which the hele_shaw unknown is stored
  const unsigned p_nodal_index = p_index_hele_shaw();
  
  // Set the value of n_intpt
  const unsigned n_intpt = integral_pt()->nweight();
  
  // Integers to store the local equation and unknown numbers
  int local_eqn = 0, local_unknown = 0;
  
  // Loop over the integration points
  for(unsigned ipt=0; ipt<n_intpt; ipt++)
   {
    // Get the integral weight
    double w = integral_pt()->weight(ipt);

    // Get local coordinate
    Vector<double> s(2);
    for(unsigned i=0;i<2;i++)
     {
      s[i] = this->integral_pt()->knot(ipt,i);
     }
    
    // Call the derivatives of the shape and test functions
    double J = dshape_and_dtest_eulerian_at_knot_hele_shaw(ipt,psi,dpsidx,
                                                           test,dtestdx);
    
    // Premultiply the weights and the Jacobian
    double W = w*J;
    
    // Calculate local values of unknown
    // Allocate and initialise to zero
    double interpolated_p = 0.0;
    Vector<double> interpolated_x(2,0.0);
    Vector<double> interpolated_dpdx(2,0.0);

    // // hierher start shite
    // {
    //  oomph_info << "hierher shite 2" << std::endl;
    //  hierher_shite();
    //  oomph_info << "hierher back from shite 2" << std::endl;
     
    //  // Loop over nodes
    //  for(unsigned l=0; l<n_node; l++)
    //   {
    //    // // Get the nodal value of the hele_shaw unknown
    //    // double p_value = raw_nodal_value(l,p_nodal_index);
    //    // interpolated_p += p_value*psi(l);
    //    // Loop over directions
    //    for(unsigned j=0; j<2; j++)
    //     {
    //      oomph_info << "inline: l,j " << l << " " << j << std::endl;
    //      double tmp=raw_nodal_position(l,j);
    //      oomph_info << "raw_nodal_position(l,j)" << tmp << std::endl;
    //     }
    //   }
    // }
    // // hierher end shite

    // Calculate function value and derivatives:
    //-----------------------------------------
    // Loop over nodes
    for(unsigned l=0; l<n_node; l++)
     {
      // Get the nodal value of the hele_shaw unknown
      double p_value = raw_nodal_value(l,p_nodal_index);
      interpolated_p += p_value*psi(l);
      // Loop over directions
      for(unsigned j=0; j<2; j++)
       {
        interpolated_x[j] += raw_nodal_position(l,j)*psi(l);
        interpolated_dpdx[j] += p_value*dpsidx(l,j);
       }
     }
    
    
    // Get gap width and (Eulerian!) wall velocity
    // (Pass in psi and dpsidx to avoid their recomputation
    // in fsi context.
    double h = 1.0;
    double dhdt = 0.0;
    get_upper_wall_data(s,interpolated_x,psi,dpsidx,h,dhdt);
    
    
    // Assemble residuals and Jacobian
    //--------------------------------
    
    
    // Loop over the test functions
    for(unsigned l=0; l<n_node; l++)
     {
      // Get the local equation
      local_eqn = nodal_local_eqn(l,p_nodal_index);
      
      // If it's not a boundary condition
      if(local_eqn >= 0)
       {
        // Wall velocity (RHS). Note that dhdt will already
        // (have to!) include the ALE term in a multi-physics
        // context.
        residuals[local_eqn] += dhdt*test(l)*W;
        
        
        // The HeleShaw bit itself
        for(unsigned k=0; k<2; k++)
         {
          residuals[local_eqn] += pow(h,3)*interpolated_dpdx[k]*dtestdx(l,k)*W;
         }
        
        // Calculate the jacobian
        //-----------------------
        if(flag)
         {
          // Loop over the velocity shape functions again
          for(unsigned l2=0; l2<n_node; l2++)
           {
            local_unknown = nodal_local_eqn(l2,p_nodal_index);
            // If at a non-zero degree of freedom add in the entry
            if(local_unknown >= 0)
             {
              // Add contribution to Elemental Matrix
              for(unsigned i=0; i<2; i++)
               {
                jacobian(local_eqn,local_unknown)
                 += pow(h,3)*dpsidx(l2,i)*dtestdx(l,i)*W;
               }
             }
           }
         }
       }
     }    
   } // End of loop over integration points
 }
 
 
 
 
 
//======================================================================
/// Self-test:  Return 0 for OK
//======================================================================
 unsigned HeleShawEquations::self_test()
 {
  
  bool passed = true;
  
// Check lower-level stuff
  if (FiniteElement::self_test()!=0)
   {
    passed = false;
   }
  
// hierher: fill in missing self-tests
  
// Return verdict
  if (passed)
   {
    return 0;
   }
  else
   {
    return 1;
   }  
 }
 
 
 
//======================================================================
/// Output function:
///
///   x,y,u   or    x,y,z,u
///
/// nplot points in each coordinate direction
//======================================================================
 void  HeleShawEquations::output(std::ostream &outfile,
                                 const unsigned &nplot)
 {
  
  // Vector of local coordinates and velocity
  Vector<double> s(2);
  Vector<double> velocity(2);
  
  // Tecplot header info
  outfile << tecplot_zone_string(nplot);
  
  // Loop over plot points
  unsigned num_plot_points = nplot_points(nplot);
  for (unsigned iplot=0; iplot<num_plot_points; iplot++)
   {    
    // Get local coordinates and velocity at plot point
    get_s_plot(iplot,nplot,s);
    get_velocity(s,velocity);
    
    for(unsigned i=0; i<2; i++)
     {
      outfile << interpolated_x(s,i) << " ";
     }
    outfile << velocity[0] << " "
            << velocity[1] << " "
            << interpolated_p_hele_shaw(s) << "\n";
   }
  
  // Write tecplot footer (e.g. FE connectivity lists)
  write_tecplot_zone_footer(outfile,nplot);  
 }
 
 
//======================================================================
/// C-style output function:
///
///   x,y,u   or    x,y,z,u
///
/// nplot points in each coordinate direction
//======================================================================
 void  HeleShawEquations::output(FILE* file_pt,
                                 const unsigned &nplot)
 {
  
  // hierher make consistent with c++ output
  
  // Vector of local coordinates
  Vector<double> s(2);
  Vector<double> velocity(2);
  
  // Tecplot header info
  fprintf(file_pt,"%s",tecplot_zone_string(nplot).c_str());
  
  // Loop over plot points
  unsigned num_plot_points = nplot_points(nplot);
  for (unsigned iplot=0; iplot<num_plot_points; iplot++)
   {
    // Get local coordinates of plot point
    get_s_plot(iplot,nplot,s);
    get_velocity(s,velocity);
    
    for(unsigned i=0; i<2; i++)
     {
      fprintf(file_pt,"%g ",interpolated_x(s,i));
     }
    fprintf(file_pt,"%g \n",velocity[0] );
    fprintf(file_pt,"%g \n",velocity[1] );
    fprintf(file_pt,"%g \n",interpolated_p_hele_shaw(s));
   }
  
  // Write tecplot footer (e.g. FE connectivity lists)
  write_tecplot_zone_footer(file_pt,nplot);
 }
 
 
 
//======================================================================
/// Output exact solution
///
/// Solution is provided via function pointer.
/// Plot at a given number of plot points.
///
///   x,y,u_exact    or    x,y,z,u_exact
//======================================================================
 void HeleShawEquations::output_fct(std::ostream &outfile,
                                    const unsigned &nplot,
                                    FiniteElement::SteadyExactSolutionFctPt exact_soln_pt)
 {
  
  
  // Vector of local coordinates
  Vector<double> s(2);
  
  // hierher make consistent with c++ output
    
  // Vector for coordintes
  Vector<double> x(2);

  // Tecplot header info
  outfile << tecplot_zone_string(nplot);
  
  // Exact solution vector: u,v,p
  Vector<double> exact_soln(3);
  
  // Loop over plot points
  unsigned num_plot_points = nplot_points(nplot);
  for (unsigned iplot=0; iplot<num_plot_points; iplot++)
   {    
    // Get local coordinates of plot point
    get_s_plot(iplot,nplot,s);

    // Get x position as Vector
    interpolated_x(s,x);
    
    // Get exact solution at this point
    (*exact_soln_pt)(x,exact_soln);
    
    // Output x,y,...,u_exact
    for(unsigned i=0; i<2; i++)
     {
      outfile << x[i] << " ";
     }
    
    // Output "exact solution"
    for(unsigned i=0; i<3; i++)
     {
      outfile << exact_soln[i] << " ";
     }
    
    outfile << std::endl;
    
   }
  
  // Write tecplot footer (e.g. FE connectivity lists)
  write_tecplot_zone_footer(outfile,nplot);
 }
 
 
//======================================================================
/// Validate against exact solution
///
/// Solution is provided via function pointer.
/// Plot error at a given number of plot points.
///
//======================================================================
 void HeleShawEquations::compute_error(std::ostream &outfile,
                                       FiniteElement::SteadyExactSolutionFctPt 
                                       exact_soln_pt,
                                       double& error, double& norm)
 {
  
  // hierher make consistent with c++ output
  
  // Initialise
  norm = 0.0;
  error = 0.0;
  
  // Vector of local coordinates
  Vector<double> s(2);
  Vector<double> velocity(2);
  
  // Vector for coordintes
  Vector<double> x(2);
  
  // Find out how many nodes there are in the element
  unsigned n_node = nnode();
  
  Shape psi(n_node);
  
  // Set the value of n_intpt
  unsigned n_intpt = integral_pt()->nweight();
  
  // Tecplot
  outfile << "ZONE" << std::endl;
  
  // Exact solution Vector (here a scalar)
  Vector<double> exact_soln(3);
  
  // Loop over the integration points
  for(unsigned ipt=0; ipt<n_intpt; ipt++)
   {    
    // Assign values of s
    for(unsigned i=0; i<2; i++)
     {
      s[i] = integral_pt()->knot(ipt,i);
     }

    // Get the integral weight
    double w = integral_pt()->weight(ipt);
    
    // Get jacobian of mapping
    double J = J_eulerian(s);
    
    // Premultiply the weights and the Jacobian
    double W = w*J;
    
    // Get x position as Vector
    interpolated_x(s,x);
    
    // Get FE function values
    Vector<double> p_fe(3);
    
    // Velocities go first (a bit naught but OK -- we're
    // only filling in the first two values)
    get_velocity(s,p_fe);
    
    // pressure is last
    p_fe[2] = interpolated_p_hele_shaw(s);
    
    // Get exact solution at this point
    (*exact_soln_pt)(x,exact_soln);
    
    // Output x,y,...,error
    for(unsigned i=0; i<2; i++)
     {
      outfile << x[i] << " ";
     }
    outfile
     << exact_soln[0] - p_fe[0] << " "
     << exact_soln[1] - p_fe[1] << " "
     << exact_soln[2] - p_fe[2] << "\n";
    
    // Add to error and norm
    for(unsigned i=0; i<3; i++)
     {
      norm += exact_soln[i]*exact_soln[i]*W;
      error += (exact_soln[i]-p_fe[i])*(exact_soln[i] - p_fe[i])*W;
     }
   }
 }
 
 
 
 
 
//====================================================================
// Force build of templates
//====================================================================
 template class QHeleShawElement<2>;
 template class QHeleShawElement<3>;
 template class QHeleShawElement<4>;
 
 
}




//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC//           Version 0.90. August 3, 2009.
//LIC// 
//LIC// Copyright (C) 2006-2009 Matthias Heil and Andrew Hazel
//LIC// 
//LIC// This library is free software; you can redistribte it and/or
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

/// Non-inline functions for HartreeFock elements

// Include the appropriate headers
#include "hartree_fock_elements.h"


namespace oomph
{


//======================================================================
 /// \short Static member data indicating which nodal value is to be plotted
 /// (a bit hacky but cheaper than having to loop over all 
 /// elements to tell them what to plot...)
//======================================================================
 template <unsigned DIM>
 unsigned HartreeFockEquations<DIM>::Value_to_be_plotted=0;

//======================================================================
/// Set the data for the number of Variables at each node.
//======================================================================
 template<unsigned DIM, unsigned NNODE_1D>
 const unsigned QHartreeFockElement<DIM,NNODE_1D>::Initial_Nvalue = 2;


//======================================================================
/// Create storage for n_stored_eigenvectors eigenvectors
//======================================================================
 template<unsigned DIM>
 void HartreeFockEquations<DIM>::create_storage_for_eigenvectors(
  const unsigned& n_stored_eigenvectors)
 { 
  /// Find out how many nodes there are in the element
  const unsigned n_node = nnode();
  
  /// Loop over nodes
  for(unsigned l=0;l<n_node;l++) 
   {
    node_pt(l)->resize(n_stored_eigenvectors+1);
    unsigned n_value=node_pt(l)->nvalue();
    for (unsigned i=1;i<n_value;i++)
     {
      node_pt(l)->pin(i);
     }
   }
  
  /// Remember number of stored eigenvectors
  Nstored_eigenvectors=n_stored_eigenvectors;

 }

 
//======================================================================
/// Copy values currently stored in nodal value 0 into 
/// storage for e-th eigenvector [e=0....n_stored_eigenvectors-1]. 
/// This requires a two-step process: First the
/// eigenvectors must be assigned to nodal values, using 
/// Problem::assign_eigenvector_to_dofs(...), then this function
/// should be used to store the eigenvectors permanently.
//======================================================================
template <unsigned DIM>
void HartreeFockEquations<DIM>::store_eigenvector(const unsigned& e)
{
#ifdef PARANOID
 if (e>Nstored_eigenvectors)
  {
   throw OomphLibError(
    "Trying to store more eigenvectors than there's storage for!",
    "HartreeFockEquations<DIM>::store_eigenvector()",
    OOMPH_EXCEPTION_LOCATION);
  }
#endif

 /// Loop over nodes
 const unsigned n_node=nnode();
 for(unsigned l=0;l<n_node;l++) 
  {
   node_pt(l)->set_value(e+1,node_pt(l)->value(0));
  }
}

//======================================================================
/// Compute  hierher
//======================================================================
template <unsigned DIM>
void  HartreeFockEquations<DIM>::
fill_in_contribution_to_jacobian_and_mass_matrix(Vector<double> &residuals, 
                                              DenseMatrix<double> &jacobian, 
                                              DenseMatrix<double> &mass_matrix)
{
 /// Retrieve the index of this element
 unsigned element_number= Element_number;

 /// Find out how many nodes there are
 const unsigned n_node = nnode();

 /// Set up memory for the shape and test functions
 Shape psi(n_node), test(n_node);
 DShape dpsidx(n_node,DIM), dtestdx(n_node,DIM);

 /// Index at which the hartree_fock unknown is stored
 const unsigned u_nodal_index = u_index_hartree_fock();
 
 /// Set the value of n_intpt
 const unsigned n_intpt = integral_pt()->nweight();

 /// Integers to store the local and global equation and unknown numbers
 int local_eqn=0, local_unknown=0, global_unknown=0;

 /// Doubles to store source terms
 double nuclear_source=0.0, repulsion_source=0.0, exchange_source=0.0;

 /// Loop over the integration points
 for(unsigned ept=0;ept<n_intpt;ept++)
  {
   /// Get the integral weight
   double w = integral_pt()->weight(ept);

   /// Call the derivatives of the shape and test functions
   double J = dshape_and_dtest_eulerian_at_knot_hartree_fock(ept,psi,dpsidx,
                                                        test,dtestdx);

   /// Premultiply the weights and the Jacobian
   double W = w*J;

   /// Allocate coords and initialise to zero
   Vector<double> interpolated_x(DIM,0.0);

   /// Calculate the integration point coordinates
   
   /// Loop over nodes
   for(unsigned l=0;l<n_node;l++) 
    {
     /// Loop over directions
     for(unsigned j=0;j<DIM;j++)
      {
       interpolated_x[j] += raw_nodal_position(l,j)*psi(l);
      }
    }


   /// Get Nuclear source term
   get_nuclear_source(interpolated_x,nuclear_source);

   /// Get Electron repulsion source term
   get_electron_repulsion(element_number,ept,repulsion_source);


   /// Assemble Jacobian and Mass Matrix
       
   /// Loop over the test functions
   for(unsigned l=0;l<n_node;l++)
    {
     /// Get the local equation
     local_eqn = nodal_local_eqn(l,u_nodal_index);

     /// If it's not a boundary condition
     if(local_eqn >= 0)
      {
       /// Loop over the shape functions again
       for(unsigned l2=0;l2<n_node;l2++)
        { 
         /// Set local unknown
         local_unknown = nodal_local_eqn(l2,u_nodal_index);

         /// Get global equation from node
         global_unknown = this->node_pt(l2)->eqn_number(u_nodal_index);
         
         /// If at a non-zero degree of freedom add in the entry
         if(local_unknown >= 0)
          {
           ///  Get local electron exchange source
           get_electron_exchange(element_number,
                                 ept,
                                 global_unknown,
                                 exchange_source);
           
           /// Add contribution to mass matrix
           mass_matrix(local_eqn,local_unknown)
            += W*psi(l)*psi(l2);

           /// Add contribution to stiffness matrix from nuclear potential
           jacobian(local_eqn,local_unknown) 
            -= W*psi(l)*nuclear_source*psi(l2);

           /// Add contribution to stiffness matrix from electron repulsion
           jacobian(local_eqn,local_unknown) 
            +=W*psi(l)*repulsion_source*psi(l2);

           /// Add contribution to stiffness matrix from local 
           /// electron exchange
           jacobian(local_eqn,local_unknown)
            -=W*psi(l)*exchange_source;

           for(unsigned i=0;i<DIM;i++)
            {
             /// Add contribution to stiffness matrix from kinetic energy
             jacobian(local_eqn,local_unknown) 
              += 0.5*dpsidx(l2,i)*dtestdx(l,i)*W;
            }

          }

        } /// End of loop over unknown shape functions
       
       //ashflag Loop over external data
       unsigned n_external_data=this->nexternal_data();
       for(unsigned ex=0;ex<n_external_data;ex++)
        {
         /// Return the local equation number
         local_unknown = this->external_local_eqn(ex,u_nodal_index);
         
         /// Return the global equation number
         global_unknown = this->external_data_pt(ex)
          ->eqn_number(u_nodal_index);


         if(local_unknown >= 0)
           {
            ///  Get non-local electron exchange source
            get_electron_exchange(element_number,
                                  ept,
                                  global_unknown,
                                  exchange_source);
            
            /// Add contribution to stiffness matrix from 
            /// non-local electron exchange
            jacobian(local_eqn,local_unknown)
             -=W*psi(l)*exchange_source;
           }

        } /// end of loop over external data
       
      }
    } /// end of loop over test functions

  } /// end of loop over integration points

} /// end of fill_in_contribution_to_jacobian_and_mass_matrix 


//======================================================================
/// Self-test:  Return 0 for OK
//======================================================================
template <unsigned DIM>
unsigned  HartreeFockEquations<DIM>::self_test()
{

 bool passed=true;

 /// Check lower-level stuff
 if (FiniteElement::self_test()!=0)
  {
   passed=false;
  }

 /// Return verdict
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
template <unsigned DIM>
void  HartreeFockEquations<DIM>::output(std::ostream &outfile, 
                                        const unsigned &nplot)
{
 //Vector of local coordinates
 Vector<double> s(DIM);
 
 // Tecplot header info
 outfile << tecplot_zone_string(nplot);
 
 // Loop over plot points
 unsigned num_plot_points=nplot_points(nplot);

for (unsigned iplot=0;iplot<num_plot_points;iplot++)
  {
   
   // Get local coordinates of plot point
   get_s_plot(iplot,nplot,s);

   for(unsigned i=0;i<DIM;i++) 
    {
     outfile << interpolated_x(s,i) << " ";
    }
   outfile << std::pow(interpolated_u_hartree_fock(s,Value_to_be_plotted),2.0)
           << std::endl;

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
template <unsigned DIM>
void  HartreeFockEquations<DIM>::output(FILE* file_pt,
                                        const unsigned &nplot)
{
 //Vector of local coordinates
 Vector<double> s(DIM);
 
 // Tecplot header info
 fprintf(file_pt,"%s",tecplot_zone_string(nplot).c_str());

 // Loop over plot points
 unsigned num_plot_points=nplot_points(nplot);
 for (unsigned iplot=0;iplot<num_plot_points;iplot++)
  {
   // Get local coordinates of plot point
   get_s_plot(iplot,nplot,s);
   
   for(unsigned i=0;i<DIM;i++) 
    {
     fprintf(file_pt,"%g ",interpolated_x(s,i));
    }
   fprintf(file_pt,"%g \n",interpolated_u_hartree_fock(s,Value_to_be_plotted));
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
template <unsigned DIM>
void HartreeFockEquations<DIM>::output_fct(std::ostream &outfile, 
                                       const unsigned &nplot, 
                  FiniteElement::SteadyExactSolutionFctPt exact_soln_pt)
{

 
 //Vector of local coordinates
 Vector<double> s(DIM);
  
  // Vector for coordintes
  Vector<double> x(DIM);
  
 // Tecplot header info
 outfile << tecplot_zone_string(nplot);
 
 // Exact solution Vector (here a scalar)
 Vector<double> exact_soln(1);
 
 // Loop over plot points
 unsigned num_plot_points=nplot_points(nplot);
 

 for (unsigned iplot=0;iplot<num_plot_points;iplot++)
  {
   
   // Get local coordinates of plot point
   get_s_plot(iplot,nplot,s);
   
   // Get x position as Vector
   interpolated_x(s,x);
   
   // Get exact solution at this point
   (*exact_soln_pt)(x,exact_soln);
   
   //Output x,y,...,u_exact
   for(unsigned i=0;i<DIM;i++)
    {
     outfile << x[i] << " ";
    }
   outfile << exact_soln[0] << std::endl;  
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
template <unsigned DIM>
void HartreeFockEquations<DIM>::compute_error(std::ostream &outfile, 
                                          FiniteElement::SteadyExactSolutionFctPt exact_soln_pt,
                                          double& error, double& norm)
{ 
 
 // Initialise
 error=0.0;
 norm=0.0;
 
 //Vector of local coordinates
 Vector<double> s(DIM);
 
 // Vector for coordintes
 Vector<double> x(DIM);
 
 //Find out how many nodes there are in the element
 unsigned n_node = nnode();
 
 Shape psi(n_node);
 
 //Set the value of n_intpt
 unsigned n_intpt = integral_pt()->nweight();
  
 // Tecplot 
 outfile << "ZONE" << std::endl;
 
 // Exact solution Vector (here a scalar)
 Vector<double> exact_soln(1);
 
 //Loop over the integration points
 for(unsigned ipt=0;ipt<n_intpt;ipt++)
  {
   
   //Assign values of s
   for(unsigned i=0;i<DIM;i++)
    {
     s[i] = integral_pt()->knot(ipt,i);
    }
   
   //Get the integral weight
   double w = integral_pt()->weight(ipt);
   
   // Get jacobian of mapping
   double J=J_eulerian(s);
   
   //Premultiply the weights and the Jacobian
   double W = w*J;
   
   // Get x position as Vector
   interpolated_x(s,x);
   
   // Get FE function value
   double u_fe=interpolated_u_hartree_fock(s,0);
   
   // Get exact solution at this point
   (*exact_soln_pt)(x,exact_soln);
   
   //Output x,y,...,error
   for(unsigned i=0;i<DIM;i++)
    {
     outfile << x[i] << " ";
    }
   outfile << exact_soln[0] << " " << exact_soln[0]-u_fe << std::endl;  
   
   // Add to error and norm
   norm+=exact_soln[0]*exact_soln[0]*W;
   error+=(exact_soln[0]-u_fe)*(exact_soln[0]-u_fe)*W;
   
  }
}





//====================================================================
// Force build of templates
//====================================================================
template class HartreeFockEquations<1>;
template class QHartreeFockElement<1,2>;
template class QHartreeFockElement<1,3>;
template class QHartreeFockElement<1,4>;

template class HartreeFockEquations<2>;
template class QHartreeFockElement<2,2>;
template class QHartreeFockElement<2,3>;
template class QHartreeFockElement<2,4>;

template class HartreeFockEquations<3>;
template class QHartreeFockElement<3,2>;
template class QHartreeFockElement<3,3>;
template class QHartreeFockElement<3,4>;

}

//Non-inline functions for Helmholtz elements
#include "helmholtz_elements.h"

namespace oomph
{

//======================================================================
/// Set the data for the number of Variables at each node, always one
/// in every case
//======================================================================
 template<unsigned DIM, unsigned NNODE_1D>
 const unsigned QHelmholtzElement<DIM,NNODE_1D>::Initial_Nvalue = 1;


//======================================================================
/// Compute element residual Vector and/or element Jacobian matrix 
/// 
/// flag=1: compute both
/// flag=0: compute only residual Vector
///
/// Pure version without hanging nodes
//======================================================================
template <unsigned DIM>
void HelmholtzEquations<DIM>::
fill_in_generic_residual_contribution_helmholtz(Vector<double> &residuals, 
                                              DenseMatrix<double> &jacobian, 
                                              unsigned flag) 
{
 //Find out how many nodes there are
 const unsigned n_node = nnode();

 //Set up memory for the shape and test functions
 Shape psi(n_node), test(n_node);
 DShape dpsidx(n_node,DIM), dtestdx(n_node,DIM);

 //Index at which the helmholtz unknown is stored
 const unsigned u_nodal_index = u_index_helmholtz();
 
 //Set the value of n_intpt
 const unsigned n_intpt = integral_pt()->nweight();

 //Integers to store the local equation and unknown numbers
 int local_eqn=0, local_unknown=0;

 //Loop over the integration points
 for(unsigned ipt=0;ipt<n_intpt;ipt++)
  {
   //Get the integral weight
   double w = integral_pt()->weight(ipt);

   //Call the derivatives of the shape and test functions
   double J = dshape_and_dtest_eulerian_at_knot_helmholtz(ipt,psi,dpsidx,
                                                        test,dtestdx);
       
   //Premultiply the weights and the Jacobian
   double W = w*J;

   //Calculate local values of the pressure and velocity components
   //Allocate and initialise to zero
//   double interpolated_u=0.0;
   Vector<double> interpolated_x(DIM,0.0);
   Vector<double> interpolated_dudx(DIM,0.0);
   
   //Calculate function value and derivatives:
   //-----------------------------------------
   // Loop over nodes
   for(unsigned l=0;l<n_node;l++) 
    {
     //Get the nodal value of the helmholtz unknown
     double u_value = raw_nodal_value(l,u_nodal_index);

     // Loop over directions
     for(unsigned j=0;j<DIM;j++)
      {
       interpolated_x[j] += raw_nodal_position(l,j)*psi(l);
       interpolated_dudx[j] += u_value*dpsidx(l,j);
      }
    }

   //Get source function
   //-------------------
   double source;
   get_source_helmholtz(interpolated_x,source);

   // Access the local coordinates for this integration point
   // (these will have been assigned earlier when the other element
   //  was assigned (currently in the driver two_mesh_helmholtz.cc))
   // Note: interaction parameter default to zero for now
   Vector<double> s_other(DIM);
   s_other=external_element_local_coord(0,ipt);

   // Add contribution to source from "other" element
   HelmholtzEquations<DIM>* source_el_pt=dynamic_cast<HelmholtzEquations<DIM>*>
    (external_element_pt(0,ipt));
   source+=source_el_pt->interpolated_u_helmholtz(s_other);

//    oomph_info << "ipt=" << ipt << ", interpolated_x=(" << interpolated_x[0]
//               << "," << interpolated_x[1] << "), s_other=("
//               << s_other[0] << "," << s_other[1] << "), source="
//               << source << std::endl;

   // Assemble residuals and Jacobian
   //--------------------------------
       
   // Loop over the test functions
   for(unsigned l=0;l<n_node;l++)
    {
     //Get the local equation
     local_eqn = nodal_local_eqn(l,u_nodal_index);
     /*IF it's not a boundary condition*/
     if(local_eqn >= 0)
      {
       // Add body force/source term here and the Helmholtz bit
       residuals[local_eqn] += source*test(l)*W;
//       residuals[local_eqn] += (interpolated_u+source)*test(l)*W;
             
       // The operator bit
       for(unsigned k=0;k<DIM;k++)
        {
         residuals[local_eqn] += interpolated_dudx[k]*dtestdx(l,k)*W;
        }

       // Calculate the jacobian
       //-----------------------
       if(flag)
        {
         //Loop over the velocity shape functions again
         for(unsigned l2=0;l2<n_node;l2++)
          { 
           local_unknown = nodal_local_eqn(l2,u_nodal_index);
           //If at a non-zero degree of freedom add in the entry
           if(local_unknown >= 0)
            {
             oomph_info << "check ";
             //Add contribution to Elemental Matrix
             for(unsigned i=0;i<DIM;i++)
              {
               jacobian(local_eqn,local_unknown) 
                += (
//                 psi(l2)*test(l)+
                 dpsidx(l2,i)*dtestdx(l,i)
                 )*W;
              }
            }
          }
        }
      }
    }

  } // End of loop over integration points
}   

// //======================================================================
// /// Flush external data, repopulate uniquely, and add equation numbers?
// //======================================================================
// template <unsigned DIM>
// void HelmholtzEquations<DIM>::assign_unique_external_data_helper()
// {
//  //Does external data exist?
//  const unsigned n_external_data = nexternal_data();
//  const unsigned n_internal_data = ninternal_data();
// //  oomph_info << "n_external_data=" << n_external_data << std::endl;
//  //If it does, back it up
//  if (n_external_data > 0)
//   {
//    Vector<Data*> backed_up_external_data_pt(n_external_data);
// //   Data** backed_up_external_data_pt=new Data*[n_external_data];
//    for (unsigned i=0;i<n_external_data;i++)
//     {
//      backed_up_external_data_pt[i] = external_data_pt(i);
//     }
//    // Now flush the external data
//    flush_external_data();
//    // Loop over backed up external data, checking it against nodal 
//    // data and internal data for duplicates; add it if it is unique
//    for (unsigned i=0;i<n_external_data;i++)
//     {
//      bool is_unique=true;
//      // Check against nodal data
//      unsigned n_node=nnode();
//      for (unsigned j=0;j<n_node;j++)
//       {
//        if (backed_up_external_data_pt[i]==node_pt(j))
//         {
// //          oomph_info << "backed_up_external_data_pt[" << i << "]="
// //                     << backed_up_external_data_pt[i] << ", node_pt(" << j
// //                     << ")=" << node_pt(j) << std::endl;
//          is_unique=false;
//         }
//       }

//      // Check against internal data (if it hasn't already been found)
//      if (is_unique)
//       {
//        for (unsigned k=0;k<n_internal_data;k++)
//         {
//          if (backed_up_external_data_pt[i]==internal_data_pt(k))
//           {
//            is_unique=false;
//           }
//         }
//       }

//      // If both of these checks do not find a duplicate...
//      unsigned extern_index;
//      if (is_unique)
//       {
//        // .. then add it back in again
//        extern_index=add_external_data(backed_up_external_data_pt[i]);
//       }
//     }
// //   oomph_info << "after uniqueness check, nexternal_data()=" 
// //              << nexternal_data() << std::endl;
//   }
// }


//======================================================================
/// Self-test:  Return 0 for OK
//======================================================================
template <unsigned DIM>
unsigned  HelmholtzEquations<DIM>::self_test()
{

 bool passed=true;

 // Check lower-level stuff
 if (FiniteElement::self_test()!=0)
  {
   passed=false;
  }

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
template <unsigned DIM>
void  HelmholtzEquations<DIM>::output(std::ostream &outfile, 
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
   outfile << interpolated_u_helmholtz(s) << std::endl;   
   
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
void  HelmholtzEquations<DIM>::output(FILE* file_pt,
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
   fprintf(file_pt,"%g \n",interpolated_u_helmholtz(s));
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
void HelmholtzEquations<DIM>::output_fct(std::ostream &outfile, 
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
void HelmholtzEquations<DIM>::compute_error(std::ostream &outfile, 
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
   double u_fe=interpolated_u_helmholtz(s);
   
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
template class QHelmholtzElement<1,2>;
template class QHelmholtzElement<1,3>;
template class QHelmholtzElement<1,4>;

template class QHelmholtzElement<2,2>;
template class QHelmholtzElement<2,3>;
template class QHelmholtzElement<2,4>;

template class QHelmholtzElement<3,2>;
template class QHelmholtzElement<3,3>;
template class QHelmholtzElement<3,4>;

}


//Non-inline functions for Fisher elements
#include "fisher_elements.h"


namespace oomph
{


//2D Fisher elements



/// Fisher equations static data
template<unsigned DIM>
double FisherEquations<DIM>::Default_Fisher_Factor_Value = 1.0;


//======================================================================
// Set the data for the number of Variables at each node
//======================================================================
template<>
const unsigned QFisherElement<1,4>::Initial_Nvalue[4]={1,1,1,1};
template<>
const unsigned QFisherElement<1,3>::Initial_Nvalue[3]={1,1,1};
template<>
const unsigned QFisherElement<1,2>::Initial_Nvalue[2]={1,1};

template<>
const unsigned QFisherElement<2,4>::Initial_Nvalue[16]
={1,1,1,1,1,1,1,1,
  1,1,1,1,1,1,1,1};
template<>
const unsigned QFisherElement<2,3>::Initial_Nvalue[9]
={1,1,1,1,1,1,1,1,1};
template<>
const unsigned QFisherElement<2,2>::Initial_Nvalue[4]
={1,1,1,1};

template<>
const unsigned QFisherElement<3,2>::Initial_Nvalue[8]
={1,1,1,1,1,1,1,1};
template<>
const unsigned QFisherElement<3,3>::Initial_Nvalue[27]
={1,1,1,1,1,1,1,1,1,
  1,1,1,1,1,1,1,1,1,
  1,1,1,1,1,1,1,1,1};
template<>
const unsigned QFisherElement<3,4>::Initial_Nvalue[64]
={1,1,1,1,1,1,1,1,
  1,1,1,1,1,1,1,1,
  1,1,1,1,1,1,1,1,
  1,1,1,1,1,1,1,1,
  1,1,1,1,1,1,1,1,
  1,1,1,1,1,1,1,1,
  1,1,1,1,1,1,1,1,
  1,1,1,1,1,1,1,1};



//======================================================================
/// Compute element residual Vector and/or element Jacobian matrix 
/// 
/// flag=1: compute both
/// flag=0: compute only residual Vector
///
/// Pure version without hanging nodes
//======================================================================
template <unsigned DIM>
void  FisherEquations<DIM>::
fill_in_generic_residual_contribution(Vector<double> &residuals, 
                                  DenseMatrix<double> &jacobian, 
                                  unsigned flag) 
{
 //Find out how many nodes there are
 unsigned n_node = nnode();
  
 //Set up memory for the shape and test functions
 Shape psi(n_node), test(n_node);
 DShape dpsidx(n_node,DIM), dtestdx(n_node,DIM);
 
 //Set the value of n_intpt
 unsigned n_intpt = integral_pt()->nweight();
   
 //Set the Vector to hold local coordinates
 Vector<double> s(DIM);

 //Integers to hold the local equation and unknowns
 int local_eqn=0, local_unknown=0;

 //Loop over the integration points
 for(unsigned ipt=0;ipt<n_intpt;ipt++)
  {

   //Assign values of s
   for(unsigned i=0;i<DIM;i++) s[i] = integral_pt()->knot(ipt,i);

   //Get the integral weight
   double w = integral_pt()->weight(ipt);

   //Call the derivatives of the shape and test functions
   double J = dshape_and_dtest_eulerian_at_knot(ipt,psi,dpsidx,test,dtestdx);

   //Premultiply the weights and the Jacobian
   double W = w*J;

   //Allocate memory for local quantities and initialise to zero
   double interpolated_u=0.0;
   double dudt=0.0;
   Vector<double> interpolated_x(DIM,0.0);
   Vector<double> interpolated_dudx(DIM,0.0);
   Vector<double> mesh_velocity(DIM,0.0);

   //Calculate function value and derivatives:
   // Loop over nodes
   for(unsigned l=0;l<n_node;l++) 
    {
     interpolated_u+=u(l)*psi(l);
     dudt+=du_dt(l)*psi(l);
     // Loop over directions
     for(unsigned j=0;j<DIM;j++)
      {
       mesh_velocity[j] += dnodal_position_dt(l,j)*psi(l);
       interpolated_x[j] += nodal_position(l,j)*psi(l);
       interpolated_dudx[j] += u(l)*dpsidx(l,j);
      }
    }


   //Get source function
   //-------------------
   double source;
   get_source(time(),interpolated_x,source);

   // Assemble residuals and Jacobian
   //--------------------------------
       
   // Loop over the test functions
   for(unsigned l=0;l<n_node;l++)
    {
     local_eqn = u_local_eqn(l);
     /*IF it's not a boundary condition*/
     if(local_eqn >= 0)
      {

       // Add body force/source term and time derivative
       residuals[local_eqn] += (dudt-
                                fisher_factor()*
                                 fisher_factor()*
                                interpolated_u*(1.0-interpolated_u)
                                +source)*test(l)*W;
           
       // The mesh velocity bit
       for(unsigned k=0;k<DIM;k++)
        {
         residuals[local_eqn] -= 
          mesh_velocity[k]*interpolated_dudx[k]*test(l)*W;
        }
       
       // Laplace operator
       for(unsigned k=0;k<DIM;k++)
        {
         residuals[local_eqn] += 
          interpolated_dudx[k]*dtestdx(l,k)*W;
        }


       // Calculate the jacobian
       //-----------------------
       if(flag)
        {
         //Loop over the velocity shape functions again
         for(unsigned l2=0;l2<n_node;l2++)
          { 
           local_unknown = u_local_eqn(l2);
           //If at a non-zero degree of freedom add in the entry
           if(local_unknown >= 0)
            {
             // Mass matrix
             jacobian(local_eqn,local_unknown) 
              += test(l)*psi(l2)*node_pt(l2)->time_stepper_pt()->weight(1,0)*W;

             // Fisher bits
             jacobian(local_eqn,local_unknown) 
              -= test(l)*psi(l2)*fisher_factor()*(1.0-2.0*interpolated_u)*W;

             // Laplace operator
             for(unsigned i=0;i<DIM;i++)
              {
               jacobian(local_eqn,local_unknown) 
                += dpsidx(l2,i)*(dtestdx(l,i)-mesh_velocity[i]*test(l))*W;
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
template <unsigned DIM>
unsigned  FisherEquations<DIM>::self_test()
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
void  FisherEquations<DIM>::output(std::ostream &outfile, 
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
   outfile << interpolated_u(s) << std::endl;   
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
void  FisherEquations<DIM>::output(FILE* file_pt,
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
    fprintf(file_pt,"%g \n",interpolated_u(s));
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
void FisherEquations<DIM>::output_fct(std::ostream &outfile, 
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
/// Output exact solution at time t
/// 
/// Solution is provided via function pointer.
/// Plot at a given number of plot points.
///
///   x,y,u_exact    or    x,y,z,u_exact
//======================================================================
template <unsigned DIM>
void FisherEquations<DIM>::output_fct(std::ostream &outfile, 
                                      const unsigned &nplot,
                                      const double& time, 
                         FiniteElement::UnsteadyExactSolutionFctPt exact_soln_pt)

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
   (*exact_soln_pt)(time,x,exact_soln);
   
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
void FisherEquations<DIM>::compute_error(std::ostream &outfile, 
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
  
 // Tecplot header info
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
   double u_fe=interpolated_u(s);
   
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




//======================================================================
/// Validate against exact solution at time t.
/// 
/// Solution is provided via function pointer.
/// Plot error at a given number of plot points.
///
//======================================================================
template<unsigned DIM>
void FisherEquations<DIM>::compute_error(std::ostream &outfile, 
                    FiniteElement::UnsteadyExactSolutionFctPt exact_soln_pt,
                    const double& time, double& error, double& norm)

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
   double u_fe=interpolated_u(s);

   // Get exact solution at this point
   (*exact_soln_pt)(time,x,exact_soln);

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
template class QFisherElement<1,2>;
template class QFisherElement<1,3>;
template class QFisherElement<1,4>;

template class QFisherElement<2,2>;
template class QFisherElement<2,3>;
template class QFisherElement<2,4>;

template class QFisherElement<3,2>;
template class QFisherElement<3,3>;
template class QFisherElement<3,4>;


}

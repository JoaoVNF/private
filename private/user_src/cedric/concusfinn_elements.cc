
//Non-inline functions for YoungLaplace elements
#include "concusfinn_elements.h"

namespace oomph
{


// Default output format
bool YoungLaplaceEquations::Output_meniscus_and_spines=false;


//======================================================================
// Set the data for the number of Variables at each node
//======================================================================
template<>
const unsigned QYoungLaplaceElement<4>::Initial_Nvalue[16]
={1,1,1,1,1,1,1,1,
  1,1,1,1,1,1,1,1};
template<>
const unsigned QYoungLaplaceElement<3>::Initial_Nvalue[9]={1,1,1,1,1,1,1,1,1};
template<>
const unsigned QYoungLaplaceElement<2>::Initial_Nvalue[4]={1,1,1,1};


//======================================================================
/// Get exact position vector to meniscus
//======================================================================
void  YoungLaplaceEquations::exact_position(const Vector<double>& s, 
                                               Vector<double>& r,
                       FiniteElement::SteadyExactSolutionFctPt exact_soln_pt)
{
 // Get global coordinates
 Vector<double> x(2);
 interpolated_x(s,x);

 // Exact solution Vector (here a scalar)
 Vector<double> exact_soln(1);
 
 // Get exact solution at this point
 (*exact_soln_pt)(x,exact_soln);
 
 if (!use_spines())
  {
   r[0]=x[0]; r[1]=x[1]; r[2]=exact_soln[0];
  }
 else
  {
   /// Get spines values
   Vector<double> spine_base(3,0.0);
   Vector<double> spine(3,0.0);
   Vector< Vector<double> > dspine_base;
   allocate_vector_of_vectors(2,3,dspine_base);
   Vector< Vector<double> > dspine;
   allocate_vector_of_vectors(2,3,dspine);
   
   get_spine_base(x,spine_base,dspine_base);
   get_spine(x,spine,dspine);
   
   /// Global Eulerian cooordinates
   for (unsigned j=0; j<3; j++)
    {
     r[j]=spine_base[j]+exact_soln[0]*spine[j];
    }
  }
 
}
 
//======================================================================
/// Get position vector to meniscus
//======================================================================
 void  YoungLaplaceEquations::position(const Vector<double>& s, 
                                       Vector<double>& r)
 {
  
  // Get global coordinates
  Vector<double> x(2);
  interpolated_x(s,x);
  
  //Displacement along spine (or cartesian displacement)
  double u=interpolated_u(s);
  
  // cartesian calculation case
  if (!use_spines())
   { 
    r[0]=x[0]; r[1]=x[1]; r[2]=u;
   }
  // spine case
  else
   {
    /// Get spines values
    Vector<double> spine_base(3,0.0);
    Vector<double> spine(3,0.0);
    Vector< Vector<double> > dspine_base;
    allocate_vector_of_vectors(2,3,dspine_base);
    Vector< Vector<double> > dspine;
    allocate_vector_of_vectors(2,3,dspine);
    
    get_spine_base(x,spine_base,dspine_base);
    get_spine(x,spine,dspine);
    
    /// Global Eulerian cooordinates
    for (unsigned j=0; j<3; j++)
     {
      r[j]=spine_base[j]+u*spine[j];
     }
   }
  
 }
 
 
 
//======================================================================
/// Compute element residual vector. Pure version without hanging nodes
//======================================================================
 void  YoungLaplaceEquations::fill_in_contribution_to_residuals(
  Vector<double> &residuals)
 {
  
  //Find out how many nodes there are
  unsigned n_node = nnode();
 
  //Set up memory for the shape functions
  Shape psi(n_node);
  DShape dpsidzeta(n_node,2);
 
  //Set the value of n_intpt
  unsigned n_intpt = integral_pt()->nweight();
 
  //Integers to store the local equation numbers
  int local_eqn=0; 
 
  //Loop over the integration points
  for(unsigned ipt=0;ipt<n_intpt;ipt++)
   {
    //Get the integral weight
    double w = integral_pt()->weight(ipt);
   
    //Call the derivatives of the shape and test functions
    double J = dshape_eulerian_at_knot(ipt,psi,dpsidzeta);
   
    //Premultiply the weights and the Jacobian
    double W = w*J;
   
    //Calculate local values of displacement along spine and its derivatives
    //Allocate and initialise to zero
    double interpolated_u=0.0;
    Vector<double> interpolated_zeta(2,0.0);
    Vector<double> interpolated_dudzeta(2,0.0);
   
    //Calculate function value and derivatives:
    //-----------------------------------------
    // Loop over nodes
    for(unsigned l=0;l<n_node;l++) 
     {
      interpolated_u += u(l)*psi(l);
      // Loop over directions
      for(unsigned j=0;j<2;j++)
       {
        interpolated_zeta[j] += nodal_position(l,j)*psi(l);
        interpolated_dudzeta[j] += u(l)*dpsidzeta(l,j);
       }
     }
   
 

    // Allocation and definition of variables necessary for
    // further calculations
   
    /// "Simple" case
    ///--------------
    double nonlinearterm=1.0;
    double sqnorm=0.0;
    double normplusone=1.0;
   
    /// Spine case 
    ///-----------

    // Derivs of position vector w.r.t. global intrinsic coords 
    Vector<Vector<double> > dRdzeta; 
    allocate_vector_of_vectors(2,3,dRdzeta);
   
    // Unnormalised normal
    Vector<double> N_unnormalised(3,0.0);
   
    // Spine and spine basis vectors, entries initialised to zero
    Vector<double> spine_base(3,0.0), spine(3,0.0);

    // Derivative of spine basis vector w.r.t to the intrinsic
    // coordinates: dspine_base[i,j] = j-th component of the deriv.
    // of the spine basis vector w.r.t. to the i-th global intrinsic
    // coordinate
    Vector< Vector<double> > dspine_base; 
    allocate_vector_of_vectors(2,3,dspine_base);

    // Derivative of spine vector w.r.t to the intrinsic
    // coordinates: dspine[i,j] = j-th component of the deriv.
    // of the spine vector w.r.t. to the i-th global intrinsic
    // coordinate
    Vector< Vector<double> > dspine; 
    allocate_vector_of_vectors(2,3,dspine);
   
    // Vector v_\alpha contains the numerator of the variations of the 
    // area element {\cal A}^{1/2} w.r.t. the components of dR/d\zeta_\alpha
    Vector<double> area_variation_numerator_0(3,0.0);
    Vector<double> area_variation_numerator_1(3,0.0);
 
    //No spines
    //---------
    if (!use_spines())
     {
      for (unsigned j=0; j<2; j++)
       { 
        sqnorm += interpolated_dudzeta[j]*interpolated_dudzeta[j];
       }
      nonlinearterm=1.0/sqrt(1.0+sqnorm);
      normplusone=1.0+sqnorm;
     }
   
    //Spines
    //------
    else
     {
      /// Get the spines
      get_spine_base(interpolated_zeta, spine_base, dspine_base);
      get_spine(interpolated_zeta, spine, dspine);

      /// calculation of dR/d\zeta_\alpha
      for (unsigned alpha=0;alpha<2;alpha++)
       {
        // Product rule for d(u {\bf S} ) / d \zeta_\alpha
        Vector<double> dudzeta_times_spine(3,0.0);
        scalar_times_vector(interpolated_dudzeta[alpha],
                            spine,dudzeta_times_spine);

        Vector<double> u_times_dspinedzeta(3,0.0);
        scalar_times_vector(interpolated_u,dspine[alpha],u_times_dspinedzeta);

        Vector<double> d_u_times_spine_dzeta(3,0.0);
        vector_sum(dudzeta_times_spine,
                   u_times_dspinedzeta,
                   d_u_times_spine_dzeta);

        // Add derivative of spine base
        vector_sum(d_u_times_spine_dzeta,dspine_base[alpha],dRdzeta[alpha]); 
       }
     
      /// Get the unnormalized normal
      cross_product(dRdzeta[0],dRdzeta[1],N_unnormalised);
     
      // Tmp storage
      Vector<double> v_tmp_1(3,0.0);
      Vector<double> v_tmp_2(3,0.0);
     
      // Calculation of 
      // |dR/d\zeta_1|^2 dR/d\zeta_0 - <dR/d\zeta_0,dR/d\zeta_1>dR/d\zeta_1
      scalar_times_vector(pow(two_norm(dRdzeta[1]),2), dRdzeta[0], v_tmp_1);
      scalar_times_vector(-1*scalar_product(dRdzeta[0],dRdzeta[1]), 
                          dRdzeta[1], v_tmp_2);
      vector_sum(v_tmp_1,v_tmp_2,area_variation_numerator_0);

      // Calculation of 
      // |dR/d\zeta_0|^2 dR/d\zeta_1 - <dR/d\zeta_0,dR/d\zeta_1>dR/d\zeta_0
      scalar_times_vector(pow(two_norm(dRdzeta[0]),2), dRdzeta[1], v_tmp_1);
      scalar_times_vector(-1*scalar_product(dRdzeta[0],dRdzeta[1]), 
                          dRdzeta[0], v_tmp_2);
      vector_sum(v_tmp_1,v_tmp_2,area_variation_numerator_1);
     
     }
   

    // Assemble residuals
    //-------------------
       
    // Loop over the test (shape) functions
    for(unsigned l=0;l<n_node;l++)
     {
      //Get the local equation
      local_eqn = u_local_eqn(l);
     
      /*IF it's not a boundary condition*/
      if(local_eqn >= 0)
       {
       
        // "simple" calculation case
        if (!use_spines())
         { 
          // Add source term: The curvature
          residuals[local_eqn] += get_kappa()*psi(l)*W; 
         
          // The YoungLaplace bit itself
          for(unsigned k=0;k<2;k++)
           {
            residuals[local_eqn] += nonlinearterm*
             interpolated_dudzeta[k]*dpsidzeta(l,k)*W;
           }
         }
       
        // Spine calculation case
        else 
         {
          // Calculation of d(u S)/d\zeta_0
          //-------------------------------
          Vector<double> v_tmp_1(3,0.0);
          scalar_times_vector(dpsidzeta(l,0), spine, v_tmp_1);

          Vector<double> v_tmp_2(3,0.0);
          scalar_times_vector(psi(l),dspine[0], v_tmp_2);

          Vector<double> d_uS_dzeta0(3,0.0);
          vector_sum(v_tmp_1,v_tmp_2,d_uS_dzeta0);

          // Add contribution to residual 
          residuals[local_eqn] += 
           W*scalar_product(area_variation_numerator_0,d_uS_dzeta0)
           /two_norm(N_unnormalised);
         
          // Calculation of d(u S)/d\zeta_1
          scalar_times_vector(dpsidzeta(l,1), spine, v_tmp_1);
          scalar_times_vector(psi(l),dspine[1], v_tmp_2);
          Vector<double> d_uS_dzeta1(3,0.0);
          vector_sum(v_tmp_1,v_tmp_2,d_uS_dzeta1);

          // Add contribution to residual
          residuals[local_eqn] += 
           W*scalar_product(area_variation_numerator_1,d_uS_dzeta1)
           /two_norm(N_unnormalised);

          // Curvature contribution to the residual : kappa N S test
          residuals[local_eqn] += W*get_kappa()*
           scalar_product(N_unnormalised,spine)*psi(l);
         }
       }
     }
   
   } // End of loop over integration points
 
 }



//======================================================================
/// Self-test:  Return 0 for OK
//======================================================================
 unsigned  YoungLaplaceEquations::self_test() 
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
 void  YoungLaplaceEquations::output(std::ostream &outfile, 
                                     const unsigned &nplot)
 {


  //Vector of local coordinates
  Vector<double> s(2);
 
  // Tecplot header info
  outfile << tecplot_zone_string(nplot);
 
  // Loop over plot points
  unsigned num_plot_points=nplot_points(nplot);
  for (unsigned iplot=0;iplot<num_plot_points;iplot++)
   {
   
    // Get local coordinates of plot point
    get_s_plot(iplot,nplot,s);

    // Compute intrinsic coordinates
    Vector<double> xx(2,0.0);
    for(unsigned i=0;i<2;i++) 
     {
      xx[i]=interpolated_x(s,i) ;
     }
   

    // Choose output format
    if (Output_meniscus_and_spines)
     {
      // Output intrinsic coordinates
      for(unsigned i=0;i<2;i++) 
       {
        outfile << xx[i] << " ";
       }
      
      //Output displacment along spine
      outfile << interpolated_u(s) << " ";
     }
    else
     {
     
      /// Calculate the cartesian coordinates of point on meniscus
      Vector<double> r(3,0.0);
      position(s,r);
     
      // Output positon on meniscus
      for(unsigned i=0;i<3;i++) 
       {
        outfile << r[i] << " ";
       }
     
     
      // Select output format
      if (use_spines())
       {
        // Get spine stuff
        Vector<double> spine_base(3,0.0), spine(3,0.0);
        Vector< Vector<double> > dspine_base; 
        allocate_vector_of_vectors(2,3,dspine_base);
        Vector< Vector<double> > dspine; 
        allocate_vector_of_vectors(2,3,dspine);
       
        if (use_spines())
         { 
          /// Get the spines
          get_spine_base(xx, spine_base, dspine_base);
          get_spine(xx, spine, dspine);
         }
       
       
        // Output spine base
        for(unsigned i=0;i<3;i++) 
         {
          outfile << spine_base[i] << " ";
         }  
       
        // Output spines
        for(unsigned i=0;i<3;i++) 
         {
          outfile << spine[i] << " ";
         }
       }
     }

    // Done
    outfile << std::endl;
   
   }
 
  // Write tecplot footer (e.g. FE connectivity lists)
  write_tecplot_zone_footer(outfile,nplot);
 
 }

//*******************************************************
// CEDRIC: Temporarily taken out -- can be reinstated once
// the standard version has reached its final form
//*******************************************************

// //======================================================================
// /// C-style output function:
// ///
// ///   x,y,u   or    x,y,z,u
// ///
// /// nplot points in each coordinate direction
// //======================================================================
// template <2>
// void  YoungLaplaceEquations<2>::output(FILE* file_pt,
//                                     const unsigned &nplot)
// {
//  //Vector of local coordinates
//  Vector<double> s(2);
 
//  // Tecplot header info
//  fprintf(file_pt,"%s",tecplot_zone_string(nplot).c_str());

//  // Loop over plot points
//  unsigned num_plot_points=nplot_points(nplot);
//  for (unsigned iplot=0;iplot<num_plot_points;iplot++)
//   {
//    // Get local coordinates of plot point
//    get_s_plot(iplot,nplot,s);
   
//    for(unsigned i=0;i<2;i++) 
//     {
//      fprintf(file_pt,"%g ",interpolated_x(s,i));
//     }
//    fprintf(file_pt,"%g \n",interpolated_u(s));
//   }

//  // Write tecplot footer (e.g. FE connectivity lists)
//  write_tecplot_zone_footer(file_pt,nplot);
// }



 //======================================================================
 /// Output exact solution
 /// 
 /// Solution is provided via function pointer.
 /// Plot at a given number of plot points.
 ///
 ///   x,y,u_exact    or    x,y,z,u_exact
 //======================================================================
  void YoungLaplaceEquations::output_fct(std::ostream &outfile, 
                                          const unsigned &nplot, 
                  FiniteElement::SteadyExactSolutionFctPt exact_soln_pt)
 {

  //Vector of local coordinates
  Vector<double> s(2);
  
  // Vector for coordinates
  Vector<double> x(2);
  
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

   /// Calculate the cartesian coordinates of point on meniscus
   Vector<double> r_exact(3,0.0);
   exact_position(s,r_exact,exact_soln_pt);
   
   //Output x_exact,y_exact,z_exact
   for(unsigned i=0;i<3;i++)
    {
     outfile << r_exact[i] << " ";
    }

   // Done
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
 void YoungLaplaceEquations::compute_error(std::ostream &outfile,  
                    FiniteElement::SteadyExactSolutionFctPt exact_soln_pt, 
                    double& error, double& norm)
 {

  // Initialise
  error=0.0;
  norm=0.0;
  
  //Vector of local coordinates
  Vector<double> s(2);
  
  // Vector for coordinates
  Vector<double> x(2);
  
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
    for(unsigned i=0;i<2;i++)
     {
      s[i] = integral_pt()->knot(ipt,i);
     }
    
    //Get the integral weight
    double w = integral_pt()->weight(ipt);
    
    // Get jacobian of mapping
    double J=J_eulerian(s);
    
    //Premultiply the weights and the Jacobian
    double W = w*J;
    
    /// Calculate the cartesian coordinates of point on meniscus
    Vector<double> r(3,0.0);
    position(s,r);
    
    /// Calculate the exact position
    Vector<double> r_exact(3,0.0);
    exact_position(s,r_exact,exact_soln_pt);
    
    //Output x,y,...,error
    for(unsigned i=0;i<2;i++)
     {
      outfile << r[i] << " ";
     }
    
    for(unsigned i=0;i<2;i++)
     {
      outfile << r_exact[i] << " ";
     }
    
    outfile << std::endl;  
    
    // Add to error and norm
    norm+=0.0;
    for(unsigned i=0;i<2;i++)
     {
      error+=(r[i]-r_exact[i])*(r[i]-r_exact[i])*W;
     }
   }
  
 }



//====================================================================
// Force build of templates
//====================================================================
template class QYoungLaplaceElement<2>;
template class QYoungLaplaceElement<3>;
template class QYoungLaplaceElement<4>;

}





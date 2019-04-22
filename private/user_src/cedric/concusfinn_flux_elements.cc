#include "concusfinn_flux_elements.h"
#include "concusfinn_elements.h"
#include "refineable_concusfinn_elements.h"


namespace oomph
{

//===========================================================================
/// Constructor, takes the pointer to the "bulk" element, the 
/// index of the fixed local coordinate and its value represented
/// by an integer (+/- 1), indicating that the face is located
/// at the max. or min. value of the "fixed" local coordinate
/// in the bulk element.
//===========================================================================
template<class ELEMENT>
YoungLaplaceFluxElement<ELEMENT>::
YoungLaplaceFluxElement(FiniteElement* const &bulk_el_pt, 
                        const int &face_index) : 
 FaceGeometry<ELEMENT>(), FaceElement()
  { 

   // Let the bulk element build the FaceElement, i.e. setup the pointers 
   // to its nodes (by referring to the appropriate nodes in the bulk
   // element), etc.
   bulk_el_pt->build_face_element(face_index,this);
 
   // Initialise the prescribed contact angle pointer to zero
   Prescribed_cos_gamma_pt = 0;
 
#ifdef PARANOID
   // Extract the dimension of the problem from the dimension of 
   // the first node
   unsigned dim = node_pt(0)->ndim();
   if (dim!=2) 
    {
     throw OomphLibError(
      "YoungLaplaceFluxElement only work with 2D nodes",
      "YoungLaplaceFluxElement<ELEMENT>::YoungLaplaceFluxElement<ELEMENT>",
      OOMPH_EXCEPTION_LOCATION);
    }
#endif
  
  }
 


//===========================================================================
/// Compute the element's contribution to the residual vector 
/// (old version)
//===========================================================================
/*template<class ELEMENT>
void YoungLaplaceFluxElement<ELEMENT>::
old_fill_in_contribution_to_residuals(Vector<double> &residuals)

{
 //Find out how many nodes there are
 unsigned n_node = nnode();

 //Set up memory for the shape (test) functions,
 // and derivatives
 Shape psif(n_node);
 DShape dpsifds(n_node,1);
  
 //Set the value of Nintpt
 unsigned n_intpt = integral_pt()->nweight();
 
 //Set the Vector to hold local coordinates
  Vector<double> s(1);
 
 //Integers to hold the local equation and unknown numbers
 int local_eqn=0;

 //Loop over the integration points
 //--------------------------------
 for(unsigned ipt=0;ipt<n_intpt;ipt++)
  {

   //Assign values of s
   for(unsigned i=0;i<1;i++) {s[i] = integral_pt()->knot(ipt,i);}
   
   //Get the integral weight
   double w = integral_pt()->weight(ipt);

   //Find the shape (test) functions and their deivatives w.r.t. 
   // the local coordinates
   dshape_local(s,psif,dpsifds);

   //Need to find position to feed into flux function, initialise to zero
   //Allocate and initialise to zero
   double interpolated_u=0.0;
   double interpolated_duds=0.0;
   Vector<double> interpolated_x(2,0.0);
   Vector<double> interpolated_dxds(2,0.0);
   
   //Calculate function value and derivatives:
   //-----------------------------------------
   // Loop over nodes
   for(unsigned l=0;l<n_node;l++) 
    {
     interpolated_u += u(l)*psif(l);
     interpolated_duds += u(l)*dpsifds(l,0);

     // Loop over directions
     for(unsigned j=0;j<2;j++)
      {
       interpolated_x[j] += nodal_position(l,j)*psif(l);
       interpolated_dxds[j] += nodal_position(l,j)*dpsifds(l,0);
      }
    }

   //Get the prescribed cos gamma
   double cos_gamma=prescribed_cos_gamma();

   // Some initialisation
   double nm=1.0; 
   double alpha_n=1.0;

   //Precalculation spine
   //---------------------

   ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>(this->bulk_element_pt());
   
   if (bulk_elem_pt->use_spines())
    {
     ///Get the spines at position x1 x2
     ///We are in the 2D case !!!!!!!!!!!
     ///---------------------------------
          
     /// Vector allocation tangent, normal, spines
     Vector<double> dRds(3,0.0);
     Vector<double> N(3,0.0); N[0]=0.0; N[1]=0.0; N[2]=1.0;
     Vector<double> spine_B(3,0.0), spine_S(3,0.0);

     /// and spines derivatives
     Vector< Vector<double> >  dspine_B, dspine_S;
     bulk_elem_pt->allocate_vector_of_vectors(2,3,dspine_B);
     bulk_elem_pt->allocate_vector_of_vectors(2,3,dspine_S);

     /// Calculate the spines at position x using bulk element functions
     bulk_elem_pt->get_spine_base(interpolated_x,spine_B,dspine_B);
     bulk_elem_pt->get_spine(interpolated_x,spine_S,dspine_S);

     // Initialisation of vectors
     Vector<double> v_1(3,0.0), v_2(3,0.0);
     Vector<double> v_3(3,0.0), v_4(3,0.0);
     
     // Pseudo jacobian
     double J=interpolated_dxds[bulk_s_index(0)];

     // Calculation of dRds
     //--------------------

     // Calculation of v_1=J(u S,phi)
     bulk_elem_pt->scalar_times_vector(interpolated_u*J,
                                       dspine_S[bulk_s_index(0)],v_1);

     // Calculation of v_2=w,s S
     bulk_elem_pt->scalar_times_vector(interpolated_duds,spine_S,v_2);

     // Calculation of v_4=J B,phi
     bulk_elem_pt->scalar_times_vector(J,dspine_B[bulk_s_index(0)],v_4);

     // Get dRds as the sum of v_1, v_2 and v_4
     bulk_elem_pt->vector_sum(v_1,v_2,v_3);
     bulk_elem_pt->vector_sum(v_3,v_4,dRds);
     
     // Norm of the tangent
     nm=bulk_elem_pt->two_norm(dRds);

     //By construction 
     //!!! the spines, I think, Nboundary and S are in the same plane)!!!!
     //the component in the fixed direction of the dRdx is nil. 
     if ( bulk_s_index(0) == 1 )
      {
       N[0]=0.0;
       N[1]=-dRds[2]/nm;
       N[2]=dRds[1]/nm;
      }
     else
      {
       N[1]=0.0;
       N[0]=-dRds[2]/nm;
       N[2]=dRds[0]/nm;
      }

     //Now get the N-component of Spine S in the (dRdx,N) plane
     alpha_n= bulk_elem_pt->scalar_product(N,spine_S);

     // Initialisation of vectors
     Vector<double> n_wall(3,0.0);
     Vector<double> axe_ez(3,0.0);
     axe_ez[2]=1.0;

     ELEMENT::cross_product(dRds,N,n_wall);

     // Normalise
     double norm=ELEMENT::two_norm(n_wall);
     for (unsigned i=0;i<3;i++) n_wall[i]/=norm;

     //Now get the TxN-component of Spine S in the (dRdx,N) plane
     double alp = bulk_elem_pt->scalar_product(n_wall,spine_S);

     if ( alpha_n <0 ){ alpha_n=-1.0*alpha_n;}
       if ( alp <0 ) { alp=-1.0*alp;}
     //Now add to the appropriate equations
   
     //Loop over the test functions
     for(unsigned l=0;l<n_node;l++)
      {
       local_eqn = u_local_eqn(l);
       //IF it's not a boundary condition
       if(local_eqn >= 0)
        {
         //Add the prescribed flux terms
         residuals[local_eqn] -= (alpha_n*cos_gamma + alp*sqrt(1-cos_gamma*cos_gamma))*psif[l]*w*nm;
         
        }
      }
    } // End of if use spines 
   else  // cartesian case
     {

       // Get the jacobian of the mapping of the contact line
       double J=J_eulerian(s);
       
       // Check for the length of the element
       double length_of_face_element=0.0;
       
       //Loop over the test functions
       for(unsigned l=0;l<n_node;l++)
	 {
	   local_eqn = u_local_eqn(l);
	   //IF it's not a boundary condition
	   if(local_eqn >= 0)
	     {
	       //Add the prescribed flux terms
	       residuals[local_eqn] -= cos_gamma*psif[l]*w*J;
	       
	       // Imposed traction doesn't depend upon the solution, 
	       // --> the Jacobian is always zero, so no Jacobian
	       // terms are required
	       
	       length_of_face_element+=J*w;
	     }
	 }
       // std::cout << "longueur de l'élément :" << length_of_face_element << std::endl;
     } // End of cartesian case 
  }
}*/




//================================================================
/// Compute the element's contribution to the residual vector 
/// (new version)
//================================================================
template<class ELEMENT>
void YoungLaplaceFluxElement<ELEMENT>::
new_fill_in_contribution_to_residuals(Vector<double> &residuals)
{
 //Find out how many nodes there are
 unsigned n_node = nnode();

 //Set up memory for the shape functions
 Shape psi(n_node);
  
 //Number of integration points
 unsigned n_intpt = integral_pt()->nweight();
 
 //Set the Vector to hold local coordinates
 // Note: We need the coordinate itself below for the evaluation
 // of the contact line vectors even though we use the *at_knot
 // version for the various shape-function-related functions
 Vector<double> s(1); 
 
 //Integers to hold the local equation and unknown numbers
 int local_eqn=0;

 //Loop over the integration points
 for(unsigned ipt=0;ipt<n_intpt;ipt++)
  {

   //Assign value of s
   s[0] = integral_pt()->knot(ipt,0);
   
   //Get the integral weight
   double w = integral_pt()->weight(ipt); 

   //Find the shape functions
   shape_at_knot(ipt,psi); 

   //Get the prescribed cos_gamma
   double cos_gamma=prescribed_cos_gamma();

   // Get the various contact line vectors
   Vector<double> tangent(3);
   Vector<double> normal(3);
   Vector<double> spine(3);
   double norm_of_drds;
   contact_line_vectors(s,tangent,normal,spine,norm_of_drds);

   // Get beta factor:

   // Cross product of spine and tangent to contact line is
   // the wall normal
   Vector<double> wall_normal(3);
   ELEMENT::cross_product(spine,tangent,wall_normal);
 
   // Normalise
   double norm=ELEMENT::two_norm(wall_normal);
   for (unsigned i=0;i<3;i++) wall_normal[i]/=norm;
   
   // Take cross product with tangent to get the normal to the
   // contact line parallel to wall
   Vector<double> normal_to_contact_line_parallel_to_wall(3);
   ELEMENT::cross_product(tangent,wall_normal,
                          normal_to_contact_line_parallel_to_wall);

   // Normalise spine
   Vector<double> normalised_spine(3);
   double spine_norm=ELEMENT::two_norm(spine);
   for (unsigned i=0;i<3;i++) normalised_spine[i]=spine[i]/spine_norm;

   double beta=ELEMENT::scalar_product(normalised_spine,
                                   normal_to_contact_line_parallel_to_wall);

                                    
   //Now add to the appropriate equations
   
   //Loop over the test functions
   for(unsigned l=0;l<n_node;l++)
    {
     local_eqn = u_local_eqn(l);
     /*IF it's not a boundary condition*/
     if(local_eqn >= 0)
      {
       //Add to residual:
       residuals[local_eqn] -= beta*cos_gamma*psi[l]*norm_of_drds*w;
      }
    }
  }

//  //Find out how many nodes there are
//  unsigned n_node = nnode();

//  //Set up memory for the shape functions
//  Shape psif(n_node);
  
//  //Number of integration points
//  unsigned n_intpt = integral_pt()->nweight();
 
//  //Set the Vector to hold local coordinates
//  Vector<double> s(1);
 
//  //Integers to hold the local equation and unknown numbers
//  int local_eqn=0;

//  ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>(this->bulk_element_pt());

//  //Loop over the integration points
//  for(unsigned ipt=0;ipt<n_intpt;ipt++)
//   {

//    //Assign value of s
//    s[0] = integral_pt()->knot(ipt,0);
   
//    //Get the integral weight
//    double w = integral_pt()->weight(ipt);

//    //Find the shape functions
//    shape(s,psif);

//    //Get the prescribed cos_gamma
//    double cos_gamma=prescribed_cos_gamma();


//    if ( bulk_elem_pt->use_spines() )
//      {
//        // Get the various contact line vectors
//        Vector<double> tangent(3);
//        Vector<double> normal(3);
//        Vector<double> spine(3);
//        double norm_of_drds;
//        contact_line_vectors(s,tangent,normal,spine,norm_of_drds);
       
//        // Get beta and lambda factors:
       
//        // Cross product of spine and tangent to contact line is
//        // the wall normal
//        Vector<double> wall_normal(3);
//        Vector<double> axe_ez(3,0.0);
//        axe_ez[2]=1.0;

//        // ELEMENT::cross_product(spine,tangent,wall_normal);
//        ELEMENT::cross_product(axe_ez,tangent,wall_normal);

//        // Normalise
//        double norm=ELEMENT::two_norm(wall_normal);
//        for (unsigned i=0;i<3;i++) wall_normal[i]/=norm;
       
//        // Take cross product with tangent to get the normal to the
//        // contact line parallel to wall
//        Vector<double> normal_to_contact_line_parallel_to_wall(3);
//        ELEMENT::cross_product(tangent,wall_normal,
// 			      normal_to_contact_line_parallel_to_wall);
       
//        // Normalise spine
//        Vector<double> normalised_spine(3);
//        double spine_norm=ELEMENT::two_norm(spine);
//        for (unsigned i=0;i<3;i++) normalised_spine[i]=spine[i]/norm;
       
//        double beta=ELEMENT::scalar_product(normalised_spine,
//       					   normal_to_contact_line_parallel_to_wall);
//        double lambda=ELEMENT::scalar_product(normalised_spine,
//       					   wall_normal);   
       
//        // Temporary check that the new axis are a regular system of axis
//        if ( beta<0 ){ beta=-1.0*beta;}
//        if ( lambda<0 ) { lambda=-1.0*lambda;}

//        //Now add to the appropriate equations
   
//        //Loop over the test functions
//        for(unsigned l=0;l<n_node;l++)
// 	 {
// 	   local_eqn = u_local_eqn(l);
// 	   /*IF it's not a boundary condition*/
// 	   if(local_eqn >= 0)
// 	     {
// 	       //Add to residual:
// 	       residuals[local_eqn] -= ( lambda*sqrt(1-cos_gamma*cos_gamma)
// 					+ beta*cos_gamma                  )
// 		                               *psif[l]*norm_of_drds*w;
// 	     }
// 	 }
//      }
//    else
//      {

//        // Get the jacobian of the mapping of the contact line
//        double J=J_eulerian(s);
       
//        //Loop over the test functions
//        for(unsigned l=0;l<n_node;l++)
// 	 {
// 	   local_eqn = u_local_eqn(l);
// 	   /*IF it's not a boundary condition*/
// 	   if(local_eqn >= 0)
// 	     {
// 	       //Add the prescribed flux terms
// 	       residuals[local_eqn] -= cos_gamma*psif[l]*w*J;
	       
// 	       // Imposed traction doesn't depend upon the solution, 
// 	       // --> the Jacobian is always zero, so no Jacobian
// 	       // terms are required
	       
// 	     }
// 	 }
//      }

//   }


}




//========================================================================
/// Get the actual contact angle
//========================================================================
template<class ELEMENT>
double YoungLaplaceFluxElement<ELEMENT>::actual_cos_contact_angle(
 const Vector<double>& s)
{  
  // Get pointer to bulk element
  ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>(this->bulk_element_pt());

  double cos_gamma=0.0;
  
  // Spine case
  if ( bulk_elem_pt->use_spines() )
    {
      // Get various contact line vectors
      Vector<double> tangent(3);
      Vector<double> normal(3);
      Vector<double> spine(3);
      double norm_of_drds;
      contact_line_vectors(s,tangent,normal,spine,norm_of_drds);
      
      // Get the wall normal: Both the tangent to the contact line
      // and the spine vector are tangential to the wall:
      Vector<double> wall_normal(3); 
      Vector<double> axe_ez(3,0.0);
      axe_ez[2]=1.0;
      //        ELEMENT::cross_product(spine,tangent,wall_normal);
      ELEMENT::cross_product(axe_ez,tangent,wall_normal);
      
      // Normalise
      double norm=0.0;
      for (unsigned i=0;i<3;i++) norm+=wall_normal[i]*wall_normal[i];
      for (unsigned i=0;i<3;i++)
	{ 
	  wall_normal[i]/=sqrt(norm);
	}
      
      // Get the auxiliary unit vector that's normal to the contact line and tangent
      // to the wall
      Vector<double> aux(3);
      ELEMENT::cross_product(tangent,wall_normal,aux);
      
      // Cos of contact angle is dot product with wall normal
      cos_gamma=ELEMENT::scalar_product(aux,normal);

    }

  // Cartesian case
  else
    {
     // Which local coordinate is const in bulk element?
     //unsigned s_fixed_index_in_bulk=s_fixed_index();
     
     // Assign the other value (we'll only ever be in 2D here!)
     //unsigned s_variable_index_in_bulk=0;
     //if (0==s_fixed_index_in_bulk) s_variable_index_in_bulk=1;
     
     // Get local coordinates in bulk element by copy construction
      Vector<double> s_bulk(local_coordinate_in_bulk(s));
      
      // Number of nodes in bulk element
      unsigned nnode_bulk=bulk_elem_pt->nnode();

      // Dimension of (= number of local coordinates in) bulk element 
      //unsigned dim_bulk=bulk_elem_pt->dim();
      
// #ifdef PARANOID
//       if (dim_bulk!=2) 
//        {
//         throw OomphLibError(
//          "YoungLaplaceFluxElements only work with 2D bulk elements",
//          "YoungLaplaceFluxElement::actual_cos_contact_angle",
//          OOMPH_EXCEPTION_LOCATION);
//        }
// #endif

      //Set up memory for the shape functions
      Shape psi(nnode_bulk);
      DShape dpsidzeta(nnode_bulk,2);

      //Call the derivatives of the shape and test functions
      //in the bulk -- must do this via s because this point
      //is not an integration point the bulk element!
      (void)bulk_elem_pt->dshape_eulerian(s_bulk,psi,dpsidzeta);

      // Get the gradient at s
      Vector<double> gradient_u(2,0.0);
 
       //Calculate function value and derivatives:
      //-----------------------------------------
      // Loop over nodes
      for(unsigned l=0;l<nnode_bulk;l++) 
	{
	  // Loop over directions
	  for(unsigned j=0;j<2;j++)
	    { gradient_u[j]+= bulk_elem_pt->u(l)*dpsidzeta(l,j); }
	}

      // Get the outer unit normal to boundary
      Vector<double> outer_normal(2,0.0);
      outer_unit_normal(s,outer_normal);

      // Compute the cosinus of the angle
      double gradient_norm_2=ELEMENT::two_norm(gradient_u)*ELEMENT::two_norm(gradient_u);
      cos_gamma=ELEMENT::scalar_product(gradient_u,outer_normal)/sqrt(1+gradient_norm_2);

    }                                                     
    
 return cos_gamma;

}


//====================================================
/// Get the length of the element
//====================================================
template<class ELEMENT>
double YoungLaplaceFluxElement<ELEMENT>::get_projected_length()
{
  //Length of the face element
  double projected_length=0.0;
  double length=0.0;

  ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>(this->bulk_element_pt());

  //Find out how many nodes there are
  unsigned n_node = nnode();
  
  //Set up memory for the shape functions
  Shape psif(n_node);
  
  //Number of integration points
  unsigned n_intpt = integral_pt()->nweight();
      
  //Set the Vector to hold local coordinates
  Vector<double> s(1);
  
  //Loop over the integration points
  for(unsigned ipt=0;ipt<n_intpt;ipt++)
    {
      
      //Assign value of s
      s[0] = integral_pt()->knot(ipt,0);
      
      //Get the integral weight
      double w = integral_pt()->weight(ipt);
      
      // Get various contact line vectors
      Vector<double> tangent(3);
      Vector<double> normal(3);
   
      double norm_of_drds=0.0;
      
      // Spine case
      if ( bulk_elem_pt->use_spines() )
	{
	  Vector<double> spine(3);
	  contact_line_vectors(s,tangent,normal,spine,norm_of_drds);
	}
      // Cartesian case
      else
	{
	  contact_line_vectors(s,tangent,normal,norm_of_drds);
	}

      // Increment projected length in the x-y plane
      projected_length+=w*norm_of_drds*sqrt(1.0-tangent[2]*tangent[2]);
      length+=w*norm_of_drds;
	  
    }
 
  return projected_length;

}


//========================================================================
/// Get tangent and normal to contact line and the spine itself (this allows
/// the wall normal to be worked out by a cross product). Final
/// argument gives the norm of dR/ds where R is the vector to the
/// contact line and s the local coordinate in the element
//========================================================================
template<class ELEMENT>
void YoungLaplaceFluxElement<ELEMENT>::contact_line_vectors(
 const Vector<double>& s, 
 Vector<double>& tangent, 
 Vector<double>& normal,
 Vector<double>& spine,
 double& norm_of_drds)
{

 // Get pointer to bulk element
// ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>(this->bulk_element_pt());
                   
 // Dimension of (= number of local coordinates in) bulk element 
 //unsigned dim_bulk=bulk_elem_pt->dim();

// #ifdef PARANOID
//  if (dim_bulk!=2) 
//   {
//    throw OomphLibError(
//     "YoungLaplaceFluxElements only work with 2D bulk elements",
//     "YoungLaplaceFluxElement::contact_line_vectors",
//          OOMPH_EXCEPTION_LOCATION);
//   }
// #endif

 throw OomphLibError(
  "This is not consistent with the new FaceElement interfaces",
  "YoungLaplaceFluxElement::contact_line_vectors",
  OOMPH_EXCEPTION_LOCATION);
 
 /*
 // Which local coordinate is const in bulk element?
 unsigned s_fixed_index_in_bulk=s_fixed_index();

 // Assign the other value (we'll only ever be in 2D here!)
 unsigned s_variable_index_in_bulk=0;
 if (0==s_fixed_index_in_bulk) s_variable_index_in_bulk=1;

 // Get local coordinates in bulk element by copy construction
 Vector<double> s_bulk(local_coordinate_in_bulk(s));
 
 // Number of nodes in bulk element
 unsigned nnode_bulk=bulk_elem_pt->nnode();

 // Get the shape functions and their derivatives w.r.t. 
 // to the local coordinates in the bulk element
 Shape psi_bulk(nnode_bulk);
 DShape dpsi_bulk(nnode_bulk,dim_bulk);
 bulk_elem_pt->dshape_local(s_bulk,psi_bulk,dpsi_bulk);

 // Displacement along spine
 double interpolated_u=0.0;

 // Derivative of u w.r.t. bulk local coordinate that varies along 
 // the face
 double interpolated_duds_variable=0.0;
 double interpolated_duds_fixed=0.0;

 // Global intrinsic coordinates of current point for evaluation of
 // spine and spine base
 Vector<double> interpolated_zeta(dim_bulk,0.0);

 // Derivative of global intrinsic coordinates with respect to
 // bulk local coordinates
 Vector<Vector<double> > interpolated_dzetads;
 ELEMENT::allocate_vector_of_vectors(dim_bulk,dim_bulk,interpolated_dzetads);

 // Loop over nodes                                                           
 for(unsigned l=0;l<nnode_bulk;l++)
  {                    
   interpolated_u+=bulk_elem_pt->u(l)*psi_bulk(l);
                                             
   interpolated_duds_variable += 
    bulk_elem_pt->u(l)*dpsi_bulk(l,s_variable_index_in_bulk);

   interpolated_duds_fixed += 
    bulk_elem_pt->u(l)*dpsi_bulk(l,s_fixed_index_in_bulk);

   // Loop over directions                                                    
   for(unsigned j=0;j<dim_bulk;j++)
    {                    
     interpolated_zeta[j]+=      
      bulk_elem_pt->nodal_position(l,j)*psi_bulk(l);
     // Loop over derivatives w.r.t. to bulk coordinates
     for (unsigned i=0;i<dim_bulk;i++)
      {
       interpolated_dzetads[j][i] +=
        bulk_elem_pt->nodal_position(l,j)*dpsi_bulk(l,i);
      }
    }                                                                         
  }          
             
 // Auxiliary vector (tangent to non-fixed bulk coordinate but
 // not necessarily normal to contact line)
 Vector<double> aux_vector(3);
 double tang_norm=0.0;
 double aux_norm=0.0;
                                                    
 if ( bulk_elem_pt->use_spines() )
   {
     // Get the spine and spine base vector at this point from the bulk element
     Vector<double> spine_base(3,0.0);
     Vector< Vector<double> > dspine;
     ELEMENT::allocate_vector_of_vectors(2,3,dspine);
     Vector< Vector<double> > dspine_base;
     ELEMENT::allocate_vector_of_vectors(2,3,dspine_base);
     bulk_elem_pt->get_spine(interpolated_zeta,spine,dspine);
     bulk_elem_pt->get_spine_base(interpolated_zeta,spine_base,dspine_base);

     // Derivative of spine and spine base w.r.t. local coordinate in 
     // FaceElement:
     Vector<double> dspine_ds_variable(3,0.0);
     Vector<double> dspine_base_ds_variable(3,0.0);
     Vector<double> dspine_ds_fixed(3,0.0);
     Vector<double> dspine_base_ds_fixed(3,0.0);
     for (unsigned i=0;i<3;i++)
       {
	 
	 dspine_ds_variable[i]+=
	   dspine[0][i]*interpolated_dzetads[0][s_variable_index_in_bulk]+
	   dspine[1][i]*interpolated_dzetads[1][s_variable_index_in_bulk];
	 
	 dspine_base_ds_variable[i]+=
	   dspine_base[0][i]*interpolated_dzetads[0][s_variable_index_in_bulk]+
	   dspine_base[1][i]*interpolated_dzetads[1][s_variable_index_in_bulk];
	 
	 
	 dspine_ds_fixed[i]+=
	   dspine[0][i]*interpolated_dzetads[0][s_fixed_index_in_bulk]+
	   dspine[1][i]*interpolated_dzetads[1][s_fixed_index_in_bulk];
	 
	 dspine_base_ds_fixed[i]+=
	   dspine_base[0][i]*interpolated_dzetads[0][s_fixed_index_in_bulk]+
	   dspine_base[1][i]*interpolated_dzetads[1][s_fixed_index_in_bulk];
	 
       }
     
     // Auxiliary vector (tangent to non-fixed bulk coordinate but
     // not necessarily normal to contact line)
     for (unsigned i=0;i<3;i++)
       {
	 tangent[i]=
	   dspine_base_ds_variable[i]+
	   interpolated_duds_variable*spine[i]+
	   interpolated_u*dspine_ds_variable[i];
	 tang_norm+=tangent[i]*tangent[i];
	 
	 
	 aux_vector[i]=
	   dspine_base_ds_fixed[i]+
	   interpolated_duds_fixed*spine[i]+
	   interpolated_u*dspine_ds_fixed[i];
	 aux_norm+=aux_vector[i]*aux_vector[i];
       }
     
   }

 // Cartesian case
 else
   {
     for(unsigned i=0;i<2;i++)
       {
	 tangent[i]= interpolated_dzetads[s_variable_index_in_bulk][i];
	 tang_norm+=tangent[i]*tangent[i];
	 aux_vector[i]= interpolated_dzetads[s_fixed_index_in_bulk][i];
	 aux_norm+=aux_vector[i]*aux_vector[i];
	 
       }

     tangent[2]=interpolated_duds_variable; tang_norm+=tangent[2]*tangent[2];
     aux_vector[2]=interpolated_duds_fixed; aux_norm+=aux_vector[2]*aux_vector[2];
     
   }
    
 // Norm of the rate of change
 norm_of_drds=sqrt(tang_norm);
 
 // Normalise
 double tang_norm_fact=1.0/sqrt(tang_norm);
 double aux_norm_fact=1.0/sqrt(aux_norm);
 for (unsigned i=0;i<3;i++)
   {
     tangent[i]*=tang_norm_fact;
     aux_vector[i]*=aux_norm_fact;
   }
 
 
 // Normal to meniscus is the cross product between the
 // two contact line vectors:
 Vector<double> meniscus_normal(3);
 ELEMENT::cross_product(tangent,aux_vector,meniscus_normal);
 
 // Calculate the norm
 double men_norm_fact=0.0;
 for (unsigned i=0;i<3;i++)
   {
     men_norm_fact+=meniscus_normal[i]*meniscus_normal[i];
   }
 
 // Normalise and adjust direction
 double sign=-double(normal_sign());
 for (unsigned i=0;i<3;i++)
   {
//     meniscus_normal[i]/=sqrt(men_norm_fact);
     meniscus_normal[i]*=sign/sqrt(men_norm_fact);
   }
 
 // The actual (bi) normal to the contact line is given
 // by another cross product
 ELEMENT::cross_product(meniscus_normal,tangent,normal);
 */ 
}


//============================================================
// Build the required elements
//============================================================
template class YoungLaplaceFluxElement<QYoungLaplaceElement<2> >;
template class YoungLaplaceFluxElement<QYoungLaplaceElement<3> >;
template class YoungLaplaceFluxElement<QYoungLaplaceElement<4> >;

template class YoungLaplaceFluxElement<RefineableQYoungLaplaceElement<2> >;
template class YoungLaplaceFluxElement<RefineableQYoungLaplaceElement<3> >;
template class YoungLaplaceFluxElement<RefineableQYoungLaplaceElement<4> >;

}

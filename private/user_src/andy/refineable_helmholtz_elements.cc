#include "refineable_helmholtz_elements.h"


namespace oomph
{


//========================================================================
/// Add element's contribution to the elemental 
/// residual vector and/or Jacobian matrix.
/// flag=1: compute both
/// flag=0: compute only residual vector
//========================================================================
template<unsigned DIM>
void RefineableHelmholtzEquations<DIM>::
fill_in_generic_residual_contribution_helmholtz(Vector<double> &residuals, 
                                              DenseMatrix<double> &jacobian, 
                                              unsigned flag)
{
//Find out how many nodes there are in the element
 unsigned n_node = nnode();

//Set up memory for the shape and test functions
 Shape psi(n_node), test(n_node);
 DShape dpsidx(n_node,DIM), dtestdx(n_node,DIM);

//Set the value of n_intpt
 unsigned n_intpt = integral_pt()->nweight();

//The local index at which the helmholtz variable is stored
 unsigned u_nodal_index = this->u_index_helmholtz(); 

//Integers to store the local equation and unknown numbers
 int local_eqn=0, local_unknown=0;

// Local storage for pointers to hang_info objects
 HangInfo *hang_info_pt=0, *hang_info2_pt=0;

//Loop over the integration points
 for(unsigned ipt=0;ipt<n_intpt;ipt++)
  {
   //Get the integral weight
   double w = integral_pt()->weight(ipt);
 
   //Call the derivatives of the shape and test functions
   double J = 
    this->dshape_and_dtest_eulerian_at_knot_helmholtz(ipt,psi,dpsidx,test,dtestdx);
 
   //Premultiply the weights and the Jacobian
   double W = w*J;
 
   //Calculate local values of the function
//   double interpolated_u=0.0;
 
   //This needs to be a Vector to be ANSI C++ (Initialise to zero) 
   Vector<double> interpolated_x(DIM,0.0);
   Vector<double> interpolated_dudx(DIM,0.0);
 
   //Calculate function value and derivatives:
   //-----------------------------------------
 
   // Loop over nodes
   for(unsigned l=0;l<n_node;l++) 
    {
     //Get the helmholtz value from the node 
     //(hanging-ness will be taken into account)
     double u_value = this->nodal_value(l,u_nodal_index);

     // Loop over directions
     for(unsigned j=0;j<DIM;j++)
      {
       interpolated_x[j] += nodal_position(l,j)*psi(l);
       interpolated_dudx[j] += u_value*dpsidx(l,j);
      }
    }
 
   //Get body force
   double source;
   this->get_source_helmholtz(interpolated_x,source);
 
   // Access the local coordinates for this integration point
   // (these will have been assigned earlier when the other element
   //  was assigned (currently in the driver two_mesh_helmholtz.cc))
   // Note: interaction parameter default to zero for now
   Vector<double> s_other(DIM);
   s_other=this->external_element_local_coord(0,ipt);

   // Add contribution to source from "other" element
   // source+=this->interpolated_u_helmholtz(s);
   HelmholtzEquations<DIM>* source_el_pt=dynamic_cast<HelmholtzEquations<DIM>*>
    (this->external_element_pt(0,ipt));
   source+=source_el_pt->interpolated_u_helmholtz(s_other);

   // Assemble residuals and Jacobian
 
   // Loop over the nodes for the test functions 
   for(unsigned l=0;l<n_node;l++)
    {
     //Local variables used to store the number of master nodes and the
     //weight associated with the shape function if the node is hanging
     unsigned n_master=1; double hang_weight=1.0;
     //Local bool (is the node hanging)
     bool is_node_hanging = this->node_pt(l)->is_hanging();

     //If the node is hanging, get the number of master nodes
     if(is_node_hanging)
      {
       hang_info_pt = this->node_pt(l)->hanging_pt();
       n_master = hang_info_pt->nmaster();
      }
     //Otherwise there is just one master node, the node itself
     else
      {
       n_master = 1;
      }
   
     //Loop over the master nodes
     for(unsigned m=0;m<n_master;m++)
      {
       //Get the local equation number and hang_weight
       //If the node is hanging
       if(is_node_hanging)
        {
         //Read out the local equation number from the m-th master node
         local_eqn =  this->local_hang_eqn(hang_info_pt->master_node_pt(m),
                                           u_nodal_index);
         //Read out the weight from the master node
         hang_weight = hang_info_pt->master_weight(m);
        }
       //If the node is not hanging
       else
        {
         //The local equation number comes from the node itself
         local_eqn = this->nodal_local_eqn(l,u_nodal_index);
         //The hang weight is one
         hang_weight = 1.0;
        }
     
       //If the nodal equation is not a boundary condition
       if(local_eqn >= 0)
        {
//         oomph_info << "source=" << source << ", hang_weight="
//                    << hang_weight << std::endl;

         // Add body force/source term here - the Helmholtz bit 
         residuals[local_eqn] += source*test(l)*W*hang_weight;
         // (Has the hangingness already been taken into account??)
//         residuals[local_eqn] += source*test(l)*W;

       
         // The operator bit [again, has hangingness already been done?]
         for(unsigned k=0;k<DIM;k++)
          {
           residuals[local_eqn] += 
            interpolated_dudx[k]*dtestdx(l,k)*W*hang_weight;
          }
       
         // Calculate the Jacobian
         if(flag)
          {
           //Local variables to store the number of master nodes
           //and the weights associated with each hanging node
           unsigned n_master2=1; double hang_weight2=1.0;
           //Loop over the nodes for the variables
           for(unsigned l2=0;l2<n_node;l2++)
            { 
             //Local bool (is the node hanging)
             bool is_node2_hanging = this->node_pt(l2)->is_hanging();
             //If the node is hanging, get the number of master nodes
             if(is_node2_hanging)
              {
               hang_info2_pt = this->node_pt(l2)->hanging_pt();
               n_master2 = hang_info2_pt->nmaster();
              }
             //Otherwise there is one master node, the node itself
             else
              {
               n_master2 = 1;
              }
           
             //Loop over the master nodes
             for(unsigned m2=0;m2<n_master2;m2++)
              {
               //Get the local unknown and weight
               //If the node is hanging
               if(is_node2_hanging)
                {
                 //Read out the local unknown from the master node
                 local_unknown = 
                  this->local_hang_eqn(hang_info2_pt->master_node_pt(m2),
                                       u_nodal_index);
                 //Read out the hanging weight from the master node
                 hang_weight2 = hang_info2_pt->master_weight(m2);
                }
               //If the node is not hanging
               else
                {
                 //The local unknown number comes from the node
                 local_unknown = this->nodal_local_eqn(l2,u_nodal_index);
                 //The hang weight is one
                 hang_weight2 = 1.0;
                }

               //If the unknown is not pinned
               if(local_unknown >= 0)
                {
//                 oomph_info << "add to jacobian" << std::endl;
                 //Add contribution to Elemental Matrix
                 for(unsigned i=0;i<DIM;i++)
                  {
                   jacobian(local_eqn,local_unknown) += 
                    (
                     dpsidx(l2,i)*dtestdx(l,i)
//                    missing the source term here?
                     )*W*hang_weight*hang_weight2;
                  }
                }
              } //End of loop over master nodes
            } //End of loop over nodes
          } //End of Jacobian calculation
       
        } //End of case when residual equation is not pinned
      } //End of loop over master nodes for residual
    } //End of loop over nodes
 
  } // End of loop over integration points
}


//====================================================================
// Force build of templates
//====================================================================
template class RefineableQHelmholtzElement<2,2>;
template class RefineableQHelmholtzElement<2,3>;
template class RefineableQHelmholtzElement<2,4>;

template class RefineableQHelmholtzElement<3,2>;
template class RefineableQHelmholtzElement<3,3>;
template class RefineableQHelmholtzElement<3,4>;

}

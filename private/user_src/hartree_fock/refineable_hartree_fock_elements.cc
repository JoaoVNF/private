#include "refineable_hartree_fock_elements.h"


namespace oomph
{


//========================================================================
/// Add element's contribution to the elemental  hierher
//========================================================================
template<unsigned DIM>
void RefineableHartreeFockEquations<DIM>::
fill_in_contribution_to_jacobian_and_mass_matrix(Vector<double> &residuals, 
                                                 DenseMatrix<double> &jacobian, 
                                                 DenseMatrix<double> &mass_matrix)
{
 /// Retrieve the index of this element
 unsigned e=this->get_element_number();

 /// Find out how many nodes there are in the element
 unsigned n_node = nnode();

 /// Set up memory for the shape and test functions
 Shape psi(n_node), test(n_node);
 DShape dpsidx(n_node,DIM), dtestdx(n_node,DIM);

 /// Set the value of n_intpt
 unsigned n_intpt = integral_pt()->nweight();

 /// The local index at which the hartree_fock variable is stored
 unsigned u_nodal_index = this->u_index_hartree_fock(); 

 ///Integers to store the local equation and unknown numbers
 int local_eqn=0, local_unknown=0, global_unknown=0; //ashflag

 /// Local storage for pointers to hang_info objects
 HangInfo *hang_info_pt=0, *hang_info2_pt=0;

 ///Loop over the integration points
 for(unsigned ept=0;ept<n_intpt;ept++)
  {
   /// Get the integral weight
   double w = integral_pt()->weight(ept);
 
   /// Call the derivatives of the shape and test functions
   double J = 
    this->dshape_and_dtest_eulerian_at_knot_hartree_fock(ept,psi,dpsidx,test,dtestdx);
   
   /// Premultiply the weights and the Jacobian
   double W = w*J;
 
   /// Calculate local values of the function
   /// double interpolated_u=0.0;
 
   /// This needs to be a Vector to be ANSI C++ (Initialise to zero) 
   Vector<double> interpolated_x(DIM,0.0);
   /// Vector<double> interpolated_dudx(DIM,0.0);
 
   //Calculate function value and derivatives:
   //-----------------------------------------
 
   /// Loop over nodes
   for(unsigned l=0;l<n_node;l++) 
    {
     //Get the hartree_fock value from the node 
     //(hanging-ness will be taken into account
     //double u_value = this->nodal_value(l,u_nodal_index);
     //The hartree_fock value is stored at the nodes
     //interpolated_u += u_value*psi(l);
     // Loop over directions
     for(unsigned j=0;j<DIM;j++)
      {
       interpolated_x[j] += nodal_position(l,j)*psi(l);
       //interpolated_dudx[j] += u_value*dpsidx(l,j);
      }
    }
   
   /// Get nuclear source
   double nuclear_source=0.0;
   this->get_nuclear_source(interpolated_x,nuclear_source);
   
   /// Get electron repulsion source
   double repulsion_source = 0.0;
   this->get_electron_repulsion(e,ept,repulsion_source);
   
   
   /// Assemble residuals and Jacobian
   
   /// Loop over the nodes for the test functions 
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
//        // Add body force/source term here 
//        residuals[local_eqn] += source*test(l)*W*hang_weight;
         
//        // The HartreeFock bit itself
//        for(unsigned k=0;k<DIM;k++)
//         {
//          residuals[local_eqn] += 
//           interpolated_dudx[k]*dtestdx(l,k)*W*hang_weight;
//         }
         
         // Calculate the Jacobian
         //if(flag)
         {
          //Local variables to store the number of master nodes
          //and the weights associated with each hanging node
          unsigned n_master2=1; double hang_weight2=1.0;
          
          
          //Loop over the nodes for the variables
          for(unsigned l2=0;l2<n_node;l2++)
           {           
            //ashflag
            /// Get global equation from node
            global_unknown = this->node_pt(l2)->eqn_number(u_nodal_index);

            /// Get electron exchange source
            double exchange_source=0.0;
            this->get_electron_exchange(e,ept,global_unknown,exchange_source);
            
            
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
                //Add contribution to elemental mass matrix
                mass_matrix(local_eqn,local_unknown)
                 += psi(l)*psi(l2)*W*hang_weight*hang_weight2;
                
                //Add contribution to stiffness matrix from nuclear potential
                jacobian(local_eqn,local_unknown)
                 -= W*psi(l)*nuclear_source*psi(l2)*hang_weight*hang_weight2;
                
                //Add contribution to stiffness matrix from electron repulsion
                jacobian(local_eqn,local_unknown)
                 += W*psi(l)*repulsion_source*psi(l2)*hang_weight*hang_weight2;
                
                //Add contribution to stiffness matrix from electron exchange
                jacobian(local_eqn,local_unknown)
                 -= W*psi(l)*exchange_source*hang_weight*hang_weight2;
                
                for(unsigned i=0;i<DIM;i++)
                 {
                  jacobian(local_eqn,local_unknown)
                   += 0.5*dpsidx(l2,i)*dtestdx(l,i)*W*hang_weight*hang_weight2;
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
template class RefineableQHartreeFockElement<2,2>;
template class RefineableQHartreeFockElement<2,3>;
template class RefineableQHartreeFockElement<2,4>;

template class RefineableQHartreeFockElement<3,2>;
template class RefineableQHartreeFockElement<3,3>;
template class RefineableQHartreeFockElement<3,4>;

}

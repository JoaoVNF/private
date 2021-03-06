using namespace std;
using namespace oomph;

namespace oomph
{
 
 
//=======================================================================
/// 1D Free surface elements for Hele Shaw problems with pseudo-elastic
/// mesh motion
//======================================================================
 template <class ELEMENT>
  class HeleShawInterfaceElement :  public virtual FaceGeometry<ELEMENT>,
  public virtual FaceElement
  {
    public:
   
   /// \short Function pointer to function which provides bubble pressure as a
   /// function of vector x. This should usually be a constant - should default
   /// to zero. Any perturbations to the system can be supplied through a 
   /// spatially varying bubble pressure function.
   typedef void (*BubblePressureFctPt)(const Vector<double>& x, 
                                       double& p_bubble);
   
   /// \short Function pointer to function which provides wall speed as a 
   /// function of x. This allows us to solve for bubble motion in a moving 
   /// frame. A constant wall speed does not affect the mass conservation 
   /// equations, but does feature in the kinematic equation for interface 
   /// motion.
   typedef void (*WallSpeedFctPt)(const Vector<double>& x, 
                                  Vector<double>& U_wall);
   
   /// Access function: Pointer to source function
   BubblePressureFctPt& bubble_pressure_fct_pt() 
    {
     return Bubble_pressure_fct_pt;
    }
   
   /// Access function: Pointer to source function. Const version
   BubblePressureFctPt bubble_pressure_fct_pt() const 
   {
    return Bubble_pressure_fct_pt;
   }
   
   WallSpeedFctPt& wall_speed_fct_pt() {return Wall_speed_fct_pt;}
   
   WallSpeedFctPt wall_speed_fct_pt() const {return Wall_speed_fct_pt;}
   
    private:
   
   /// Pointer to function that specifies the gap width and wall velocity
   BubblePressureFctPt Bubble_pressure_fct_pt;
   
   WallSpeedFctPt Wall_speed_fct_pt;
   
   /// \short Pointer to the aspect_ratio ratio: reference gap width / 
   /// in-plane lengthscale
   double* Aspect_ratio_pt;
   
   /// Pointer to the inverse Capillary number
   double *Ca_inv_pt;
   
   /// Pointer to the Strouhal number
   double *St_pt;

   /// Default value for physical constants
   static double Default_Physical_Constant_Value;

   /// \short ID of additional unknowns introduced by this face element
   /// (smoothed components of derivative of tangent vector, and
   /// Lagrange multiplier)
   unsigned Id;

   /// Index of the nodal value at which the pressure is stored
   unsigned P_index_interface;

   /// Equation number of the equation for the Lagrange multiplier
   int lagrange_local_eqn(const unsigned &j)
   {
    // Get the index of the nodal value associated with Lagrange multiplier
    // NOTE: It's the first new one
    unsigned lagr_index=dynamic_cast<BoundaryNodeBase*>(node_pt(j))->
     index_of_first_value_assigned_by_face_element(Id)+0;
    
    // Return nodal value
    return this->nodal_local_eqn(j,lagr_index);
   }
   
   /// Return the Lagrange multiplier at local node j
   double &lagrange(const unsigned &j)
    {
     // Get the index of the nodal value associated with Lagrange multiplier
     // NOTE: It's the first new one
     unsigned lagr_index=dynamic_cast<BoundaryNodeBase*>(node_pt(j))->
      index_of_first_value_assigned_by_face_element(Id)+0;
     
     // Return (ref) to value
     return *node_pt(j)->value_pt(lagr_index);
    }
   
   /// \short Equation number of equation that does the projection
   /// for the derivative of the i-th component of the tangent vector at node j
   int projected_tangent_deriv_local_eqn(const unsigned& j,
                                         const unsigned &i)
   {
    // Get the index of the nodal value associated with projected tang vector
    unsigned tang_index=dynamic_cast<BoundaryNodeBase*>(node_pt(j))->
     index_of_first_value_assigned_by_face_element(Id)+1+i;
    
    // Return nodal value
    return this->nodal_local_eqn(j,tang_index);
   }
      
   /// \short Return i-th component of projected deriv of tangent vector
   /// at local node j
   double &projected_tangent_deriv(const unsigned& j, const unsigned &i)
    {
     // Get the index of the nodal value associated with tang deriv
     unsigned veloc_index=dynamic_cast<BoundaryNodeBase*>(node_pt(j))->
      index_of_first_value_assigned_by_face_element(Id)+1+i;
     
     // Return ref to value
     return *node_pt(j)->value_pt(veloc_index);
    }
   
   /// \short Helper function to calculate the residuals and
   /// (if flag==1) the Jacobian of the equations.
   void fill_in_generic_residual_contribution_hele_shaw_interface(
    Vector<double> &residuals,
    DenseMatrix<double> &jacobian,
    unsigned flag);
   
   
   
    public:
   
   /// \short Constructor, pass a pointer to the bulk element and the face 
   /// index of the bulk element to which the element is to be attached to.
   /// The optional identifier can be used to distinguish the additional nodal 
   /// values created by this element (in order: 
   /// 0: Lagrange multiplier; 
   /// 1: x-component of smoothed derivative of tangent vector; 
   /// 2: y-component of smoothed derivative of tangent vector)
   /// from those created by other FaceElements.
   HeleShawInterfaceElement(FiniteElement* const &element_pt,
                            const int &face_index,
                            const unsigned &id=0) :
    
    FaceGeometry<ELEMENT>(), Id(id)
    {
     
     // Attach the geometrical information to the element
     // This function also assigned nbulk_value from required_nvalue of the
     // bulk element
     element_pt->build_face_element(face_index,this);
     
     Wall_speed_fct_pt = 0;
     Bubble_pressure_fct_pt = 0;
          
     // Set the Strouhal number to the default value
     St_pt = &Default_Physical_Constant_Value;
     
     // Find the index at which the pressure stored
     // from the bulk element
     ELEMENT* bulk_element_pt = dynamic_cast<ELEMENT*>(element_pt);
     this->P_index_interface = bulk_element_pt->p_index_hele_shaw();
     
     // Read out the number of nodes on the face
     unsigned n_node_face = this->nnode();
     
     // Set the additional data values in the face
     // There are three additional values at each node -- the
     // Lagrange multiplier and the two components of the smoothed
     // derivative of the tangent vector
     Vector<unsigned> additional_data_values(n_node_face);
     for(unsigned i=0; i<n_node_face; i++)
      {
       additional_data_values[i] = 3;
      }
     
     // Now add storage and set the map containing
     // the position of the first entry of this face element's
     // additional values.
     add_additional_values(additional_data_values,id);
    }
    
    // rmodif
    /// Return the contribution to the integral of the square of the deviation
    /// from the mean radius
    void compute_deviation_from_mean_radius(
     const double& mean_radius, double& deviation_contribution)
    {
     
     // Initialise
     deviation_contribution = 0.0;
     
     // Find out how many nodes there are
     unsigned n_node = this->nnode();
     
     // Set up memory for the shape functions
     Shape psif(n_node);
     DShape dpsifds(n_node,1);
     
     // Set the value of n_intpt
     unsigned n_intpt = this->integral_pt()->nweight();
     
     // Storage for the local coordinate
     Vector<double> s(1);
     
     // Storage for the global coordinate
     Vector<double> x(2);
     
     // Loop over the integration points
     for(unsigned ipt=0; ipt<n_intpt; ipt++)
      {
       
       // Get the local coordinate at the integration point
       s[0] = integral_pt()->knot(ipt,0);
       
       // Get the integral weight
       double W = this->integral_pt()->weight(ipt);
       
       // Call the derivatives of the shape function at the knot point
       this->dshape_local_at_knot(ipt,psif,dpsifds);
       
       // Compute what we need...
       Vector<double> interpolated_tangent(2,0.0);
       Vector<double> interpolated_x(2,0.0);
       
       // Loop over the shape functions
       for(unsigned l=0; l<n_node; l++)
        {
         // Loop over directional components
         for(unsigned i=0; i<2; i++)
          {         
           // Spatial bits
           interpolated_x[i] += this->nodal_position(l,i)*psif(l);
           interpolated_tangent[i] += this->nodal_position(l,i)*dpsifds(l,0);
          }
        }
       
       // Calculate the length of the tangent Vector
       double tlength = interpolated_tangent[0]*interpolated_tangent[0] +
        interpolated_tangent[1]*interpolated_tangent[1];
       
       // Set the Jacobian of the line element
       double J = sqrt(tlength);
       
       // Add contribution
       deviation_contribution += W*J*
        (sqrt(interpolated_x[0]*interpolated_x[0] +
              interpolated_x[1]*interpolated_x[1]) - mean_radius)*
        (sqrt(interpolated_x[0]*interpolated_x[0] +
              interpolated_x[1]*interpolated_x[1]) - mean_radius);
      }
    }
    
    // rmodif
    /// Return the contribution to the mean radius
    void compute_mean_radius(double& arc_length_contribution, 
                             double& mean_radius_nod_contribution)
    {
     
     // Initialise
     mean_radius_nod_contribution = 0.0;
     arc_length_contribution = 0.0;

     // Find out how many nodes there are
     unsigned n_node = this->nnode();

     // Set up memory for the shape functions
     Shape psif(n_node);
     DShape dpsifds(n_node,1);
     
     // Set the value of n_intpt
     unsigned n_intpt = this->integral_pt()->nweight();
     
     // Storage for the local coordinate
     Vector<double> s(1);
     
     // Storage for the global coordinate
     Vector<double> x(2);
     
     // Loop over the integration points
     for(unsigned ipt=0; ipt<n_intpt; ipt++)
      {
       
       // Get the local coordinate at the integration point
       s[0] = integral_pt()->knot(ipt,0);
       
       // Get the integral weight
       double W = this->integral_pt()->weight(ipt);
       
       // Call the derivatives of the shape function at the knot point
       this->dshape_local_at_knot(ipt,psif,dpsifds);
     
       // Compute what we need...
       Vector<double> interpolated_tangent(2,0.0);
       Vector<double> interpolated_x(2,0.0);
     
       // Loop over the shape functions
       for(unsigned l=0; l<n_node; l++)
        {
         // Loop over directional components
         for(unsigned i=0; i<2; i++)
          {         
           // Spatial bits
           interpolated_x[i] += this->nodal_position(l,i)*psif(l);
           interpolated_tangent[i] += this->nodal_position(l,i)*dpsifds(l,0);
          }
        }
       
       // Calculate the length of the tangent Vector
       double tlength = interpolated_tangent[0]*interpolated_tangent[0] +
        interpolated_tangent[1]*interpolated_tangent[1];
       
       // Set the Jacobian of the line element
       double J = sqrt(tlength);
       
       // Add contribution
       arc_length_contribution += W*J;
       mean_radius_nod_contribution += W*J*
        sqrt(interpolated_x[0]*interpolated_x[0]+
             interpolated_x[1]*interpolated_x[1]);
      }
    }
    
    /// Fill in contribution to residuals and Jacobian
    void fill_in_contribution_to_jacobian(Vector<double> &residuals,
                                          DenseMatrix<double> &jacobian)
    {
     fill_in_generic_residual_contribution_hele_shaw_interface(residuals,
                                                               jacobian,1);
     
     /// ALICE: Need this for elastic version?
     this->fill_in_jacobian_from_solid_position_by_fd(jacobian);
     
     /// Alice: For my sanity later on. External dependencies must still be 
     /// added in driver code. May be useful for finding steady solutions where
     /// frame speed or capillary number are free parameters. 
     /// The residuals depend linearly on these parameters, so fd is not 
     /// completely unreasonable.
     /// NB. Similar hack in Poisson flux elements?
     this->fill_in_jacobian_from_external_by_fd(jacobian);
    }
    
    /// Calculate the residuals by calling the generic residual contribution.
    void fill_in_contribution_to_residuals(Vector<double> &residuals)
    {
     // Add the residual contributions
     fill_in_generic_residual_contribution_hele_shaw_interface(
      residuals,GeneralisedElement::Dummy_matrix,0);
    }
    
    /// \short The "global" intrinsic coordinate of the element when
    /// viewed as part of a geometric object should be given by
    /// the FaceElement representation, by default
    double zeta_nodal(const unsigned &n, const unsigned &k,
                      const unsigned &i) const
    {return FaceElement::zeta_nodal(n,k,i);}
    
    /// \short Virtual function that specifies the non-dimensional
    /// surface tension as a function of local position within the element.
    /// The default behaviour is a constant surface tension of value 1.0
    /// This function can be overloaded in more
    /// specialised elements to incorporate variations in surface tension.
    virtual double sigma(const Vector<double> &s_local)
    {return 1.0;}
    
    /// The value of the inverse Capillary number
    const double& ca_inv() const {return *Ca_inv_pt;}
    
    /// Pointer to the inverse Capillary number
    double*& ca_inv_pt() {return Ca_inv_pt;}
    
    /// The value of the Strouhal number
    const double &st() const {return *St_pt;}
    
    /// The pointer to the Strouhal number
    double* &st_pt() {return St_pt;}
    
    /// Aspect ratio: Reference gap width / in-plane lengthscale
    double aspect_ratio() const
    {
#ifdef PARANOID
     if (Aspect_ratio_pt==0)
      {
       throw OomphLibError(
        "Aspect_ratio_pt has not been set yet for HeleShaw elements",
        "HeleShawInterfaceElement::aspect_ratio()",
        OOMPH_EXCEPTION_LOCATION);
      }
#endif
     return *Aspect_ratio_pt;
    }
    
    /// \short Pointer to aspect ratio
    double*& aspect_ratio_pt()
     {
      return Aspect_ratio_pt;
     }    
    
    /// \short Pointer to aspect ratio. Const version.
    double* aspect_ratio_pt() const
     {
      return Aspect_ratio_pt;
     }
    
    /// Return the value of the external pressure
    double get_p_bubble(Vector<double>& x) const
    {
     if (Bubble_pressure_fct_pt==0)
      {
       return 0.0;
      }
     else
      {
       double pressure;
       (*Bubble_pressure_fct_pt)(x,pressure);
       return pressure;
      }
    }
    
    /// Return the value of the horizontal speed of the plates
    void get_wall_velocity(Vector<double>&x, Vector<double>& U)
    {
     if (Wall_speed_fct_pt==0)
      {
       U[0]=0.0;
       U[1]=0.0;
      }
     else
      {
       (*Wall_speed_fct_pt)(x,U);
      }
    }
    
    /// Overload the output functions
    void output(std::ostream &outfile) {FiniteElement::output(outfile);}
    
    /// Output the element
    void output(std::ostream &outfile, const unsigned &n_plot)
    {
     // Local coordinate
     Vector<double> s(1);
     
     // Find out how many nodes there are
     unsigned n_node = this->nnode();
     
     // Set up memory for the shape functions
     Shape psif(n_node);
     DShape dpsifds(n_node,1);
     
     // Get the value of the inverse Capillary number
     double Ca_inv = ca_inv();
     
     outfile << "ZONE\n";
     
     // Loop over plot points
     for (unsigned i_plot=0; i_plot<n_plot; i_plot++)
      {
       // Get coordinate
       get_s_plot(i_plot,n_plot,s);
       
       
       // Call the derivatives of the shape function at the knot point
       this->dshape_local(s,psif,dpsifds);
       
       // Compute what we need...
       Vector<double> interpolated_x(2,0.0);
       Vector<double> interpolated_dx_dt(2,0.0);
       Vector<double> tau(2,0.0);
       
       // Loop over the shape functions
       for(unsigned l=0; l<n_node; l++)
        {
         // Loop over directional components
         for(unsigned i=0; i<2; i++)
          {
           // Smoothed tangent vector
           tau[i] += projected_tangent_deriv(l,i)*psif(l);
           
           // Spatial bits
           interpolated_x[i] += this->nodal_position(l,i)*psif(l);
           interpolated_dx_dt[i] += this->dnodal_position_dt(l,i)*psif(l);
          }
        }
       
       // Now calculate the unit normal vector
       Vector<double> interpolated_n(2);
       outer_unit_normal(s,interpolated_n);
       
       // Also get the (possibly variable) surface tension
       double sigma_local = this->sigma(s);
       
       // Assemble the surface tension and normal speed terms
       double sigma_kappa=0.0;
       for(unsigned k=0; k<2; k++)
        {
         sigma_kappa+=sigma_local*Ca_inv*tau[k]*interpolated_n[k];
        }
       
       outfile << interpolated_x[0] << " "
               << interpolated_x[1] << " "
               << tau[0] << " "
               << tau[1] << " "
               << interpolated_dx_dt[0] << " "
               << interpolated_dx_dt[1] << " "
               << sigma_kappa << " "
               << std::endl;
      }
    }
    
    /// Overload the C-style output function
    void output(FILE* file_pt) {FiniteElement::output(file_pt);}
    
    /// C-style Output function
    void output(FILE* file_pt, const unsigned &n_plot)
    {
     // hierher fill this in when we've finalised what lives in output
    }
    
  };

 
 
 
//============================================================
/// Default value for physical constant (static)
//============================================================
 template<class ELEMENT>
  double HeleShawInterfaceElement<ELEMENT>::Default_Physical_Constant_Value = 
  1.0;
 
 
 
 
//=======================================================================
/// Calculate the residual contribution (kinematic and dynamic BC and
/// Lagrange multiplier contribution from pseudo-elastic node updates
//=======================================================================
 template<class ELEMENT>
  void HeleShawInterfaceElement<ELEMENT>::
  fill_in_generic_residual_contribution_hele_shaw_interface(
   Vector<double> &residuals,
   DenseMatrix<double> &jacobian,
   unsigned flag)
  {
   
   // Find out how many nodes there are
   unsigned n_node = this->nnode();
   
   // Set up memory for the shape functions
   Shape psif(n_node);
   DShape dpsifds(n_node,1);
   
   // Set the value of n_intpt
   unsigned n_intpt = this->integral_pt()->nweight();

   // Get the value of the inverse Capillary number
   double Ca_inv = ca_inv();
   
   // Get the value of the Strouhal numer
   double St = st();

   // Integers to store the local equation numbers
   int local_eqn=0;
   int local_unknown;
   
   // Storage for the local cooridinate
   Vector<double> s(1);
   
   // Loop over the integration points
   for(unsigned ipt=0; ipt<n_intpt; ipt++)
    {
     // Get the local coordinate at the integration point
     s[0] = integral_pt()->knot(ipt,0);
     
     // Get the integral weight
     double W = this->integral_pt()->weight(ipt);
     
     // Call the derivatives of the shape function at the knot point
     this->dshape_local_at_knot(ipt,psif,dpsifds);
     
     // Compute what we need...
     double interpolated_p=0.0;
     double interpolated_lagrange = 0.0;
     Vector<double> interpolated_tangent(2,0.0);
     Vector<double> interpolated_x(2,0.0);
     Vector<double> interpolated_dx_dt(2,0.0);
     Vector<double> tau(2,0.0);
     
     // Loop over the shape functions
     for(unsigned l=0; l<n_node; l++)
      {
       interpolated_lagrange += lagrange(l)*psif[l];
       interpolated_p += node_pt(l)->value(P_index_interface)*psif[l];
       
       // Loop over directional components
       for(unsigned i=0; i<2; i++)
        {
         // Smoothed tangent vector
         tau[i] += projected_tangent_deriv(l,i)*psif(l);
         
         // Spatial bits
         interpolated_x[i] += this->nodal_position(l,i)*psif(l);
         interpolated_dx_dt[i] += this->dnodal_position_dt(l,i)*psif(l);
         interpolated_tangent[i] += this->nodal_position(l,i)*dpsifds(l,0);
        }
      }
     
     // Get the value of the bubble pressure
     double p_bubble = get_p_bubble(interpolated_x);
     
     // Get the velocity of the plates
     Vector<double> U_wall(2,0.0);
     get_wall_velocity(interpolated_x,U_wall);
     
     // Calculate the length of the tangent Vector
     double tlength = interpolated_tangent[0]*interpolated_tangent[0] +
      interpolated_tangent[1]*interpolated_tangent[1];
     
     // Set the Jacobian of the line element
     double J = sqrt(tlength);
     
     // Normalise the tangent Vector
     interpolated_tangent[0] /= J;
     interpolated_tangent[1] /= J;
     
     // Now calculate the unit normal vector
     Vector<double> interpolated_n(2);
     outer_unit_normal(ipt,interpolated_n);
     
     // Also get the (possibly variable) surface tension
     // double sigma_local = this->sigma(s);
     
     // Assemble the surface tension and normal speed terms
     // double sigma_kappa = 0.0;

     // Curvature of the bubble
     double kappa = 0.0;

     // Speed of the bubble boundary
     double normal_speed = 0.0;

     // Horizontal speed of the plates
     double normal_speed_wall = 0.0;

     for(unsigned k=0; k<2; k++)
      {
       // sigma_kappa += sigma_local*Ca_inv*tau[k]*interpolated_n[k];
       kappa += tau[k]*interpolated_n[k];
       normal_speed += interpolated_dx_dt[k]*interpolated_n[k];
       normal_speed_wall += U_wall[k]*interpolated_n[k];
      }
     
     // Get pointer to bulk element
     ELEMENT* bulk_el_pt = dynamic_cast<ELEMENT*>(
      this->bulk_element_pt());
     
     // Get h from bulk!
     Vector<double> s_bulk(2);
     get_local_coordinate_in_bulk(s,s_bulk);
     double h = 0.0;
     double dhdt = 0.0;

     //---------------------------------------------------------
     // NOTE: This is dangerous/wrong. We're really supposed to
     // pass in dpsi/dx (2D!) which we don't have/need here.
     // However it's only used for the computation of the ALE bits
     // in dh/dt and dh/dt isn't used here! Phew...
     //---------------------------------------------------------
     // Number of nodes in bulk element
     unsigned n_nod_bulk=bulk_el_pt->nnode();
     DShape dpsidx_bulk(n_nod_bulk,2);
     Shape psi_bulk(n_nod_bulk);
     bulk_el_pt->shape(s_bulk,psi_bulk);
     bulk_el_pt->get_upper_wall_data(s_bulk,interpolated_x,
                                     psi_bulk,dpsidx_bulk,h,dhdt);
     
     // Non-dim gap-width
     double local_aspect_ratio = aspect_ratio();
          
     // Loop over the shape functions
     for(unsigned l=0; l<n_node; l++)
      {
       
       // Eqns to determine the smoothed derivatives of the tangent vector
       // ----------------------------------------------------------------
       for(unsigned i=0; i<2; i++)
        {
         local_eqn=projected_tangent_deriv_local_eqn(l,i);
         
         // If it's not a boundary condition
         if(local_eqn >= 0)
          {
           residuals[local_eqn] +=
            (tau[i]*psif(l)*J+interpolated_tangent[i]*dpsifds(l,0))*W;
          }
         
         // Do the Jacobian calculation
         if(flag)
          {
           // Loop over the nodes
           for(unsigned l2=0; l2<n_node; l2++)
            {
             // Derivatives w.r.t. solid positions will be handled by FDing.
             // Only 'tau' above is not solid position.
             local_unknown = projected_tangent_deriv_local_eqn(l2,i);
             if(local_unknown >= 0)
              {
               jacobian(local_eqn,local_unknown) +=
                psif[l2]*psif(l)*W*J;
              }
            }
          } // End of Jacobian calculation
        }
       
       
       // Eqn for Lagrange multiplier (dynamic free surface condition)
       //-------------------------------------------------------------
       local_eqn = lagrange_local_eqn(l);
       if(local_eqn >= 0)
        {
         residuals[local_eqn] +=
          (interpolated_p-p_bubble
           + local_aspect_ratio*Ca_inv*
           (2.0/h+local_aspect_ratio*kappa)/12.0)*psif(l)*W*J;
         
         if (flag)
          {
           // Loop over the nodes
           for(unsigned l2=0; l2<n_node; l2++)
            {
             // Derivatives w.r.t. solid positions will be handled by FDing
             /// Lagrange residual depends on interpolated pressure
             local_unknown = nodal_local_eqn(l2,P_index_interface);
             if(local_unknown >= 0)
              {
               // This is for interpolated pressure in residual equation
               jacobian(local_eqn,local_unknown) +=
                psif[l2]*psif(l)*W*J;
              }
             
             /// Lagrange residual also depends on smoothed derivatives through
             /// curvature term
             for (unsigned i2=0; i2<2; i2++)
              {
               local_unknown = projected_tangent_deriv_local_eqn(l2,i2);
               if (local_unknown>=0)
                {
                 // This bit relates to d_kappa /d tau.
                 jacobian(local_eqn, local_unknown) +=
                  local_aspect_ratio*Ca_inv*
                  local_aspect_ratio/12.0*psif(l)*W*J*
                  psif(l2)*interpolated_n[i2];
                }
              }
            }
          }
        }
       
       // Contribution to bulk equation (kinematic bc)
       //---------------------------------------------
       local_eqn=nodal_local_eqn(l,P_index_interface);
       if (local_eqn>=0)
        {
         residuals[local_eqn] += h*St*normal_speed*psif(l)*W*J;
         // residuals[local_eqn] += h*normal_speed_wall*psif(l)*W*J;
         
         /// These residuals depend only on solid position, and
         /// possibly external data.
        }
       
       // Lagrange multiplier contributions to pseudo_solid equations
       //------------------------------------------------------------
       for(unsigned i=0; i<2; i++)
        {
         local_eqn = this->position_local_eqn(l,0,i);
         
         if(local_eqn >= 0)
          {
           // Add in a "Lagrange multiplier"
           residuals[local_eqn] +=
            -interpolated_lagrange*interpolated_n[i]*psif[l]*W*J;
           
           // Do the Jacobian calculation
           if(flag)
            {
             // Loop over the nodes
             for(unsigned l2=0; l2<n_node; l2++)
              {
               // Derivatives w.r.t. solid positions will be handled by FDing
               // That leaves the "lagrange multipliers" only
               local_unknown = lagrange_local_eqn(l2);
               if(local_unknown >= 0)
                {
                 jacobian(local_eqn,local_unknown) +=
                  -psif[l2]*interpolated_n[i]*psif[l]*W*J;
                }
              }
            } // End of Jacobian calculation
          }
        }       
      } // End of loop over shape functions
    } // End of loop over integration points
  }
 
 
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
 
 
 
 
}


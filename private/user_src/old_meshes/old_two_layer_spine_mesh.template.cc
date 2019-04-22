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
#ifndef OOMPH_TWO_LAYER_SPINE_MESH_TEMPLATE_CC
#define OOMPH_TWO_LAYER_SPINE_MESH_TEMPLATE_CC

#include "two_layer_spine_mesh.template.h"
#include "rectangular_quadmesh.template.cc"


namespace oomph
{


///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

// //=========================================================================
// /// \short Constructor for a non-uniform two layer spine mesh, with 
// /// element layout in the lower fluid reflected in the upper. Three
// /// distinct y regions need numbers of element specified and two 
// /// x regions.
// //=========================================================================
// template <class ELEMENT, class INTERFACE_ELEMENT >
// Axisym2x6TwoLayerSpineMesh<ELEMENT,INTERFACE_ELEMENT>::
// Axisym2x6TwoLayerSpineMesh(const unsigned &nxa, const unsigned &nxb, 
//                            const unsigned &nya1, const unsigned &nyb1,
//                            const unsigned &nyc1, const unsigned &nya2, 
//                            const unsigned &nyb2, const unsigned &nyc2,
//                            const double &lx, 
//                            const double &h1, const double &h2,
//                            TimeStepper* time_stepper_pt) :
//  TwoLayerSpineMesh<ELEMENT, INTERFACE_ELEMENT >()
// {
//  // We've called the "generic" constructor for the RectangularQuadMesh
//  // which doesn't do much...
//  // Now set up the parameters that characterise the mesh geometry
//  // then call the constructor

//  Nxa = nxa; Nxb = nxb;
//  Nya1 = nya1; Nyb1 = nyb1; Nyc1 = nyc1;
//  Nya2 = nya2; Nyb2 = nyb2; Nyc2 = nyc2;
//  Xfraction = 0.8;
//  Yfraction1 = 0.2;
//  Yfraction2 = 0.8;
 
//  // Check validity of Xfraction, Yfraction1 and Yfraction2.
//  if (Yfraction1 < 0.0 || Yfraction1>Yfraction2)
//   {
//    throw OomphLibError("Invalid Yfraction1",
//                        "AxiSym2x6TwoLayerSpineMesh",
//                        OOMPH_EXCEPTION_LOCATION);
//   }
//  if (Yfraction2 < Yfraction1 || Yfraction2>1.0)
//   {
//    throw OomphLibError("Invalid Yfraction2",
//                        "AxiSym2x6TwoLayerSpineMesh",
//                        OOMPH_EXCEPTION_LOCATION);
//   }
//  if (Xfraction<0.0 || Xfraction>1.0)
//   {
//    throw OomphLibError("Invalid xfraction",
//                        "AxiSym2x6TwoLayerSpineMesh",
//                        OOMPH_EXCEPTION_LOCATION);
//   }


//  // Number of elements in x direction
//  RectangularQuadMesh<ELEMENT >::Nx = Nxa+Nxb; 
 
//  // Number of elements in bottom and top layers
//  this->Ny1 = Nya1+Nyb1+Nyc1;
//  this->Ny2 = Nya2+Nyb2+Nyc2;
 
//  // Number of elements in y direction
//  RectangularQuadMesh<ELEMENT >::Ny = this->Ny1 + this->Ny2;
 
//  // Min. x coordinate
//  RectangularQuadMesh<ELEMENT >::Xmin = 0.0; 
 
//  // Max. x coordinate
//  RectangularQuadMesh<ELEMENT >::Xmax = lx; 
 
//  // Min. y coordinate
//  RectangularQuadMesh<ELEMENT >::Ymin = 0.0; 
 
//  // Max. y coordinate
//  RectangularQuadMesh<ELEMENT >::Ymax = h1+h2;
 
//  // Periodic?
//  RectangularQuadMesh<ELEMENT >::Xperiodic = false;
 
//  // Set height of upper and lower layers
//  this->H1 = h1;
//  this->H2 = h2;
 
//  // Now build the mesh: 
//  this->build_two_layer_mesh(time_stepper_pt);
 
// }


//=========================================================================
/// \short Constructor for a non-uniform two layer spine mesh, with 
/// element layout in the lower fluid reflected in the upper. Three
/// distinct y regions need numbers of element specified and two 
/// x regions. The fractions of these regions are specified in this 
/// constructor.
// PATRICKFLAG this is the one I've been editing...
//=========================================================================
template <class ELEMENT, class INTERFACE_ELEMENT >
Axisym2x6TwoLayerSpineMesh<ELEMENT,INTERFACE_ELEMENT>::
Axisym2x6TwoLayerSpineMesh(const unsigned &nxa, const unsigned &nxb, 
                           const double &x_frac, 
                           const unsigned &nya1, const unsigned &nyb1, 
                           const unsigned &nyc1, const unsigned &nya2,
                           const unsigned &nyb2, const unsigned &nyc2, 
                           const double &y1_frac1, const double &y1_frac2,
                           const double &y2_frac1, const double &y2_frac2,
                           const double &lx, 
                           const double &h1, const double &h2,
                           TimeStepper* time_stepper_pt) :
 TwoLayerSpineMesh<ELEMENT,INTERFACE_ELEMENT>(nxa+nxb,nya1+nyb1+nyc1,
                                              nya2+nyb2+nyc2,lx,h1,h2,
                                              false,false,time_stepper_pt)
{
 // Mesh can only be built with 2D Qelements.
 MeshChecker::assert_geometric_element<QElementGeometricBase,ELEMENT>(2);

 // We've called the "generic" constructor for the TwoLayerSpineMesh and
 // set the "build_mesh" flag to false. This is done so that we do not
 // call "build_two_layer_mesh(...)" prematurely. We now set up the
 // parameters that characterise this particular mesh's geometry before
 // calling "build_two_layer_mesh(...)".
 
 Nxa = nxa; Nxb = nxb;
 Nya1 = nya1; Nyb1 = nyb1; Nyc1 = nyc1;
 Nya2 = nya2; Nyb2 = nyb2; Nyc2 = nyc2;
 Xfraction = x_frac;
 Y1fraction1 = y1_frac1; Y1fraction2 = y1_frac2;
 Y2fraction1 = y2_frac1; Y2fraction2 = y2_frac2;

 // Check validaty of Xfraction
 if (Xfraction<0.0 || Xfraction>1.0)
  {
   throw OomphLibError("Invalid Xfraction",
                       "Axisym2x6TwoLayerSpineMesh",
                       OOMPH_EXCEPTION_LOCATION);
  }

 // Check validaty of Y1fraction1 and Y2fraction1
 if (Y1fraction1<0.0 || Y1fraction1>Y1fraction2 || Y1fraction1>1.0)
  {
   throw OomphLibError("Invalid Y1fraction1",
                       "Axisym2x6TwoLayerSpineMesh",
                       OOMPH_EXCEPTION_LOCATION);
  }
 if (Y2fraction1<0.0 || Y2fraction1>Y2fraction2 || Y2fraction1>1.0)
  {
   throw OomphLibError("Invalid Y2fraction1",
                       "Axisym2x6TwoLayerSpineMesh",
                       OOMPH_EXCEPTION_LOCATION);
  }

 // Check validaty of Y1fraction2 and Y2fraction2
 if (Y1fraction2<0.0 || Y1fraction2<Y1fraction1 || Y1fraction2>1.0)
  {
   throw OomphLibError("Invalid Y1fraction2",
                       "Axisym2x6TwoLayerSpineMesh",
                       OOMPH_EXCEPTION_LOCATION);
  }
 if (Y2fraction2<0.0 || Y2fraction2<Y2fraction1 || Y2fraction2>1.0)
  {
   throw OomphLibError("Invalid Y2fraction2",
                       "Axisym2x6TwoLayerSpineMesh",
                       OOMPH_EXCEPTION_LOCATION);
  }

 
 // Number of elements in x direction
// RectangularQuadMesh<ELEMENT >::Nx = Nxa+Nxb; 
 
 // // Store number of elements in bottom and top layers
//  this->Ny1 = Nya1+Nyb1+Nyc1;
//  this->Ny2 = Nya2+Nyb2+Nyc2;
 
// //  // Number of elements in y direction
// //  RectangularQuadMesh<ELEMENT >::Ny = this->Ny1 + this->Ny2;
 
// //  // Min. x coordinate
// //  RectangularQuadMesh<ELEMENT >::Xmin = 0.0; 
 
// //  // Max. x coordinate
// //  RectangularQuadMesh<ELEMENT >::Xmax = lx; 
 
// //  // Min. y coordinate
// //  RectangularQuadMesh<ELEMENT >::Ymin = 0.0; 
 
// //  // Max. y coordinate
// //  RectangularQuadMesh<ELEMENT >::Ymax = h1+h2;
 
// //  // Periodic?
// //  RectangularQuadMesh<ELEMENT >::Xperiodic = false;
 
//  // Store height of upper and lower layers
//  this->H1 = h1;
//  this->H2 = h2;

 // Now build the mesh: 
 this->build_two_layer_mesh(time_stepper_pt);
 
}


/// \short The spacing function for the x co-ordinates with two 
/// regions.
template <class ELEMENT, class INTERFACE_ELEMENT >
double Axisym2x6TwoLayerSpineMesh<ELEMENT,INTERFACE_ELEMENT>::
x_spacing_function(unsigned xelement, unsigned xnode,
                   unsigned yelement, unsigned ynode)
 {
  // Set up some spacing parameters

  // Left region starts at Xmin
  double Xmin =  RectangularQuadMesh<ELEMENT >::Xmin;
  
  // Right region ends at Xmax
  double Xmax =  RectangularQuadMesh<ELEMENT >::Xmax;

  // Number of nodes per element (in one direction)
  unsigned n_p =  RectangularQuadMesh<ELEMENT >::Np;

  // Left region starts at Xmin
  double x1init = Xmin;

  // Right region starts at Xmin + Xfraction(Xmax-Xmin)
  double x2init = Xmin + Xfraction*(Xmax-Xmin);

  // Assuming uniform spacing, calculate the spacing between
  // the nodes in each region

  // Left region has a length Xfraction*(Xmax-Xmin)
  double x1step = Xfraction*(Xmax-Xmin)/((n_p-1)*Nxa);

  // Right region has a length (1.0-Xfraction)*(Xmax-Xmin)
  double x2step = (1.0-Xfraction)*(Xmax-Xmin)/((n_p-1)*Nxb);
  
  // Now set up the particular spacing in each region
  if(xelement < Nxa)
   {
    return (x1init + x1step*((n_p-1)*xelement + xnode));
   }
  else
   {
    return (x2init + x2step*((n_p-1)*(xelement-Nxa) + xnode));
   }
 }

/// \short The spacing function for the y co-ordinates with three
/// regions in each fluid.
template <class ELEMENT, class INTERFACE_ELEMENT >
double Axisym2x6TwoLayerSpineMesh<ELEMENT,INTERFACE_ELEMENT>::
y_spacing_function(unsigned xelement, unsigned xnode,
                   unsigned yelement, unsigned ynode)
{
 //Set up some spacing parameters
 //The lower region a starts at Ymin
 double Ymin = RectangularQuadMesh<ELEMENT >::Ymin;
 //The interface is at Ymid
 double Ymid = this->H1;
 //The upper region a ends at Ymax
 double Ymax = RectangularQuadMesh<ELEMENT >::Ymax;
 //Number of nodes per element
 unsigned n_p = RectangularQuadMesh<ELEMENT >::Np;
 //The lower region a starts at Ymin
 double y1init = Ymin;
 //The lower region b starts at Ymin + Y1fraction1*(Ymid-Ymin)
 double y2init = Ymin + Y1fraction1*(Ymid-Ymin);
 //The lower region c starts at Ymin + Y1fraction2*(Ymid-Ymin)
 double y3init = Ymin + Y1fraction2*(Ymid-Ymin);
 //The upper region c starts at Ymid
 double y4init = Ymid;
 //The upper region b starts at Ymax - Y2fraction2*(Ymax-Ymid)
 double y5init = Ymax - Y2fraction2*(Ymax-Ymid);
 //The upper region a starts at Ymax - Y2fraction1*(Ymax-Ymid)
 double y6init = Ymax - Y2fraction1*(Ymax-Ymid);
 //Calculate the space between each node in each region,
 //Assumming uniform spacing
 //Lower region a has a length Y1fraction1(Ymid-Ymin)
 double y1step = Y1fraction1*(Ymid-Ymin)/((n_p-1)*Nya1);
 //Lower region b has a length (Y1fraction2-Y1fraction1)*(Ymid-Ymin)
 double y2step = (Y1fraction2-Y1fraction1)*(Ymid-Ymin)/((n_p-1)*Nyb1);
 //Lower region c has a length (1.0-Y1fraction2)*(Ymid-Ymin)
 double y3step = (1.0-Y1fraction2)*(Ymid-Ymin)/((n_p-1)*Nyc1);
 //Upper region c has a length (1.0-Y1fraction2)*(Ymax-Ymid)
 double y4step = (1.0-Y2fraction2)*(Ymax-Ymid)/((n_p-1)*Nyc2);
 //Upper region b has a length (Y1fraction2-Y1fraction1)*(Ymax-Ymid)
 double y5step = (Y2fraction2-Y2fraction1)*(Ymax-Ymid)/((n_p-1)*Nyb2);
 //Upper region a has a length Y1fraction1(Ymax-Ymid)
 double y6step = Y2fraction1*(Ymax-Ymid)/((n_p-1)*Nya2);

 //Now return the actual node position, it's different in the two
 //regions, of course
 if(yelement < Nya1) 
  {
   return (y1init + y1step*((n_p-1)*yelement + ynode));
  }
 else if (yelement < Nya1+Nyb1)
  {
   return (y2init + y2step*((n_p-1)*(yelement-Nya1) + ynode));
  }
 else if (yelement < Nya1+Nyb1+Nyc1)
  {
   return (y3init + y3step*((n_p-1)*(yelement-Nya1-Nyb1) + ynode));
  }
 else if (yelement < Nya1+Nyb1+Nyc1+Nyc2)
  {
   return (y4init + y4step*((n_p-1)*(yelement-Nya1-Nyb1-Nyc1) + ynode));
  }
 else if (yelement < Nya1+Nyb1+Nyc1+Nyc2+Nyb2)
  {
   return (y5init + y5step*((n_p-1)*(yelement-Nya1-Nyb1-Nyc1-Nyc2) + ynode));
  }
 else
  {
   return (y6init + 
           y6step*((n_p-1)*(yelement-Nya1-Nyb1-Nyc1-Nyc2-Nyb2) + ynode));
  }
}

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////


// //=========================================================================
// /// \short Constructor for a non-uniform two layer spine mesh, with 
// /// element layout in the lower fluid reflected in the upper. Three
// /// distinct y regions need numbers of element specified and two 
// /// x regions.
// //=========================================================================
// template <class ELEMENT, class INTERFACE_ELEMENT >
// Axisym3x6TwoLayerSpineMesh<ELEMENT,INTERFACE_ELEMENT>::
// Axisym3x6TwoLayerSpineMesh(const unsigned &nxa, const unsigned &nxb, 
//                            const unsigned &nxc,
//                            const unsigned &nya1, const unsigned &nyb1, 
//                            const unsigned &nyc1, const unsigned &nya2, 
//                            const unsigned &nyb2, const unsigned &nyc2, 
//                            const double &lx, 
//                            const double &h1, const double &h2,
//                            TimeStepper* time_stepper_pt) :
//  TwoLayerSpineMesh<ELEMENT, INTERFACE_ELEMENT >()
// {
//  // We've called the "generic" constructor for the RectangularQuadMesh
//  // which doesn't do much...
//  // Now set up the parameters that characterise the mesh geometry
//  // then call the constructor

//  Nxa = nxa; Nxb = nxb; Nxc = nxc; 
//  Nya1 = nya1; Nyb1 = nyb1; Nyc1 = nyc1;
//  Nya2 = nya2; Nyb2 = nyb2; Nyc2 = nyc2;
//  Xfraction1 = 0.2;
//  Xfraction2 = 0.8;
//  Yfraction1 = 0.2;
//  Yfraction2 = 0.8;
 
//  // Check validity of Xfraction, Yfraction1 and Yfraction2.
//  if (Yfraction1 < 0.0 || Yfraction1>Yfraction2)
//   {
//    throw OomphLibError("Invalid Yfraction1",
//                        "AxiSymm3x6TwoLayerSpineMesh",
//                        OOMPH_EXCEPTION_LOCATION);
//   }
//  if (Yfraction2 < Yfraction1 || Yfraction2>1.0)
//   {
//    throw OomphLibError("Invalid Yfraction2",
//                        "AxiSymm3x6TwoLayerSpineMesh",
//                        OOMPH_EXCEPTION_LOCATION);
//   }
//  if (Xfraction1<0.0 || Xfraction1>Xfraction2)
//   {
//    throw OomphLibError("Invalid Xfraction1",
//                        "AxiSymm3x6TwoLayerSpineMesh",
//                        OOMPH_EXCEPTION_LOCATION);
//   }
//  if (Xfraction2 < Xfraction1 || Xfraction2>1.0)
//   {
//     throw OomphLibError("Invalid Xfraction2",
//                         "AxiSymm3x6TwoLayerSpineMesh",
//                         OOMPH_EXCEPTION_LOCATION);
//   }

//  // Number of elements in x direction
//  RectangularQuadMesh<ELEMENT >::Nx = Nxa+Nxb+Nxc; 
 
//  // Number of elements in bottom and top layers
//  this->Ny1 = Nya1+Nyb1+Nyc1;
//  this->Ny2 = Nya2+Nyb2+Nyc2;
 
//  // Number of elements in y direction
//  RectangularQuadMesh<ELEMENT >::Ny = this->Ny1 + this->Ny2;
 
//  // Min. x coordinate
//  RectangularQuadMesh<ELEMENT >::Xmin = 0.0; 
 
//  // Max. x coordinate
//  RectangularQuadMesh<ELEMENT >::Xmax = lx; 
 
//  // Min. y coordinate
//  RectangularQuadMesh<ELEMENT >::Ymin = 0.0; 
 
//  // Max. y coordinate
//  RectangularQuadMesh<ELEMENT >::Ymax = h1;
 
//  // Periodic?
//  RectangularQuadMesh<ELEMENT >::Xperiodic = false;
 
//  // Set height of upper and lower layers
//  this->H1 = h1;
//  this->H2 = h2;
 
//  // Now build the mesh: 
//  this->build_two_layer_mesh(time_stepper_pt);
 
// }


//=========================================================================
/// \short Constructor for a non-uniform two layer spine mesh, with 
/// element layout in the lower fluid reflected in the upper. Three
/// distinct y regions need numbers of element specified and two 
/// x regions. The fractions of these regions are specified in this 
/// constructor.
// PATRICKFLAG I've also been editing this one...
//=========================================================================
template <class ELEMENT, class INTERFACE_ELEMENT >
Axisym3x6TwoLayerSpineMesh<ELEMENT,INTERFACE_ELEMENT>::
Axisym3x6TwoLayerSpineMesh(const unsigned &nxa, const unsigned &nxb, 
                           const unsigned &nxc, 
                           const double &x_frac1, const double &x_frac2,
                           const unsigned &nya1, const unsigned &nyb1, 
                           const unsigned &nyc1, const unsigned &nya2, 
                           const unsigned &nyb2, const unsigned &nyc2,
                           const double &y1_frac1, const double &y1_frac2,
                           const double &y2_frac1, const double &y2_frac2,
                           const double &lx, 
                           const double &h1, const double &h2,
                           TimeStepper* time_stepper_pt) :
 TwoLayerSpineMesh<ELEMENT,INTERFACE_ELEMENT>(nxa+nxb+nxc,nya1+nyb1+nyc1,
                                              nya2+nyb2+nyc2,lx,h1,h2,
                                              false,false,time_stepper_pt)
{
 // Mesh can only be built with 2D Qelements.
 MeshChecker::assert_geometric_element<QElementGeometricBase,ELEMENT>(2);

 // We've called the "generic" constructor for the TwoLayerSpineMesh and
 // set the "build_mesh" flag to false. This is done so that we do not
 // call "build_two_layer_mesh(...)" prematurely. We now set up the
 // parameters that characterise this particular mesh's geometry before
 // calling "build_two_layer_mesh(...)".

 Nxa = nxa; Nxb = nxb; Nxc = nxc; 
 Nya1 = nya1; Nyb1 = nyb1; Nyc1 = nyc1;
 Nya2 = nya2; Nyb2 = nyb2; Nyc2 = nyc2;
 Xfraction1 = x_frac1;  Xfraction2 = x_frac2;
 Y1fraction1 = y1_frac1; Y1fraction2 = y1_frac2;
 Y2fraction1 = y2_frac1; Y2fraction2 = y2_frac2;

 // Check validaty of Xfraction1
 if (Xfraction1<0.0 || Xfraction1>Xfraction2 || Xfraction1>1.0)
  {
   throw OomphLibError("Invalid Xfraction1",
                       "Axisym3x6TwoLayerSpineMesh",
                       OOMPH_EXCEPTION_LOCATION);
  }

 // Check validaty of Xfraction2
 if (Xfraction2<0.0 || Xfraction2<Xfraction1 || Xfraction2>1.0)
  {
   throw OomphLibError("Invalid Xfraction2",
                       "Axisym3x6TwoLayerSpineMesh",
                       OOMPH_EXCEPTION_LOCATION);
  }

 // Check validaty of Y1fraction1 and Y2fraction1
 if (Y1fraction1<0.0 || Y1fraction1>Y1fraction2 || Y1fraction1>1.0)
  {
   throw OomphLibError("Invalid Y1fraction1",
                       "Axisym3x6TwoLayerSpineMesh",
                       OOMPH_EXCEPTION_LOCATION);
  }
 if (Y2fraction1<0.0 || Y2fraction1>Y2fraction2 || Y2fraction1>1.0)
  {
   throw OomphLibError("Invalid Y2fraction1",
                       "Axisym3x6TwoLayerSpineMesh",
                       OOMPH_EXCEPTION_LOCATION);
  }

 // Check validaty of Y1fraction2 and Y2fraction2
 if (Y1fraction2<0.0 || Y1fraction2<Y1fraction1 || Y1fraction2>1.0)
  {
   throw OomphLibError("Invalid Y1fraction2",
                       "Axisym3x6TwoLayerSpineMesh",
                       OOMPH_EXCEPTION_LOCATION);
  }
 if (Y2fraction2<0.0 || Y2fraction2<Y2fraction1 || Y2fraction2>1.0)
  {
   throw OomphLibError("Invalid Y2fraction2",
                       "Axisym3x6TwoLayerSpineMesh",
                       OOMPH_EXCEPTION_LOCATION);
  }
 
//  // Number of elements in x direction
//  RectangularQuadMesh<ELEMENT >::Nx = Nxa+Nxb+Nxc; 
 
//  // Number of elements in bottom and top layers
//  this->Ny1 = Nya1+Nyb1+Nyc1;
//  this->Ny2 = Nya2+Nyb2+Nyc2;
 
//  // Number of elements in y direction
//  RectangularQuadMesh<ELEMENT >::Ny = this->Ny1 + this->Ny2;
 
//  // Min. x coordinate
//  RectangularQuadMesh<ELEMENT >::Xmin = 0.0; 
 
//  // Max. x coordinate
//  RectangularQuadMesh<ELEMENT >::Xmax = lx; 
 
//  // Min. y coordinate
//  RectangularQuadMesh<ELEMENT >::Ymin = 0.0; 
 
//  // Max. y coordinate
//  RectangularQuadMesh<ELEMENT >::Ymax = h1+h2;
 
//  // Periodic?
//  RectangularQuadMesh<ELEMENT >::Xperiodic = false;
 
//  // Set height of upper and lower layers
//  this->H1 = h1;
//  this->H2 = h2;

 // Now build the mesh: 
 this->build_two_layer_mesh(time_stepper_pt);
 
}


/// \short The spacing function for the x co-ordinates with two 
/// regions.
template <class ELEMENT, class INTERFACE_ELEMENT >
double Axisym3x6TwoLayerSpineMesh<ELEMENT,INTERFACE_ELEMENT>::
x_spacing_function(unsigned xelement, unsigned xnode,
                   unsigned yelement, unsigned ynode)
 {
  //Set up some spacing parameters
  //region a starts at Xmin
  double Xmin =  RectangularQuadMesh<ELEMENT >::Xmin;
  //region c ends at Xmax
  double Xmax =  RectangularQuadMesh<ELEMENT >::Xmax;
  //Number of nodes per element
  unsigned n_p =  RectangularQuadMesh<ELEMENT >::Np;
  //region a starts at Xmin
  double x1init = Xmin;
  //region b starts at Xmin + Xfraction1*(Xmax-Xmin)
  double x2init = Xmin + Xfraction1*(Xmax-Xmin);
  //region c starts at Xmin + Xfraction2*(Xmax-Xmin)
  double x3init = Xmin + Xfraction2*(Xmax-Xmin);
  //Calculate the spacing between the nodes in each region
  //Assuming uniform spacing
  //region a has a length Xfraction1*(Xmax-Xmin)
  double x1step = Xfraction1*(Xmax-Xmin)/((n_p-1)*Nxa);
  //region b has a length (Xfraction2-Xfraction1)*(Xmax-Xmin)
  double x2step = (Xfraction2-Xfraction1)*(Xmax-Xmin)/((n_p-1)*Nxb);
  //region c has a length (1.0-Xfraction2)*(Xmax-Xmin)
  double x3step = (1.0-Xfraction2)*(Xmax-Xmin)/((n_p-1)*Nxc);
  
  //Now set up the particular spacing 
  //(it's different in the two different regions)
  if(xelement < Nxa)
   {
    return (x1init + x1step*((n_p-1)*xelement + xnode));
   }
  else
   {
    if(xelement < Nxa+Nxb)
     {
      return (x2init + x2step*((n_p-1)*(xelement-Nxa) + xnode));
     }
    else
     {
      return (x3init + x3step*((n_p-1)*(xelement-Nxa-Nxb) + xnode));
     }
   }
 }

/// \short The spacing function for the y co-ordinates with three
/// regions in each fluid.
template <class ELEMENT, class INTERFACE_ELEMENT >
double Axisym3x6TwoLayerSpineMesh<ELEMENT,INTERFACE_ELEMENT>::
y_spacing_function(unsigned xelement, unsigned xnode,
                   unsigned yelement, unsigned ynode)
{
 //Set up some spacing parameters
 //The lower region a starts at Ymin
 double Ymin = RectangularQuadMesh<ELEMENT >::Ymin;
 //The interface is at Ymid
 double Ymid = this->H1;
 //The upper region a ends at Ymax
 double Ymax = RectangularQuadMesh<ELEMENT >::Ymax;
 //Number of nodes per element
 unsigned n_p = RectangularQuadMesh<ELEMENT >::Np;
 //The lower region a starts at Ymin
 double y1init = Ymin;
 //The lower region b starts at Ymin + Y1fraction1*(Ymid-Ymin)
 double y2init = Ymin + Y1fraction1*(Ymid-Ymin);
 //The lower region c starts at Ymin + Y1fraction2*(Ymid-Ymin)
 double y3init = Ymin + Y1fraction2*(Ymid-Ymin);
 //The upper region c starts at Ymid
 double y4init = Ymid;
 //The upper region b starts at Ymax - Y2fraction2*(Ymax-Ymid)
 double y5init = Ymax - Y2fraction2*(Ymax-Ymid);
 //The upper region a starts at Ymax - Y2fraction1*(Ymax-Ymid)
 double y6init = Ymax - Y2fraction1*(Ymax-Ymid);
 //Calculate the space between each node in each region,
 //Assumming uniform spacing
 //Lower region a has a length Y1fraction1(Ymid-Ymin)
 double y1step = Y1fraction1*(Ymid-Ymin)/((n_p-1)*Nya1);
 //Lower region b has a length (Y1fraction2-Y1fraction1)*(Ymid-Ymin)
 double y2step = (Y1fraction2-Y1fraction1)*(Ymid-Ymin)/((n_p-1)*Nyb1);
 //Lower region c has a length (1.0-Y1fraction2)*(Ymid-Ymin)
 double y3step = (1.0-Y1fraction2)*(Ymid-Ymin)/((n_p-1)*Nyc1);
 //Upper region c has a length (1.0-Y1fraction2)*(Ymax-Ymid)
 double y4step = (1.0-Y2fraction2)*(Ymax-Ymid)/((n_p-1)*Nyc2);
 //Upper region b has a length (Y1fraction2-Y1fraction1)*(Ymax-Ymid)
 double y5step = (Y2fraction2-Y2fraction1)*(Ymax-Ymid)/((n_p-1)*Nyb2);
 //Upper region a has a length Y1fraction1(Ymax-Ymid)
 double y6step = Y2fraction1*(Ymax-Ymid)/((n_p-1)*Nya2);

 //Now return the actual node position, it's different in the two
 //regions, of course
 if(yelement < Nya1) 
  {
   return (y1init + y1step*((n_p-1)*yelement + ynode));
  }
 else if (yelement < Nya1+Nyb1)
  {
   return (y2init + y2step*((n_p-1)*(yelement-Nya1) + ynode));
  }
 else if (yelement < Nya1+Nyb1+Nyc1)
  {
   return (y3init + y3step*((n_p-1)*(yelement-Nya1-Nyb1) + ynode));
  }
 else if (yelement < Nya1+Nyb1+Nyc1+Nyc2)
  {
   return (y4init + y4step*((n_p-1)*(yelement-Nya1-Nyb1-Nyc1) + ynode));
  }
 else if (yelement < Nya1+Nyb1+Nyc1+Nyc2+Nyb2)
  {
   return (y5init + y5step*((n_p-1)*(yelement-Nya1-Nyb1-Nyc1-Nyc2) + ynode));
  }
 else
  {
   return (y6init + 
           y6step*((n_p-1)*(yelement-Nya1-Nyb1-Nyc1-Nyc2-Nyb2) + ynode));
  }
}

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////



//=========================================================================
/// \short Constructor for a non-uniform two layer spine mesh, with 
/// element layout in the lower fluid reflected in the upper. Three
/// distinct y regions need numbers of element specified and two 
/// x regions.
//=========================================================================
template <class ELEMENT, class INTERFACE_ELEMENT >
Axisym3x8TwoLayerSpineMesh<ELEMENT,INTERFACE_ELEMENT>::
Axisym3x8TwoLayerSpineMesh(const unsigned &nxa, const unsigned &nxb, 
                           const unsigned &nxc, 
                           const unsigned &nya1, const unsigned &nyb1, 
                           const unsigned &nyc1, const unsigned &nyd1, 
                           const unsigned &nya2, const unsigned &nyb2, 
                           const unsigned &nyc2, const unsigned &nyd2, 
                           const double &lx, 
                           const double &h1, const double &h2,
                           TimeStepper* time_stepper_pt) :
 TwoLayerSpineMesh<ELEMENT, INTERFACE_ELEMENT >()
{
 // Mesh can only be built with 2D Qelements.
 MeshChecker::assert_geometric_element<QElementGeometricBase,ELEMENT>(2);

 // We've called the "generic" constructor for the RectangularQuadMesh
 // which doesn't do much...
 // Now set up the parameters that characterise the mesh geometry
 // then call the constructor

 Nxa = nxa; Nxb = nxb; Nxc = nxc; 
 Nya1 = nya1; Nyb1 = nyb1; Nyc1 = nyc1; Nyd1 = nyd1;
 Nya2 = nya2; Nyb2 = nyb2; Nyc2 = nyc2; Nyd2 = nyd2;
 Xfraction1 = 0.2;
 Xfraction2 = 0.8;
 Yfraction1 = 0.2;
 Yfraction2 = 0.8;
 Yfraction3 = 0.93;
 
 // Check validity of Xfraction, Yfraction1 and Yfraction2.
 if (Yfraction1 < 0.0 || Yfraction1>Yfraction2)
  {
   throw OomphLibError("Invalid Yfraction1",
                       "AxiSymm3x8TwoLayerSpineMesh",
                       OOMPH_EXCEPTION_LOCATION);
  }
 if (Yfraction2 < Yfraction1 || Yfraction2>Yfraction3)
  {
    throw OomphLibError("Invalid Yfraction2",
                       "AxiSymm3x8TwoLayerSpineMesh",
                       OOMPH_EXCEPTION_LOCATION);
  }
 if (Yfraction3 < Yfraction2 || Yfraction3>1.0)
  {
   throw OomphLibError("Invalid Yfraction3",
                       "AxiSymm3x8TwoLayerSpineMesh",
                       OOMPH_EXCEPTION_LOCATION);
  }
 if (Xfraction1<0.0 || Xfraction1>Xfraction2)
  {
   throw OomphLibError("Invalid Xfraction1",
                       "AxiSymm3x8TwoLayerSpineMesh",
                       OOMPH_EXCEPTION_LOCATION);
  }
 if (Xfraction2 < Xfraction1 || Xfraction2>1.0)
  {
   throw OomphLibError("Invalid Xfraction2",
                       "AxiSymm3x8TwoLayerSpineMesh",
                       OOMPH_EXCEPTION_LOCATION);
  }

 // Number of elements in x direction
 RectangularQuadMesh<ELEMENT >::Nx = Nxa+Nxb+Nxc; 
 
 // Number of elements in bottom and top layers
 this->Ny1 = Nya1+Nyb1+Nyc1+Nyd1;
 this->Ny2 = Nya2+Nyb2+Nyc2+Nyd2;
 
 // Number of elements in y direction
 RectangularQuadMesh<ELEMENT >::Ny = this->Ny1 + this->Ny2;
 
 // Min. x coordinate
 RectangularQuadMesh<ELEMENT >::Xmin = 0.0; 
 
 // Max. x coordinate
 RectangularQuadMesh<ELEMENT >::Xmax = lx; 
 
 // Min. y coordinate
 RectangularQuadMesh<ELEMENT >::Ymin = 0.0; 
 
 // Max. y coordinate
 RectangularQuadMesh<ELEMENT >::Ymax = h1+h2;
 
 // Periodic?
 RectangularQuadMesh<ELEMENT >::Xperiodic = false;
 
 // Set height of upper and lower layers
 this->H1 = h1;
 this->H2 = h2;
 
 // Now build the mesh: 
 this->build_two_layer_mesh(time_stepper_pt);
 
}


//=========================================================================
/// \short Constructor for a non-uniform two layer spine mesh, with 
/// element layout in the lower fluid reflected in the upper. Three
/// distinct y regions need numbers of element specified and two 
/// x regions. The fractions of these regions are specified in this 
/// constructor.
//=========================================================================
template <class ELEMENT, class INTERFACE_ELEMENT >
Axisym3x8TwoLayerSpineMesh<ELEMENT,INTERFACE_ELEMENT>::
Axisym3x8TwoLayerSpineMesh(const unsigned &nxa, const unsigned &nxb,
                           const unsigned &nxc,
                           const double &Xfrac1, const double &Xfrac2,
                           const unsigned &nya1, const unsigned &nyb1,
                           const unsigned &nyc1, const unsigned &nyd1, 
                           const unsigned &nya2, const unsigned &nyb2,
                           const unsigned &nyc2, const unsigned &nyd2, 
                           const double &Yfrac1, const double &Yfrac2,
                           const double &Yfrac3, 
                           const double &lx, 
                           const double &h1, const double &h2,
                           TimeStepper* time_stepper_pt) :
 TwoLayerSpineMesh<ELEMENT,INTERFACE_ELEMENT>()
{
 // We've called the "generic" constructor for the TwoLayerSpineMesh
 // which doesn't do much...
 // Now set up the parameters that characterise the mesh geometry
 // then call the constructor

 Nxa = nxa; Nxb = nxb; Nxc = nxc; 
 Nya1 = nya1; Nyb1 = nyb1; Nyc1 = nyc1; Nyd1 = nyd1;
 Nya2 = nya2; Nyb2 = nyb2; Nyc2 = nyc2; Nyd2 = nyd2;
 Xfraction1 = Xfrac1;
 Xfraction2 = Xfrac2;
 Yfraction1 = Yfrac1;
 Yfraction2 = Yfrac2;
 Yfraction3 = Yfrac3;

 // Check validity of Xfraction, Yfraction1 and Yfraction2.
 if (Yfraction1 < 0.0 || Yfraction1>Yfraction2)
  {
   throw OomphLibError("Invalid Yfraction1",
                       "AxiSymm3x8TwoLayerSpineMesh",
                       OOMPH_EXCEPTION_LOCATION);
  }
 if (Yfraction2 < Yfraction1 || Yfraction2>Yfraction3)
  {
   throw OomphLibError("Invalid Yfraction2",
                       "AxiSymm3x8TwoLayerSpineMesh",
                       OOMPH_EXCEPTION_LOCATION);
  }
 if (Yfraction3 < Yfraction2 || Yfraction3>1.0)
  {
   throw OomphLibError("Invalid Yfraction3",
                       "AxiSymm3x8TwoLayerSpineMesh",
                       OOMPH_EXCEPTION_LOCATION);
  }
 if (Xfraction1<0.0 || Xfraction1>Xfraction2)
  {
    throw OomphLibError("Invalid Xfraction1",
                       "AxiSymm3x8TwoLayerSpineMesh",
                       OOMPH_EXCEPTION_LOCATION);
  }
 if (Xfraction2 < Xfraction1 || Xfraction2>1.0)
  {
   throw OomphLibError("Invalid Xfraction2",
                       "AxiSymm3x8TwoLayerSpineMesh",
                       OOMPH_EXCEPTION_LOCATION);
  }
 
 // Number of elements in x direction
 RectangularQuadMesh<ELEMENT >::Nx = Nxa+Nxb+Nxc; 
 
 // Number of elements in bottom and top layers
 this->Ny1 = Nya1+Nyb1+Nyc1+Nyd1;
 this->Ny2 = Nya1+Nyb1+Nyc1+Nyd1;
 
 // Number of elements in y direction
 RectangularQuadMesh<ELEMENT >::Ny = this->Ny1 + this->Ny2;
 
 // Min. x coordinate
 RectangularQuadMesh<ELEMENT >::Xmin = 0.0; 
 
 // Max. x coordinate
 RectangularQuadMesh<ELEMENT >::Xmax = lx; 
 
 // Min. y coordinate
 RectangularQuadMesh<ELEMENT >::Ymin = 0.0; 
 
 // Max. y coordinate
 RectangularQuadMesh<ELEMENT >::Ymax = h1+h2;
 
 // Periodic?
 RectangularQuadMesh<ELEMENT >::Xperiodic = false;
 
 // Set height of upper and lower layers
 this->H1 = h1;
 this->H2 = h2;

 // Now build the mesh: 
 this->build_two_layer_mesh(time_stepper_pt);
 
}


/// \short The spacing function for the x co-ordinates with two 
/// regions.
template <class ELEMENT, class INTERFACE_ELEMENT >
double Axisym3x8TwoLayerSpineMesh<ELEMENT,INTERFACE_ELEMENT>::
x_spacing_function(unsigned xelement, unsigned xnode,
                   unsigned yelement, unsigned ynode)
 {
  //Set up some spacing parameters
  //region a starts at Xmin
  double Xmin =  RectangularQuadMesh<ELEMENT >::Xmin;
  //region c ends at Xmax
  double Xmax =  RectangularQuadMesh<ELEMENT >::Xmax;
  //Number of nodes per element
  unsigned n_p =  RectangularQuadMesh<ELEMENT >::Np;
  //region a starts at Xmin
  double x1init = Xmin;
  //region b starts at Xmin + Xfraction1*(Xmax-Xmin)
  double x2init = Xmin + Xfraction1*(Xmax-Xmin);
  //region c starts at Xmin + Xfraction2*(Xmax-Xmin)
  double x3init = Xmin + Xfraction2*(Xmax-Xmin);
  //Calculate the spacing between the nodes in each region
  //Assuming uniform spacing
  //region a has a length Xfraction1*(Xmax-Xmin)
  double x1step = Xfraction1*(Xmax-Xmin)/((n_p-1)*Nxa);
  //region b has a length (Xfraction2-Xfraction1)*(Xmax-Xmin)
  double x2step = (Xfraction2-Xfraction1)*(Xmax-Xmin)/((n_p-1)*Nxb);
  //region c has a length (1.0-Xfraction2)*(Xmax-Xmin)
  double x3step = (1.0-Xfraction2)*(Xmax-Xmin)/((n_p-1)*Nxc);
  
  //Now set up the particular spacing 
  //(it's different in the two different regions)
  if(xelement < Nxa)
   {
    return (x1init + x1step*((n_p-1)*xelement + xnode));
   }
  else
   {
    if(xelement < Nxa+Nxb)
     {
      return (x2init + x2step*((n_p-1)*(xelement-Nxa) + xnode));
     }
    else
     {
      return (x3init + x3step*((n_p-1)*(xelement-Nxa-Nxb) + xnode));
     }
   }
 }

/// \short The spacing function for the y co-ordinates with three
/// regions in each fluid.
template <class ELEMENT, class INTERFACE_ELEMENT >
double Axisym3x8TwoLayerSpineMesh<ELEMENT,INTERFACE_ELEMENT>::
y_spacing_function(unsigned xelement, unsigned xnode,
                   unsigned yelement, unsigned ynode)
{
 //Set up some spacing parameters
 //The lower region a starts at Ymin
 double Ymin = RectangularQuadMesh<ELEMENT >::Ymin;
 //The interface is at Ymid
 double Ymid = this->H1;
 //The upper region a ends at Ymax
 double Ymax = RectangularQuadMesh<ELEMENT >::Ymax;
 //Number of nodes per element
 unsigned n_p = RectangularQuadMesh<ELEMENT >::Np;
 //The lower region a starts at Ymin
 double y1init = Ymin;
 //The lower region b starts at Ymin + Yfraction1*(Ymid-Ymin)
 double y2init = Ymin + Yfraction1*(Ymid-Ymin);
 //The lower region c starts at Ymin + Yfraction2*(Ymid-Ymin)
 double y3init = Ymin + Yfraction2*(Ymid-Ymin);
 //The lower region d starts at Ymin + Yfraction3*(Ymid-Ymin)
 double y4init = Ymin + Yfraction3*(Ymid-Ymin);
 //The upper region d starts at Ymid
 double y5init = Ymid;
 //The upper region c starts at Ymax - Yfraction3*(Ymax-Ymid)
 double y6init = Ymax - Yfraction3*(Ymax-Ymid);
 //The upper region b starts at Ymax - Yfraction2*(Ymax-Ymid)
 double y7init = Ymax - Yfraction2*(Ymax-Ymid);
 //The upper region a starts at Ymax - Yfraction1*(Ymax-Ymid)
 double y8init = Ymax - Yfraction1*(Ymax-Ymid);
 //Calculate the space between each node in each region,
 //Assumming uniform spacing
 //Lower region a has a length Yfraction1*(Ymid-Ymin)
 double y1step = Yfraction1*(Ymid-Ymin)/((n_p-1)*Nya1);
 //Lower region b has a length (Yfraction2-Yfraction1)*(Ymid-Ymin)
 double y2step = (Yfraction2-Yfraction1)*(Ymid-Ymin)/((n_p-1)*Nyb1);
 //Lower region c has a length (Yfraction3-Yfraction2)*(Ymin-Ymin)
 double y3step = (Yfraction3-Yfraction2)*(Ymid-Ymin)/((n_p-1)*Nyc1);
 //Lower region d has length (1.0-Yfraction3)*(Ymin-Ymin)
 double y4step = (1.0-Yfraction3)*(Ymid-Ymin)/((n_p-1)*Nyd1);
 //Upper region d has a length = Lower region d
 double y5step = (1.0-Yfraction3)*(Ymax-Ymid)/((n_p-1)*Nyd2);;
 //Upper region c has a length = lower region c
 double y6step = (Yfraction3-Yfraction2)*(Ymax-Ymid)/((n_p-1)*Nyc2);
 //Upper region b has a length = lower region b
 double y7step = (Yfraction2-Yfraction1)*(Ymax-Ymid)/((n_p-1)*Nyb2);
 //Upper region a has a length = lower region a
 double y8step = Yfraction1*(Ymax-Ymid)/((n_p-1)*Nya2);
 
 //Now return the actual node position, it's different in the two
 //regions, of course
 if(yelement < Nya1) 
  {
   return (y1init + y1step*((n_p-1)*yelement + ynode));
  }
 else if (yelement < Nya1+Nyb1)
  {
   return (y2init + y2step*((n_p-1)*(yelement-Nya1) + ynode));
  }
 else if (yelement < Nya1+Nyb1+Nyc1)
  {
   return (y3init + y3step*((n_p-1)*(yelement-Nya1-Nyb1) + ynode));
  }
 else if (yelement < Nya1+Nyb1+Nyc1+Nyd1)
  {
   return (y4init + y4step*((n_p-1)*(yelement-Nya1-Nyb1-Nyc1) + ynode));
  }
 else if (yelement < Nya1+Nyb1+Nyc1+Nyd1+Nyd2)
  {
   return (y5init + y5step*((n_p-1)*(yelement-Nya1-Nyb1-Nyc1-Nyd1) + ynode));
  }
 else if (yelement < Nya1+Nyb1+Nyc1+Nyd1+Nyd2+Nyc2)
  {
   return (y6init + 
           y6step*((n_p-1)*(yelement-Nya1-Nyb1-Nyc1-Nyd1-Nyd2) + ynode));
  }
 else if (yelement < Nya1+Nyb1+Nyc1+Nyd1+Nyd2+Nyc2+Nyb2)
  {
   return (y7init + 
           y7step*((n_p-1)*(yelement-Nya1-Nyb1-Nyc1-Nyd1-Nyd2-Nyc2) + ynode));
  }
 else
  {
   return (y8init + 
           y8step*((n_p-1)*(yelement-Nya1-Nyb1-Nyc1-Nyd1-Nyd2-Nyc2-Nyb2) + 
                   ynode));
  }
}

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

}
#endif

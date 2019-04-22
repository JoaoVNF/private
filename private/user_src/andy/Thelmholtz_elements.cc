//Non-inline functions for Helmholtz elements
#include "Thelmholtz_elements.h"


namespace oomph
{



/////////////////////////////////////////////////////////////////////////
// THelmholtzElement
/////////////////////////////////////////////////////////////////////////



//======================================================================
// Set the data for the number of Variables at each node, always 1
//======================================================================
template<unsigned DIM, unsigned NNODE_1D>
const unsigned THelmholtzElement<DIM,NNODE_1D>::Initial_Nvalue = 1;

//====================================================================
// Force build of templates
//====================================================================
template class THelmholtzElement<1,2>;
template class THelmholtzElement<1,3>;
template class THelmholtzElement<1,4>;

template class THelmholtzElement<2,2>;
template class THelmholtzElement<2,3>;
template class THelmholtzElement<2,4>;

template class THelmholtzElement<3,2>;
template class THelmholtzElement<3,3>;

}

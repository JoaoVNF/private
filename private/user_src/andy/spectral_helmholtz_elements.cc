//Non-inline functions and static data for spectral helmholtz elements
#include "spectral_helmholtz_elements.h"


namespace oomph
{

//==========================================================================
// Static value that returns the number of variables at each node, always 1
//==========================================================================
template<unsigned DIM, unsigned NNODE_1D>
const unsigned QSpectralHelmholtzElement<DIM,NNODE_1D>::Initial_Nvalue = 1;

template class QSpectralHelmholtzElement<1,2>;
template class QSpectralHelmholtzElement<1,3>;
template class QSpectralHelmholtzElement<1,4>;
template class QSpectralHelmholtzElement<1,5>;
template class QSpectralHelmholtzElement<1,6>;
template class QSpectralHelmholtzElement<1,7>;
template class QSpectralHelmholtzElement<1,8>;
template class QSpectralHelmholtzElement<1,9>;

template class QSpectralHelmholtzElement<2,2>;
template class QSpectralHelmholtzElement<2,3>;
template class QSpectralHelmholtzElement<2,4>;
template class QSpectralHelmholtzElement<2,5>;
template class QSpectralHelmholtzElement<2,6>;
template class QSpectralHelmholtzElement<2,7>;


template class QSpectralHelmholtzElement<3,2>;
template class QSpectralHelmholtzElement<3,3>;
template class QSpectralHelmholtzElement<3,4>;
template class QSpectralHelmholtzElement<3,5>;
template class QSpectralHelmholtzElement<3,6>;
template class QSpectralHelmholtzElement<3,7>;

}

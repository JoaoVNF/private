#ifndef OOMPH_HELLOWORLD_HEADER
#define OOMPH_HELLOWORLD_HEADER

// Config header generated by autoconfig
#ifdef HAVE_CONFIG_H
  #include <oomph-lib-config.h>
#endif

// oomph-lib headers
#include "../generic/matrices.h"


namespace oomph
{




//=================================================================
/// Hello world class -- uses object from generic library
//=================================================================
class HelloWorld
{

public:


 /// Constructor
 HelloWorld()
  {
   oomph_info << "In HelloWorld constructor" << std::endl;
   Some_matrix.resize(2);
  }

 /// Say it inline
 void say_it_inline()
  {
   oomph_info << "Inline: Hello world! " << Some_matrix.nrow() << std::endl;
   
  }

 /// Say it externally
 void say_it_external();

  private:

 /// Private full matrix
 DenseMatrix<double> Some_matrix;


};


}

#endif
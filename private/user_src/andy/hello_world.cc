#include "hello_world.h"


namespace oomph
{

 /// Say it externally
 void AndysHelloWorld::say_it_external()
  {
   oomph_info << "AndysHelloWorld object says this: Hello world! " 
             << Some_matrix.nrow() << " from a compiled function! " 
             << std::endl;
  }

}

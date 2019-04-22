#include "hello_world.h"


namespace oomph
{

 /// Say it externally
 void JacksHelloWorld::say_it_external()
  {
   oomph_info << "JacksHelloWorld object says: Hello world! " 
             << Some_matrix.nrow() << " from a compiled function! " 
             << std::endl;
  }

}

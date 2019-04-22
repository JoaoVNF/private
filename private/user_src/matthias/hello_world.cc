#include "hello_world.h"


namespace oomph
{

 /// Say it externally
 void HelloWorld::say_it_external()
  {
   oomph_info << "External: HELLO WORLD! " << Some_matrix.nrow() << std::endl;
  }

}

#ifndef LMESH2OOMPH_D_HEADER
#define LMESH2OOMPH_D_HEADER

#include "line_mesh.h"
#include "../generic/refineable_mesh.h"
#include "refineable_line_element.h"

namespace oomph
{

 //===========================================================================
 /// Intermediate mesh class that implements the mesh adaptation functions
 /// specified in the RefineableMesh class for meshes that contain the
 /// refineable variant of QElement s [The class ELEMENT provided as the
 /// template parameter must be of type RefineableQElement<1>].
 /// 
 /// Mesh adaptation/refinement is implemented by BinaryTree procedures
 /// and any concrete implementation of this class needs to provide a
 /// BinaryTreeForest representation of the initial (coarse) mesh.
 //===========================================================================
 template <class ELEMENT>
  class RefineableLineMesh : public virtual RefineableMesh<ELEMENT>, 
  public virtual LineMeshBase
  {
   
   public:
   
   /// Constructor: Set up static binary tree data
   RefineableLineMesh()
    {
     // BinaryTree static data needs to be setup before binary tree-based
     // mesh refinement works
     BinaryTree::setup_static_data();
    }
   
   /// Broken copy constructor
   RefineableLineMesh(const RefineableLineMesh& dummy) 
    { 
     BrokenCopy::broken_copy("RefineableLineMesh");
    } 
   
   /// Broken assignment operator
   void operator=(const RefineableLineMesh&) 
    {
     BrokenCopy::broken_assign("RefineableLineMesh");
    }
   
   /// Destructor:
   virtual ~RefineableLineMesh() {}
   
   /// \short Set up the tree forest associated with the Mesh. 
   /// Forwards call to setup_binary_tree_forest().
   virtual void setup_tree_forest()
    {
     setup_binary_tree_forest();
    }
   
   /// Set up BinaryTreeForest. Wipes any existing tree structure and
   /// regards the currently active elements as the root trees in the forest.
   void setup_binary_tree_forest()
    {
     // This wipes all elements/binary trees in the tree representation
     // but leaves the leaf elements alone
     if (this->Forest_pt!=0) delete this->Forest_pt;
     
     // Each finite element in the coarse base mesh gets associated with
     // (the root of) a BinaryTree. Store BinaryTreeRoots in vector:
     Vector<TreeRoot*> trees_pt;
     
     // Determine number of elements in mesh
     const unsigned n_element=this->nelement();
     
     // Loop over all elements, build corresponding BinaryTree and store
     // BinaryTreeRoots in vector:
     for (unsigned e=0;e<n_element;e++)
      {
       // Get pointer to full element type 
       ELEMENT* el_pt = dynamic_cast<ELEMENT*>(this->element_pt(e));
       
       // Build associated binary tree(root) -- pass pointer to corresponding
       // finite element and add the pointer to vector of binary tree (roots):
       trees_pt.push_back(new BinaryTreeRoot(el_pt));
      } 
     
     // Plant BinaryTreeRoots in BinaryTreeForest
     this->Forest_pt = new BinaryTreeForest(trees_pt);
    }
   
  };
 
} // End of namespace

#endif

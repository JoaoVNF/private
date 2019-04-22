// Domain partitioning

#ifndef OOMPH_PARTITIONING_HEADER 
#define OOMPH_PARTITIONING_HEADER 


// Config header generated by autoconfig
#ifdef HAVE_CONFIG_H
  #include <oomph-lib-config.h>
#endif


// ooomph-lib includes
#include "../generic/Vector.h"
#include "../generic/problem.h"


namespace oomph
{

//==================================================================
// Interfaces to METIS functions
//==================================================================
extern "C" 
{
/// \short Metis graph partitioning function -- decomposes 
/// nodal graph based on minimum edgecut
 void METIS_PartGraphKway(int *, int *, int *, int *, int *, 
                          int *, int *, int *, int *, int *, int *);

/// \short Metis graph partitioning function -- decomposes 
/// nodal graph based on minimum communication volume
void METIS_PartGraphVKway(int *, int *, int *, int *, int *,
                          int *, int *, int *, int *, int *, int *);
}








//====================================================================
/// Namespace for METIS graph partitioning routines
//====================================================================
namespace METIS
{
 
 /// \short Partition mesh uniformly by dividing elements
 /// equally over the partitions, in the order
 /// in which they are returned by problem.
 /// On return, element_domain[ielem] contains the number
 /// of the domain [0,1,...,ndomain-1] to which 
 /// element ielem has been assigned.
 extern void uniform_partition_mesh(Problem* problem_pt,
                                    const unsigned& ndomain,
                                    Vector<unsigned>& element_domain);

 
 /// \short Use METIS to assign each element to a domain.
 /// On return, element_domain[ielem] contains the number
 /// of the domain [0,1,...,ndomain-1] to which 
 /// element ielem has been assigned.
 /// - objective=0: minimise edgecut.
 /// - objective=1: minimise total communications volume.
 /// .
 /// Partioning is based on nodal graph of mesh.
  extern void partition_mesh(Problem* problem_pt,
                             const unsigned& ndomain,
                             const unsigned& objective,
                             Vector<unsigned>& element_domain);


 /// \short Use METIS to assign each element to a domain.
 /// On return, element_domain[ielem] contains the number
 /// of the domain [0,1,...,ndomain-1] to which 
 /// element ielem has been assigned.
 /// - objective=0: minimise edgecut.
 /// - objective=1: minimise total communications volume.
 /// .
 /// Partioning is based on nodal graph of mesh.
  extern void partition_mesh(Mesh* mesh_pt,
                             const unsigned& ndomain,
                             const unsigned& objective,
                             Vector<unsigned>& element_domain);

//  /// \short Use METIS to assign each element to a domain.
//  /// On return, element_domain[ielem] contains the number
//  /// of the domain [0,1,...,ndomain-1] to which 
//  /// element ielem has been assigned.
//  /// - objective=0: minimise edgecut.
//  /// - objective=1: minimise total communications volume.
//  /// .
//  /// Partioning is based on "Data" graph of mesh.
//  extern void partition_mesh_data(Problem* problem_pt,
//                                  const unsigned& ndomain,
//                                  const unsigned& objective,
//                                 Vector<unsigned>& element_domain);


}


}

#endif

#include "partitioning.h"
#include "mesh.h"


namespace oomph
{

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
 void uniform_partition_mesh(Problem* problem_pt,
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
 void partition_mesh(Problem* problem_pt,
                     const unsigned& ndomain,
                     const unsigned& objective,
                     Vector<unsigned>& element_domain)
 {
  partition_mesh(problem_pt->mesh_pt(),
                 ndomain,
                 objective,
                 element_domain);
 }
 

 /// \short Use METIS to assign each element to a domain.
 /// On return, element_domain[ielem] contains the number
 /// of the domain [0,1,...,ndomain-1] to which 
 /// element ielem has been assigned.
 /// - objective=0: minimise edgecut.
 /// - objective=1: minimise total communications volume.
 /// .
 /// Partioning is based on nodal graph of mesh.
 void partition_mesh(Mesh* mesh_pt,
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
//  void partition_mesh_data(Problem* problem_pt,
//                           const unsigned& ndomain,
//                           const unsigned& objective,
//                           Vector<unsigned>& element_domain);

}





//==================================================================
/// Partition mesh uniformly by dividing elements
/// equally over the partitions, in the order
/// in which they are returned by problem.
/// On return, element_domain[ielem] contains the number
/// of the domain [0,1,...,ndomain-1] to which 
/// element ielem has been assigned.
//==================================================================
void METIS::uniform_partition_mesh(Problem* problem_pt,
                                   const unsigned& ndomain,
                                   Vector<unsigned>& element_domain)
{
 
 // Number of elements
 unsigned nelem=problem_pt->mesh_pt()->nelement();
 
#ifdef PARANOID
 if (nelem!=element_domain.size())
  {
   std::ostringstream error_stream;
   error_stream << "element_domain Vector has wrong length " 
                << nelem << " " << element_domain.size() << std::endl;

   throw OomphLibError(error_stream.str(),
                       "METIS::uniform_partition_mesh()",
                       OOMPH_EXCEPTION_LOCATION);
  }
#endif
 
 
 
 // Uniform partitioning
 unsigned nel_per_domain=int(float(nelem)/float(ndomain));
 for (unsigned ielem=0;ielem<nelem;ielem++)
  {
   unsigned idomain=unsigned(float(ielem)/float(nel_per_domain));
   element_domain[ielem]=idomain;
  }
 
}





//==================================================================
/// Use METIS to assign each element to a domain.
/// On return, element_domain[ielem] contains the number
/// of the domain [0,1,...,ndomain-1] to which 
/// element ielem has been assigned.
/// - objective=0: minimise edgecut.
/// - objective=1: minimise total communications volume.
/// .
/// Partioning is based on nodal graph of mesh.
//==================================================================
void METIS::partition_mesh(Mesh* mesh_pt, const unsigned& ndomain,
                           const unsigned& objective,
                           Vector<unsigned>& element_domain)
{

 // Number of elements
 unsigned nelem=mesh_pt->nelement();

#ifdef PARANOID
 if (nelem!=element_domain.size())
  {
   std::ostringstream error_stream;
   error_stream << "element_domain Vector has wrong length " 
                << nelem << " " << element_domain.size() << std::endl;
   
   throw OomphLibError(error_stream.str(),
                       "METIS::partition_mesh()",
                       OOMPH_EXCEPTION_LOCATION);
  }
#endif

 // Setup nodal graph
 //------------------

 // Start timer
 clock_t cpu_start=clock();

 // Use map to associate node with number
 std::map<Node*,unsigned> node_number;

 unsigned nnode=mesh_pt->nnode();
 for (unsigned inode=0;inode<nnode;inode++)
  {
   Node* node_pt=mesh_pt->node_pt(inode);
   node_number[node_pt]=inode;
  }

  
 // Vector of sets which store the nodes that are connected with a given node
 Vector<std::set<unsigned> > connected_nodes;
 connected_nodes.resize(nnode);


 // Loop over all elements
 for (unsigned ielem=0;ielem<nelem;ielem++)
  {
   // Loop over all element-related nodes 
   unsigned nnode=mesh_pt->finite_element_pt(ielem)->nnode();
   for (unsigned inode=0;inode<nnode;inode++)
    {
     // Get pointer to node
     Node* inode_pt=mesh_pt->
      finite_element_pt(ielem)->node_pt(inode);

     // Get global node number:
     unsigned global_inode=node_number[inode_pt];
      
     // Now loop over all other nodes in the element -- they are
     // connected
     for (unsigned jnode=0;jnode<nnode;jnode++)
      {
       // Get pointer to node
       Node* jnode_pt=mesh_pt->
        finite_element_pt(ielem)->node_pt(jnode);

       // Get global node number:
       unsigned global_jnode=node_number[jnode_pt];

       // Don't store self-references
       if (inode!=jnode)
        {
         // Node jnode is connected with node inode:
         connected_nodes[global_inode].insert(global_jnode);

         //...and vice versa (might not be but we need a symmetric
         // graph in METIS. Nobody (?) seems to know any better
         // ways for dealing with structurally unsymmetric matrices.
         connected_nodes[global_jnode].insert(global_inode);

        }
      }
    }
  }    


 // Now convert into C-style packed array for interface with METIS
 int* xadj = new int[nnode+1];
 Vector<int> adjacency_vector;

 // Initialise counters
 unsigned ientry=0;

 // Loop over all nodes
 for (unsigned inode=0;inode<nnode;inode++)
  {

   // First entry for current node
   xadj[inode]=ientry;

   // Loop over nodes that are connected to current node
   typedef std::set<unsigned>::iterator IT;
   for (IT it=connected_nodes[inode].begin();
        it!=connected_nodes[inode].end();it++)
    {
     // Copy into adjacency array
     adjacency_vector.push_back(*it);

     // We've just made another entry
     ientry++;
    }

   // Entry after last entry for current node:
   xadj[inode+1]=ientry;
    
  }

 // Number of edges in graph 
 unsigned nedges=ientry;


 // Move into C-style arrayfor interface with METIS
 int* adjcncy = new int[nedges];

 for (unsigned i=0;i<nedges;i++)
  {
   adjcncy[i]=adjacency_vector[i];
  }
  

 // End timer
 clock_t cpu_end=clock();

 // Doc
 double cpu0=double(cpu_end-cpu_start)/CLOCKS_PER_SEC;
 oomph_info  
  << "CPU time for setup of METIS data structures            [nelem="
  << nelem <<"]: " 
  << cpu0
  << " sec" << std::endl;


 // Call METIS graph partitioner
 //-----------------------------

 // Start timer
 cpu_start=clock();

 // Number of vertices in graph
 int nvertex=nnode;

 // No vertex weights 
 int* vwgt=0;

 // No edge weights
 int* adjwgt=0;

 // Flag indicating that graph isn't weighted: 0; vertex weights only: 2
 int wgtflag=0;

 // Use C-style numbering (first array entry is zero)
 int numflag=0; 

 // Number of desired partitions
 int nparts=ndomain;

 // Use default options
 int* options= new int[10];
 options[0]=0;

 // Number of cut edges in graph
 int* edgecut = new int[nnode];

 // Array containing the partition information
 int* part = new int[nnode];


 // Call partitioner
 if (objective==0)
  {
   // Partition with the objective of minimising the edge cut
   METIS_PartGraphKway(&nvertex, xadj,adjcncy,vwgt,adjwgt, 
                       &wgtflag, &numflag, &nparts, options, edgecut, part);
  }
 else if (objective==1)
  {
   // Partition with the objective of minimising the total communication 
   // volume  
   METIS_PartGraphVKway(&nvertex, xadj, adjcncy,vwgt, adjwgt, 
                        &wgtflag, &numflag, &nparts, options, edgecut, part);
  }
 else
  {
   std::ostringstream error_stream;
   error_stream 
    << "Wrong objective for METIS. objective = " << objective << std::endl;

   throw OomphLibError(error_stream.str(),"METIS::partition_mesh()",
                       OOMPH_EXCEPTION_LOCATION);
  }
   
 
 // End timer
 cpu_end=clock();

 // Doc
 double cpu1=double(cpu_end-cpu_start)/CLOCKS_PER_SEC;
 oomph_info  
  << "CPU time for METIS mesh partitioning                   [nelem="
  << nelem <<"]: " 
  << cpu1
  << " sec" << std::endl;

 // Now turn into mesh partitioning
 //--------------------------------

 // Start timer
 cpu_start=clock();

 // Vector holding the number of elements assigned to a given partition
 Vector<unsigned> nelements_in_partition(ndomain,0);
 
 // Boolean STL vector to check if element has been assigned to domain
 std::vector<bool> element_is_assigned(nelem);


 // Vector of maps holding the number of nodes in the elemnt
 //that are contained in a certain
 // partition
 Vector<std::map<unsigned,unsigned> > partition_count(nelem);
  

 // First loop over all elements to find out which partitions their
 // domains are part of
 for (unsigned ielem=0;ielem<nelem;ielem++)
  {

   element_is_assigned[ielem]=false;

   // Loop over all element-related nodes
   unsigned nnode=mesh_pt->finite_element_pt(ielem)->nnode();
   for (unsigned inode=0;inode<nnode;inode++)
    {
     // Get pointer to node
     Node* node_pt=mesh_pt->finite_element_pt(ielem)->
      node_pt(inode);

     // This partition has just received another node
     partition_count[ielem][part[node_number[node_pt]]]++;

    }

   // Check if all nodes of the element are in the same partition
   if (partition_count[ielem].size()==1)
    {
     // Assign it
     element_domain[ielem]=partition_count[ielem].begin()->first;

     // Increment counter for elements in current partition
     nelements_in_partition[partition_count[ielem].begin()->first]++;

     // Label it as assigned
     element_is_assigned[ielem]=true;

    }

  }

 // Loop over all elements again to deal with the ones
 // that stradle two partitions
 for (unsigned ielem=0;ielem<nelem;ielem++)
  {

   if (!element_is_assigned[ielem])
    {

     // Find out which partition has most nodes in common with present element
     unsigned max_part=0;
     unsigned part_max=0;
     typedef std::map<unsigned,unsigned>::iterator IT;
     for (IT it=partition_count[ielem].begin();
          it!=partition_count[ielem].end();it++)
      {
       if (it->second>max_part)
        {
         part_max=it->first;
         max_part=it->second;
        }
      }

     // Check for dead heat:
     unsigned count_max=0;
     for (IT it=partition_count[ielem].begin();
          it!=partition_count[ielem].end();it++)
      {
       // How many partitions have the max. count?
       if (it->second==max_part)
        {
         count_max++;
        }
      }
      
     unsigned partition=0;
     if (count_max==1)
      {
       partition=part_max;
      }
     else
      {
       unsigned fewest_elements=1000000;
       for (IT it=partition_count[ielem].begin();
            it!=partition_count[ielem].end();it++)
        {
         // If this this partition is (one of) the most frequently
         // shared partitions, check if it's got the fewest elements
         // assigned
         if (it->second==max_part)
          {
           if (nelements_in_partition[it->first]<fewest_elements)
            {
             fewest_elements=nelements_in_partition[it->first];
             partition=it->first;
            }
          }
        }
      }

     // Assign to the partition that's shared by most nodes
     // and (if there are multiple ones) that has the fewest elements
     // assigned to it:
     element_domain[ielem]=partition; 

     // Increment counter for elements in current partition
     nelements_in_partition[partition]++;

    }
  }
 
 // End timer
 cpu_end=clock();

 // Doc
 double cpu2=double(cpu_end-cpu_start)/CLOCKS_PER_SEC;
 oomph_info  
  << "CPU time for conversion to oomph-lib mesh partitioning [nelem="
  << nelem <<"]: " 
  << cpu2
  << " sec" << std::endl;


  oomph_info  
  << "CPU time ratio (METIS/oomph-lib)                       [nelem="
  << nelem <<"]: " 
  << cpu1/(cpu0+cpu2)*100.0
  << "% " << std::endl;
 
 

 // Cleanup
 delete [] xadj;
 delete [] adjcncy;
 delete [] part;
 delete [] edgecut;
 delete [] options;

}



//  //=================================================================
//  /// \short Use METIS to assign each element to a domain.
//  /// On return, element_domain[ielem] contains the number
//  /// of the domain [0,1,...,ndomain-1] to which 
//  /// element ielem has been assigned.
//  /// - objective=0: minimise edgecut.
//  /// - objective=1: minimise total communications volume.
//  /// .
//  /// Partioning is based on "Data" graph of mesh.
//  //=================================================================
//  void METIS::partition_mesh_data(Problem* problem_pt,
//                                  const unsigned& ndomain,
//                                  const unsigned& objective,
//                                  Vector<unsigned>& element_domain)
//  {

//   // Number of elements
//   unsigned nelem=problem_pt->mesh_pt()->nelement();

// #ifdef PARANOID
//   if (nelem!=element_domain.size())
//    {
//     cout << "element_domain Vector has wrong length " 
//          << nelem << " " << element_domain.size() << std::endl;
//     assert(false);
//    }
// #endif

//   // Setup "Data" graph
//   //-------------------

//   // Loop over all elements and collect the items of Data they depend on.
//   // Stick them into a set
//   std::set<Data*> all_data_pt;
//   for (unsigned ielem=0;ielem<nelem;ielem++)
//    {
//     // Loop over all element-related items of data 
//     unsigned ndata=problem_pt->mesh_pt()->element_pt(ielem)->ndata();
//     for (unsigned idata=0;idata<ndata;idata++)
//      {
//       // Get pointer to data
//       Data* data_pt=problem_pt->mesh_pt()->element_pt(ielem)->data_pt(idata);
//       all_data_pt.insert(data_pt);
//      }
//    }


//   // Use map to associate Data with number
//   std::map<Data*,unsigned> data_number;

//   unsigned count=0;
//   unsigned ndata=all_data_pt.size();
//   typedef std::set<Data*>::iterator IT;
//   for (IT it=all_data_pt.begin();it!=all_data_pt.end();it++)
//    {
//     data_number[*it]=count;
//     count++;
//    }

//   // Wipe...
//   all_data_pt.clear();

//   // Vector of sets which store the items of data that are connected with 
//   // a given data
//   Vector<std::set<unsigned> > connected_data;
//   connected_data.resize(ndata);


//   // Loop over all elements
//   for (unsigned ielem=0;ielem<nelem;ielem++)
//    {
//     // Loop over all element-related items of data 
//     unsigned ndata=problem_pt->mesh_pt()->element_pt(ielem)->ndata();
//     for (unsigned idata=0;idata<ndata;idata++)
//      {
//       // Get pointer to data
//       Data* idata_pt=problem_pt->mesh_pt()->element_pt(ielem)->data_pt(idata);

//       // Get global data number:
//       unsigned global_idata=data_number[idata_pt];
      
//       // Now loop over all other items of data in the element -- they are
//       // connected
//       for (unsigned jdata=0;jdata<ndata;jdata++)
//        {
//         // Get pointer to data
//         Data* jdata_pt=problem_pt->mesh_pt()->
//          element_pt(ielem)->data_pt(jdata);

//         // Get global data number:
//         unsigned global_jdata=data_number[jdata_pt];

//         // Don't store self-references
//         if (idata!=jdata)
//          {
//           // Data jdata is connected with data idata:
//           connected_data[global_idata].insert(global_jdata);

//           //...and vice versa (might not be but we need a symmetric
//           // graph in METIS. Nobody (?) seems to know any better
//           // ways for dealing with structurally unsymmetric matrices.
//           connected_data[global_jdata].insert(global_idata);

//          }
//        }
//      }
//    }    


//   // Now convert into C-style packed array for interface with METIS
//   int* xadj = new int[ndata+1];
//   Vector<int> adjacency_vector;

//   // Initialise counters
//   unsigned ientry=0;

//   // Loop over all items of data
//   for (unsigned idata=0;idata<ndata;idata++)
//    {

//     // First entry for current data
//     xadj[idata]=ientry;

//     // Loop over items of data that are connected to current data
//     typedef std::set<unsigned>::iterator IT;
//     for (IT it=connected_data[idata].begin();
//          it!=connected_data[idata].end();it++)
//      {
//       // Copy into adjacency array
//       adjacency_vector.push_back(*it);

//       // We've just made another entry
//       ientry++;
//      }

//     // Entry after last entry for current data:
//     xadj[idata+1]=ientry;
    
//    }

//   // Number of edges in graph 
//   unsigned nedges=ientry;


//   // Move into C-style arrayfor interface with METIS
//   int* adjcncy = new int[nedges];

//   for (unsigned i=0;i<nedges;i++)
//    {
//     adjcncy[i]=adjacency_vector[i];
//    }
  
  
//   // Call METIS graph partitioner
//   //-----------------------------

//   // Number of vertices in graph
//   int nvertex=ndata;

//   // No vertex weights 
//   int* vwgt=0;

//   // No edge weights
//   int* adjwgt=0;

//   // Flag indicating that graph isn't weighted: 0; vertex weights only: 2
//   int wgtflag=0;

//   // Use C-style numbering (first array entry is zero)
//   int numflag=0; 

//   // Number of desired partitions
//   int nparts=ndomain;

//   // Use default options
//   int* options= new int[10];
//   options[0]=0;

//   // Number of cut edges in graph
//   int* edgecut = new int[ndata];

//   // Array containing the partition information
//   int* part = new int[ndata];


//   // Call partitioner
//   if (objective==0)
//    {
//     // Partition with the objective of minimising the edge cut
//     METIS_PartGraphKway(&nvertex, xadj,adjcncy,vwgt,adjwgt, 
//                         &wgtflag, &numflag, &nparts, options, edgecut, part);
//    }
//   else if (objective==1)
//    {
//     // Partition with the objective of minimising the total communication 
//     // volume  
//     METIS_PartGraphVKway(&nvertex, xadj, adjcncy,vwgt, adjwgt, 
//                       &wgtflag, &numflag, &nparts, options, edgecut, part);
//    }
//   else
//    {
//     cout << "Wrong objective for METIS. objective = " << objective << std::endl;
//     assert(false);
//    }
   
 

//   // Now turn into mesh partitioning
//   //--------------------------------

//   // Vector holding the number of elements assigned to a given partition
//   Vector<unsigned> nelements_in_partition(ndomain,0);
 
//   // Boolean STL ector to check if element has been assigned to domain
//   vector<bool> element_is_assigned(nelem);


//   // Vector of maps holding the number of items of data in the element
//   // that are contained in a certain partition
//   Vector<std::map<unsigned,unsigned> > partition_count(nelem);
  

//   // First loop over all elements to find out which partitions their
//   // domains are part of
//   for (unsigned ielem=0;ielem<nelem;ielem++)
//    {

//     element_is_assigned[ielem]=false;

//     // Loop over all element-related items of data
//     unsigned ndata=problem_pt->mesh_pt()->element_pt(ielem)->ndata();
//     for (unsigned idata=0;idata<ndata;idata++)
//      {
//       // Get pointer to data
//       Data* data_pt=problem_pt->mesh_pt()->element_pt(ielem)->
//        data_pt(idata);

//       // This partition has just received another data
//       partition_count[ielem][part[data_number[data_pt]]]++;

//      }

//     // Check if all items of data of the element are in the same partition
//     if (partition_count[ielem].size()==1)
//      {
//       // Assign it
//       element_domain[ielem]=partition_count[ielem].begin()->first;

//       // Increment counter for elements in current partition
//       nelements_in_partition[partition_count[ielem].begin()->first]++;

//       // Label it as assigned
//       element_is_assigned[ielem]=true;

//      }

//    }

//   // Loop over all elements again to deal with the ones
//   // that stradle two partitions
//   for (unsigned ielem=0;ielem<nelem;ielem++)
//    {

//     if (!element_is_assigned[ielem])
//      {

//       // Find out which partition has most items of data in common 
//       // with present element
//       unsigned max_part=0;
//       unsigned part_max;
//       typedef std::map<unsigned,unsigned>::iterator IT;
//       for (IT it=partition_count[ielem].begin();
//            it!=partition_count[ielem].end();it++)
//        {
//         if (it->second>max_part)
//          {
//           part_max=it->first;
//           max_part=it->second;
//          }
//        }

//       // Check for dead heat:
//       unsigned count_max=0;
//       for (IT it=partition_count[ielem].begin();
//            it!=partition_count[ielem].end();it++)
//        {
//         // How many partitions have the max. count?
//         if (it->second==max_part)
//          {
//           count_max++;
//          }
//        }
      
//       unsigned partition;
//       if (count_max==1)
//        {
//         partition=part_max;
//        }
//       else
//        {
//         unsigned fewest_elements=1000000;
//         for (IT it=partition_count[ielem].begin();
//              it!=partition_count[ielem].end();it++)
//          {
//           // If this this partition is (one of) the most frequently
//           // shared partitions, check if it's got the fewest elements
//           // assigned
//           if (it->second==max_part)
//            {
//             if (nelements_in_partition[it->first]<fewest_elements)
//              {
//               fewest_elements=nelements_in_partition[it->first];
//               partition=it->first;
//              }
//            }
//          }
//        }

//       // Assign to the partition that's shared by most items of data
//       // and (if there are multiple ones) that has the fewest elements
//       // assigned to it:
//       element_domain[ielem]=partition; 

//       // Increment counter for elements in current partition
//       nelements_in_partition[partition]++;

//      }
//    }


//   // Cleanup
//   delete [] xadj;
//   delete [] adjcncy;
//   delete [] part;
//   delete [] edgecut;
//   delete [] options;

//  }



}

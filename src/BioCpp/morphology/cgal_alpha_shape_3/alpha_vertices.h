#ifndef ALPHA_VERTICES_H
#define ALPHA_VERTICES_H

#include <vector>
#include <list>

#include <BioCpp/pdb/pdb_MODEL_format.h>
namespace BioCpp{
namespace morphology{

template<typename alpha_shape_t, typename atom_t>
void alpha_vertices( typename BioCpp::pdb::model<atom_t>::type& all_info, 
                     double radius, 
                     std::vector<typename alpha_shape_t::Vertex_handle>& surface,
                     typename alpha_shape_t::Mode mode= alpha_shape_t::GENERAL,
                     std::list<typename alpha_shape_t::Classification_type> classification = {alpha_shape_t::REGULAR,alpha_shape_t::EXTERIOR,alpha_shape_t::SINGULAR} ){

  alpha_shape_t as( all_info.begin(), all_info.end(), radius, mode );

  for(typename std::list< typename alpha_shape_t::Classification_type >::iterator it = classification.begin(); it!=classification.end(); ++it){
    as.get_alpha_shape_vertices( std::back_inserter(surface), *it );
  }
  
}

template<typename alpha_shape_t, typename atom_t>
void alpha_vertices( std::vector<atom_t>& all_info, 
                     double radius, 
                     std::vector<typename alpha_shape_t::Vertex_handle>& surface,
                     typename alpha_shape_t::Mode mode= alpha_shape_t::GENERAL,
                     std::list<typename alpha_shape_t::Classification_type> classification = {alpha_shape_t::REGULAR,alpha_shape_t::EXTERIOR,alpha_shape_t::SINGULAR} ){

  alpha_shape_t as( all_info.begin(), all_info.end(), radius, mode );

  for(typename std::list< typename alpha_shape_t::Classification_type >::iterator it = classification.begin(); it!=classification.end(); ++it){
    as.get_alpha_shape_vertices( std::back_inserter(surface), *it );
  }
}

}
}
#endif

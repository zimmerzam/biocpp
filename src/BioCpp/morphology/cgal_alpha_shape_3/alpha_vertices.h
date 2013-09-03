/* ************************************************************************** */
/*                                                                            */
/*    Copyright 2013 Stefano Zamuner                                          */
/*                                                                            */
/*    This file is part of BioCpp.                                            */
/*                                                                            */
/*    BioCpp is free software: you can redistribute it and/or modify          */
/*    it under the terms of the GNU General Public License as published by    */
/*    the Free Software Foundation, either version 3 of the License, or       */
/*    (at your option) any later version.                                     */
/*                                                                            */
/*    BioCpp is distributed in the hope that it will be useful,               */
/*    but WITHOUT ANY WARRANTY; without even the implied warranty of          */
/*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           */
/*    GNU General Public License for more details.                            */
/*                                                                            */
/*    You should have received a copy of the GNU General Public License       */
/*    along with BioCpp.  If not, see <http://www.gnu.org/licenses/>.         */
/*                                                                            */
/* ************************************************************************** */

#ifndef ALPHA_VERTICES_H
#define ALPHA_VERTICES_H

#include <vector>
#include <list>

#include <BioCpp/pdb/pdb_MODEL_format.h>
namespace BioCpp{
namespace morphology{

template<typename alpha_shape_t, typename atom_t>
void alpha_surface( typename BioCpp::pdb::model<atom_t>::type& all_info, 
                     double radius, 
                     std::vector<typename alpha_shape_t::Vertex_handle>& surface,
                     typename alpha_shape_t::Mode mode= alpha_shape_t::GENERAL,
                     std::list<typename alpha_shape_t::Classification_type> classification = 
                     std::list<typename alpha_shape_t::Classification_type>{alpha_shape_t::REGULAR,alpha_shape_t::EXTERIOR,alpha_shape_t::SINGULAR} ){

  alpha_shape_t as( all_info.begin(), all_info.end(), radius, mode );

  for(typename std::list< typename alpha_shape_t::Classification_type >::iterator it = classification.begin(); it!=classification.end(); ++it){
    as.get_alpha_shape_vertices( std::back_inserter(surface), *it );
  }
  
}

template<typename alpha_shape_t, typename atom_t>
void alpha_surface( std::vector<atom_t>& all_info, 
                     double radius, 
                     std::vector<typename alpha_shape_t::Vertex_handle>& surface,
                     typename alpha_shape_t::Mode mode= alpha_shape_t::GENERAL,
                     std::list<typename alpha_shape_t::Classification_type> classification = 
                     std::list<typename alpha_shape_t::Classification_type>{alpha_shape_t::REGULAR,alpha_shape_t::EXTERIOR,alpha_shape_t::SINGULAR} ){

  alpha_shape_t as( all_info.begin(), all_info.end(), radius, mode );

  for(typename std::list< typename alpha_shape_t::Classification_type >::iterator it = classification.begin(); it!=classification.end(); ++it){
    as.get_alpha_shape_vertices( std::back_inserter(surface), *it );
  }
}

}
}
#endif

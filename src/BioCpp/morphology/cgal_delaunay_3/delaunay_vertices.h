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

#ifndef DELAUNAY_VERTICES_H
#define DELAUNAY_VERTICES_H

#include <vector>
#include <set>

#include <BioCpp/pdb/pdb_MODEL_format.h>
namespace BioCpp{
namespace morphology{

template<typename delaunay_t, typename atom_t>
void delaunay_surface( typename BioCpp::pdb::model<atom_t>::type& all_info, 
                        std::vector<typename delaunay_t::Vertex_handle>& surface ){

  delaunay_t DT( all_info.begin(), all_info.end() );

  std::set< std::array<typename delaunay_t::Vertex_handle,3> > faces;
  for (typename delaunay_t::Finite_cells_iterator it = DT.finite_cells_begin(); it != DT.finite_cells_end(); ++it){
    for( int excl = 0; excl < 4; ++excl ){
      std::map<int, int> order;
      for(int i = 0; i < 4; ++i){
        if(i!=excl){
          order[ it->vertex(i)->point().serial ] = i;
        }
      }  
      
      std::map<int, int>::iterator ord0 = order.begin();
      std::map<int, int>::iterator ord1 = ord0; ++ord1;
      std::map<int, int>::iterator ord2 = ord1; ++ord2;
      
      std::array<typename delaunay_t::Vertex_handle,3> arr;
      arr[0] = it->vertex(ord0->second);
      arr[1] = it->vertex(ord1->second); 
      arr[2] = it->vertex(ord2->second);
      
      typename std::set< std::array<typename delaunay_t::Vertex_handle,3> >::iterator ele = faces.find( arr );
      if( ele==faces.end() ){
        faces.insert( arr );
      }
      else{
        if( ele!=faces.end() ){
          faces.erase( ele );
        }
      }
    }
  }
  std::set< int > ins;
  for( typename std::set< std::array<typename delaunay_t::Vertex_handle,3> >::iterator it = faces.begin(); it!=faces.end(); ++it ){
    for(int i = 0; i < 3; ++i){
      int serial = (*it)[i]->point().serial;
      if( ins.find(serial) == ins.end() ){
        surface.push_back( (*it)[i] );
        ins.insert( serial );
      }
    }
  }
}

}
}
#endif

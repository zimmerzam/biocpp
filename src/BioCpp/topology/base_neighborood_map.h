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

#ifndef NEIGHBOROOD_MAP_H
#define NEIGHBOROOD_MAP_H

#include <map>
#include <set>

namespace BioCpp{

/*!
    \brief Neighborood map
*/
template <typename obj>
class base_neighborood_map{
    typedef std::map<obj, std::set<obj> > neighbor_list;
  private:
    neighbor_list neighbors;
  public:
		typedef typename neighbor_list::iterator iterator;
		typedef typename neighbor_list::const_iterator const_iterator;
		typedef typename neighbor_list::reverse_iterator reverse_iterator;
		
		iterator begin(){return neighbors.begin();};
		iterator end(){return neighbors.end();};
		reverse_iterator rbegin(){return neighbors.rbegin();};
		reverse_iterator rend(){return neighbors.rend();};
		
		std::set< obj >& operator [] (obj key){
		  return neighbors[key];
		};
		
		unsigned int size(){
		  return neighbors.size();
		}
		
		unsigned int size(obj key){
		  return neighbors[key].size();
		}
};

} // end namespace
#endif

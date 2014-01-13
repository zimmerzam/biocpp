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

#ifndef BIOCPP_BASE_CONTACT_MAP_H
#define BIOCPP_BASE_CONTACT_MAP_H

#include <unordered_map>

namespace BioCpp{
namespace base{

/*! \brief Describes a table
		
		This simple `struct` can be used in order to store all the H-bridges of
		a particular system.
		@tparam T usually is the resSeq (`int`) of a residue, or an `iterator`.  
*/
template <typename T, typename value_t>
class contact_map{
	protected:
	  typedef std::unordered_map< std::pair< T, T >, value_t > map;
		map data;
	public:
	  /*! \brief Void constructor */
		contact_map(){};
		
		/*! \brief iterators */
		typedef typename map::iterator iterator; /*!< Iterator over the children */
		typedef typename map::const_iterator const_iterator; /*!< Const iterator over the children */
//		typedef typename map::reverse_iterator reverse_iterator; /*!< Reverse iterator */
		
		iterator begin(){return data.begin();}; /*!< Iterator to the first child */
		iterator end(){return data.end();}; /*!< \note This is not the last child item */
//		reverse_iterator rbegin(){return data.rbegin();}; /*!< Iterator to the last child */
//		reverse_iterator rend(){return data.rend();} /*!< \note This is not the first child */
		
		/*! \brief Get an element */
		value_t& operator[]( std::pair< T, T > key ){return data[key];}
		
		/*! \brief Get an element */
		value_t& operator()( T key1, T key2 ){
		  return (*this)[std::make_pair(key1, key2)];
		}
};

} // end namsepace
} // end namespace

#endif

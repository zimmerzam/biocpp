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

#ifndef BASE_H_BRIDGE_MAP_H
#define BASE_H_BRIDGE_MAP_H

#include <map>

namespace BioCpp{

/*! \brief Describes the H-bridge network
		
		This simple `struct` can be used in order to store all the H-bridges of
		a particular system.
		@tparam T usually is the resSeq (`int`) of a residue, or an `iterator`.  
*/
template <typename T>
class base_h_bridge_map{
		typedef std::map< std::pair< T, T >, bool > bridge_map;
	private:
		bridge_map h_bridge;
	public:
	  /*! \brief Void constructor */
		base_h_bridge_map(){};
		/*! \brief Generic constructor, to be specialized */
		template <typename structure>
		base_h_bridge_map(structure& cmp);
		
		/*! \brief iterators */
		typedef typename bridge_map::iterator iterator; /*!< Iterator over the children */
		typedef typename bridge_map::const_iterator const_iterator; /*!< Const iterator over the children */
		typedef typename bridge_map::reverse_iterator reverse_iterator; /*!< Reverse iterator */
		
		iterator begin(){return h_bridge.begin();}; /*!< Iterator to the first child */
		iterator end(){return h_bridge.end();}; /*!< \note This is not the last child item */
		reverse_iterator rbegin(){return h_bridge.rbegin();}; /*!< Iterator to the last child */
		reverse_iterator rend(){return h_bridge.rend();} /*!< \note This is not the first child */
		
		/*! \brief Get an element */
		bool& operator[]( std::pair< T, T > key ){return h_bridge[key];}
};

} // end namespace

#endif

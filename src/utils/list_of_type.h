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

#ifndef LIST_OF_TYPE_H
#define LIST_OF_TYPE_H

#include <map>

/*! \brief Simple `struct` that allows the initialization of `std::map` object
		with an array of elements */
template<class K, class V>
struct map_list_of_type {
  typedef std::map<K, V> Map;
  Map data;
  map_list_of_type(K k, V v) { data[k] = v; }
  map_list_of_type& operator()(K k, V v) { data[k] = v; return *this; }
  operator Map const&() const { return data; }
};

#endif

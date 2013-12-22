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

#ifndef BIOCPP_DPSS_H_BRIDGE_MAP_CONSTRUCTOR
#define BIOCPP_DPSS_H_BRIDGE_MAP_CONSTRUCTOR

#include "base_h_bridge_map.hxx"

namespace BioCpp{
namespace dpss{

template <typename T, typename complex_t>
class base_h_bridge_map_constructor{
  public:
    typedef typename base_h_bridge_map<T>::type h_bridge_map_type;
    h_bridge_map_type& operator()(complex_t& cmp);
};

} // end dpss
} // end BioCpp
#endif

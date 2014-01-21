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

#ifndef BIOCPP_STANDARD_BASE_DPSS_H_BRIDGE_MAP_CONSTRUCTOR
#define BIOCPP_STANDARD_BASE_DPSS_H_BRIDGE_MAP_CONSTRUCTOR

#include "../base.hpp"
#include "h_bridge_map.hpp"
#include <BioCpp/proteins/dpss/base_h_bridge_map_constructor.hxx>

namespace BioCpp{
namespace standard{
namespace base{

typedef BioCpp::dpss::base_h_bridge_map_constructor<BioCpp::standard::base::chain::iterator, 
                                                    BioCpp::standard::base::complex>
  h_bridge_map_constructor;

}
}
}
#endif

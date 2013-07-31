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

#ifndef BIOCPP_STANDARD_DPSS_H_BRIDGE_MAP
#define BIOCPP_STANDARD_DPSS_H_BRIDGE_MAP

#include <BioCpp/standard/base/chain.h>
#include <BioCpp/polimers/proteins/dpss/base_h_bridge_map.h>

namespace BioCpp{
namespace standard{
namespace base{
namespace dpss{

typedef typename BioCpp::dpss::base_h_bridge_map<BioCpp::standard::base::chain::iterator>::type h_bridge_map;

}
}
}
}
#endif

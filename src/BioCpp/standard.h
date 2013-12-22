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

#ifndef BIOCPP_STANDARD_H
#define BIOCPP_STANDARD_H

// The following are template classes. The only template parameter is the atom type.
#include "standard/proteins/model.hxx"
#include "standard/proteins/residue.hxx"
#include "standard/proteins/chain.hxx"
#include "standard/proteins/complex.hxx"
#include "standard/proteins/complex_constructor.hxx"

// The following are fully specialized classes. The atom type is "BioCpp::base_atom"
#include "standard/proteins/base/model.hpp"
#include "standard/proteins/base/residue.hpp"
#include "standard/proteins/base/chain.hpp"
#include "standard/proteins/base/complex.hpp"
#include "standard/proteins/base/complex_constructor.hpp"

#include "standard/proteins/base/dpss/h_bridge_map.hpp"
#include "standard/proteins/base/dpss/h_bridge_map_constructor.hpp"

// The following are fully specialized classes. The atom type is
// "BioCpp::morphology::cgal_extensible_kernel::point<BioCpp::base_atom>"
#if defined BIOCPP_INCLUDE_CGAL
  #include "standard/proteins/morphology/atom.hpp"
  #include "standard/proteins/morphology/model.hpp"
  #include "standard/proteins/morphology/residue.hpp"
  #include "standard/proteins/morphology/chain.hpp"
  #include "standard/proteins/morphology/complex.hpp"
  #include "standard/proteins/morphology/complex_constructor.hpp"
  #include "standard/proteins/morphology/kernel.hpp"
  #include "standard/proteins/morphology/vertex_base.hpp"
  #include "standard/proteins/morphology/cell_base.hpp"
  #include "standard/proteins/morphology/triangulation_data_structure.hpp"
  #include "standard/proteins/morphology/triangulation_3.hpp"
  #include "standard/proteins/morphology/alpha_shape_3.hpp"
  #include "standard/proteins/morphology/delaunay_3.hpp"
#endif

//namespace BioCpp{
//namespace standard{
//typedef BioCpp::base_h_bridge_map<chain::iterator> h_bridge_map;
//}
//}
//  #include <CGAL/basic.h>
//  #include <CGAL/Filtered_kernel.h>
//  #include <CGAL/Delaunay_triangulation_3.h>
//  #include <CGAL/Alpha_shape_3.h>
//}
#endif

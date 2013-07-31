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

#include "standard/ids/atom_id.h"
#include "standard/ids/element_id.h"
#include "standard/ids/amino_acid_id.h"
#include "standard/ids/moiety_id.h"

// The following are template classes. The only template parameter is the atom type.
#include "standard/model.h"
#include "standard/residue.h"
#include "standard/chain.h"
#include "standard/complex.h"
#include "standard/complex_constructor.h"

// The following are fully specialized classes. The atom type is "BioCpp::base_atom"
#include "standard/base/model.h"
#include "standard/base/residue.h"
#include "standard/base/chain.h"
#include "standard/base/complex.h"
#include "standard/base/complex_constructor.h"

#include "standard/base/dpss/h_bridge_map.h"
#include "standard/base/dpss/h_bridge_map_constructor.h"

// The following are fully specialized classes. The atom type is
// "BioCpp::morphology::cgal_extensible_kernel::point<BioCpp::base_atom>"
#if defined BIOCPP_INCLUDE_CGAL
  #include "standard/morphology/atom.h"
  #include "standard/morphology/model.h"
  #include "standard/morphology/residue.h"
  #include "standard/morphology/chain.h"
  #include "standard/morphology/complex.h"
  #include "standard/morphology/complex_constructor.h"
  #include "standard/morphology/kernel.h"
  #include "standard/morphology/vertex_base.h"
  #include "standard/morphology/cell_base.h"
  #include "standard/morphology/triangulation_data_structure.h"
  #include "standard/morphology/triangulation_3.h"
  #include "standard/morphology/alpha_shape_3.h"
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

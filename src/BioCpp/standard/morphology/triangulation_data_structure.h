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

#ifndef BIOCPP_STANDARD_MORPHOLOGY_TRIANGULATION_DATA_STRUCTURE_H
#define BIOCPP_STANDARD_MORPHOLOGY_TRIANGULATION_DATA_STRUCTURE_H

#include "vertex_base.h"
#include "cell_base.h"
#include <CGAL/Delaunay_triangulation_3.h>

namespace BioCpp{
namespace standard{
namespace morphology{

typedef CGAL::Triangulation_data_structure_3<vertex_base,cell_base> triangulation_data_structure;

}
}
}

#endif

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

#ifndef MORPHOLOGY_KERNEL_BBOX_H
#define MORPHOLOGY_KERNEL_BBOX_H

#include <CGAL/Origin.h>
#include <CGAL/Bbox_3.h>
#include "point.h"

namespace BioCpp{
namespace morphology{
namespace cgal_extensible_kernel{

template <typename ConstructBbox>
class bbox : public ConstructBbox {
public:
  using ConstructBbox::operator();

  template < typename point_t >
  CGAL::Bbox_3 operator()(const point<point_t>& p) const {
    return CGAL::Bbox_3(p.x(), p.y(), p.x(), p.y());
  }
};

}
}
}

#endif 

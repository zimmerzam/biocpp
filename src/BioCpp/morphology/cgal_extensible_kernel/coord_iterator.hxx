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

#ifndef MORPHOLOGY_KERNEL_COORD_ITERATOR_H
#define MORPHOLOGY_KERNEL_COORD_ITERATOR_H

#include "point.hxx"

namespace BioCpp{
namespace morphology{
namespace cgal_extensible_kernel{

class coord_iterator {
public:
  template < typename point_t >
  const double* operator()(const point<point_t>& p){
    return &p.x();
  }

  template < typename point_t >
  const double* operator()(const point<point_t>& p, int k){
    if(k==0) return &p.x();
    if(k==1) return &p.y();
    if(k==2) return &p.z();
    const double* pzptr = &p.z();
    pzptr++;
    return pzptr;
  }
};

}
}
}

#endif

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

#ifndef SURFACE_AREA_LCPO_H
#define SURFACE_AREA_LCPO_H

#include "../base_atom/has_coordinate.h"
#include "../topology/base_neighborood_map.h"
#include "two_sphere_buried_surface_area.h"

namespace BioCpp{

/*! \brief Estimate the surface area exposed to solvent by using the LCPO method 

    According to LCPO algorithm the exposed surface area \f$ A_i \f$ of a sphere \f$ i \f$ is given by:
    \f[
      A_{i} = p_1S_1 + p_2\sum_j A_{ij} + p_3\sum_{jk}A_{jk} + p_{3_1}\sum_jA_{ij}\left(\sum_kA_{jk}\right),
    \f]
    where \f$ S_1 \f$ is the total surface area of the sphere and \f$ A_{ij} \f$ is the surface of sphere \f$ i \f$
    buried by sphere \f$ j \f$.
*/
template< typename T,
          typename = typename std::enable_if< has_coordinate<T>::value, T >::type >
double SurfaceAreaLCPO( T& center, base_neighborood_map<T*>& map, const double& radius){
  double sum2 = 0., sum3 = 0., adj3 = 0.;
  for( typename std::set< T* >::iterator k = map[&center].begin(); k != map[&center].end(); ++k ){
    double tmp_sum2 = TwoSphereBuriedSurfaceArea(radius, radius, (center.coordinate - (*k)->coordinate).norm() );
    double tmp_sum3 = 0.;
    for( typename std::set< T* >::iterator j = map[(*k)].begin(); j != map[(*k)].end(); ++j ){
      if( map[&center].find( (*j) ) != map[&center].end() and (*j)!=(&center) ){
        tmp_sum3 += TwoSphereBuriedSurfaceArea(radius, radius, ( (*j)->coordinate - (*k)->coordinate).norm() );
			}
    }
    sum2 += tmp_sum2;           // contribution from two-spheres 
    sum3 += tmp_sum3;           // contribution from three-spheres 
    adj3 += tmp_sum2*tmp_sum3;  // adjust-term
  }
	return 0.23348*4*M_PI*radius*radius - 0.072627*sum2 - 0.00020079*sum3 + 0.000079670*adj3;

}

template< typename T,
          typename container,
          typename = typename std::enable_if< is_container_of< container, T >::value, T >::type,
          typename = typename std::enable_if< has_coordinate<T>::value, T >::type >
double SurfaceAreaLCPO( container& center, base_neighborood_map<T*>& map, const double& radius){
  double area = 0;
  for(typename container::iterator it = center.begin(); it!= center.end(); ++it){
    area += SurfaceAreaLCPO( *it, map, radius);
  }
  return area;
}


} //end namespace

#endif

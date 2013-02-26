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

#include "base_neighborood_map.h"
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
template< typename my_residue >
class BaseSurfaceAreaLCPO{
  private:
    double param1  ; /*!< Parameter that multiplies \f$ S_1 \f$ */
    double param2  ; /*!< Parameter for the two-body term */
    double param3  ; /*!< Parameter for the three-body term */
    double param3_1; /*!< Parameter used to adjust the three-body term */

  public:
    double value;  /*!< The final surface area */
    double radius; /*!< The radius of each atom */
    
    /*! \brief Neighborood map */
    base_neighborood_map< typename my_residue::iterator > map;

    /*! \brief Void constructor */
    BaseSurfaceAreaLCPO();
    
    /*! \brief Apply the LCPO method */
    void operator() (typename my_residue::iterator& center );    
};

template< typename my_residue >
BaseSurfaceAreaLCPO<my_residue>::BaseSurfaceAreaLCPO(){
  param1   =  0.23348;
  param2   = -0.072627;
  param3   = -0.00020079;
  param3_1 =  0.000079670;
  value = 0.;
};

template <typename my_residue>
void BaseSurfaceAreaLCPO<my_residue>::operator()( typename my_residue::iterator& center ){
  double sum2 = 0., sum3 = 0., adj3 = 0.;
  for( typename std::set< typename my_residue::iterator >::iterator k = map[center].begin(); k != map[center].end(); ++k ){
    double tmp_sum2 = TwoSphereBuriedSurfaceArea(radius, radius, (center->coordinate - (*k)->coordinate).norm() );
    double tmp_sum3 = 0.;
    for( typename std::set< typename my_residue::iterator >::iterator j = map[(*k)].begin(); j != map[(*k)].end(); ++j ){
      if( map[center].find( (*j) ) != map[center].end() and (*j)!=center ){
        tmp_sum3 += TwoSphereBuriedSurfaceArea(radius, radius, ( (*j)->coordinate - (*k)->coordinate).norm() );
			}
    }
    sum2 += tmp_sum2;           // two-spheres contribution
    sum3 += tmp_sum3;           // three-spheres contribution
    adj3 += tmp_sum2*tmp_sum3;  // three-spheres adjust-term
  }
	value += param1*4*M_PI*radius*radius + param2*sum2 + param3*sum3 + param3_1*adj3;
}

} //end namespace

#endif

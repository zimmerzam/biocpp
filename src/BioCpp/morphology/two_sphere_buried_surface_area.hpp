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

#ifndef TWO_SPHERE_BURIED_SURFACE_AREA
#define TWO_SPHERE_BURIED_SURFACE_AREA

namespace BioCpp{

/*! \brief Compute the surface of a sphere buried by another sphere

    Surface \f$ A_i \f$ of a sphere \f$ i \f$ of radius \f$ r_i \f$ that is buried by a sphere
    \f$ j \f$ of radius \f$ r_j \f$ is given by
    \f[
      A_i = 2\pi r_{i}\left( r_{i} - \frac{ d_{ij} }{ 2 } - \frac{r_{i}^2 - r_{j}^2}{ 2d_{ij} } \right),
    \f]
    where \f$ d_{ij} \f$ is the center-to-center distance between the two spheres.

    \param radius1 Radius of the first sphere
    \param radius2 Radius of the second sphere
    \param distance the spheres center-to-center distance
    \return The surface of a sphere of radius1 buried from a sphere of radius2
*/
inline double TwoSphereBuriedSurfaceArea(double radius1, double radius2, double distance){
  if( radius1+radius2 < distance ){
    return 0.;
  }
	if(radius1 == radius2){
    return 2*M_PI*radius1*( radius1 - distance/2. );
  }
  return 2*M_PI*radius1*( radius1 - distance/2. - ( radius1*radius1 - radius2*radius2 )/(2.*distance) );
}

} //end namespace
#endif

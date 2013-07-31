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

#ifndef H_BRIDGE_ENERGY
#define H_BRIDGE_ENERGY

#include <Eigen/Core>

namespace BioCpp{
namespace dpss{

/*! \brief Estimates the H-bridge energy.
		H-bridge energy is here estimated by using the DPSS formula
		\f[ E=27.88\left( \frac{1}{r(on)} + \frac{1}{r(ch)} - \frac{1}{r(oh)} - \frac{1}{r(cn)} \right) \f]
		where \f$ r(xy) \f$ is the distance between atom \f$x\f$ and \f$y\f$ and is the same as `dist_xy` in the code.
		\return the H-bridge energy, according to DPSS definition.
		@param dist_on distance between the oxygen-like atom and the nitrogen-like one
		@param dist_ch distance between the carbon-like atom and the hydrogen-like one
		@param dist_oh distance between the oxygen-like atom and the hydrogen-like one
		@param dist_cn distance between the carbon-like atom and the nitrogen-like one
		\note Usually an H-bridge is assigned if \f$E\f$ is less than a cutoff value, more often \f$E<-0.5kcal/mol\f$.

*/
inline double h_bridge_energy( double dist_on, double dist_ch, 
												double dist_oh, double dist_cn );

/*! \brief Estimates the H-bridge energy.

		H-bridge energy is here estimated by using the DPSS formula
		\return the H-bridge energy, according to DPSS definition.
		@param c is the carbon-like atom coordinate
		@param o is the oxygen-like atom coordinate
		@param n is the nitrogen-like atom coordinate
		@param h is the hydrogen-like atom coordinate
		\note Usually an H-bridge is assigned if \f$E\f$ is less than a cutoff value, more often \f$E<-0.5kcal/mol\f$.
		\see h_bridge_energy
*/
inline double h_bridge_energy( Eigen::Vector3d c, Eigen::Vector3d o, 
												Eigen::Vector3d n, Eigen::Vector3d h );



inline double h_bridge_energy( Eigen::Vector3d c, Eigen::Vector3d o, 
												Eigen::Vector3d n, Eigen::Vector3d h ){
	double dist_on = ( o - n ).norm();
	double dist_ch = ( c - h ).norm();
	double dist_oh = ( o - h ).norm();
	double dist_cn = ( c - n ).norm();
	return h_bridge_energy( dist_on, dist_ch, dist_oh, dist_cn );
}

inline double h_bridge_energy( double dist_on, double dist_ch, double dist_oh, double dist_cn ){
	return 27.888*( 1./dist_on + 1./dist_ch - 1./dist_oh - 1./dist_cn );
}

/*! \example dpss.cpp */

} // end dpss
} // end BioCpp
#endif

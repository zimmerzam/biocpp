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

#include "bb_oxygen.hpp"

namespace BioCpp{
namespace protein{

inline Eigen::Vector3d BbOxygenPosition( Eigen::Vector3d ca, Eigen::Vector3d c, Eigen::Vector3d n ){
	Eigen::Vector3d a = ca-c;
	Eigen::Vector3d b = n-c;
	a.normalize();
	b.normalize();
	double den = a(0)*b(1) - a(1)*b(0);
	Eigen::Vector3d r;
	r(0) = ( b(1)*cos_CA_C_O_angle - a(1)*cos_O_C_N_angle )/den;
	r(1) = ( a(0)*cos_O_C_N_angle - b(0)*cos_CA_C_O_angle )/den;
	r(2) = 1. - r(0)*r(0) - r(1)*r(1);
	return c+C_O_bond_length*r;
}

} // end namespace
}

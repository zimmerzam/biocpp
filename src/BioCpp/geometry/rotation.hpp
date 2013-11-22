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

#ifndef BIOCPP_GEOMETRY_ROTATION_T_H
#define BIOCPP_GEOMETRY_ROTATION_T_H

#include <Eigen/Dense>

template < typename atom_t >
struct rotation_t{
	Eigen::Vector3d center;
	Eigen::Matrix3d rotation_m;
	
	rotation_t( Eigen::Vector3d& cent, Eigen::Matrix3d& rot_m ): center(cent), rotation_m(rot_m) {};
	rotation_t( Eigen::Vector3d  cent = Eigen::Vector3d(0.,0.,0.), Eigen::Matrix3d rot_m = Eigen::Matrix3d::Identity() ): center(cent), rotation_m(rot_m){};
	
	void operator()(atom_t& atm, Eigen::Vector3d& cent, Eigen::Matrix3d& rot_m){
		atm.coordinate=rot_m*(atm.coordinate-cent) + cent;
	}
	
	void operator()(atom_t& atm){
		(*this)(atm, center, rotation_m);
	}
	
};

#endif

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

#ifndef BIOCPP_GEOMETRY_TRASLATION_T_H
#define BIOCPP_GEOMETRY_TRASLATION_T_H

#include <Eigen/Dense>

template < typename atom_t >
struct traslation_t{
	Eigen::Vector3d displacement;
	
	traslation_t( Eigen::Vector3d  displ = Eigen::Vector3d(0.,0.,0.) ): displacement(displ){};
	
	void operator()(atom_t& atm, Eigen::Vector3d& displ){
		atm.coordinate+=displ;
	}
	
	void operator()(atom_t& atm){
		(*this)(atm, displacement);
	}
};

#endif

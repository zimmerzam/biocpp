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

#ifndef BETA_CARBON_RECONSTRUCTION_H
#define BETA_CARBOM_RECONSTRUCTION_H

#include <Eigen/Core>

namespace BioCpp{
namespace protein{
/*! \brief Estimates the position of the beta_carbon atom.

		\return the coordinate of the beta carbon atom given the position of the alpha carbon of near residues. 
		@param cap coordinate of the carbon CA of the previous residue
		@param ca  coordinate of the carbon CA of the same residue
		@param caf coordinate of the carbon CA of the following residue 
		
		\note TODO */
inline Eigen::Vector3d BetaCarbonPosition( Eigen::Vector3d cap, Eigen::Vector3d ca, Eigen::Vector3d caf );

} // end namespace
}
#endif

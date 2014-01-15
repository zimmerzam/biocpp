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

#ifndef BACKBONE_OXYGEN_RECONSTRUCTION_H
#define BACKBONE_OXYGEN_RECONSTRUCTION_H

#include <Eigen/Core>

namespace BioCpp{
namespace protein{

/*! \brief Estimates the position of the hydrogen backbone atom.

		In order to determine the secondary structures of a protein, the positions of C', O, N and H
		backbone atoms are required. When the hydrogen atom is missing in the pdb, its coordinate can 
		be estimated by using this function.
		\return the coordinate of the hydrogen backbone atom given the position of the
		nitrogen atom of the same residue as well as the position of carbon and oxygen atoms
		of the previous residue. 
		@param ca coordinate of the carbon CA of the same residue
		@param c  coordinate of the carbon C' of the same residue
		@param n  coordinate of the nitrogen of the following residue */
inline Eigen::Vector3d BbOxygenPosition( Eigen::Vector3d ca, Eigen::Vector3d c, Eigen::Vector3d n );

} // end namespace
}
#endif

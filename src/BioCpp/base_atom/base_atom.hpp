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

#ifndef ATOM_H
#define ATOM_H

#include <Eigen/Core>

namespace BioCpp{
namespace base{

/*! \brief ATOM line in pdb file

		According to [pdb specification 3.3 for ATOM record ](http://www.wwpdb.org/documentation/format33/sect9.html#ATOM ).
		This is just a simple atom: you can define your own 'atom' with other properties/methods by inheriting from this class. 
*/
class atom{
  public:
    int serial;                     /*!< serial (progressive) number */
 	  int id;                         /*!< atom identifier (CA, CB, ...) */
 	  char altLoc;                    /*!< alternative location */
 	  int resName;                    /*!< residue name */
 	  char chainId;                   /*!< chain identifier */
 	  int resSeq;                     /*!< residue number */
 	  char iCode;                     /*!< code for insertion of residues */
 	  Eigen::Vector3d coordinate;     /*!< coordinate */
 	  double occupancy;               /*!< occupancy */
 	  double tempFactor;              /*!< temperature factor */
 	  int element;                    /*!< element type (carbon, oxygen, ...) */
 	  double charge;                  /*!< charge */
};

}
}

#endif

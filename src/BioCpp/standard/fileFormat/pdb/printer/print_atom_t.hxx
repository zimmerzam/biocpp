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

#ifndef BIOCPP_PDB_PRINT_ATOM_LINE_H
#define BIOCPP_PDB_PRINT_ATOM_LINE_H

#include <string>
#include <ostream>
#include <iomanip>
#include <cmath>
#include <Eigen/Dense>

namespace BioCpp{
namespace io{
namespace pdb{

/*! \brief Print an ATOM line from a atom_info
		
		@param out the output stream (i.e. `std::cout`)
		@param info the atom_info to be printed 
*/
template <typename eleDict, typename atmDict, typename resDict>
struct print_atom_t{
  std::ostream& out;
  eleDict& element_dictionary;
  atmDict& atom_dictionary;
  resDict& residue_dictionary;

  print_atom_t(std::ostream& dev, eleDict& e, atmDict& a, resDict& r): 
       out(dev), element_dictionary(e), atom_dictionary(a), residue_dictionary(r) {};

  /*! \brief Print an ATOM line
		@param out the output stream (i.e. `std::cout`)
		@param serial atom serial number
		@param id atom identifier
		@param altLoc alternative location
		@param resName residue name
		@param chainId chain identifier
		@param resSeq residue number
		@param iCode code for insertion of residues
		@param coordinate atom coordinate
		@param occupancy occupancy
		@param tempFactor temperature factor
		@param element element type
		@param charge atom charge 
  */
  std::ostream& operator()( int& serial, int& id, char& altLoc, 
														int& resName, char& chainId, int& resSeq, char& iCode, 
														Eigen::Vector3d& coordinate, double& occupancy, 
														double& tempFactor, int& element, double& charge ){
	  out << "ATOM  ";
	  out << std::setw(5) << serial;
	  out << " ";
	  out << atom_dictionary.id_to_string[id];
	  out << altLoc;
	  out << residue_dictionary.id_to_string[resName];
	  out << " ";
	  out << chainId;
	  out << std::setw(4) << resSeq;
	  out << iCode;
	  out << "   ";
	  out << std::fixed << std::setprecision(3);
	  out << std::setw(8) << coordinate(0);
	  out << std::setw(8) << coordinate(1);
	  out << std::setw(8) << coordinate(2);
	  out << std::fixed << std::setprecision(2);
	  out << std::setw(6) << occupancy;
	  out << std::setw(6) << tempFactor;
	  out << "          "; // blank
	  out << element_dictionary.id_to_string[element];
	  if( charge == 0. )
		  out << "  ";
	  else{
		  out << std::setprecision(0) << fabs(charge);
		  if( charge > 0 )
			  out << '+';
		  else if(charge < 0){
			  out << '-';
		  }
	  }
	  out << std::endl;
	  return out;
  };
  
  template <typename atom_t>
  std::ostream& operator()( atom_t& info){
	  return (*this)( info.serial, info.id, info.altLoc, 
									  info.resName, info.chainId, info.resSeq, 
									  info.iCode, info.coordinate, info.occupancy, 
									  info.tempFactor, info.element, info.charge );
  };

};

} // end namespace
} // end namespace
} // end namespace
#endif

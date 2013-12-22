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

#ifndef BIOCPP_PDB_PARSE_ATOM_H
#define BIOCPP_PDB_PARSE_ATOM_H

#include <string>
#include <stdexcept>
#include <Eigen/Dense>
#include "../sections_and_records/sections_and_records.hpp"

#include <BioCpp/standard/chemicalComponentDictionary/element_dictionary_standard.hpp>
#include <BioCpp/standard/chemicalComponentDictionary/atom_dictionary_standard.hpp>
#include <BioCpp/standard/chemicalComponentDictionary/residue_dictionary_standard.hpp>

namespace BioCpp{
namespace io{
namespace pdb{

/*! \brief Get atom line info from a pdb ATOM line 
		
		\return a atom_info containing all informations from the ATOM line
		@param line the ATOM line
*/
template <typename atom_t, typename eleDict, typename atmDict, typename resDict>
atom_t parseAtom(std::string& line, eleDict& e, atmDict& a, resDict& r){
	atom_t info;
	if(get_record(line) == ATOM){
	  info.serial = atoi(line.substr(6, 5).c_str()); 
  	info.id = a.string_to_id.find(line.substr(12, 4))!=a.string_to_id.end() ? a.string_to_id[line.substr(12, 4)] : -1;
  	info.altLoc = line.substr(16, 1).c_str()[0]=='A' ? ' ' : line.substr(16, 1).c_str()[0];
  	info.resName = r.string_to_id.find( line.substr(17, 3) )!=r.string_to_id.end() ? r.string_to_id[line.substr(17, 3)] : -1;
  	info.chainId = isdigit( line.substr(21, 1).c_str()[0] ) ? char( 'A'+atoi(line.substr(21, 1).c_str()) ) : line.substr(21, 1).c_str()[0];
  	if(info.chainId == ' ') info.chainId = 'A';
  	info.resSeq = atoi( line.substr(22, 4).c_str() );
  	info.iCode = line.substr(26, 1).c_str()[0];
  	info.coordinate = Eigen::Vector3d( atof(line.substr(30, 8).c_str()), atof(line.substr(38, 8).c_str()),
  	                           atof(line.substr(46, 8).c_str()));
	
  	try{
  	  info.occupancy = atof(line.substr(54, 6).c_str()) ? atof(line.substr(54, 6).c_str()) : 0.;
  	}
  	catch(std::out_of_range& oor){
  	  info.occupancy = 0.0;
  	}
	
  	try{
  	  info.tempFactor = atof(line.substr(60, 6).c_str()) ? atof(line.substr(60, 6).c_str()) : 0.;
  	}
  	catch (std::out_of_range& oor){
  	  info.tempFactor = 0.0;
  	}
  	info.element = a.definition[info.id].element;
  	info.charge = 0;
  	try {
  	  if( isdigit( line.substr(78, 2).c_str()[0] ) and ( line.substr(78, 2).c_str()[1]=='-' or line.substr(78, 2).c_str()[1]=='+') ){
  	    info.charge = line.substr(78, 2).c_str()[1]=='-' ? -atof(line.substr(78, 1).c_str()) : atof(line.substr(78, 1).c_str());
  	  }
  	}
  	catch (std::out_of_range& oor){
  	  info.charge = 0;
  	}
  }
  return info;
};

template <typename atom_t>
atom_t parseAtom(std::string& line){
  return parseAtom<atom_t, BioCpp::element::dictionary_t, 
                   BioCpp::atom::dictionary_t, BioCpp::residue::dictionary_t >(
                     line, BioCpp::element::dictionary, 
                     BioCpp::atom::dictionary, BioCpp::residue::dictionary
                   ); 
}

} // end namespace
} // end namespace
} // end namespace

#endif

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

#ifndef BIOCPP_PDB_PARSEMODEL_H
#define BIOCPP_PDB_PARSEMODEL_H

#include <BioCpp/io_files/model/model.hxx>
#include "parseAtom.hxx"
#include "../sections_and_records/sections_and_records.hpp"

namespace BioCpp{
namespace io{
namespace pdb{

/*! \brief read a buffer of `char` describing a model and get structured informations.

		@param buffer a string containing a pdb model
		\return a model containing all the info read from the buffer
*/  
template <typename atom_t, typename eleDict, typename atmDict, typename resDict>
typename BioCpp::io::model<atom_t>::type parseModel( char buffer[], eleDict& e, atmDict& a, resDict& r ){
	std::vector< std::pair<int,atom_t> > all_info;
	char* c_line = strtok(buffer, "\n");
	while(c_line){
		std::string line(c_line, std::find(c_line, c_line + 70, '\0'));
		if(get_record(line) == ATOM){
			atom_t info = parseAtom<atom_t, eleDict, atmDict, resDict>(line, e, a, r);
			all_info.push_back( std::make_pair( info.serial, info ) );
		}
		c_line = strtok(NULL, "\n");
	}
	return typename BioCpp::io::model<atom_t>::type(all_info);
}

template <typename atom_t>
typename BioCpp::io::model<atom_t>::type parseModel( char buffer[]){
  return parseModel< atom_t, BioCpp::element::dictionary_t, 
                     BioCpp::atom::dictionary_t, BioCpp::residue::dictionary_t >(
                         buffer, BioCpp::element::dictionary, 
                         BioCpp::atom::dictionary, BioCpp::residue::dictionary
                     );
}
} //end namespace
} //end namespace
} //end namespace
#endif

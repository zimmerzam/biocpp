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

#ifndef PDB_MODEL_FORMAT_H
#define PDB_MODEL_FORMAT_H

#include "model.h"
#include "pdb_ATOM_format.h"
#include "pdb_sections_and_records.h"

namespace BioCpp{
namespace pdb{

/*! \brief read a buffer of `char` describing a model and get structured informations.

		@param buffer a string containing a pdb model
		\return a model containing all the info read from the buffer
*/  
template <typename atom_t>
typename model<atom_t>::type read_model_record( char buffer[], int option_flag ){
	std::vector< std::pair<int,atom_t> > all_info;
	char* c_line = strtok(buffer, "\n");
	while(c_line){
		std::string line(c_line, std::find(c_line, c_line + 70, '\0'));
		if(get_record(line) == ATOM){
			atom_t info = read_atom_line<atom_t>(line);
			if(info.element == BioCpp::element::H and option_flag == 1){
			  continue;
			}
			all_info.push_back( std::make_pair( info.serial, info ) );
		}
		c_line = strtok(NULL, "\n");
	}
	return typename model<atom_t>::type(all_info);
}

} // end namespace
} //end namespace
#endif

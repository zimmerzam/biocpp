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

#include "../polimers/base_container.h"
#include "pdb_ATOM_format.h"
#include "pdb_sections_and_records.h"

namespace BioCpp{
namespace pdb{

/*! \brief A container with all the atoms of a model in the pdb file

		@tparam int is the serial number of the `atom`
		@tparam atom_info stores the info relative to a specific atom
*/
typedef base_container<int, atom_info> model;

/*! \brief read a buffer of `char` describing a model and get structured informations.

		@param buffer a string containing a pdb model
		\return a model containing all the info read from the buffer
*/
model read_model_record( char buffer[] ){
	std::vector< std::pair<int,atom_info> > all_info;
	char* c_line = strtok(buffer, "\n");
	while(c_line){
		std::string line(c_line, std::find(c_line, c_line + 70, '\0'));
		if(get_record(line) == ATOM){
			atom_info info = read_atom_line(line);
			all_info.push_back( std::make_pair( info.serial, info ) );
		}
		c_line = strtok(NULL, "\n");
	}
	return model(all_info);
}

} // end namespace
} //end namespace
#endif

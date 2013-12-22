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

#ifndef BIOCPP_PDB_PARSE_SEQRES_H
#define BIOCPP_PDB_PARSE_SEQRES_H

#include <BioCpp/io_files/seqres/seqres_record.hpp>
#include "../sections_and_records/sections_and_records.hpp"
#include <string>
#include <cstring>
#include <algorithm>
#include <sstream>

namespace BioCpp{
namespace io{
namespace pdb{

/*! \brief Read a SEQRES record from a buffer
		
		\return a seqres_record object
		@param buffer a buffer containing a pdb SEQRES record
*/
template <typename resDict>
BioCpp::io::seqres_record parseSeqres( char buffer[], resDict& res_dict ){
	BioCpp::io::seqres_record record;
	char* c_line = strtok(buffer, "\n");
	while(c_line){
		std::string line(c_line, std::find(c_line, c_line + 80, '\0'));
		char chain_id;
		std::string s_res = "   ";
		if(get_record(line) == SEQRES){
			chain_id = isdigit( line.substr(11, 1).c_str()[0] ) ? char( 'A'+atoi(line.substr(11, 1).c_str()) ) : line.substr(11, 1).c_str()[0];
			chain_id = (chain_id==' ') ? 'A' : chain_id;
			for(int col = 19; col < 70; col+=4){
        try {
          s_res = line.substr(col,3);
        }
        catch (std::out_of_range& oor){
          break;
        }
        if( s_res != "   " ){
          int res = res_dict.string_to_id[ s_res ];
          if( record.find( chain_id )!=record.end() )
            record[ chain_id ] += res_dict.definition[res].one_letter_name;
          else{
            std::stringstream ss;
            ss << res_dict.definition[res].one_letter_name;
            record.insert( std::make_pair(chain_id, ss.str() ) );
          }
        }
      }
		}
		c_line = strtok(NULL, "\n");
	}
	return record;
}

BioCpp::io::seqres_record parseSeqres( char buffer[] ){
  return parseSeqres<residue::dictionary_t>(buffer, residue::dictionary);
};


} // end namespace
} // end namespace
} // end namespace
#endif

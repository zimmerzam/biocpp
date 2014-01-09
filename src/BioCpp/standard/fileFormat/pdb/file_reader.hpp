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

#ifndef STD_PDB_FILE_READER_HPP
#define STD_PDB_FILE_READER_HPP

#include <BioCpp/io_files/base_file_reader.hpp>
#include <BioCpp/io_files/model/model.hxx>
#include "parser/parseAtom.hxx"
#include "parser/parseModel.hxx"
#include "parser/parseSeqres.hxx"

namespace BioCpp{
namespace io{
namespace pdb{

/*! \brief This is a `pdb file` object.

    Use this to read a pdb file and get a list of atom_info.
    \example pdb.cpp
*/
class file: public BioCpp::io::file{
  public:
    file(const char* pdb_name, int init_flag);
    ~file();
    
    template<typename atom_t, typename eleDict, typename atmDict, typename resDict>
    typename BioCpp::io::model<atom_t>::type readModel(int mdl, 
                                           eleDict& ele_dict, atmDict& atm_dict,
                                           resDict& res_dict){
			if(mdl > this->n_models){
    		typename BioCpp::io::model<atom_t>::type record;
    		return record;
  		}
  		int length = this->model_end_pos[mdl] - this->model_beg_pos[mdl];
  		char mdl_buf[ length + 1 ];
 		 	std::memcpy(mdl_buf, this->buffer + this->model_beg_pos[mdl], length);
	  	return parseModel<atom_t, eleDict, atmDict, resDict>
  		                 ( mdl_buf, ele_dict, atm_dict, res_dict );
		};

    template<typename atom_t>
    typename BioCpp::io::model<atom_t>::type readModel(int mdl){
  		return readModel<atom_t, BioCpp::element::dictionary_t, BioCpp::atom::dictionary_t, BioCpp::residue::dictionary_t>
                 		(mdl, BioCpp::element::dictionary, BioCpp::atom::dictionary, BioCpp::residue::dictionary);
		};

};

} // end namespace
} // end namespace
} // end namespace

#endif

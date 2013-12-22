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

#ifndef BIOCPP_PDB_PRINTSEQRES_H
#define BIOCPP_PDB_PRINTSEQRES_H

#include <string>
#include <iomanip>
#include <ostream>

#include "seqres_record.hpp"

namespace BioCpp{
namespace pdb{

/*! \brief Print a SEQRES record
		
		Print a SEQRES record according to [pdb specification 3.3 for ATOM record](http://www.wwpdb.org/documentation/format33/sect3.html#SEQRES)	
		@param out the output stream (i.e. `std::cout`)
		@param chainId the chain identifier
		@param sequence the FASTA sequence */
template <typename resDict>
std::ostream& printSeqres( std::ostream& out, char chainId, std::string sequence, resDict& res_dict){
	int numRes = sequence.size();
	int line = numRes/13+1;
	for(int serNum = 0; serNum < line; ++serNum){
		out << "SEQRES ";
		out << std::setw(3) << serNum+1;
		out << " ";
		out << chainId;
		out << " ";
		out << std::setw(4) << numRes;
		out << "  ";
		std::string res;
		for(int i =0; i < std::min(13, numRes-13*serNum); ++i){
			try{
				res = sequence.substr(13*serNum + i,1);
			}
			catch(int e){
				break;
			}
			out << res_dict.id_to_string[ res_dict.string_to_id[res] ];
			out << " ";
		}
		out << std::endl;
	}
	return out;
}

/*! \brief Print a SEQRES record 
		
		@param out the output stream (i.e. `std::cout`)
		@param seqres a seqres_record (i.e. BioCpp::pdb.TseqRes)
		\see print_seqres_record( std::ostream&, char, std::string), pdb
*/
template <typename resDict>
std::ostream& printSeqres( std::ostream& out, const seqres_record& seqres, resDict& res_dict){
	for(seqres_record::const_iterator seq=seqres.begin(); seq!=seqres.end(); ++seq){
		printSeqres<resDict>(out, seq->first, seq->second, res_dict);
	}
	return out;
}

template<>
std::ostream& printSeqres( std::ostream& out, const seqres_record& seqres){
  return printSeqres<residue::dictionary_t>(out, seqres, residue::dictionary);
}

} // end namespace
} // end namespace
#endif

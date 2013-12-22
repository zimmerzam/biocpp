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

#ifndef BIOCPP_PDB_PARSE_SEQRES_FROM_FASTA_H
#define BIOCPP_PDB_PARSE_SEQRES_FROM_FASTA_H

#include <string>
#include <algorithm>
#include <cctype>

#include "seqres_record.hpp"

namespace BioCpp{
namespace pdb{

/*! \brief Read a SEQRES record from a FASTA string 
		
		\return a seqres_record object
		@param fasta a string containing info about the sequences. Chain identifier has to be
		">chain_id_1", i.e. ">A"
*/
seqres_record parseSeqresFromFasta(std::string& fasta);

}
}

#endif

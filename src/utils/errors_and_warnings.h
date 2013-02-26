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

#ifndef ERRORS_AND_WARNINGS_H
#define ERRORS_AND_WARNINGS_H

namespace BioCpp{

enum warning{
  WAR_NONE                     = 0,      /*!< No warnings */
	PDB_SEQRES_NOT_FOUND         = (1<<0), /*!< SeqRes section is missing */
	PDB_BACKBONE_HOLE            = (1<<1), /*!< Some residue is missing */
	ALIGN_SEQRES_NEQ_FASTA       = (1<<2)  /*!< Theoretical sequence is different from the structure one */
};

enum error{
  ERR_NONE                       = 0,      /*!< No errors */
	PDB_COORDINATE_NOT_FOUND       = (1<<0), /*!< Coordinate section is missing */
	ALIGN_MISSING_TSEQRES          = (1<<1), /*!< Theoretical sequence is missing */
	ALIGN_FAILED                   = (1<<2), /*!< Alignment failed */
  ALIGN_NOT_A_VALID_TSEQRES      = (1<<3), /*!< Theoretical sequence is not valid */
  ALIGN_NOT_A_VALID_RSEQRES      = (1<<4), /*!< Structure sequence is not valid */
  UNSPECIFIED                    = (1<<5)  /*!< Unspecified error */
};

/*!  */
std::map< BioCpp::warning, std::string > warning_to_string = 
        map_list_of_type< BioCpp::warning, std::string >
                (WAR_NONE, "Success ")
                (PDB_SEQRES_NOT_FOUND, "No SEQRES section found ")
                (PDB_BACKBONE_HOLE, "One or more residues are missing ")
                (ALIGN_SEQRES_NEQ_FASTA, "SEQRES section is different from the input fasta sequence ")
                ;

std::map< BioCpp::error, std::string > error_to_string = 
        map_list_of_type< BioCpp::error, std::string >
                (ERR_NONE, "Success ")
                (PDB_COORDINATE_NOT_FOUND, "No ATOM section found ")
                (ALIGN_MISSING_TSEQRES, "Alignment is not possible due to missing target sequence ")
                (ALIGN_FAILED, "Alignment failed ")
                (ALIGN_NOT_A_VALID_TSEQRES, "Not a valid TseqRes ")
                (ALIGN_NOT_A_VALID_RSEQRES, "Not a valid RseqRes ")
                (UNSPECIFIED, "Unspecified problem. Sorry ")
                ;
} //end namespace
#endif

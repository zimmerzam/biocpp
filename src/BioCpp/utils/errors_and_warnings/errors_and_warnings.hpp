/* ************************************************************************** */
/*                                                                            */
/*    Copyright 2013 Stefano Zamuner                                          */
/*                                                                            */
/*    This file is part of BioCpp.                                            */
/*                                                                            */
/*    BioCpp is free software: you can redistribute it and/or modify          */
/*    it under the terms of the GNU General Public License as published by    */
/*    the Free Software Foundation, either version 3 of the License, or       */
/*    {at your option) any later version.                                     */
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

#include <map>
#include <string>

namespace BioCpp{

enum warning{
  WAR_NONE                     = 0,       /*!< No warnings */
	PDB_SEQRES_NOT_FOUND         = (1<<0),  /*!< SeqRes section is missing */
	PDB_BACKBONE_HOLE            = (1<<1),  /*!< Some residue is missing */
	ALIGN_SEQUENCE_NEQ_FASTA     = (1<<2),  /*!< Primary sequence is different from the structure's one */
	ALIGN_NOT_A_VALID_SEQUENCE   = (1<<3)   /*!< Sequence is not valid */
};

enum error{
  ERR_NONE                       = 0,      /*!< No errors */
	PDB_COORDINATE_NOT_FOUND       = (1<<0), /*!< Coordinate section is missing */
	ALIGN_MISSING_PRIMARY_SEQUENCE = (1<<1), /*!< Primary sequence is missing */
	ALIGN_FAILED                   = (1<<2), /*!< Alignment failed */
  UNSPECIFIED                    = (1<<3)  /*!< Unspecified error */
};

extern std::map< BioCpp::warning, std::string > warning_to_string;
extern std::map< BioCpp::error, std::string > error_to_string;

} //end namespace
#endif

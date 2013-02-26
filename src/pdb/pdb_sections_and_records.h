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

#ifndef PDB_SECTIONS_AND_RECORDS_H
#define PDB_SECTIONS_AND_RECORDS_H

#include <map>
#include "../utils/list_of_type.h"

namespace BioCpp{

/*! \brief A list of pdb sections
*/
enum pdb_section_id{
  PDB_UNKNOWN_SECTION, /*!< Unknown section */
	PDB_TITLE_SEC, /*!< Title */
	PDB_PRIMARY_STRUCTURE, /*!< Everything about the primary structure */
	PDB_HETEROGEN,
	PDB_SECONDARY_STRUCTURE, /*!< Everything about the secondary structure */
	PDB_CONNECTIVITY_ANNOTATION,
	PDB_MISCELLANEOUS_FEATURES, 
	PDB_CRYSTALLOGRAPHIC,
	PDB_COORDINATE,  /*!< In this section there is a description of each atom belonging to the structure */
	PDB_CONNECTIVITY,
	PDB_BOOKKEEPING
};

/*! \brief A list of pdb records

		Each line of pdb file begins with a six-char code describing the line content.
*/
enum pdb_record_id{
	PDB_UNKNOWN_RECORD, /*!< Unknown record */
	PDB_HEADER, /*!< Header record */
	PDB_OBSLTE,
	PDB_TITLE,  /*!< Title record */
	PDB_SPLT,
	PDB_CAVEAT,
	PDB_COMPND,
	PDB_SOURCE,
	PDB_KEYWDS,
	PDB_EXPDTA,
	PDB_NUMMDL,
	PDB_MDLTYP,
	PDB_AUTHOR,
	PDB_REVDAT,
	PDB_SPRSDE,
	PDB_JRNL,
	PDB_REMARK,
	PDB_DBREF,
	PDB_DBREF1,
	PDB_DBREF2,
	PDB_SEQADV,
	PDB_SEQRES, /*!< This record describes the primary sequence of the protein */
	PDB_MODRES,
	PDB_HET,
	PDB_FORMUL,
	PDB_HETNAM,
	PDB_HETSYN,
	PDB_HELIX,
	PDB_SHEET,
	PDB_SSBOND,
	PDB_LINK,
	PDB_CISPEP,
	PDB_SITE,
	PDB_CRYST1,
	PDB_MTRIXn,
	PDB_ORIGXn,
	PDB_SCALEn,
	PDB_MODEL,  /*!< Model record. Specifies the beginning of a new model */
	PDB_ATOM,  /*!< An atom entry in coordinate section */
	PDB_ANISOU,
	PDB_TER, /*!< Specifies the end of a chain */
	PDB_HETATM,
	PDB_ENDMDL, /*!< Specifies the end of a model */
	PDB_CONECT,
	PDB_MASTER,
	PDB_END /*!< Specifies the end of the pdb */
};

/*! Map a six-char string to its corresponding pdb_record_id

*/
std::map< std::string, pdb_record_id > string_to_pdb_record_id = 
				map_list_of_type< std::string, pdb_record_id >
        	("HEADER",PDB_HEADER)("OBSLTE",PDB_OBSLTE)("TITLE ",PDB_TITLE)("SPLT  ",PDB_SPLT)
					("CAVEAT",PDB_CAVEAT)("COMPND",PDB_COMPND)("SOURCE",PDB_SOURCE)("KEYWDS",PDB_KEYWDS)
					("EXPDTA",PDB_EXPDTA)("NUMMDL",PDB_NUMMDL)("MDLTYP",PDB_MDLTYP)("AUTHOR",PDB_AUTHOR)
					("REVDAT",PDB_REVDAT)("SPRSDE",PDB_SPRSDE)("JRNL  ",PDB_JRNL)("REMARK",PDB_REMARK)
					("DBREF ",PDB_DBREF)("DBREF1",PDB_DBREF1)("DBREF2",PDB_DBREF2)("SEQADV",PDB_SEQADV)
					("SEQRES",PDB_SEQRES)("MODRES",PDB_MODRES)("HET   ",PDB_HET)("FORMUL",PDB_FORMUL)
					("HETNAM",PDB_HETNAM)("HETSYN",PDB_HETSYN)("HELIX ",PDB_HELIX)("SHEET ",PDB_SHEET)
					("SSBOND",PDB_SSBOND)("LINK  ",PDB_LINK)("CISPEP",PDB_CISPEP)("SITE  ",PDB_SITE)
					("CRYST1",PDB_CRYST1)("MTRIXn",PDB_MTRIXn)("ORIGXn",PDB_ORIGXn)("SCALEn",PDB_SCALEn)
					("MODEL ",PDB_MODEL)("ATOM  ",PDB_ATOM)("ANISOU",PDB_ANISOU)("TER   ",PDB_TER)
					("HETATM",PDB_HETATM)("ENDMDL",PDB_ENDMDL)("CONECT",PDB_CONECT)("MASTER",PDB_MASTER)
					("END   ",PDB_END);

/*! \brief Get the pdb record described by a pdb line

		\return the corresponding pdb_record_id
		@param line a pdb line
*/
inline pdb_record_id get_pdb_record(std::string& line){
	if(line.size() < 6)
		line+="      ";
	return string_to_pdb_record_id[line.substr(0,6)];
};

} //end namespace

#endif

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
namespace pdb{

/*! \brief A list of pdb sections
*/
enum section_id{
  UNKNOWN_SECTION, /*!< Unknown section */
	TITLE_SEC, /*!< Title */
	PRIMARY_STRUCTURE, /*!< Everything about the primary structure */
	HETEROGEN,
	SECONDARY_STRUCTURE, /*!< Everything about the secondary structure */
	CONNECTIVITY_ANNOTATION,
	MISCELLANEOUS_FEATURES, 
	CRYSTALLOGRAPHIC,
	COORDINATE,  /*!< In this section there is a description of each atom belonging to the structure */
	CONNECTIVITY,
	BOOKKEEPING
};

/*! \brief A list of pdb records

		Each line of pdb file begins with a six-char code describing the line content.
*/
enum record_id{
	UNKNOWN_RECORD, /*!< Unknown record */
	HEADER, /*!< Header record */
	OBSLTE,
	TITLE,  /*!< Title record */
	SPLT,
	CAVEAT,
	COMPND,
	SOURCE,
	KEYWDS,
	EXPDTA,
	NUMMDL,
	MDLTYP,
	AUTHOR,
	REVDAT,
	SPRSDE,
	JRNL,
	REMARK,
	DBREF,
	DBREF1,
	DBREF2,
	SEQADV,
	SEQRES, /*!< This record describes the primary sequence of the protein */
	MODRES,
	HET,
	FORMUL,
	HETNAM,
	HETSYN,
	HELIX,
	SHEET,
	SSBOND,
	LINK,
	CISPEP,
	SITE,
	CRYST1,
	MTRIXn,
	ORIGXn,
	SCALEn,
	MODEL,  /*!< Model record. Specifies the beginning of a new model */
	ATOM,  /*!< An atom entry in coordinate section */
	ANISOU,
	TER, /*!< Specifies the end of a chain */
	HETATM,
	ENDMDL, /*!< Specifies the end of a model */
	CONECT,
	MASTER,
	END /*!< Specifies the end of the pdb */
};

/*! Map a six-char string to its corresponding record_id

*/
std::map< std::string, record_id > string_to_record_id = 
				map_list_of_type< std::string, record_id >
        	("HEADER",HEADER)("OBSLTE",OBSLTE)("TITLE ",TITLE)("SPLT  ",SPLT)
					("CAVEAT",CAVEAT)("COMPND",COMPND)("SOURCE",SOURCE)("KEYWDS",KEYWDS)
					("EXPDTA",EXPDTA)("NUMMDL",NUMMDL)("MDLTYP",MDLTYP)("AUTHOR",AUTHOR)
					("REVDAT",REVDAT)("SPRSDE",SPRSDE)("JRNL  ",JRNL)("REMARK",REMARK)
					("DBREF ",DBREF)("DBREF1",DBREF1)("DBREF2",DBREF2)("SEQADV",SEQADV)
					("SEQRES",SEQRES)("MODRES",MODRES)("HET   ",HET)("FORMUL",FORMUL)
					("HETNAM",HETNAM)("HETSYN",HETSYN)("HELIX ",HELIX)("SHEET ",SHEET)
					("SSBOND",SSBOND)("LINK  ",LINK)("CISPEP",CISPEP)("SITE  ",SITE)
					("CRYST1",CRYST1)("MTRIXn",MTRIXn)("ORIGXn",ORIGXn)("SCALEn",SCALEn)
					("MODEL ",MODEL)("ATOM  ",ATOM)("ANISOU",ANISOU)("TER   ",TER)
					("HETATM",HETATM)("ENDMDL",ENDMDL)("CONECT",CONECT)("MASTER",MASTER)
					("END   ",END);

/*! \brief Get the pdb record described by a pdb line

		\return the corresponding record_id
		@param line a pdb line
*/
inline record_id get_record(std::string& line){
	if(line.size() < 6)
		line+="      ";
	return string_to_record_id[line.substr(0,6)];
};

} //end namespace
} //end namespace

#endif

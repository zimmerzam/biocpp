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

#include "sections_and_records.hpp"

namespace BioCpp{
namespace io{
namespace pdb{

std::map< std::string, record_id > string_to_record_id = {
        	{"HEADER",HEADER},{"OBSLTE",OBSLTE},{"TITLE ",TITLE },{"SPLT  ",SPLT  },
					{"CAVEAT",CAVEAT},{"COMPND",COMPND},{"SOURCE",SOURCE},{"KEYWDS",KEYWDS},
					{"EXPDTA",EXPDTA},{"NUMMDL",NUMMDL},{"MDLTYP",MDLTYP},{"AUTHOR",AUTHOR},
					{"REVDAT",REVDAT},{"SPRSDE",SPRSDE},{"JRNL  ",JRNL  },{"REMARK",REMARK},
					{"DBREF ",DBREF },{"DBREF1",DBREF1},{"DBREF2",DBREF2},{"SEQADV",SEQADV},
					{"SEQRES",SEQRES},{"MODRES",MODRES},{"HET   ",HET   },{"FORMUL",FORMUL},
					{"HETNAM",HETNAM},{"HETSYN",HETSYN},{"HELIX ",HELIX },{"SHEET ",SHEET },
					{"SSBOND",SSBOND},{"LINK  ",LINK  },{"CISPEP",CISPEP},{"SITE  ",SITE  },
					{"CRYST1",CRYST1},{"MTRIXn",MTRIXn},{"ORIGXn",ORIGXn},{"SCALEn",SCALEn},
					{"MODEL ",MODEL },{"ATOM  ",ATOM  },{"ANISOU",ANISOU},{"TER   ",TER   },
					{"HETATM",HETATM},{"ENDMDL",ENDMDL},{"CONECT",CONECT},{"MASTER",MASTER},
					{"END   ",END}
};

/*! \brief Get the pdb record described by a pdb line

		\return the corresponding record_id
		@param line a pdb line
*/
record_id get_record(std::string& line){
	if(line.size() < 6)
		line+="      ";
	return string_to_record_id[line.substr(0,6)];
};

} //end namespace
} //end namespace
} //end namespace

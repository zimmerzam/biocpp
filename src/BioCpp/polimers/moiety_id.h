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

#ifndef MOIETY_ID
#define MOIETY_ID

#include <map>

namespace BioCpp{
namespace moiety{

enum id{ UNK,       /*!< Unknown moiety */ 
         PRI_AMINE, /*!< Primary Amine NH2 */
         SEC_AMINE, /*!< Secondary Amine NH */
         DISULFIDE, /*!< Disulfide SS */
         THIOL,     /*!< Thiol SH */
         ALCOHOL,   /*!< Alcohol OH */
         CH_GROUP,  /*!< CH */
         CARBONYL,  /*< Carbonil CO */
         CARBOXYL,  /*!< Carboxyl COOH*/
         CHARGED_CARBONYL, /*< Carbonil CO- */
         PHENYL,    /*!< Phenil C6H5*/
         PYRROLE,   /*!< Pyrrole C4H4NH*/
         INDOLE     /*!< Indole C8H6NH*/
        };

std::map< std::string, id > string_to_id = {
    {"PAM", PRI_AMINE},{"SAM", SEC_AMINE},
    {"DSU", DISULFIDE},{"THI", THIOL},
    {"ALC", ALCOHOL},{"CXY", CARBOXYL},{"CNY", CARBONYL},
    {"CHg", CH_GROUP},{"CCY", CHARGED_CARBONYL},
    {"PHE", PHENYL},{"PYR", PYRROLE},{"IND",INDOLE},
    {"UNK", UNK}
  };

std::map< id, std::string > id_to_string = {
    {PRI_AMINE, "PAM"},{SEC_AMINE, "SAM"},
    {DISULFIDE, "DSU"},{THIOL, "THI"},
    {ALCOHOL, "ALC"},{CARBOXYL, "CXY"},{CARBONYL, "CNY"},
    {CH_GROUP, "CHg"},{CHARGED_CARBONYL, "CCY"},
    {PHENYL, "PHE"},{PYRROLE, "PYR"},{INDOLE, "IND"},
    {UNK, "UNK"}
  };

} // end moiety
} // end BioCpp

std::ostream& operator << (std::ostream& out, const BioCpp::moiety::id moiety_id ) {
  out << BioCpp::moiety::id_to_string[moiety_id];
  return out;
}

#endif

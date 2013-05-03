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

#ifndef MOIETY_ID
#define MOIETY_ID

#include <map>
#include "../utils/list_of_type.h"

namespace BioCpp{
namespace moiety{

enum id{ UNK,       /*!< Unknown moiety */ 
         PRI_AMINE, /*!< Primary Amine NH2 */
         SEC_AMINE, /*!< Secondary Amine NH */
         DISULFIDE, /*!< Disulfide SS */
         THIOL,     /*!< Thiol SH */
         ALCOHOL,   /*!< Alcohol OH */
         CARBONYL,  /*< Carbonil CO */
         CARBOXYL   /*!< Carboxyl COOH*/
        };

std::map< std::string, id > string_to_id = 
        map_list_of_type< std::string, id >
                ("PRI_AM  ", PRI_AMINE)("SEC_AM  ", SEC_AMINE)
                ("DISULF  ", DISULFIDE)("THIOL   ", THIOL)
                ("ALCOHOL ", ALCOHOL)  ("CARBOXYL", CARBOXYL)("CARBONYL", CARBONYL)
                ("UNK     ", UNK);

std::map< id, std::string > id_to_string = 
        map_list_of_type< id, std::string >
                (PRI_AMINE, "PRI_AM  ")(SEC_AMINE, "SEC_AM  ")
                (DISULFIDE, "DISULF  ")(THIOL, "THIOL   ")
                (ALCOHOL, "ALCOHOL ")(CARBOXYL, "CARBOXYL")(CARBONYL, "CARBONYL")
                (UNK, "UNK     ");

} // end moiety
} // end BioCpp

std::ostream& operator << (std::ostream& out, BioCpp::moiety::id moiety_id ) {
  out << BioCpp::moiety::id_to_string[moiety_id];
  return out;
}

#endif

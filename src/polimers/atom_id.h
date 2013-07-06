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

#ifndef ATOM_ID
#define ATOM_ID

#include <map>
#include <cstring>

#include "../utils/list_of_type.h"

namespace BioCpp{
namespace atom{

/*! \brief id list of the standard atoms that can be found in proteins

    All other atoms are identified with X_ id.
*/
enum id { X_, C_, CA, CB, 
          CD, CD1, CD2, CD3, CE, CE1, CE2, CE3, 
          CG, CG1, CG2, CG3, CZ, CZ1, CZ2, CZ3,
          CH1, CH2, 
          H_, HA, HA1, HA2, HA3,
          HB, HB1, HB2, HB3, 
          HD, HD1, HD2, HD3, HD11, HD12, HD13, HD21, HD22, HD23,
          HE, HE1, HE2, HE3, HE11, HE12, HE13, HE21, HE22, HE23,
          HG, HG1, HG2, HG3, HG11, HG12, HG13, HG21, HG22, HG23,
          HH, HH1, HH2, HH3, HH11, HH12, HH13, HH21, HH22, HH23,
          HN, HN1, HN2, HN3, HN11, HN12, HN13, HN21, HN22, HN23, 
          HZ, HZ1, HZ2, HZ3,
          H1, H2, H3,
          N_, ND, ND1, ND2, NE, NE1, NE2, NH, NH1, NH2, NZ,
          O_, OD, OD1, OD2, OE, OE1, OE2, OH, OH1, OH2, OXT, OG, OG1, OG2,
          SD, SG
        };
/*! \brief Map a string to an atom id

    All strings that do not correspond to a valid atom nomenclature are 
    mapped onto X_ by default

    \code
      #include <iostream>
      #include "../src/BioCpp.h"

      int main(){

        std::cout << BioCpp::atom::string_to_id["2HE2"] << std::endl; // output: "HE22"

        return 0;
      }
    \endcode
*/
std::map< std::string, id > string_to_id = {
    { " X  ", X_},
    { " C  ", C_},{ " CA ", CA},{ " CB ", CB},{ " CD ", CD},
    { " CD1", CD1},{ " CD2", CD2},{ " CD3", CD3},{ " CE ", CE},
    { " CE1", CE1},{ " CE2", CE2},{ " CE3", CE3},{ " CG ", CG},
    { " CG1", CG1},{ " CG2", CG2},{ " CG3", CG3},{ " CZ ", CZ},
    { " CZ1", CZ1},{ " CZ2", CZ2},{ " CZ3", CZ3},{ " CH1", CH1},
    { " CH2", CH2},{ " H  ", H_},
    { " HA ", HA},{ " HA1", HA1},{ "1HA ", HA1},{ " HA2", HA2},{ "2HA ", HA2},{ " HA3", HA3},{ "3HA ", HA3},
    { " HB ", HB},{ " HB1", HB1},{ "1HB ", HB1},{ " HB2", HB2},{ "2HB ", HB2},{ " HB3", HB3},{ "3HB ", HB3},
    { " HD ", HD},{ " HD1", HD1},{ "1HD ", HD1},{ " HD2", HD2},{ "2HD ", HD2},{ " HD3", HD3},{ "3HD ", HD3},
    { "HD11", HD11},{ "1HD1", HD11},{ "HD12", HD12},{ "2HD1", HD12},{ "HD13", HD13},{ "3HD1", HD13},
    { "HD21", HD21},{ "1HD2", HD21},{ "HD22", HD22},{ "2HD2", HD22},{ "HD23", HD23},{ "3HD2", HD23},
    { " HE ", HE},{ " HE1", HE1},{ "1HE ", HE1},{ " HE2", HE2},{ "2HE ", HE2},{ " HE3", HE3},{ "3HE ", HE3},
    { "HE11", HE11},{ "1HE1", HE11},{ "HE12", HE12},{ "2HE1", HE12},{ "HE13", HE13},{ "3HE1", HE13},
    { "HE21", HE21},{ "1HE2", HE21},{ "HE22", HE22},{ "2HE2", HE22},{ "HE23", HE23},{ "3HE2", HE23},
    { " HG ", HG},{ " HG1", HG1},{ "1HG ", HG1},{ " HG2", HG2},{ "2HG ", HG2},{ " HG3", HG3},{ "3HG ", HG3},
    { "HG11", HG11},{ "1HG1", HG11},{ "HG12", HG12},{ "2HG1", HG12},{ "HG13", HG13},{ "3HG1", HG13},
    { "HG21", HG21},{ "1HG2", HG21},{ "HG22", HG22},{ "2HG2", HG22},{ "HG23", HG23},{ "3HG2", HG23},
    { " HH ", HH},{ " HH1", HH1},{ "1HH ", HH1},{ " HH2", HH2},{ "2HH ", HH2},{ " HH3", HH3},{ "3HH ", HH3},
    { "HH11", HH11},{ "1HH1", HH11},{ "HH12", HH12},{ "2HH1", HH12},{ "HH13", HH13},{ "3HH1", HH13},
    { "HH21", HH21},{ "1HH2", HH21},{ "HH22", HH22},{ "2HH2", HH22},{ "HH23", HH23},{ "3HH2", HH23},
    { " HN ", HN},{ " HN1", HN1},{ "1HN ", HN1},{ " HN2", HN2},{ "2HN ", HN2},{ " HN3", HN3},{ "3HN", HN3},
    { "HN11", HN11},{ "1HN1", HN11},{ "HN12", HN12},{ "2HN1", HN12},{ "HN13", HN13},{ "3HN1", HN13},
    { "HN21", HN21},{ "1HN2", HN21},{ "HN22", HN22},{ "2HN2", HN22},{ "HN23", HN23},{ "3HN2", HN23},
    { " HZ ", HZ},{ " HZ1", HZ1},{ "1HZ ", HZ1},{ " HZ2", HZ2},{ "2HZ ", HZ2},{ " HZ3", HZ3},{ "3HZ ", HZ3},
    { " H1 ", H1},{ " 1H ", H1 },{ " H2 ", H2},{ " 2H ", H2 },{ " H3 ", H3},{ " 3H ", H3 },
    { " N  ", N_},
    { " ND ", ND},{ " ND1", ND1},{ " ND2", ND2},
    { " NE ", NE},{ " NE1", NE1},
    { " NE2", NE2},{ " NH ", NH},{ " NH1", NH1},{ " NH2", NH2},
    { " NZ ", NZ},{ " O  ", O_},{ " OD ", OD},{ " OD1", OD1},
    { " OD2", OD2},{ " OE ", OE},{ " OE1", OE1},{ " OE2", OE2},
    { " OH ", OH},{ " OH1", OH1},{ " OH2", OH2},{ " OXT", OXT},
    { " OG ", OG},{ " OG1", OG1},{ " OG2", OG2},{ " SD ", SD},{ " SG ", SG}
  };

/*! \brief Map atom id to a four-letter string

    \code
      #include <iostream>
      #include "../src/BioCpp.h"

      int main(){

        std::cout << BioCpp::atom::id_to_string[BioCpp::atom::N_] << std::endl; // output: " N  

        return 0;
      }
    \endcode
*/
std::map< id, std::string > id_to_string = {
    { X_, " X  "},
    { C_, " C  "},{ CA, " CA "},{ CB, " CB "},{ CD, " CD "},
    { CD1, " CD1"},{ CD2, " CD2"},{ CD3, " CD3"},{ CE, " CE "},
    { CE1, " CE1"},{ CE2, " CE2"},{ CE3, " CE3"},{ CG, " CG "},
    { CG1, " CG1"},{ CG2, " CG2"},{ CG3, " CG3"},{ CZ, " CZ "},
    { CZ1, " CZ1"},{ CZ2, " CZ2"},{ CZ3, " CZ3"},{ CH1, " CH1"},
    { CH2, " CH2"},{ H_, " H  "},
    { HA, " HA "},{ HA1, " HA1"},{ HA2, " HA2"},{ HA3, " HA3"},
    { HB, " HB "},{ HB1, " HB1"},{ HB2, " HB2"},{ HB3, " HB3"},
    { HD, " HD "},{ HD1, " HD1"},{ HD2, " HD2"},{ HD3, " HD3"},
    { HD11, "HD11"},{ HD12, "HD12"},{ HD13, "HD13"},
    { HD21, "HD21"},{ HD22, "HD22"},{ HD23, "HD23"},
    { HE, " HE "},{ HE1, " HE1"},{ HE2, " HE2"},{ HE3, " HE3"},
    { HE11, "HE11"},{ HE12, "HE12"},{ HE13, "HE13"},
    { HE21, "HE21"},{ HE22, "HE22"},{ HE23, "HE23"},
    { HG, " HG "},{ HG1, " HG1"},{ HG2, " HG2"},{ HG3, " HG3"},
    { HG11, "HG11"},{ HG12, "HG12"},{ HG13, "HG13"},
    { HG21, "HG21"},{ HG22, "HG22"},{ HG23, "HG23"},
    { HH, " HH "},{ HH1, " HH1"},{ HH2, " HH2"},{ HH3, " HH3"},
    { HH11, "HH11"},{ HH12, "HH12"},{ HH13, "HH13"},
    { HH21, "HH21"},{ HH22, "HH22"},{ HH23, "HH23"},
    { HN, " HN "},{ HN1, " HN1"},{ HN2, " HN2"},{ HN3, " HN3"},
    { HN11, "HN11"},{ HN12, "HN12"},{ HN13, "HN13"},
    { HN21, "HN21"},{ HN22, "HN22"},{ HN23, "HN23"},
    { HZ, " HZ "},{ HZ1, " HZ1"},{ HZ2, " HZ2"},{ HZ3, " HZ3"},
    { H1, " H1 "},{ H2, " H2 "},{ H3, " H3 "},
    { N_, " N  "},
    { ND, " ND "},{ ND1, " ND1"},{ ND2, " ND2"},
    { NE, " NE "},{ NE1, " NE1"},
    { NE2, " NE2"},{ NH, " NH "},{ NH1, " NH1"},{ NH2, " NH2"},
    { NZ, " NZ "},{ O_, " O  "},{ OD, " OD "},{ OD1, " OD1"},
    { OD2, " OD2"},{ OE, " OE "},{ OE1, " OE1"},{ OE2, " OE2"},
    { OH, " OH "},{ OH1, " OH1"},{ OH2, " OH2"},{ OXT, " OXT"},
    { OG, " OG "},{ OG1, " OG1"},{ OG2, " OG2"},
    { SD, " SD "},{ SG, " SG "}
  };

} //end atom
} //end BioCpp

/*! \brief Print an id as a four-letter string

    \code
      #include <iostream>
      #include "../src/BioCpp.h"

      int main(){

        std::cout << BioCpp::atom::CA << std::endl; // output: " CA "
        std::cout << BioCpp::atom::string_to_id["2HE2"] << std::endl; // output: "HE22"

        return 0;
      }
    \endcode

*/
std::ostream& operator << (std::ostream& out, BioCpp::atom::id atom_id ){
  out << BioCpp::atom::id_to_string[atom_id];
  return out;
}

#endif

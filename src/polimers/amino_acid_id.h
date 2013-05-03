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

#ifndef AMINO_ACID_ID
#define AMINO_ACID_ID

#include <map>
#include "../utils/list_of_type.h"

namespace BioCpp{
namespace amino_acid{

/*! \brief id list of the twenty standard amino-acids

    All non-standard amino-acids are identified with the UNK id.
*/
enum id{ UNK, /*!< Unknown residue */ 
         ALA, /*!< Alanine */
         ARG, /*!< Arginine */
         ASN, /*!< Asparagine */
         ASP, /*!< Aspartic acid */
         CYS, /*!< Cysteine */
         GLN, /*!< Glutamine */
         GLU, /*!< Glutamic acid */
         GLY, /*!< Glycine */
         HIS, /*!< Histidine */
         ILE, /*!< Isoleucine */
         LEU, /*!< Leucine */
         LYS, /*!< Lysine */
         MET, /*!< Methionine */
         PHE, /*!< Phenylalanine */
         PRO, /*!< Proline */
         SER, /*!< Serine */
         THR, /*!< Threonine */
         TRP, /*!< Tryptophan */
         TYR, /*!< Tyrosine */
         VAL  /*!< Valine */
        };

/*! \brief Map a string to an amino-acid id

    All strings that do not correspond to a valid amino-acid nomenclature are 
    mapped onto UNK by default

    \code
      #include <iostream>
      #include "../src/BioCpp.h"

      int main(){

        std::cout << BioCpp::amino_acid::string_to_id["GLY"] << std::endl; // output: "GLY"

        return 0;
      }
    \endcode
*/
std::map< std::string, id > string_to_id = 
        map_list_of_type< std::string, id >
                ("ALA", ALA)("ARG", ARG)("ASN", ASN)("ASP", ASP)
                ("CYS", CYS)("GLN", GLN)("GLU", GLU)("GLY", GLY)
                ("HIS", HIS)("HIE", HIS)("ILE", ILE)("LEU", LEU)
                ("LYS", LYS)("MET", MET)("PHE", PHE)("PRO", PRO)
                ("SER", SER)("THR", THR)("TRP", TRP)("TYR", TYR)
                ("VAL", VAL)("UNK", UNK)
                ("A", ALA)("R", ARG)("N", ASN)("D", ASP)
                ("C", CYS)("Q", GLN)("E", GLU)("G", GLY)
                ("H", HIS)("I", ILE)("L", LEU)("K", LYS)
                ("M", MET)("F", PHE)("P", PRO)("S", SER)
                ("T", THR)("W", TRP)("Y", TYR)("V", VAL)
                ("X", UNK);

/*! \brief Map an id to its corresponding three-letter string

    \code
      #include <iostream>
      #include "../src/BioCpp.h"

      int main(){

        std::cout << BioCpp::amino_acid::string_to_id["GLY"] << std::endl; // output: "GLY"

        return 0;
      }
    \endcode
*/
std::map< id, std::string > id_to_string = 
        map_list_of_type< id, std::string >
                (ALA, "ALA")(ARG, "ARG")(ASN, "ASN")(ASP, "ASP")
                (CYS, "CYS")(GLN, "GLN")(GLU, "GLU")(GLY, "GLY")
                (HIS, "HIS")(ILE, "ILE")(LEU, "LEU")(LYS, "LYS")
                (MET, "MET")(PHE, "PHE")(PRO, "PRO")(SER, "SER")
                (THR, "THR")(TRP, "TRP")(TYR, "TYR")(VAL, "VAL")
                (UNK, "UNK");

/*! \brief Map an id to its corresponding one-letter string

    \code
      #include <iostream>
      #include "../src/BioCpp.h"

      int main(){

        std::cout << BioCpp::amino_acid::id_to_1_letter[BioCpp::amino_acid::PRO] << std::endl; // output: "P"

        return 0;
      }
    \endcode
*/
std::map< id, std::string > id_to_1_letter = 
        map_list_of_type< id, std::string >
                (ALA, "A")(ARG, "R")(ASN, "N")(ASP, "D")
                (CYS, "C")(GLN, "Q")(GLU, "E")(GLY, "G")
                (HIS, "H")(ILE, "I")(LEU, "L")(LYS, "K")
                (MET, "M")(PHE, "F")(PRO, "P")(SER, "S")
                (THR, "T")(TRP, "W")(TYR, "Y")(VAL, "V")
                (UNK, "X");
} // end amino_acid
} //end BioCpp

/*! \brief Print an id as a three-letter string

    \code
      #include <iostream>
      #include "../src/BioCpp.h"

      int main(){

        std::cout << BioCpp::amino_acid::ALA << std::endl; // output: "ALA"
        std::cout << BioCpp::amino_acid::string_to_id["GLY"] << std::endl; // output: "GLY"

        return 0;
      }
    \endcode
*/
std::ostream& operator << (std::ostream& out, BioCpp::amino_acid::id amino_acid_id ) {
  out << BioCpp::amino_acid::id_to_string[amino_acid_id];
  return out;
}

#endif

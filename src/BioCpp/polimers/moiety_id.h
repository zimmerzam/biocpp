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
#include <list>
#include <string>
#include <BioCpp/polimers/atom_id.h>
#include <BioCpp/polimers/amino_acid_id.h>

namespace BioCpp{
namespace moiety{

enum id{ UNK,       /*!< Unknown moiety */
         ALCOHOL,   /*!< Alcohol OH */
         CH2,       /*!< CH2 */
         CH3,       /*!< CH3 */
         CARBONYL,  /*!< Carbonil CO */
         CARBOXYL,  /*!< Carboxylate CO2- */
         PEPTIDE,   /*!< Peptide group CONH */
         SCH3,      /*!< CH3 Sulfide SCH3 */
         THIOL,     /*!< Thiol SH*/
         N3H5,      /*!< N3H5+ */
         SEC_AMINE, /*!< Secondary_amine NH */
         PRI_AMINE, /*!< Primary amine NH2 */
         PRT_AMINE, /*!< Protonated amine NH3+ */
         PHENYL,    /*!< C6H5 */
         C3N2H3,    /*!< C3N2H3 */
         C8NH6      /*!< C8NH6 */
        };

std::map< std::string, id > string_to_id = {
    {"ALC", ALCOHOL}, {"CH2", CH2},
    {"CH3", CH3}, {"CNY", CARBONYL},
    {"CXY", CARBOXYL}, {"SC3", SCH3},
    {"THI", THIOL}, {"N3H", N3H5},
    {"SAM", SEC_AMINE}, {"PAM", PRI_AMINE},
    {"RAM", PRT_AMINE}, {"PHE", PHENYL},
    {"PNT", C3N2H3},{"TRP", C8NH6},
    {"PEP", PEPTIDE},{"UNK", UNK}
  };

std::map< id, std::string > id_to_string = {
    {ALCOHOL, "ALC"}, {CH2, "CH2"},
    {CH3, "CH3"}, {CARBONYL, "CNY"},
    {CARBOXYL, "CXY"}, {SCH3, "SC3"},
    {THIOL, "THI"}, {N3H5, "N3H"},
    {SEC_AMINE, "SAM"}, {PRI_AMINE, "PAM"},
    {PRT_AMINE, "RAM"}, {PHENYL, "PHE"},
    {C3N2H3, "PNT"},{C8NH6, "TRP"},
    {PEPTIDE, "PEP"},{UNK, "UNK"}
  };

typedef std::pair< id, std::list<BioCpp::atom::id> > info;

std::list< info > ResidueToMoietyInfo( BioCpp::amino_acid::id resId, bool first, bool last ){
  std::list< info > info;

  // first residue along the chain...  
  if( first and resId!= BioCpp::amino_acid::PRO ){
    info.push_back( std::make_pair( BioCpp::moiety::PRT_AMINE, std::list<BioCpp::atom::id>{{ BioCpp::atom::N_ , BioCpp::atom::H1  }} ) );
  }
  else if( first and resId == BioCpp::amino_acid::PRO ){
    info.push_back( std::make_pair( BioCpp::moiety::PRI_AMINE, std::list<BioCpp::atom::id>{{ BioCpp::atom::N_ , BioCpp::atom::H1, BioCpp::atom::H2  }} ) );
  }

  // last residue 
  if(last){
    info.push_back( std::make_pair( BioCpp::moiety::CARBOXYL, std::list<BioCpp::atom::id>{{ BioCpp::atom::C_ , BioCpp::atom::O_, BioCpp::atom::OXT  }} ) );
  }

  switch(resId){
    case BioCpp::amino_acid::PRO:
      info.push_back( std::make_pair( BioCpp::moiety::CH2, std::list<BioCpp::atom::id>{{ BioCpp::atom::CB, BioCpp::atom::HB1, BioCpp::atom::HB1 }} ) );
      info.push_back( std::make_pair( BioCpp::moiety::CH2, std::list<BioCpp::atom::id>{{ BioCpp::atom::CG, BioCpp::atom::HG1, BioCpp::atom::HG2 }} ) );
      info.push_back( std::make_pair( BioCpp::moiety::CH2, std::list<BioCpp::atom::id>{{ BioCpp::atom::CD, BioCpp::atom::HD1, BioCpp::atom::HD2 }} ) );
      return info;
    case BioCpp::amino_acid::GLY:
      return info;
    case BioCpp::amino_acid::ALA:
      info.push_back( std::make_pair( BioCpp::moiety::CH3, std::list<BioCpp::atom::id>{{ BioCpp::atom::CB, BioCpp::atom::HB1, BioCpp::atom::HB2, BioCpp::atom::HB3 }} ) );
      return info;
    case BioCpp::amino_acid::VAL:
      info.push_back( std::make_pair( BioCpp::moiety::CH3, std::list<BioCpp::atom::id>{{ BioCpp::atom::CG1, BioCpp::atom::HG11, BioCpp::atom::HG12, BioCpp::atom::HG13 }} ) );
      info.push_back( std::make_pair( BioCpp::moiety::CH3, std::list<BioCpp::atom::id>{{ BioCpp::atom::CG2, BioCpp::atom::HG21, BioCpp::atom::HG22, BioCpp::atom::HG23 }} ) );
      return info;
    case BioCpp::amino_acid::LEU:
      info.push_back( std::make_pair( BioCpp::moiety::CH2, std::list<BioCpp::atom::id>{{ BioCpp::atom::CB, BioCpp::atom::HB2, BioCpp::atom::HB3 }} ) );
      info.push_back( std::make_pair( BioCpp::moiety::CH3, std::list<BioCpp::atom::id>{{ BioCpp::atom::CD1, BioCpp::atom::HD11, BioCpp::atom::HD12, BioCpp::atom::HD13 }} ) );
      info.push_back( std::make_pair( BioCpp::moiety::CH3, std::list<BioCpp::atom::id>{{ BioCpp::atom::CD2, BioCpp::atom::HD21, BioCpp::atom::HD22, BioCpp::atom::HD23 }} ) );
      return info;
    case BioCpp::amino_acid::ILE:
      info.push_back( std::make_pair( BioCpp::moiety::CH2, std::list<BioCpp::atom::id>{{ BioCpp::atom::CG1, BioCpp::atom::HG12, BioCpp::atom::HG13 }} ) );
      info.push_back( std::make_pair( BioCpp::moiety::CH3, std::list<BioCpp::atom::id>{{ BioCpp::atom::CG2, BioCpp::atom::HG21, BioCpp::atom::HG22, BioCpp::atom::HG23 }} ) );
      info.push_back( std::make_pair( BioCpp::moiety::CH3, std::list<BioCpp::atom::id>{{ BioCpp::atom::CD1, BioCpp::atom::HD11, BioCpp::atom::HD12, BioCpp::atom::HD13 }} ) );
      return info;
    case BioCpp::amino_acid::MET:
      info.push_back( std::make_pair( BioCpp::moiety::CH2, std::list<BioCpp::atom::id>{{ BioCpp::atom::CB, BioCpp::atom::HB2, BioCpp::atom::HB3 }} ) );
      info.push_back( std::make_pair( BioCpp::moiety::CH2, std::list<BioCpp::atom::id>{{ BioCpp::atom::CG, BioCpp::atom::HG1, BioCpp::atom::HG2 }} ) );
      info.push_back( std::make_pair( BioCpp::moiety::SCH3, std::list<BioCpp::atom::id>{{ BioCpp::atom::SD, BioCpp::atom::CD2, BioCpp::atom::HE1, BioCpp::atom::HE2, BioCpp::atom::HE3 }} ) );
      return info;
    case BioCpp::amino_acid::CYS:
      info.push_back( std::make_pair( BioCpp::moiety::CH2, std::list<BioCpp::atom::id>{{ BioCpp::atom::CB, BioCpp::atom::HB2, BioCpp::atom::HB3 }} ) );
      info.push_back( std::make_pair( BioCpp::moiety::THIOL, std::list<BioCpp::atom::id>{{ BioCpp::atom::SG }} ) );
      return info;
    case BioCpp::amino_acid::SER:
      info.push_back( std::make_pair( BioCpp::moiety::CH2, std::list<BioCpp::atom::id>{{ BioCpp::atom::CB, BioCpp::atom::HB2, BioCpp::atom::HB3 }} ) );
      info.push_back( std::make_pair( BioCpp::moiety::ALCOHOL, std::list<BioCpp::atom::id>{{ BioCpp::atom::OG, BioCpp::atom::HG }} ) );
      return info;
    case BioCpp::amino_acid::THR:
      info.push_back( std::make_pair( BioCpp::moiety::CH3, std::list<BioCpp::atom::id>{{ BioCpp::atom::CG1, BioCpp::atom::HG1, BioCpp::atom::HG2, BioCpp::atom::HG3 }} ) );
      info.push_back( std::make_pair( BioCpp::moiety::ALCOHOL, std::list<BioCpp::atom::id>{{ BioCpp::atom::OG1, BioCpp::atom::HG1 }} ) );
      return info;
    case BioCpp::amino_acid::ASP:
      info.push_back( std::make_pair( BioCpp::moiety::CH2, std::list<BioCpp::atom::id>{{ BioCpp::atom::CB, BioCpp::atom::HB2, BioCpp::atom::HB3 }} ) );
      info.push_back( std::make_pair( BioCpp::moiety::CARBOXYL, std::list<BioCpp::atom::id>{{ BioCpp::atom::CG, BioCpp::atom::OD1, BioCpp::atom::OD2 }} ) );
      return info;
    case BioCpp::amino_acid::GLU:
      info.push_back( std::make_pair( BioCpp::moiety::CH2, std::list<BioCpp::atom::id>{{ BioCpp::atom::CB, BioCpp::atom::HB2, BioCpp::atom::HB3 }} ) );
      info.push_back( std::make_pair( BioCpp::moiety::CH2, std::list<BioCpp::atom::id>{{ BioCpp::atom::CG, BioCpp::atom::HG1, BioCpp::atom::HG2 }} ) );
      info.push_back( std::make_pair( BioCpp::moiety::CARBOXYL, std::list<BioCpp::atom::id>{{ BioCpp::atom::CD, BioCpp::atom::OE1, BioCpp::atom::OE2 }} ) );
      return info;
    case BioCpp::amino_acid::ASN:
      info.push_back( std::make_pair( BioCpp::moiety::CH2, std::list<BioCpp::atom::id>{{ BioCpp::atom::CB, BioCpp::atom::HB2, BioCpp::atom::HB3 }} ) );
      info.push_back( std::make_pair( BioCpp::moiety::PEPTIDE, std::list<BioCpp::atom::id>{{ BioCpp::atom::CG, BioCpp::atom::OD1, BioCpp::atom::ND2, BioCpp::atom::HD21, BioCpp::atom::HD22 }} ) );
      return info;
    case BioCpp::amino_acid::GLN:
      info.push_back( std::make_pair( BioCpp::moiety::CH2, std::list<BioCpp::atom::id>{{ BioCpp::atom::CB, BioCpp::atom::HB2, BioCpp::atom::HB3 }} ) );
      info.push_back( std::make_pair( BioCpp::moiety::CH2, std::list<BioCpp::atom::id>{{ BioCpp::atom::CG, BioCpp::atom::HG1, BioCpp::atom::HG2 }} ) );
      info.push_back( std::make_pair( BioCpp::moiety::PEPTIDE, std::list<BioCpp::atom::id>{{ BioCpp::atom::CD, BioCpp::atom::OE1, BioCpp::atom::NE2, BioCpp::atom::HE21, BioCpp::atom::HE22 }} ) );
      return info;
    case BioCpp::amino_acid::LYS:
      info.push_back( std::make_pair( BioCpp::moiety::CH2, std::list<BioCpp::atom::id>{{ BioCpp::atom::CB, BioCpp::atom::HB2, BioCpp::atom::HB3 }} ) );
      info.push_back( std::make_pair( BioCpp::moiety::CH2, std::list<BioCpp::atom::id>{{ BioCpp::atom::CG, BioCpp::atom::HG1, BioCpp::atom::HG2 }} ) );
      info.push_back( std::make_pair( BioCpp::moiety::CH2, std::list<BioCpp::atom::id>{{ BioCpp::atom::CD, BioCpp::atom::HD1, BioCpp::atom::HD2 }} ) );
      info.push_back( std::make_pair( BioCpp::moiety::CH2, std::list<BioCpp::atom::id>{{ BioCpp::atom::CE, BioCpp::atom::HE1, BioCpp::atom::HE2 }} ) );
      info.push_back( std::make_pair( BioCpp::moiety::PRT_AMINE, std::list<BioCpp::atom::id>{{ BioCpp::atom::NZ, BioCpp::atom::HZ1, BioCpp::atom::HZ2, BioCpp::atom::HZ3 }} ) );
      return info;
    case BioCpp::amino_acid::ARG:
      info.push_back( std::make_pair( BioCpp::moiety::CH2, std::list<BioCpp::atom::id>{{ BioCpp::atom::CB, BioCpp::atom::HB2, BioCpp::atom::HB3 }} ) );
      info.push_back( std::make_pair( BioCpp::moiety::CH2, std::list<BioCpp::atom::id>{{ BioCpp::atom::CG, BioCpp::atom::HG1, BioCpp::atom::HG2 }} ) );
      info.push_back( std::make_pair( BioCpp::moiety::CH2, std::list<BioCpp::atom::id>{{ BioCpp::atom::CD, BioCpp::atom::HD1, BioCpp::atom::HD2 }} ) );
      info.push_back( std::make_pair( BioCpp::moiety::N3H5, std::list<BioCpp::atom::id>{{ BioCpp::atom::NE, BioCpp::atom::NH1, BioCpp::atom::NH2, BioCpp::atom::HE, BioCpp::atom::HH11, BioCpp::atom::HH12, BioCpp::atom::HH21, BioCpp::atom::HH22 }} ) );
      return info;
    case BioCpp::amino_acid::HIS: 
      info.push_back( std::make_pair( BioCpp::moiety::CH2, std::list<BioCpp::atom::id>{{ BioCpp::atom::CB, BioCpp::atom::HB2, BioCpp::atom::HB3 }} ) );
      info.push_back( std::make_pair( BioCpp::moiety::C3N2H3, std::list<BioCpp::atom::id>{{ BioCpp::atom::CG, BioCpp::atom::ND1, BioCpp::atom::HD1, BioCpp::atom::CE1, BioCpp::atom::HE1, BioCpp::atom::CD2, BioCpp::atom::HD2, BioCpp::atom::NE2, BioCpp::atom::HE2}} ) );
      return info;
    case BioCpp::amino_acid::PHE:
      info.push_back( std::make_pair( BioCpp::moiety::CH2, std::list<BioCpp::atom::id>{{ BioCpp::atom::CB, BioCpp::atom::HB2, BioCpp::atom::HB3 }} ) );
      info.push_back( std::make_pair( BioCpp::moiety::PHENYL, std::list<BioCpp::atom::id>{{ BioCpp::atom::CG, BioCpp::atom::CD1, BioCpp::atom::HD1, BioCpp::atom::CE1, BioCpp::atom::HE1, BioCpp::atom::CD2, BioCpp::atom::HD2, BioCpp::atom::CE2, BioCpp::atom::HE2, BioCpp::atom::CZ, BioCpp::atom::HZ }} ) );
      return info;
    case BioCpp::amino_acid::TYR:
      info.push_back( std::make_pair( BioCpp::moiety::CH2, std::list<BioCpp::atom::id>{{ BioCpp::atom::CB, BioCpp::atom::HB2, BioCpp::atom::HB3 }} ) );
      info.push_back( std::make_pair( BioCpp::moiety::PHENYL, std::list<BioCpp::atom::id>{{ BioCpp::atom::CG, BioCpp::atom::CD1, BioCpp::atom::CE1, BioCpp::atom::CD2, BioCpp::atom::CE2, BioCpp::atom::CZ, BioCpp::atom::HD1, BioCpp::atom::HD2, BioCpp::atom::HE1, BioCpp::atom::HE2 }} ) );
      info.push_back( std::make_pair( BioCpp::moiety::ALCOHOL, std::list<BioCpp::atom::id>{{ BioCpp::atom::OH, BioCpp::atom::HH }} ) );
      return info;
    case BioCpp::amino_acid::TRP:
      info.push_back( std::make_pair( BioCpp::moiety::CH2, std::list<BioCpp::atom::id>{{ BioCpp::atom::CB, BioCpp::atom::HB2, BioCpp::atom::HB3 }} ) );
      info.push_back( std::make_pair( BioCpp::moiety::C8NH6, std::list<BioCpp::atom::id>{{ BioCpp::atom::CG, BioCpp::atom::CD1, BioCpp::atom::NE1, BioCpp::atom::CD2, BioCpp::atom::CE3, BioCpp::atom::CZ3, BioCpp::atom::CE2, BioCpp::atom::CH2, BioCpp::atom::CZ2, BioCpp::atom::HE1, BioCpp::atom::HD1, BioCpp::atom::HH2, BioCpp::atom::HE3, BioCpp::atom::HZ2, BioCpp::atom::HZ3 }} ) );
      return info;
    default:
      return info;
  }
  return info;
};

} // end moiety
} // end BioCpp

std::ostream& operator << (std::ostream& out, const BioCpp::moiety::id moiety_id ) {
  out << BioCpp::moiety::id_to_string[moiety_id];
  return out;
}



#endif

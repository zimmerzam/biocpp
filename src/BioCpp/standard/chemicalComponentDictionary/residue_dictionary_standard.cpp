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

#include "residue_dictionary_standard.hpp"
#include "atom_dictionary_standard.hpp"

namespace BioCpp{
namespace residue{

// initialize default dictionary of elements
dictionary_t dictionary(
  {
    {UNK,"UNK"},{CYS,"CYS"},{PHE,"PHE"},{LEU,"LEU"},{TRP,"TRP"},{VAL,"VAL"},
    {ILE,"ILE"},{MET,"MET"},{HIS,"HIS"},{TYR,"TYR"},{ALA,"ALA"},{GLY,"GLY"},
    {PRO,"PRO"},{ASN,"ASN"},{THR,"THR"},{SER,"SER"},{ARG,"ARG"},{GLN,"GLN"},
    {ASP,"ASP"},{LYS,"LYS"},{GLU,"GLU"}
  },
  {
    {"UNK",UNK},{"CYS",CYS},{"PHE",PHE},{"LEU",LEU},{"TRP",TRP},{"VAL",VAL},
    {"ILE",ILE},{"MET",MET},{"HIS",HIS},{"TYR",TYR},{"ALA",ALA},{"GLY",GLY},
    {"PRO",PRO},{"ASN",ASN},{"THR",THR},{"SER",SER},{"ARG",ARG},{"GLN",GLN},
    {"ASP",ASP},{"LYS",LYS},{"GLU",GLU},
    {"X",UNK},{"C",CYS},{"F",PHE},{"L",LEU},{"W",TRP},{"V",VAL},
    {"I",ILE},{"M",MET},{"H",HIS},{"Y",TYR},{"A",ALA},{"G",GLY},
    {"P",PRO},{"N",ASN},{"T",THR},{"S",SER},{"R",ARG},{"Q",GLN},
    {"D",ASP},{"K",LYS},{"E",GLU}
  },
  {
    {UNK, 
          { 'X',
            {  
              definition_t::model_t( {}, {}, -1, -1 )
            }
          } 
    },
    {CYS,
          { 'C',
            {
              definition_t::model_t(
                {
                  atom::id::N, atom::id::CA, atom::id::C, atom::id::O, atom::id::CB,
                  atom::id::SG, atom::id::H
                },
                {
                  {atom::id::N, atom::id::CA},{atom::id::CA, atom::id::C},{atom::id::C, atom::id::O},
                  {atom::id::CA, atom::id::CB},{atom::id::CB, atom::id::SG}
                },
                atom::id::N, 
                atom::id::C
              )
            }
          }
    },
    {PHE, 
          { 'F',
            {  
              definition_t::model_t(
                { 
                  atom::id::N, atom::id::CA, atom::id::C, atom::id::O, atom::id::CB, 
                  atom::id::CG, atom::id::CD1, atom::id::CD2, atom::id::CE1, atom::id::CE2,
                  atom::id::CZ
                },
                { 
                  {atom::id::N, atom::id::CA},{atom::id::CA, atom::id::C},{atom::id::C, atom::id::O},
                  {atom::id::CA, atom::id::CB},{atom::id::CB, atom::id::CG},{atom::id::CG, atom::id::CD1},
                  {atom::id::CG, atom::id::CD2},{atom::id::CD1, atom::id::CE1},{atom::id::CD2, atom::id::CE2},
                  {atom::id::CE1, atom::id::CZ},{atom::id::CE2, atom::id::CZ}
                },
                atom::id::N, 
                atom::id::C
              )
            }
          } 
    },
    {LEU, 
          { 'L',
            {  
              definition_t::model_t(
                { 
                  atom::id::N, atom::id::CA, atom::id::C, atom::id::O, atom::id::CB, 
                  atom::id::CG, atom::id::CD1, atom::id::CD2
                },
                { 
                  {atom::id::N, atom::id::CA},{atom::id::CA, atom::id::C},{atom::id::C, atom::id::O},
                  {atom::id::CA, atom::id::CB},{atom::id::CB, atom::id::CG},{atom::id::CG, atom::id::CD1},
                  {atom::id::CG, atom::id::CD2}
                },
                atom::id::N, 
                atom::id::C
              )
            }
          } 
    },
    {TRP, 
          { 'W',
            {  
              definition_t::model_t(
                { 
                  atom::id::N, atom::id::CA, atom::id::C, atom::id::O, atom::id::CB, 
                  atom::id::CG, atom::id::CD1, atom::id::CD2, atom::id::NE1, atom::id::CE2,
                  atom::id::CE3, atom::id::CZ2, atom::id::CZ3, atom::id::CH2
                },
                { 
                  {atom::id::N, atom::id::CA},{atom::id::CA, atom::id::C},{atom::id::C, atom::id::O},
                  {atom::id::CA, atom::id::CB},{atom::id::CB, atom::id::CG},{atom::id::CG, atom::id::CD1},
                  {atom::id::CG, atom::id::CD2},{atom::id::CD1, atom::id::NE1},{atom::id::CD2, atom::id::CE2},
                  {atom::id::NE1, atom::id::CE2},{atom::id::CE2, atom::id::CZ2},{atom::id::CZ2, atom::id::CH2},
                  {atom::id::CD2, atom::id::CE3},{atom::id::CE3, atom::id::CZ3},{atom::id::CH2, atom::id::CZ3}
                },
                atom::id::N, 
                atom::id::C
              )
            }
          } 
    },
    {VAL, 
          { 'V',
            {  
              definition_t::model_t(
                { 
                  atom::id::N, atom::id::CA, atom::id::C, atom::id::O, atom::id::CB, 
                  atom::id::CG1, atom::id::CG2
                },
                { 
                  {atom::id::N, atom::id::CA},{atom::id::CA, atom::id::C},{atom::id::C, atom::id::O},
                  {atom::id::CA, atom::id::CB},{atom::id::CB, atom::id::CG1},{atom::id::CB, atom::id::CG2}
                },
                atom::id::N, 
                atom::id::C
              )
            }
          } 
    },
    {ILE, 
          { 'I',
            {  
              definition_t::model_t(
                { 
                  atom::id::N, atom::id::CA, atom::id::C, atom::id::O, atom::id::CB, 
                  atom::id::CG1, atom::id::CD1, atom::id::CG2
                },
                { 
                  {atom::id::N, atom::id::CA},{atom::id::CA, atom::id::C},{atom::id::C, atom::id::O},
                  {atom::id::CA, atom::id::CB},{atom::id::CB, atom::id::CG1},{atom::id::CG1, atom::id::CD1},
                  {atom::id::CB, atom::id::CG2}
                },
                atom::id::N, 
                atom::id::C
              )
            }
          } 
    },
    {MET, 
          { 'M',
            {  
              definition_t::model_t(
                { 
                  atom::id::N, atom::id::CA, atom::id::C, atom::id::O, atom::id::CB, 
                  atom::id::CG, atom::id::SD, atom::id::CE
                },
                { 
                  {atom::id::N, atom::id::CA},{atom::id::CA, atom::id::C},{atom::id::C, atom::id::O},
                  {atom::id::CA, atom::id::CB},{atom::id::CB, atom::id::CG},{atom::id::CG, atom::id::SD},
                  {atom::id::SD, atom::id::CE}
                },
                atom::id::N, 
                atom::id::C
              )
            }
          } 
    },
    {HIS, 
          { 'H',
            {  
              definition_t::model_t(
                { 
                  atom::id::N, atom::id::CA, atom::id::C, atom::id::O, atom::id::CB, 
                  atom::id::OG1, atom::id::CG2
                },
                { 
                  {atom::id::N, atom::id::CA},{atom::id::CA, atom::id::C},{atom::id::C, atom::id::O},
                  {atom::id::CA, atom::id::CB},{atom::id::CB, atom::id::OG1},{atom::id::CB, atom::id::CG2}
                },
                atom::id::N, 
                atom::id::C
              )
            }
          } 
    },
    {TYR, 
          { 'Y',
            {  
              definition_t::model_t(
                { 
                  atom::id::N, atom::id::CA, atom::id::C, atom::id::O, atom::id::CB, 
                  atom::id::CG, atom::id::CD1, atom::id::CD2, atom::id::CE1, atom::id::CE2,
                  atom::id::CZ, atom::id::OH
                },
                { 
                  {atom::id::N, atom::id::CA},{atom::id::CA, atom::id::C},{atom::id::C, atom::id::O},
                  {atom::id::CA, atom::id::CB},{atom::id::CB, atom::id::CG},{atom::id::CG, atom::id::CD1},
                  {atom::id::CG, atom::id::CD2},{atom::id::CD1, atom::id::CE1},{atom::id::CD2, atom::id::CE2},
                  {atom::id::CE1, atom::id::CZ},{atom::id::CE2, atom::id::CZ},{atom::id::CZ, atom::id::OH}
                },
                atom::id::N, 
                atom::id::C
              )
            }
          } 
    },
    {ALA, 
          { 'A',
            {  
              definition_t::model_t(
                { 
                  atom::id::N, atom::id::CA, atom::id::C, atom::id::O, atom::id::CB 
                },
                { 
                  {atom::id::N, atom::id::CA},{atom::id::CA, atom::id::C},{atom::id::C, atom::id::O},
                  {atom::id::CA, atom::id::CB}
                },
                atom::id::N, 
                atom::id::C
              )
            }
          } 
    },
    {GLY, 
          { 'G',
            {  
              definition_t::model_t(
                { 
                  atom::id::N, atom::id::CA, atom::id::C, atom::id::O
                },
                { 
                  {atom::id::N, atom::id::CA},{atom::id::CA, atom::id::C},{atom::id::C, atom::id::O}
                },
                atom::id::N, 
                atom::id::C
              )
            }
          } 
    },
//    {PRO,{0}},
    {ASN, 
          { 'N',
            {  
              definition_t::model_t(
                { 
                  atom::id::N, atom::id::CA, atom::id::C, atom::id::O, atom::id::CB, 
                  atom::id::CG, atom::id::OD1, atom::id::ND2
                },
                { 
                  {atom::id::N, atom::id::CA},{atom::id::CA, atom::id::C},{atom::id::C, atom::id::O},
                  {atom::id::CA, atom::id::CB},{atom::id::CB, atom::id::CG},{atom::id::CG, atom::id::OD1},
                  {atom::id::CG, atom::id::ND2}
                },
                atom::id::N, 
                atom::id::C
              )
            }
          } 
    },
    {THR, 
          { 'T',
            {  
              definition_t::model_t(
                { 
                  atom::id::N, atom::id::CA, atom::id::C, atom::id::O, atom::id::CB, 
                  atom::id::OG1, atom::id::CG2
                },
                { 
                  {atom::id::N, atom::id::CA},{atom::id::CA, atom::id::C},{atom::id::C, atom::id::O},
                  {atom::id::CA, atom::id::CB},{atom::id::CB, atom::id::OG1},{atom::id::CB, atom::id::CG2}
                },
                atom::id::N, 
                atom::id::C
              )
            }
          } 
    },
    {SER, 
          { 'S',
            {  
              definition_t::model_t(
                { 
                  atom::id::N, atom::id::CA, atom::id::C, atom::id::O, atom::id::CB, 
                  atom::id::OG
                },
                { 
                  {atom::id::N, atom::id::CA},{atom::id::CA, atom::id::C},{atom::id::C, atom::id::O},
                  {atom::id::CA, atom::id::CB},{atom::id::CB, atom::id::OG}
                },
                atom::id::N, 
                atom::id::C
              )
            }
          } 
    },
    {ARG, 
          { 'R',
            {  
              definition_t::model_t(
                { 
                  atom::id::N, atom::id::CA, atom::id::C, atom::id::O, atom::id::CB, 
                  atom::id::CG, atom::id::CD, atom::id::NE, atom::id::CZ, atom::id::NH1, 
                  atom::id::NH2
                },
                { 
                  {atom::id::N, atom::id::CA},{atom::id::CA, atom::id::C},{atom::id::C, atom::id::O},
                  {atom::id::CA, atom::id::CB},{atom::id::CB, atom::id::CG},{atom::id::CG, atom::id::CD},
                  {atom::id::CD, atom::id::NE},{atom::id::NE, atom::id::CZ},{atom::id::CZ, atom::id::NH1},
                  {atom::id::CZ, atom::id::NH2}
                },
                atom::id::N, 
                atom::id::C
              )
            }
          } 
    },
    {GLN, 
          { 'Q',
            {  
              definition_t::model_t(
                { 
                  atom::id::N, atom::id::CA, atom::id::C, atom::id::O, atom::id::CB, 
                  atom::id::CG, atom::id::CD, atom::id::OE1, atom::id::NE2
                },
                { 
                  {atom::id::N, atom::id::CA},{atom::id::CA, atom::id::C},{atom::id::C, atom::id::O},
                  {atom::id::CA, atom::id::CB},{atom::id::CB, atom::id::CG},{atom::id::CG, atom::id::CD},
                  {atom::id::CD, atom::id::OE1},{atom::id::CD, atom::id::NE2}
                },
                atom::id::N, 
                atom::id::C
              )
            }
          } 
    },
    {ASP, 
          { 'D',
            {  
              definition_t::model_t(
                { 
                  atom::id::N, atom::id::CA, atom::id::C, atom::id::O, atom::id::CB, 
                  atom::id::CG, atom::id::OD1, atom::id::OD2
                },
                { 
                  {atom::id::N, atom::id::CA},{atom::id::CA, atom::id::C},{atom::id::C, atom::id::O},
                  {atom::id::CA, atom::id::CB},{atom::id::CB, atom::id::CG},{atom::id::CG, atom::id::OD1},
                  {atom::id::CG, atom::id::OD2}
                },
                atom::id::N, 
                atom::id::C
              )
            }
          } 
    },
    {LYS, 
          { 'K',
            {  
              definition_t::model_t(
                { 
                  atom::id::N, atom::id::CA, atom::id::C, atom::id::O, atom::id::CB, 
                  atom::id::CG, atom::id::CD, atom::id::CE, atom::id::NZ
                },
                { 
                  {atom::id::N, atom::id::CA},{atom::id::CA, atom::id::C},{atom::id::C, atom::id::O},
                  {atom::id::CA, atom::id::CB},{atom::id::CB, atom::id::CG},{atom::id::CG, atom::id::CD},
                  {atom::id::CD, atom::id::CE},{atom::id::CE, atom::id::NZ}
                },
                atom::id::N, 
                atom::id::C
              )
            }
          } 
    },
    {GLU, 
          { 'E',
            {  
              definition_t::model_t(
                { 
                  atom::id::N, atom::id::CA, atom::id::C, atom::id::O, atom::id::CB, 
                  atom::id::CG, atom::id::CD, atom::id::OE1, atom::id::OE2
                },
                { 
                  {atom::id::N, atom::id::CA},{atom::id::CA, atom::id::C},{atom::id::C, atom::id::O},
                  {atom::id::CA, atom::id::CB},{atom::id::CB, atom::id::CG},{atom::id::CG, atom::id::CD},
                  {atom::id::CD, atom::id::OE1},{atom::id::CD, atom::id::OE2}
                },
                atom::id::N, 
                atom::id::C
              )
            }
          } 
    }
  }
);

}; // end namespace residue
};

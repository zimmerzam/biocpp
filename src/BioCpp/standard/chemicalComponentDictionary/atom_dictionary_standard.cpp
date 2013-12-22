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

#include "atom_dictionary_standard.hpp"
#include "element_dictionary_standard.hpp"

namespace BioCpp{
namespace atom{

// initialize default dictionary of atoms
dictionary_t dictionary(
  // { ID, string }
  {
    /* Unknown Atoms */
    { id::X    , " X  "},
    /* Backbone Atoms */
    { id::N   , " N  "},{ id::CA  , " CA "},{ id::C   , " C  "},{ id::O   , " O  "},
    { id::OXT , " OXT"},{ id::H   , " H  "},{ id::HA  , " HA "},{ id::HA1 , " HA1"},
    { id::HA2 , " HA2"},{ id::HA3 , " HA3"},{ id::H1  , " H1 "},{ id::H2  , " H2 "},
    { id::H3  , " H3 "},
    /* Side Chain Atoms : Carbon */
    { id::CB  , " CB "},{ id::CD  , " CD "},{ id::CD1 , " CD1"},{ id::CD2 , " CD2"},
    { id::CD3 , " CD3"},{ id::CE  , " CE "},{ id::CE1 , " CE1"},{ id::CE2 , " CE2"},
    { id::CE3 , " CE3"},{ id::CG  , " CG "},{ id::CG1 , " CG1"},{ id::CG2 , " CG2"},
    { id::CG3 , " CG3"},{ id::CZ  , " CZ "},{ id::CZ1 , " CZ1"},{ id::CZ2 , " CZ2"},
    { id::CZ3 , " CZ3"},{ id::CH1 , " CH1"},{ id::CH2 , " CH2"},
    /* Side Chain Atoms : Nitrogen */
    { id::ND  , " ND "},{ id::ND1 , " ND1"},{ id::ND2 , " ND2"},{ id::NE  , " NE "},
    { id::NE1 , " NE1"},{ id::NE2 , " NE2"},{ id::NH  , " NH "},{ id::NH1 , " NH1"},
    { id::NH2 , " NH2"},{ id::NZ  , " NZ "},
    /* Side Chain Atoms : Oxygen */
    { id::OD  , " OD "},{ id::OD1 , " OD1"},{ id::OD2 , " OD2"},{ id::OE  , " OE "},
    { id::OE1 , " OE1"},{ id::OE2 , " OE2"},{ id::OH  , " OH "},{ id::OH1 , " OH1"},
    { id::OH2 , " OH2"},{ id::OG  , " OG "},{ id::OG1 , " OG1"},{ id::OG2 , " OG2"},
    /* Side Chain Atoms : Sulfur */
    { id::SD  , " SD "},{ id::SG  , " SG "},
    /* Side Chain Atoms : Hydrogen */
    { id::HB  , " HB "},{ id::HB1 , " HB1"},{ id::HB2 , " HB2"},{ id::HB3 , " HB3"},
    { id::HD  , " HD "},{ id::HD1 , " HD1"},{ id::HD2 , " HD2"},{ id::HD3 , " HD3"},
    { id::HD11, "HD11"},{ id::HD12, "HD12"},{ id::HD13, "HD13"},{ id::HD21, "HD21"},
    { id::HD22, "HD22"},{ id::HD23, "HD23"},{ id::HE  , " HE "},{ id::HE1 , " HE1"},
    { id::HE2 , " HE2"},{ id::HE3 , " HE3"},{ id::HE11, "HE11"},{ id::HE12, "HE12"},
    { id::HE13, "HE13"},{ id::HE21, "HE21"},{ id::HE22, "HE22"},{ id::HE23, "HE23"},
    { id::HG  , " HG "},{ id::HG1 , " HG1"},{ id::HG2 , " HG2"},{ id::HG3 , " HG3"},
    { id::HG11, "HG11"},{ id::HG12, "HG12"},{ id::HG13, "HG13"},{ id::HG21, "HG21"},
    { id::HG22, "HG22"},{ id::HG23, "HG23"},{ id::HH  , " HH "},{ id::HH1 , " HH1"},
    { id::HH2 , " HH2"},{ id::HH3 , " HH3"},{ id::HH11, "HH11"},{ id::HH12, "HH12"},
    { id::HH13, "HH13"},{ id::HH21, "HH21"},{ id::HH22, "HH22"},{ id::HH23, "HH23"},
    { id::HN  , " HN "},{ id::HN1 , " HN1"},{ id::HN2 , " HN2"},{ id::HN3 , " HN3"},
    { id::HN11, "HN11"},{ id::HN12, "HN12"},{ id::HN13, "HN13"},{ id::HN21, "HN21"},
    { id::HN22, "HN22"},{ id::HN23, "HN23"},{ id::HZ  , " HZ "},{ id::HZ1 , " HZ1"},
    { id::HZ2 , " HZ2"},{ id::HZ3 , " HZ3"}
  },
  // { string, ID }
  {
    /* Unknown Atoms */
    { " X  ", id::X   },
    /* Backbone Atoms */
    { " N  ", id::N   },{ " CA ", id::CA  },{ " C  ", id::C   },{ " O  ", id::O   },
    { " OXT", id::OXT },{ " H  ", id::H   },{ " HA ", id::HA  },{ " HA1", id::HA1 },
    { "1HA ", id::HA1 },{ " HA2", id::HA2 },{ "2HA ", id::HA2 },{ " HA3", id::HA3 },
    { "3HA ", id::HA3 },{ " H1 ", id::H1  },{ " 1H ", id::H1  },{ " H2 ", id::H2  },
    { " 2H ", id::H2  },{ " H3 ", id::H3  },{ " 3H ", id::H3  },
    /* Side Chain Atoms : Carbon */
    { " CB ", id::CB  },{ " CD ", id::CD  },{ " CD1", id::CD1 },{ " CD2", id::CD2 },
    { " CD3", id::CD3 },{ " CE ", id::CE  },{ " CE1", id::CE1 },{ " CE2", id::CE2 },
    { " CE3", id::CE3 },{ " CG ", id::CG  },{ " CG1", id::CG1 },{ " CG2", id::CG2 },
    { " CG3", id::CG3 },{ " CZ ", id::CZ  },{ " CZ1", id::CZ1 },{ " CZ2", id::CZ2 },
    { " CZ3", id::CZ3 },{ " CH1", id::CH1 },{ " CH2", id::CH2 },
    /* Side Chain Atoms : Nitrogen */
    { " ND ", id::ND  },{ " ND1", id::ND1 },{ " ND2", id::ND2 },{ " NE ", id::NE  },
    { " NE1", id::NE1 },{ " NE2", id::NE2 },{ " NH ", id::NH  },{ " NH1", id::NH1 },
    { " NH2", id::NH2 },{ " NZ ", id::NZ  },
    /* Side Chain Atoms : Oxygen */
    { " OD ", id::OD  },{ " OD1", id::OD1 },{ " OD2", id::OD2 },{ " OE ", id::OE  },
    { " OE1", id::OE1 },{ " OE2", id::OE2 },{ " OH ", id::OH  },{ " OH1", id::OH1 },
    { " OH2", id::OH2 },{ " OG ", id::OG  },{ " OG1", id::OG1 },{ " OG2", id::OG2 },
    /* Side Chain Atoms : Sulfur */
    { " SD ", id::SD  },{ " SG ", id::SG},
    /* Side Chain Atoms : Hydrogen */
    { " HB ", id::HB  },{ " HB1", id::HB1 },{ "1HB ", id::HB1 },{ " HB2", id::HB2 },
    { "2HB ", id::HB2 },{ " HB3", id::HB3 },{ "3HB ", id::HB3 },{ " HD ", id::HD  },
    { " HD1", id::HD1 },{ "1HD ", id::HD1 },{ " HD2", id::HD2 },{ "2HD ", id::HD2 },
    { " HD3", id::HD3 },{ "3HD ", id::HD3 },{ "HD11", id::HD11},{ "1HD1", id::HD11},
    { "HD12", id::HD12},{ "2HD1", id::HD12},{ "HD13", id::HD13},{ "3HD1", id::HD13},
    { "HD21", id::HD21},{ "1HD2", id::HD21},{ "HD22", id::HD22},{ "2HD2", id::HD22},
    { "HD23", id::HD23},{ "3HD2", id::HD23},{ " HE ", id::HE  },{ " HE1", id::HE1 },
    { "1HE ", id::HE1 },{ " HE2", id::HE2 },{ "2HE ", id::HE2 },{ " HE3", id::HE3 },
    { "3HE ", id::HE3 },{ "HE11", id::HE11},{ "1HE1", id::HE11},{ "HE12", id::HE12},
    { "2HE1", id::HE12},{ "HE13", id::HE13},{ "3HE1", id::HE13},{ "HE21", id::HE21},
    { "1HE2", id::HE21},{ "HE22", id::HE22},{ "2HE2", id::HE22},{ "HE23", id::HE23},
    { "3HE2", id::HE23},{ " HG ", id::HG  },{ " HG1", id::HG1 },{ "1HG ", id::HG1 },
    { " HG2", id::HG2 },{ "2HG ", id::HG2 },{ " HG3", id::HG3 },{ "3HG ", id::HG3 },
    { "HG11", id::HG11},{ "1HG1", id::HG11},{ "HG12", id::HG12},{ "2HG1", id::HG12},
    { "HG13", id::HG13},{ "3HG1", id::HG13},{ "HG21", id::HG21},{ "1HG2", id::HG21},
    { "HG22", id::HG22},{ "2HG2", id::HG22},{ "HG23", id::HG23},{ "3HG2", id::HG23},
    { " HH ", id::HH  },{ " HH1", id::HH1 },{ "1HH ", id::HH1 },{ " HH2", id::HH2 },
    { "2HH ", id::HH2 },{ " HH3", id::HH3 },{ "3HH ", id::HH3 },{ "HH11", id::HH11},
    { "1HH1", id::HH11},{ "HH12", id::HH12},{ "2HH1", id::HH12},{ "HH13", id::HH13},
    { "3HH1", id::HH13},{ "HH21", id::HH21},{ "1HH2", id::HH21},{ "HH22", id::HH22},
    { "2HH2", id::HH22},{ "HH23", id::HH23},{ "3HH2", id::HH23},{ " HN ", id::HN  },
    { " HN1", id::HN1 },{ "1HN ", id::HN1 },{ " HN2", id::HN2 },{ "2HN ", id::HN2 },
    { " HN3", id::HN3 },{ "3HN ", id::HN3 },{ "HN11", id::HN11},{ "1HN1", id::HN11},
    { "HN12", id::HN12},{ "2HN1", id::HN12},{ "HN13", id::HN13},{ "3HN1", id::HN13},
    { "HN21", id::HN21},{ "1HN2", id::HN21},{ "HN22", id::HN22},{ "2HN2", id::HN22},
    { "HN23", id::HN23},{ "3HN2", id::HN23},{ " HZ ", id::HZ  },{ " HZ1", id::HZ1 },
    { "1HZ ", id::HZ1 },{ " HZ2", id::HZ2 },{ "2HZ ", id::HZ2 },{ " HZ3", id::HZ3 },
    { "3HZ ", id::HZ3 }
  },
  // { ID, element ID }
  { 
    /* Unknown Atoms */
    { id::X   , element::id::X },
    /* Backbone Atoms */
    { id::N   , element::id::N },{ id::CA  , element::id::C },{ id::C   , element::id::C },
    { id::O   , element::id::O },{ id::OXT , element::id::O },{ id::H   , element::id::H },
    { id::HA  , element::id::H },{ id::HA1 , element::id::H },{ id::HA2 , element::id::H },
    { id::HA3 , element::id::H },{ id::H1  , element::id::H },{ id::H2  , element::id::H },
    { id::H3  , element::id::H },
    /* Side Chain Atoms : Carbon */
    { id::CB  , element::id::C },{ id::CD  , element::id::C },{ id::CD1 , element::id::C },
    { id::CD2 , element::id::C },{ id::CD3 , element::id::C },{ id::CE  , element::id::C },
    { id::CE1 , element::id::C },{ id::CE2 , element::id::C },{ id::CE3 , element::id::C },
    { id::CG  , element::id::C },{ id::CG1 , element::id::C },{ id::CG2 , element::id::C },
    { id::CG3 , element::id::C },{ id::CZ  , element::id::C },{ id::CZ1 , element::id::C },
    { id::CZ2 , element::id::C },{ id::CZ3 , element::id::C },{ id::CH1 , element::id::C },
    { id::CH2 , element::id::C },
    /* Side Chain Atoms : Nitrogen */
    { id::ND  , element::id::N },{ id::ND1 , element::id::N },{ id::ND2 , element::id::N },
    { id::NE  , element::id::N },{ id::NE1 , element::id::N },{ id::NE2 , element::id::N },
    { id::NH  , element::id::N },{ id::NH1 , element::id::N },{ id::NH2 , element::id::N },
    { id::NZ  , element::id::N },
    /* Side Chain Atoms : Oxygen */
    { id::OD  , element::id::O },{ id::OD1 , element::id::O },{ id::OD2 , element::id::O },
    { id::OE  , element::id::O },{ id::OE1 , element::id::O },{ id::OE2 , element::id::O },
    { id::OH  , element::id::O },{ id::OH1 , element::id::O },{ id::OH2 , element::id::O },
    { id::OG  , element::id::O },{ id::OG1 , element::id::O },{ id::OG2 , element::id::O },
    /* Side Chain Atoms : Sulfur */
    { id::SD  , element::id::S },{ id::SG  , element::id::S },
    /* Side Chain Atoms : Hydrogen */
    { id::HB  , element::id::H },{ id::HB1 , element::id::H },{ id::HB2 , element::id::H },
    { id::HB3 , element::id::H },{ id::HD  , element::id::H },{ id::HD1 , element::id::H },
    { id::HD2 , element::id::H },{ id::HD3 , element::id::H },{ id::HD11, element::id::H },
    { id::HD12, element::id::H },{ id::HD13, element::id::H },{ id::HD21, element::id::H },
    { id::HD22, element::id::H },{ id::HD23, element::id::H },{ id::HE  , element::id::H },
    { id::HE1 , element::id::H },{ id::HE2 , element::id::H },{ id::HE3 , element::id::H },
    { id::HE11, element::id::H },{ id::HE12, element::id::H },{ id::HE13, element::id::H },
    { id::HE21, element::id::H },{ id::HE22, element::id::H },{ id::HE23, element::id::H },
    { id::HG  , element::id::H },{ id::HG1 , element::id::H },{ id::HG2 , element::id::H },
    { id::HG3 , element::id::H },{ id::HG11, element::id::H },{ id::HG12, element::id::H },
    { id::HG13, element::id::H },{ id::HG21, element::id::H },{ id::HG22, element::id::H },
    { id::HG23, element::id::H },{ id::HH  , element::id::H },{ id::HH1 , element::id::H },
    { id::HH2 , element::id::H },{ id::HH3 , element::id::H },{ id::HH11, element::id::H },
    { id::HH12, element::id::H },{ id::HH13, element::id::H },{ id::HH21, element::id::H },
    { id::HH22, element::id::H },{ id::HH23, element::id::H },{ id::HN  , element::id::H },
    { id::HN1 , element::id::H },{ id::HN2 , element::id::H },{ id::HN3 , element::id::H },
    { id::HN11, element::id::H },{ id::HN12, element::id::H },{ id::HN13, element::id::H },
    { id::HN21, element::id::H },{ id::HN22, element::id::H },{ id::HN23, element::id::H },
    { id::HZ  , element::id::H },{ id::HZ1 , element::id::H },{ id::HZ2 , element::id::H },
    { id::HZ3 , element::id::H }
  }
);

} // end namespace atom
}

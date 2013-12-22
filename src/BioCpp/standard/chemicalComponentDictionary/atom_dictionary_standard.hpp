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

#ifndef BIOCPP_PDB_ATOM_DICTIONARY_STANDARD
#define BIOCPP_PDB_ATOM_DICTIONARY_STANDARD

#include <BioCpp/chemical_component_dictionary/atoms/atom_dictionary.hpp>

namespace BioCpp{
namespace atom{

// atom id is the order of apparance in pdb files (or casual)
enum id : signed int{
  /* Unknown Atoms */
  X  = -1  ,
  /* Backbone Atoms */
  N    =  1  , CA   =  2  , C    =  3  , O    =  4  ,
  OXT  =  7  , H    =  6  , HA   =  7  , HA1  =  8  ,
  HA2  =  9  , HA3  =  10 , H1   =  11 , H2   =  12 ,
  H3   =  13 ,
  /* Side Chain Atoms : Carbon */
  CB   =  14 , CD   =  15 , CD1  =  16 , CD2  =  17 , CD3  =  18 , CE   =  19 , 
  CE1  =  20 , CE2  =  21 , CE3  =  22 , CG   =  23 , CG1  =  24 , CG2  =  25 , 
  CG3  =  26 , CZ   =  26 , CZ1  =  28 , CZ2  =  29 , CZ3  =  30 , CH1  =  31 , 
  CH2  =  32 , 
  /* Side Chain Atoms : Nitrogen */
  ND   =  33 , ND1  =  34 , ND2  =  35 , NE   =  36 , NE1  =  37 , NE2  =  38 , 
  NH   =  39 , NH1  =  40 , NH2  =  41 , NZ   =  42 ,
  /* Side Chain Atoms : Oxygen */
  OD   =  43 , OD1  =  44 , OD2  =  45 , OE   =  46 , OE1  =  47 , OE2  =  48 , 
  OH   =  49 , OH1  =  50 , OH2  =  51 , OG   =  52 , OG1  =  53 , OG2  =  54 ,
  /* Side Chain Atoms : Sulfur */
  SD   =  56 , SG   =  57 ,
  /* Side Chain Atoms : Hydrogen */
  HB   =  58 , HB1  =  59 , HB2  =  60 , HB3  =  61 , HD   =  62 , HD1  =  63 ,
  HD2  =  64 , HD3  =  65 , HD11 =  66 , HD12 =  67 , HD13 =  68 , HD21 =  69 , 
  HD22 =  70 , HD23 =  71 , HE   =  72 , HE1  =  73 , HE2  =  74 , HE3  =  75 , 
  HE11 =  76 , HE12 =  77 , HE13 =  78 , HE21 =  79 , HE22 =  80 , HE23 =  81 ,
  HG   =  82 , HG1  =  83 , HG2  =  84 , HG3  =  85 , HG11 =  86 , HG12 =  87 , 
  HG13 =  88 , HG21 =  89 , HG22 =  90 , HG23 =  91 , HH   =  92 , HH1  =  93 , 
  HH2  =  94 , HH3  =  95 , HH11 =  96 , HH12 =  97 , HH13 =  98 , HH21 =  99 , 
  HH22 = 100 , HH23 = 101 , HN   = 102 , HN1  = 103 , HN2  = 104 , HN3  = 105 , 
  HN11 = 106 , HN12 = 107 , HN13 = 108 , HN21 = 109 , HN22 = 110 , HN23 = 111 , 
  HZ   = 114 , HZ1  = 115 , HZ2  = 116 , HZ3  = 117
};

// initialize default dictionary of atoms
extern dictionary_t dictionary;

} // end namespace atom
}

#endif

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

#ifndef BIOCPP_PDB_RESIDUE_DICTIONARY_STANDARD
#define BIOCPP_PDB_RESIDUE_DICTIONARY_STANDARD

#include <BioCpp/chemical_component_dictionary/residues/residue_dictionary.hpp>

namespace BioCpp{
namespace residue{

// in order of ??? -> ask to Antonio
enum id : signed int{
  UNK = -1 , /* Unknown residue */
  CYS =  1 , /* Cysteine */
  PHE =  2 , /* Phenylalanine */
  LEU =  3 , /* Leucine */
  TRP =  4 , /* Tryptophan */
  VAL =  5 , /* Valine */
  ILE =  6 , /* Isoleucine */
  MET =  7 , /* Methionine */
  HIS =  8 , /* Histidine */
  TYR =  9 , /* Tyrosine */
  ALA =  10, /* Alanine */
  GLY =  11, /* Glycine */
  PRO =  12, /* Proline */
  ASN =  13, /* Asparagine */
  THR =  14, /* Threonine */
  SER =  15, /* Serine */
  ARG =  16, /* Arginine */
  GLN =  17, /* Glutamine */
  ASP =  18, /* Aspartic Acid */
  LYS =  19, /* Lysine */
  GLU =  20  /* Glutamic Acid */
};

// initialize default dictionary of elements
extern dictionary_t dictionary;

}; // end namespace residue
};
#endif

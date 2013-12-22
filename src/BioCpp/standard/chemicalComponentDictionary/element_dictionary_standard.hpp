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

#ifndef BIOCPP_PDB_ELEMENT_DICTIONARY_STANDARD
#define BIOCPP_PDB_ELEMENT_DICTIONARY_STANDARD

#include <BioCpp/chemical_component_dictionary/elements/element_dictionary.hpp>

namespace BioCpp{
namespace element{

// element id is the atomic number Z
enum id : signed int {
  X = -1,  /* Unknown element */
  H =  1,  /* Hydrogen */
  C =  6,  /* Carbon */
  N =  7,  /* Nitrogen */
  O =  8,  /* Oxygen */
  S =  16  /* Sulfur */
};

// initialize default dictionary of elements
extern dictionary_t dictionary;

} // end namespace element
}

#endif

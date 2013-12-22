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

#ifndef BIOCPP_PDB_STANDARD_PRINTATOM
#define BIOCPP_PDB_STANDARD_PRINTATOM

#include "../atom/print_atom_t.hxx"
#include "element_dictionary_standard.hpp"
#include "atom_dictionary_standard.hpp"
#include "residue_dictionary_standard.hpp"
#include <iostream>

namespace BioCpp{
namespace pdb{

print_atom_t< BioCpp::element::dictionary_t,
              BioCpp::atom::dictionary_t,
              BioCpp::residue::dictionary_t 
            > print_atom( std::cout, 
                          BioCpp::element::dictionary, 
                          BioCpp::atom::dictionary, 
                          BioCpp::residue::dictionary
                        );

}
}

#endif

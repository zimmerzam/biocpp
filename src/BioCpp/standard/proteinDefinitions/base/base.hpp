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

#ifndef BIOCPP_STANDARD_BASE_DEFINITION_H
#define BIOCPP_STANDARD_BASE_DEFINITION_H

#include <BioCpp/base_atom/base_atom.hpp>
#include <BioCpp/standard/proteinDefinitions/templated.hxx>
#include <BioCpp/standard/chemicalComponentDictionary/residue_dictionary_standard.hpp>

namespace BioCpp{
namespace standard{
namespace base{

typedef BioCpp::base::atom atom;
typedef BioCpp::standard::model<atom>::type model;
typedef BioCpp::standard::residue<atom>::type residue;
typedef BioCpp::standard::chain<atom>::type chain;
typedef BioCpp::standard::complex<atom>::type complex;
typedef BioCpp::standard::complex_constructor<atom, BioCpp::residue::dictionary_t> complex_constructor;

}
}
}

#endif

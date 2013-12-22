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

#ifndef BIOCPP_PDB_STANDARD_ATOM_DICTIONARY_H
#define BIOCPP_PDB_STANDARD_ATOM_DICTIONARY_H

#include "../dictionary.hxx"
#include "../definition.hpp"

namespace BioCpp{
namespace atom{

class definition_t : public BioCpp::definition{
  public:
    int element;
    definition_t();
    definition_t( int i );
    
    void importSetting(libconfig::Setting& setting);
};

typedef dictionary< definition_t> dictionary_t;

}; // end namespace atom
};
#endif

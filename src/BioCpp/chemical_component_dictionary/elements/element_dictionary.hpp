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

#ifndef BIOCPP_PDB_STANDARD_ELEMENT_DICTIONARY_H
#define BIOCPP_PDB_STANDARD_ELEMENT_DICTIONARY_H

#include "../dictionary.hxx"
#include "../definition.hpp"
#include <array>
#include <map>

namespace BioCpp{
namespace element{

class definition_t : public BioCpp::definition{
  public:
    double vanderwaals_radius;
    double covalent_radius;
    double mass;
    definition_t(  );
    definition_t( double vdwr, double covr, double m);
    definition_t( std::array< double, 3 > info );
    
    void importSetting(libconfig::Setting& setting);
};

typedef BioCpp::dictionary< definition_t > dictionary_t;

}; // end namespace element
};
#endif

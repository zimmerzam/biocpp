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

#ifndef BIOCPP_RESIDUE_DICTIONARY_H
#define BIOCPP_RESIDUE_DICTIONARY_H

#include "../dictionary.hxx"
#include "../definition.hpp"
#include <map>
#include <set>
#include <string>

namespace BioCpp{
namespace residue{

class definition_t : public BioCpp::definition{
  public:
    class model_t{
      public:
        typedef std::set< int > atom_list_t;
        typedef std::set< std::pair<int,int> > bond_list_t;
        atom_list_t atom_list;
        bond_list_t bond_list;
        int atom_start;
        int atom_end;
        
        model_t();
        model_t( atom_list_t al, bond_list_t bl, int as, int ae );
    
        void importSetting(libconfig::Setting& setting);
    };
		
		std::string one_letter_name; 
    std::list<model_t> model;   
    definition_t();    
    definition_t(std::string n, std::list<model_t> m );
    void importSetting(libconfig::Setting& setting);
    
};

typedef BioCpp::dictionary< definition_t> dictionary_t;

}; // end namespace element
};
#endif

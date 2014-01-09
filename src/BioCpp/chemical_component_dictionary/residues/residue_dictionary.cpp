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

#include "residue_dictionary.hpp"

BioCpp::residue::definition_t::model_t::model_t(){};
BioCpp::residue::definition_t::model_t::model_t( atom_list_t al, bond_list_t bl, int as, int ae ) : 
                                      atom_list(al), 
                                      bond_list(bl), 
                                      atom_start(as), 
                                      atom_end(ae){};

void BioCpp::residue::definition_t::model_t::importSetting(libconfig::Setting& setting){
  setting.lookupValue("atom_start", atom_start);
  setting.lookupValue("atom_end", atom_end);
  libconfig::Setting& atmlist = setting["atom_list"];
  int atmlist_size = atmlist.getLength();
  for(int i =0; i!= atmlist_size; ++i){
    atom_list.insert( int(atmlist[i]) );
  }
  libconfig::Setting& bndlist = setting["bond_list"];
  int bndlist_size = bndlist.getLength();
  for(int i =0; i!= bndlist_size; ++i){
    bond_list.insert( std::make_pair( int(bndlist[i][0]), int(bndlist[i][1]) ) );
  }
};

BioCpp::residue::definition_t::definition_t(){};    
BioCpp::residue::definition_t::definition_t( char n, std::list<model_t> m ) : 
    one_letter_name(n), model(m) {};

void BioCpp::residue::definition_t::importSetting(libconfig::Setting& setting){
  setting.lookupValue("one_letter_name", one_letter_name);
  libconfig::Setting& mdllist = setting["model_list"];
  int mdl_size = mdllist.getLength();
  for(int i = 0; i != mdl_size; ++i){
    model.push_back( model_t() );
    model.back().importSetting( mdllist[i] );
  }
}

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

#include "element_dictionary.hpp"

void BioCpp::element::definition_t::importSetting(libconfig::Setting& setting){
  setting.lookupValue("vanderwaals_radius", vanderwaals_radius);
  setting.lookupValue("covalent_radius", covalent_radius);
  setting.lookupValue("mass", mass);
};

BioCpp::element::definition_t::definition_t(): vanderwaals_radius(-1.), covalent_radius(-1.), mass(-1.) {};
BioCpp::element::definition_t::definition_t( double vdwr, double covr, double m) : vanderwaals_radius(vdwr), covalent_radius(covr), mass(m) {};
BioCpp::element::definition_t::definition_t( std::array< double, 3 > info ): vanderwaals_radius(info[0]), covalent_radius(info[1]), mass(info[2]) {};

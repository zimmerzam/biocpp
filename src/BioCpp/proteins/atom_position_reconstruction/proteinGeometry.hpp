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

#ifndef PROTEIN_GEOMETRY_PARAMETERS_H
#define PROTEIN_GEOMETRY_PARAMETERS_H

namespace BioCpp{
namespace protein{

const double N_CA_bond_length = 1.449;
const double CA_C_bond_length = 1.522;
const double C_N_bond_length = 1.335;

const double N_H_bond_length = 1.0;
const double C_O_bond_length = 1.229;
const double CA_CB_bond_length = 1.526;

const double cos_C_N_CA_angle = -0.52843833472;  // angle is 121.9
const double cos_N_CA_C_angle = -0.34365969458;  // angle is 110.1
const double cos_CA_C_N_angle = -0.44775908783;  // angle is 116.6 
const double cos_CA_C_O_angle = -0.50151073715;  // angle is 120.1
const double cos_O_C_N_angle = -0.54024032047;    // angle is 122.7
const double cos_N_CA_CB_angle = -0.33380685923; // angle is 109.5
const double cos_CB_CA_C_angle = -0.35999680812; // angle is 111.1

} // end namespace
}
#endif

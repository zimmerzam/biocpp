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

#include "element_dictionary_standard.hpp"

namespace BioCpp{
namespace element{

// initialize default dictionary of elements
dictionary_t dictionary(
  // { ID, string }
  {
    {id::H, " H" },
    {id::C, " C" },
    {id::N, " N" },
    {id::O, " O" },
    {id::S, " S" }
  },
  // { string, ID }
  {
    { " H", id::H },
    { " C", id::C },
    { " N", id::N },
    { " O", id::O },
    { " S", id::S }
  },
  // { ID, {van der Waals radius, covalent radius, mass} }
  { 
    {id::H, {1.20, 0.31, 1.008 } },
    {id::C, {1.70, 0.76, 12.011} },
    {id::N, {1.55, 0.71, 14.007} },
    {id::O, {1.52, 0.66, 15.999} },
    {id::S, {1.80, 1.05, 32.06 } }
  }
);

} // end namespace element
}

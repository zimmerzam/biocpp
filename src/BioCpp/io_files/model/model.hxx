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

#ifndef BIOCPP_PDB_MODEL_H
#define BIOCPP_PDB_MODEL_H

#include <BioCpp/base_container/base_container.hxx>

namespace BioCpp{
namespace io{

/*! \brief A container with all the atoms of a model in the pdb file

		@tparam int is the serial number of the `atom`
		@tparam atom_info stores the info relative to a specific atom
*/
template <typename atom_t>
class model{
  public:
    typedef base::container<int, atom_t> type;
};

}
}

#endif

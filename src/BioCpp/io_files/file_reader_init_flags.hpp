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

#ifndef BIOCPP_PDB_INIT_FLAGS
#define BIOCPP_PDB_INIT_FLAGS

/*! \brief These flags can be passed to pdb constructor to customize its behaviour while reading a file */
namespace BioCpp{
namespace io{

enum init_flags : unsigned int {
  INIT_COMPLETE    = 0,
  INIT_FAST        = (1<<1),  /*!< Do not read the structure sequence from the first model, but only read SeqRes section */
  INIT_FIRST_MODEL = (1<<2)   /*!< Read only the first model */
};

}
}

#endif

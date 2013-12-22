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

#ifndef TOPOLOGY_CONSTRUCTOR_H
#define TOPOLOGY_CONSTRUCTOR_H

#include <BioCpp/base_container/base_container.hxx>
#include <BioCpp/io_files/model/model.hxx>
#include <BioCpp/io_files/seqres/seqres_record.hpp>

namespace BioCpp{

template< typename atom_t,
          typename vertex_t, typename edge_t
        >
class base_topology_constructor{
  protected:
    typedef topology<vertex_t, edge_t> topology_t;
  public:  
    virtual topology_t operator()( typename io::model<atom_t>::type& info, io::seqres_record& RseqRes, io::seqres_record& TseqRes, residue::dictionary_t& resdict ) = 0;
    virtual topology_t operator()( typename io::model<atom_t>::type& info, io::seqres_record& RseqRes, residue::dictionary_t& resdict ) = 0;
};

}
#endif

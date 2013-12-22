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

#ifndef COMPLEX_CONSTRUCTOR_H
#define COMPLEX_CONSTRUCTOR_H

#include <BioCpp/base_container/base_container.hxx>
#include <BioCpp/io_files/model/model.hxx>
#include <BioCpp/io_files/seqres/seqres_record.hpp>

namespace BioCpp{

template< typename atom_t, 
          typename res_child_id, typename res_id,
          typename cha_child_id, typename cha_id,
          typename cmp_child_id, typename cmp_id
        >
class base_complex_constructor{
  protected:
    typedef base::container<res_child_id, atom_t, res_id> residue_type;
    typedef base::container<cha_child_id, residue_type, cha_id> chain_type;
    typedef base::container<cmp_child_id, chain_type, cmp_id> complex_type;
  public:  

    virtual complex_type operator()( typename io::model<atom_t>::type& info, io::seqres_record& RseqRes, io::seqres_record& TseqRes ) = 0;
    virtual complex_type operator()( typename io::model<atom_t>::type& info, io::seqres_record& RseqRes ) = 0;
};

}
#endif

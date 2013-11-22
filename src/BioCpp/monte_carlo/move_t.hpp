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

#ifndef BIOCPP_MC_MOVE_T_H
#define BIOCPP_MC_MOVE_T_H

#include "move_parameter_t.hpp"
#include <BioCpp/rng/rng.hpp>

namespace BioCpp{
namespace MC{

template <typename Topology_t>
class move_t{
  protected:
    move_parameter_t<Topology_t> parameter;
    rngs& rng;
  public:
    Topology_t& target;

    move_t( move_parameter_t<Topology_t> par, Topology_t& topo, rngs& r ): parameter(par), target(topo), rng(r){};
		
		// Update parameters for the next move
    virtual void update(){};
    // Apply the move to target topology
    virtual void apply(){};
    // Apply the inverse move: undo the last apply()
    virtual void apply_inverse(){};
};

} // end namespace MC
} // end namespace BioCpp
#endif

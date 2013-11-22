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

#ifndef BIOCPP_MC_TRASLATION_T_H
#define BIOCPP_MC_TRASLATION_T_H

#include "move_t.hpp"

namespace BioCpp{
namespace MC{

template <typename Topology_t>
class traslation_t : public move_t<Topology_t>{
  protected:
  	Eigen::Vector3d displacement;
  public:
  	traslation_t( move_parameter_t<Topology_t> par, Topology_t& topo, rngs& r ):move_t<Topology_t>(par,topo,r){};
		// Update parameters for the next move
    void update(){
    	displacement = Eigen::Vector3d(move_t<Topology_t>::rng.RandomGaussianDouble( 0., 1. ), move_t<Topology_t>::rng.RandomGaussianDouble( 0., 1. ),
    	                               move_t<Topology_t>::rng.RandomGaussianDouble( 0., 1. ) ); // generate a vector on a sphere
    	displacement/=displacement.norm();                 // normalize the vector
    	displacement*=move_t<Topology_t>::rng.RandomGaussianDouble( move_t<Topology_t>::parameter.value.first, move_t<Topology_t>::parameter.value.second ); // rescale the vector
    };
    // Apply the move to target topology
    virtual void apply(){
    	typename Topology_t::vertex_iterator v_beg, v_end;
    	boost::tie(v_beg, v_end) = boost::vertices( move_t<Topology_t>::target.getGraph() );
    	for( typename Topology_t::vertex_iterator v = v_beg; v!=v_end; ++v ){
    		move_t<Topology_t>::target.getGraph()[*v].atom->coordinate += displacement;
    	}
    };
    // Apply the inverse move: undo the last apply()
    virtual void apply_inverse(){
    	typename Topology_t::vertex_iterator v_beg, v_end;
    	boost::tie(v_beg, v_end) = boost::vertices( move_t<Topology_t>::target.getGraph() );
    	for( typename Topology_t::vertex_iterator v = v_beg; v!=v_end; ++v ){
    		move_t<Topology_t>::target.getGraph()[*v].atom->coordinate -= displacement;
    	}
    };
};

} //end namespace MC
} //end namespace BioCpp
#endif

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

#ifndef BIOCPP_MC_PIVOT_T_H
#define BIOCPP_MC_PIVOT_T_H

#include "move_t.hpp"

namespace BioCpp{
namespace MC{

template <typename Topology_t>
class pivot_t : public move_t<Topology_t>{
  protected:
  	typename Topology_t::edge_iterator e_beg, e_end;
  	Eigen::Vector3d axis, center;
  	double angle;
  	Topology_t& topology_all_system;
  public:
  	pivot_t( move_parameter_t<Topology_t> par, Topology_t& topo_back, Topology_t topo_all, rngs& r ):move_t<Topology_t>(par,topo_back,r), topology_all_system(topo_all){};
		// Update parameters for the next move
    void update(){
    	boost::tie(e_beg, e_end) = boost::edges( move_t<Topology_t>::target.getGraph() );
    	int num_edges = boost::num_edges( move_t<Topology_t>::target.getGraph() );
			int edge_idx  = move_t<Topology_t>::rng.RandomUniformInteger( 0, num_edges );
			
			typename Topology_t::vertex_t u, v;
			if(edge_idx >= num_edges/2){
				for(int i=0; i!= edge_idx;++i){
					++e_beg;
				}
				u = boost::source(*e_beg, move_t<Topology_t>::target.getGraph());
      	v = boost::target(*e_beg, move_t<Topology_t>::target.getGraph());
      	axis = move_t<Topology_t>::target.getGraph()[v].atom->coordinate - move_t<Topology_t>::target.getGraph()[u].atom->coordinate;
			}
			else{
				e_end = e_beg;
				for(int i=0; i!= edge_idx;++i){
					++e_end;
				}
				u = boost::source(*e_end, move_t<Topology_t>::target.getGraph());
      	v = boost::target(*e_end, move_t<Topology_t>::target.getGraph());
      	axis = move_t<Topology_t>::target.getGraph()[v].atom->coordinate - move_t<Topology_t>::target.getGraph()[u].atom->coordinate;
			}
			center = move_t<Topology_t>::target.getGraph()[v].atom->coordinate;
			angle = move_t<Topology_t>::rng.RandomGaussianDouble( move_t<Topology_t>::parameter.value.first, move_t<Topology_t>::parameter.value.second );
    };
    // Apply the move to target topology
    virtual void apply(){
    };
    // Apply the inverse move: undo the last apply()
    virtual void apply_inverse(){
    };
};

} //end namespace MC
} //end namespace BioCpp
#endif

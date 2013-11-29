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
  	Eigen::Matrix3d rot_m, irot_m;
  	Eigen::Vector3d center;
  	typename std::vector< typename Topology_t::edge_t > edges;
  	typename Topology_t::vertex_t stop;
  	bool first_half;
  public:
  	pivot_t( move_parameter_t<Topology_t> par, Topology_t& topo, rngs& r, typename std::vector< typename Topology_t::edge_t > e_list ):move_t<Topology_t>(par,topo,r), edges(e_list){};
		// Update parameters for the next move
    void update(){
      typename Topology_t::edge_iterator e_beg, e_end;
    	boost::tie(e_beg, e_end) = boost::edges( move_t<Topology_t>::target.getGraph() );
    	unsigned int num_edges = edges.size();
			int edge_idx  = move_t<Topology_t>::rng.RandomUniformInteger( 0, num_edges-1 );
			
			typename Topology_t::vertex_t u = boost::source(edges[edge_idx], move_t<Topology_t>::target.getGraph());
      typename Topology_t::vertex_t v = boost::target(edges[edge_idx], move_t<Topology_t>::target.getGraph());
      Eigen::Vector3d axis = move_t<Topology_t>::target.getGraph()[v].atom->coordinate - move_t<Topology_t>::target.getGraph()[u].atom->coordinate;
      axis/=axis.norm();
			double angle = move_t<Topology_t>::rng.RandomGaussianDouble( move_t<Topology_t>::parameter.value.first, move_t<Topology_t>::parameter.value.second );
			center = move_t<Topology_t>::target.getGraph()[v].atom->coordinate;
			rot_m  = Eigen::AngleAxisd( angle,axis).toRotationMatrix();
			irot_m = Eigen::AngleAxisd(-angle,axis).toRotationMatrix();
			if(edge_idx < num_edges/2){
			  stop = v;
			  first_half = true;
			}
			else{
			  stop = u;
			  first_half = false;
			}
    };
    // Apply the move to target topology
    void apply(){
      typename Topology_t::vertex_iterator v_beg, v_end, v;
    	boost::tie(v_beg, v_end) = boost::vertices( move_t<Topology_t>::target.getGraph() );
    	if(first_half){
    	  for( v=v_beg; v!=v_end; ++v){
    	    if( *v==stop ){
    	      break;
    	    }
    	    move_t<Topology_t>::target.getGraph()[*v].atom->coordinate = rot_m*(move_t<Topology_t>::target.getGraph()[*v].atom->coordinate - center) + center;
    	  }
    	}
    	else{
    	  v=v_end; --v;
    	  for( ;v!=v_beg;--v){
    	    if( *v==stop ){
    	      break;
    	    }
    	    move_t<Topology_t>::target.getGraph()[*v].atom->coordinate = rot_m*(move_t<Topology_t>::target.getGraph()[*v].atom->coordinate - center) + center;
    	  }
    	}
    };
    // Apply the inverse move: undo the last apply()
    void apply_inverse(){
      typename Topology_t::vertex_iterator v_beg, v_end, v;
    	boost::tie(v_beg, v_end) = boost::vertices( move_t<Topology_t>::target.getGraph() );
    	if(first_half){
    	  for( v=v_beg; v!=v_end; ++v){
    	    if( *v==stop ){
    	      break;
    	    }
    	    move_t<Topology_t>::target.getGraph()[*v].atom->coordinate = irot_m*(move_t<Topology_t>::target.getGraph()[*v].atom->coordinate - center) + center;
    	  }
    	}
    	else{
    	  v=v_end; --v;
    	  for( ;v!=v_beg;--v){
    	    if( *v==stop ){
    	      break;
    	    }
    	    move_t<Topology_t>::target.getGraph()[*v].atom->coordinate = irot_m*(move_t<Topology_t>::target.getGraph()[*v].atom->coordinate - center) + center;
    	  }
    	}
    };
};

} //end namespace MC
} //end namespace BioCpp
#endif

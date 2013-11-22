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

#ifndef BIOCPP_MC_METROPOLIS_T_H
#define BIOCPP_MC_METROPOLIS_T_H

#include <BioCpp/rng/rng.hpp>

namespace BioCpp{
namespace MC{

template <typename CollectiveVariable, typename Topology_t>
class metropolis_t{
	public:
  	double inverseT;
	  double E_value_before;
  	rngs& rng;
	  CollectiveVariable& E;
	  metropolis_t( double iT, CollectiveVariable& e, rngs& r ): inverseT(iT), rng(r), E(e) {};
	
  	void init( Topology_t& initial ){
  		E_value_before = E(initial);
  	};
	
		bool check( Topology_t& before, Topology_t& after ){
			double E_bef = E(before);
			double E_aft = E(after );
			return check(E_bef, E_aft);
		}
	
		bool check( Topology_t& after ){
			double E_aft = E(after );
			return check(E_value_before, E_aft);
		}
	
		bool check( double E_bef, double E_aft ){
			if(E_aft < E_bef){
				E_value_before = E_aft;
				return true;
			}
			else if( log( rng.RandomUniformDouble(0.,1.) ) < inverseT*(E_bef-E_aft) ){
				E_value_before = E_aft;
				return true;
			}
			return false;
		}
	};

} // end namespace MC
} // end namespace BioCpp
#endif

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

#ifndef BIOCPP_MC_DYNAMIC_T_H
#define BIOCPP_MC_DYNAMIC_T_H

namespace BioCpp{
namespace MC{

// class for MC-dynamics
template <typename Topology_t, typename Metropolis_t>
class dynamic_t{
  protected:
  	std::vector<double>                move_prob;
	  std::vector< move_t<Topology_t>* > move_list;
	  Metropolis_t& metropolis;
	  rngs& rng;

  public:
	  dynamic_t( std::vector< move_t<Topology_t>* >& list, std::vector<double>& prob, Metropolis_t& m, rngs& r ): metropolis(m), rng(r){
		  if( prob.size()!=list.size() ){
		    std::cout << "move probabilities have to be as many as move types!" << std::endl;
		    exit(1);
		  }
		  unsigned int move_number = prob.size();
		  move_prob.reserve( move_number );
		  move_list.reserve( move_number );
		  for(std::vector<double>::iterator it = prob.begin(); it!=prob.end(); ++it){
			  if( it==prob.begin() ){
				  move_prob.push_back(*it);
			  }
			  else{
				  move_prob.push_back(*it+move_prob.back());
			  }
		  }
		  for(typename std::vector< move_t<Topology_t>* >::iterator it=list.begin(); it!=list.end(); ++it){
			  move_list.push_back( *it );
		  }
	  };

	  move_t<Topology_t>& get_move(double t){
		  unsigned int size = move_list.size();
		  for(unsigned int i=0; i!=size; ++i){
			  if(t<move_prob[i]){
				  return *move_list[i];
			  }
		  }
		  return *move_list.back();
	  }

	  void step(  ){
		  move_t<Topology_t>& move = get_move(rng.RandomUniformDouble(0.,1.));
      move.update();    // Generate new parameters for the move
		  move.apply();     // Apply the move
//		  if( excluded_volume.check() == false ){
//			  move.apply_inverse();
//		  }
		  if( metropolis.check(move.target) == false ){
			  move.apply_inverse();
		  }
		  return;
	  }	
};

} // end namespace MC
} // end namespace BioCpp
#endif

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

#ifndef BIOCPP_MD_DYNAMICS_H
#define BIOCPP_MD_DYNAMICS_H

class dynamics{
  private:
    unsigned int t_start;
    unsigned int t_end;
    double dt;
    
    bool conditional_stop;
    
    unsigned int do_every_n_timesteps;
    unsigned int equilibrate_every_n_timesteps;
  public:
    dynamics(double start, double end, double delta): t_start(start), t_end(end), dt(delta){};

    template< typename Config, typename Potential, typename ImplicitWater, typename Thermostat, typename EvolutionAlgorithm, typename Functor, typename Stop >
    void start( Config& system, Potential& V, ImplicitWater& water, Thermostat& thermostat, EvolutionAlgorithm& evolve, Functor& todo, Stop& stop );
};

template< typename Config, typename Potential, typename ImplicitWater, typename Thermostat, typename EvolutionAlgorithm, typename Functor, typename Stop >
void dynamics::start( Config& system, Potential& V, ImplicitWater& water, Thermostat& thermostat, EvolutionAlgorithm& evolve, Functor& todo, Stop& stop ){
  unsigned int step = 0;
  double t = t_start;
  while( (not conditional_stop and t < t_end) and ( conditional_stop and stop(system) ) ){
    if( step%do_every_n_timesteps == 0 ){
      todo( system );
    }
    if( step%equilibrate_every_n_timesteps == 0 ){
      thermostat( system );
    }
    
    evolve( system, V, water );
    
    t+=dt;    
    ++step;
  }
}

#endif

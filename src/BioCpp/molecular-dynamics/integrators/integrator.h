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

#ifndef BIOCPP_INTEGRATOR_H
#define BIOCPP_INTEGRATOR_H

class integrator{
  private:
    unsigned int t_start;
    unsigned int t_end;
    double dt;
    
    bool conditional_stop;
    
    unsigned int do_every_n_timesteps;
    unsigned int equilibrate_every_n_timesteps;
  public:
    integrator(double start, double end, double delta): t_start(start), t_end(end), dt(delta){};

    template< typename Config, typename Potential, typename Functor >    
    void start( Config& xvaf, Potential& V, Functor& todo );
};

template< typename Config, typename Potential, typename Functor, typename Stop& >
void integrator::start( Config& xvaf, Potential& V, Functor& todo, Stop& stop ){
  unsigned int step = 0;
  double t = t_start;
  while( (not conditional_stop and t < t_end) and ( conditional_stop and stop(xva) ) ){
    if( step%do_every_n_timesteps == 0 ){
      todo( xva, force, energy );
    }
    if( step%equilibrate_every_n_timesteps == 0 ){
      thermostat( xva, force, energy );
    }
    
    config XVAF = V( Projector.inverse() )
    
    
    t+=dt;    
    ++step;
  }
}

#endif

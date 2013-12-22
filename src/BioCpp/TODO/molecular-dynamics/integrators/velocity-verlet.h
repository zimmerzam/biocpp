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

#ifndef BIOCPP_MD_INTEGRATOR_ALGORITHM_VELOCITY_VERLET_H
#define BIOCPP_MD_INTEGRATOR_ALGORITHM_VELOCITY_VERLET_H

#include "integrator_algorithm.h"

class velocity_verlet : integrator_algorithm{
  public:
    template <typename Config, typename Potential, typename ImplicitWater>
    void step( Config& system, Potential& V, ImplicitWater& water ){
      system.gen_coordinate = system.gen_coordinate + system.gen_velocity*dt + 0.5*system.gen_acceleration*dt*dt;
      system.set( system.gen_coordinate );
      system.gen_force  = system.project( -V.gradient( system ) );
      system.gen_force += system.project( water.force( system ) );
      system.gen_velocity = system.gen_velocity + 0.5*( system.gen_acceleration + system.gen_mass.inverse()*system.gen_force )*dt;
      system.gen_acceleration = system.gen_mass.inverse()*system.gen_force;
    }
};

#endif

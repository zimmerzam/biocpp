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

#ifndef BIOCPP_MD_THERMOSTAT_VELOCITY_RESCALING_H
#define BIOCPP_MD_THERMOSTAT_VELOCITY_RESCALING_H

#include "../thermostat.h"

class velocity_rescaling : thermostat{
  
  public:
    template < typename Config >
    void equilibrate ( Config& system ){
      double kinetic_energy = system.gen_velocity.transpose()*system.gen_mass*system.gen_velocity;
      double alpha = sqrt(  system.dof*Kb*temperature/(2.*kinetic_energy) );
      system.gen_velocity = alpha*system.gen_velocity;
    }; 

};

#endif

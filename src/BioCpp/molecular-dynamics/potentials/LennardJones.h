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

#ifndef BIOCPP_MD_POTENTIAL_LJ_H
#define BIOCPP_MD_POTENTIAL_LJ_H

template< typename Parameter, unsigned int exp=6 >
class LennardJones{
  public:
    double value( Parameter& param, double distance ){
      double factor6 = pow( param.sigma/distance, exp );
      return 4.*param.epsilon*( factor6*factor6 - factor6 );
    };
    
    double first_derivative( Parameter& param, double distance ){
      double factor6 = pow( param.sigma/distance, exp );
      return 4.*exp*param.epsilon/distance*( -2.*factor6*factor6 + factor6 ); 
    }
};

#endif

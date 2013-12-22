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

#ifndef BIOCPP_MD_POTENTIAL_H
#define BIOCPP_MD_POTENTIAL_H

template < typename FunctionalForm, typename ParameterTable, typename Compute >
class potential{
  private:
    ParameterTable parameter;
    FunctionalForm equation; 
  public:
    template <typename Config>
    double value( Config& system );
    
    template< typename Config, unsigned int dof >
    Eigen::Matrix<double, dof, 1> gradient( Config& system );
    
};

#endif

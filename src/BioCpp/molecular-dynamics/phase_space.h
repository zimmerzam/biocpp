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

#ifndef BIOCPP_MD_PHASE_SPACE_H
#define BIOCPP_MD_PHASE_SPACE_H

#include <Eigen/Dense>
#include <BioCpp.h>
template< typename atom_t, unsigned int N, unsigned int projected_dof, unsigned int extended_dof = 3*N >
class phase_space{
  protected:
    typedef Eigen::Matrix< double, extended_dof, 1 >              ext_vector_t;
  public:
    typedef Eigen::Matrix< double, projected_dof, 1 >             vector_t;
    typedef Eigen::Matrix< double, projected_dof, projected_dof > matrix_t;

    typename BioCpp::pdb::model<atom_t>::type& system;

    unsigned int degrees_of_freedom;
    unsigned int size;
    
    vector_t gen_coordinate;
    vector_t gen_velocity;
    vector_t gen_acceleration;
    vector_t gen_force;
    matrix_t gen_mass;

    phase_space( typename BioCpp::pdb::model<atom_t>::type& sys ) : system(sys), degrees_of_freedom(projected_dof), size(N) {};
    
    virtual void set( vector_t& coordinate ) = 0;
    
    virtual vector_t project( ext_vector_t& vector ) = 0;
};

#endif

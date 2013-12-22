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

#ifndef MORPHOLOGY_KERNEL_ATOM_H
#define MORPHOLOGY_KERNEL_ATOM_H

#include <Eigen/Core>

namespace BioCpp{
namespace morphology{
namespace cgal_extensible_kernel{

//TODO check it point_t has coordinate
template <typename point_t>
class point : public point_t{
  public:
    const double& x() const { return this->coordinate(0); }
    const double& y() const { return this->coordinate(1); }
    const double& z() const { return this->coordinate(2); }
    double& x() { return this->coordinate(0); }
    double& y() { return this->coordinate(1); }
    double& z() { return this->coordinate(2); }
    
    point(){};
    
    point( double x, double y, double z ){
      this->coordinate = Eigen::Vector3d(x,y,z);
    }
 
};

}
}
}

#endif

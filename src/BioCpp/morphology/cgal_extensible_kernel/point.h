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

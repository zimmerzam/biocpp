#ifndef MORPHOLOGY_KERNEL_BBOX_H
#define MORPHOLOGY_KERNEL_BBOX_H

#include <CGAL/Origin.h>
#include <CGAL/Bbox_3.h>
#include "point.h"

namespace BioCpp{
namespace morphology{
namespace cgal_extensible_kernel{

template <typename ConstructBbox>
class bbox : public ConstructBbox {
public:
  using ConstructBbox::operator();

  template < typename point_t >
  CGAL::Bbox_3 operator()(const point<point_t>& p) const {
    return CGAL::Bbox_3(p.x(), p.y(), p.x(), p.y());
  }
};

}
}
}

#endif 

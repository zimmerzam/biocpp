#ifndef MORPHOLOGY_KERNEL_COORD_ITERATOR_H
#define MORPHOLOGY_KERNEL_COORD_ITERATOR_H

#include "point.h"

namespace BioCpp{
namespace morphology{
namespace kernel{

class coord_iterator {
public:
  template < typename point_t >
  const double* operator()(const point<point_t>& p){
    return &p.x();
  }

  template < typename point_t >
  const double* operator()(const point<point_t>& p, int k){
    if(k==0) return &p.x();
    if(k==1) return &p.y();
    if(k==2) return &p.z();
    const double* pzptr = &p.z();
    pzptr++;
    return pzptr;
  }
};

}
}
}

#endif

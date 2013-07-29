#ifndef MORPHOLOGY_KERNEL_KERNEL_H
#define MORPHOLOGY_KERNEL_KERNEL_H

#include <CGAL/Cartesian.h>
#include "point.h"
#include "point_constructor.h"
//#include "MySegmentC2.h"
#include "bbox.h"
#include "coord_iterator.h"

namespace BioCpp{
namespace morphology{
namespace cgal_extensible_kernel{

template < typename point_t, typename K_, typename K_Base >
class cartesian : public K_Base::template Base<K_>::Type {
  typedef typename K_Base::template Base<K_>::Type   OldK;
  
  public:
    typedef K_                                             Kernel;
    typedef point<point_t>                                 Point_3;
//    typedef MySegmentC2<Kernel>                            Segment_3;
    typedef point_constructor<point_t, Kernel, OldK>       Construct_point_3;
    typedef const double*                                  Cartesian_const_iterator_3;
    typedef coord_iterator                                 Construct_cartesian_const_iterator_3;
    typedef bbox<typename OldK::Construct_bbox_3>          Construct_bbox_3;

  Construct_point_3 construct_point_object() const { 
    return Construct_point_3(); 
  }

  Construct_bbox_3 construct_bbox_object() const { 
    return Construct_bbox_3(); 
  }

  Construct_cartesian_const_iterator_3 construct_cartesian_const_iterator_object() const { 
    return Construct_cartesian_const_iterator_3(); 
  }

  template < typename Kernel3 >
  struct Base { 
    typedef cartesian<point_t, Kernel3, K_Base>  Type; 
  };
  
};


template < typename point_t, typename FT_ >
struct kernel : public CGAL::Type_equality_wrapper< cartesian<point_t, kernel<point_t, FT_>, CGAL::Cartesian<FT_> >, kernel<point_t, FT_> >{

};

}
}
}

#endif

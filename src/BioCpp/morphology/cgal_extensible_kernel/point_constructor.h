#ifndef MORPHOLOGY_KERNEL_POINT_CONSTRUCTOR_H
#define MORPHOLOGY_KERNEL_POINT_CONSTRUCTOR_H

#include "point.h"

namespace BioCpp{
namespace morphology{
namespace cgal_extensible_kernel{

template <typename point_t, typename K, typename OldK>
class point_constructor{
  typedef typename K::RT         RT;
  typedef typename K::Point_3    Point_3;
  typedef typename K::Line_3     Line_3;
  typedef typename Point_3::Rep  Rep;
  public:
    typedef Point_3              result_type;
    
    Rep operator()(CGAL::Return_base_tag, CGAL::Origin o) const { 
      return Rep(o); 
    }

    Rep operator()(CGAL::Return_base_tag, const RT& x, const RT& y, const RT& z) const { 
      return Rep(x, y, z); 
    }

    Rep operator()(CGAL::Return_base_tag, const RT& x, const RT& y, const RT& z, const RT& w) const { 
      return Rep(x, y, z, w); 
    }

    Point_3 operator()(const CGAL::Origin&) const { 
      return point<point_t>(0, 0, 0); 
    }

    Point_3 operator()(const RT& x, const RT& y, const RT& z) const {
      return point<point_t>(x, y, z);
    }

    Point_3 operator()(const Line_3& l) const {
      typename OldK::Construct_point_3 base_operator;
      Point_3 p = base_operator(l);
      return p;
    }

    Point_3 operator()(const Line_3& l, int i) const {
      typename OldK::Construct_point_3 base_operator;
      return base_operator(l, i);
    }

    Point_3 operator()(const RT& x, const RT& y, const RT& z, const RT& w) const {
      if(w != 1){
        return point<point_t>(x/w, y/w, z/w);
      } 
      else {
        return point<point_t>(x,y,z);
      }
    }
};

}
}
}
#endif

#include "../geometry/Eigen/Core"

template<typename T>
class has_coordinate{
    struct Fallback { Eigen::Vector3d coordinate; };
    struct Derived : T, Fallback { };
 
    template<typename U, U> struct Check;
 
    typedef char ArrayOfOne[1];
    typedef char ArrayOfTwo[2];
 
    template<typename U> 
    static ArrayOfOne & func(Check<Eigen::Vector3d Fallback::*, &U::coordinate> *);
 
    template<typename U> 
    static ArrayOfTwo & func(...);
 
  public:
    typedef has_coordinate type;
    enum { value = sizeof(func<Derived>(0)) == 2 };
};

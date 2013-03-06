/*!
    \file surface_area.cpp
    \brief Compute the surface of a set of spheres that is exposed to the solvent
    
*/

#include <iostream>
#include "../src/BioCpp.h"

struct point{
  Eigen::Vector3d coordinate;
  point(double x, double y, double z){
    coordinate = Eigen::Vector3d(x,y,z);
  }
};

int main(){
  point p1(0,0,0);
  point p2(1,0,0);
  point p3(0,1,0);
  
  BioCpp::base_neighborood_map<point*> map;
  
  double radius = 1.;
  std::cout << BioCpp::SurfaceAreaLCPO(p1, map, radius) << std::endl; // output 2.934
  
  map[&p1].insert(&p2);
  map[&p2].insert(&p1);
  std::cout << BioCpp::SurfaceAreaLCPO(p1, map, radius) << std::endl; // output 2.706
  
  map[&p1].insert(&p3);map[&p2].insert(&p3);
  map[&p3].insert(&p1);map[&p3].insert(&p2);
  std::cout << BioCpp::SurfaceAreaLCPO(p1, map, radius) << std::endl; // output 2.478
  
  BioCpp::base_container<int, point, int> line;
  line.Append(1,p1);
  line.Append(2,p2);
  line.Append(3,p3);
  std::cout << BioCpp::SurfaceAreaLCPO(line, map, radius) << std::endl; // output 8.802
}

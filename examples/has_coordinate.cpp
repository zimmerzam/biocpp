#include <iostream>
#include <BioCpp/utils/has_coordinate.h>

struct yes{
  Eigen::Vector3d coordinate;
};

struct no{
  Eigen::Vector3d coordinates;
};

int main(){
  std::cout << has_coordinate<yes>::value << std::endl;
  std::cout << has_coordinate<no>::value << std::endl;
  
  return 0;
}

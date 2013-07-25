#include <BioCpp/morphology/kernel/kernel.h>

#include <CGAL/basic.h>
#include <CGAL/Filtered_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Alpha_shape_3.h>
#include <cassert>
#include <list>
#include <vector>

#include <BioCpp.h>

typedef BioCpp::morphology::kernel::kernel<BioCpp::pdb::atom_info, double>  MK;
typedef CGAL::Filtered_kernel_adaptor<MK>           K;
typedef CGAL::Alpha_shape_vertex_base_3<K>          Vb;
typedef CGAL::Alpha_shape_cell_base_3<K>            Fb;
typedef CGAL::Triangulation_data_structure_3<Vb,Fb> Tds;
typedef CGAL::Delaunay_triangulation_3<K, Tds>      Triangulation_3;
typedef CGAL::Alpha_shape_3<Triangulation_3>        Alpha_shape_3;

typedef K::Point_3  Point;

int main(){
  std::list<Point> points;
  points.push_back( Point(-1,-1,-1) );
  points.push_back( Point(0,0,1) );
  points.push_back( Point(0,1,0) );
  points.push_back( Point(1,0,0) );
  points.push_back( Point(0.01,0.01,0.01) );
  
  Alpha_shape_3 as( points.begin(), points.end() );
  
  as.set_alpha(0.5);

  std::vector<Alpha_shape_3::Vertex_handle> vertices;
  as.get_alpha_shape_vertices( std::back_inserter(vertices), Alpha_shape_3::REGULAR );

  std::size_t nbf = vertices.size();
  for(std::size_t i = 0; i < nbf; ++i){
    std::cout << "########## " << i << " ########" << std::endl;
    std::cout << vertices[i]->point() << std::endl;
  }
  
  
  return 0;
}

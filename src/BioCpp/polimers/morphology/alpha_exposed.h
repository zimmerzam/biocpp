#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Alpha_shape_3.h>

#include <fstream>
#include <list>
#include <cassert>
#include <BioCpp/geometry/Eigen/Core>

namespace BioCpp{

void AlphaExposed( std::list<Eigen::Vector3d>& points, std::list<Eigen::Vector3d>& surface, double alpha ){
  typedef CGAL::Exact_predicates_inexact_constructions_kernel Gt;

  typedef CGAL::Alpha_shape_vertex_base_3<Gt>          Vb;
  typedef CGAL::Alpha_shape_cell_base_3<Gt>            Fb;
  typedef CGAL::Triangulation_data_structure_3<Vb,Fb>  Tds;
  typedef CGAL::Delaunay_triangulation_3<Gt,Tds>       Triangulation_3;
  typedef CGAL::Alpha_shape_3<Triangulation_3>         Alpha_shape_3;

  typedef Gt::Point_3                                  Point;
  typedef Alpha_shape_3::Alpha_iterator               Alpha_iterator;
  
  std::list<Point> lp;  
  Point p;
  for( typename std::list< Eigen::Vector3d >::iterator pt = points.begin(); pt != points.end(); ++pt ){
    p = Point( (*pt)(0), (*pt)(1), (*pt)(2) );
    lp.push_back(p);
  }
  Alpha_shape_3 as(lp.begin(),lp.end());
  
  if( alpha < 0 ){
    Alpha_iterator opt = as.find_optimal_alpha(1);
    as.set_alpha(*opt);
  }
  else{
    as.set_alpha(alpha);
  }

  std::vector<Alpha_shape_3::Facet> facets;
  as.get_alpha_shape_facets( std::back_inserter(facets), Alpha_shape_3::REGULAR );

  std::size_t nbf = facets.size();
  Eigen::Vector3d pt;
  for(std::size_t i = 0; i < nbf; ++i){
    for( int j = 0; j < 3; ++j ){
      pt(0) = facets[i].first->vertex(j)->point()[0];
      pt(1) = facets[i].first->vertex(j)->point()[1];
      pt(2) = facets[i].first->vertex(j)->point()[2];
      surface.push_back( pt );
    }
  }
}

} //end namespace

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Alpha_shape_3.h>

#include <boost/iterator/zip_iterator.hpp>
#include <vector>

class BioCppVertex : public CGAL::Alpha_shape_vertex_base_3<CGAL::Exact_predicates_inexact_constructions_kernel>{
  public:
    int ciao;
    BioCppVertex(){};
};

typedef CGAL::Exact_predicates_inexact_constructions_kernel         Gt;
typedef BioCppVertex                                                Vb;
typedef CGAL::Alpha_shape_cell_base_3<Gt>                           Fb;
typedef CGAL::Triangulation_data_structure_3<Vb,Fb>                 Tds;
typedef CGAL::Delaunay_triangulation_3<Gt, Tds>                     Delaunay;
typedef Gt::Point_3                                                 Point;

typedef CGAL::Alpha_shape_3<Delaunay>                               Alpha_shape_3;
typedef Alpha_shape_3::Alpha_iterator                               Alpha_iterator;

int main()
{

  std::vector<unsigned> indices;
  indices.push_back(0);
  indices.push_back(1);
  indices.push_back(2);
  indices.push_back(3);
  indices.push_back(4);
  indices.push_back(5);  
  
  std::vector<Point> points;
  points.push_back(Point(0,0,0));
  points.push_back(Point(1,0,0));
  points.push_back(Point(0,1,0));
  points.push_back(Point(0,0,1));
  points.push_back(Point(2,2,2));
  points.push_back(Point(-1,0,1));

  
  
  Delaunay T( boost::make_zip_iterator(boost::make_tuple( points.begin(),indices.begin() )),
              boost::make_zip_iterator(boost::make_tuple( points.end(),indices.end() ) )  );

//  Alpha_shape_3 as(boost::make_zip_iterator(boost::make_tuple( points.begin(),indices.begin() )),
//              boost::make_zip_iterator(boost::make_tuple( points.end(),indices.end() ) ));

  CGAL_assertion( T.number_of_vertices() == 6 );
  
  
  // check that the info was correctly set.
  Delaunay::Finite_vertices_iterator vit;
  for (vit = T.finite_vertices_begin(); vit != T.finite_vertices_end(); ++vit)
    std::cout << vit->ciao << std::endl;  

  return 0;
}

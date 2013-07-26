#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Alpha_shape_3.h>
#include <Eigen/Core>

#include <fstream>
#include <list>
#include <cassert>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Gt;

typedef CGAL::Alpha_shape_vertex_base_3<Gt>          Vb;
typedef CGAL::Alpha_shape_cell_base_3<Gt>            Fb;
typedef CGAL::Triangulation_data_structure_3<Vb,Fb>  Tds;
typedef CGAL::Delaunay_triangulation_3<Gt,Tds>       Triangulation_3;
typedef CGAL::Alpha_shape_3<Triangulation_3>         Alpha_shape_3;

typedef Gt::Point_3                                  Point;
typedef Alpha_shape_3::Alpha_iterator               Alpha_iterator;

int main(int argc, char* argv[])
{
  const char* filename = argc==2 ? argv[1] : "./data/bunny_5000";
  std::list<Point> lp;
  
  //read input
  std::ifstream is(filename);
  int n;
  is >> n;
  std::cout << "Reading " << n << " points " << std::endl;
  Point p;
  for( ; n>0 ; n--)    {
    is >> p;
    lp.push_back(p);
  }

  // compute alpha shape
  Alpha_shape_3 as(lp.begin(),lp.end());
  std::cout << "Alpha shape computed in REGULARIZED mode by default"
	    << std::endl;

  // find optimal alpha value
  Alpha_iterator opt = as.find_optimal_alpha(1);
  std::cout << "Optimal alpha value to get one connected component is "
	    <<  *opt    << std::endl;
  as.set_alpha(*opt);
//  as.set_alpha(4);
  assert(as.number_of_solid_components() == 1);

  std::vector<Alpha_shape_3::Facet> facets;
  as.get_alpha_shape_facets( std::back_inserter(facets), Alpha_shape_3::REGULAR );

  std::size_t nbf = facets.size();
  for(std::size_t i = 0; i < nbf; ++i){
    std::cout << facets[i].first->vertex(0)->point() << std::endl
              << facets[i].first->vertex(1)->point() << std::endl
              << facets[i].first->vertex(2)->point() << std::endl;
  }
  return 0;
}

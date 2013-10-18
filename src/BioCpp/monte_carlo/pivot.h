/* ************************************************************************** */
/*                                                                            */
/*    Copyright 2013 Stefano Zamuner                                          */
/*                                                                            */
/*    This file is part of BioCpp.                                            */
/*                                                                            */
/*    BioCpp is free software: you can redistribute it and/or modify          */
/*    it under the terms of the GNU General Public License as published by    */
/*    the Free Software Foundation, either version 3 of the License, or       */
/*    (at your option) any later version.                                     */
/*                                                                            */
/*    BioCpp is distributed in the hope that it will be useful,               */
/*    but WITHOUT ANY WARRANTY; without even the implied warranty of          */
/*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           */
/*    GNU General Public License for more details.                            */
/*                                                                            */
/*    You should have received a copy of the GNU General Public License       */
/*    along with BioCpp.  If not, see <http://www.gnu.org/licenses/>.         */
/*                                                                            */
/* ************************************************************************** */

#ifndef PIVOT_H
#define PIVOT_H

#include <boost/graph/breadth_first_search.hpp>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <BioCpp/topology/topology.h>

template < typename atom_prop, typename bond_prop >
class pivot{
  protected:
  
    class bfs_rotate_visitor : public boost::default_bfs_visitor {
      protected:
        Eigen::Vector3d base;
        Eigen::Matrix3d rotation;
        
      public:
        bfs_rotate_visitor( Eigen::Vector3d& source, Eigen::Vector3d& target, double rotangle ): base(target){
          Eigen::Vector3d axis = target-source;
          axis = axis/axis.norm();
          Eigen::AngleAxisd rot(rotangle, axis);
          rotation = rot.toRotationMatrix();
        }
        
        template < typename Vertex, typename Graph >
        void discover_vertex(Vertex u, const Graph & G) const {
          if( base == G[u].atom->coordinate ){
            return;
          }
          G[u].atom->coordinate = rotation*(G[u].atom->coordinate - base) + base;
//          std::cout << G[u].atom->chainId << G[u].atom->resSeq << "  " << G[u].atom->id << "  "
//                    << G[u].atom->coordinate.transpose() << std::endl;
        }
    };
  
  public:
    void operator()( BioCpp::topology<atom_prop,bond_prop>& topo, typename BioCpp::topology<atom_prop,bond_prop>::edge_iterator& ei, double angle ){
      topo.getGraph()[(*ei)].edge.set();
      typename BioCpp::topology<atom_prop,bond_prop>::vertex_t u = boost::source(*ei, topo.getGraph());
      typename BioCpp::topology<atom_prop,bond_prop>::vertex_t v = boost::target(*ei, topo.getGraph());
      bfs_rotate_visitor rotate( topo.getGraph()[u].atom->coordinate, topo.getGraph()[v].atom->coordinate, angle );
      boost::breadth_first_search( topo.getGraph(), v, boost::visitor(rotate) );
    }

};


#endif

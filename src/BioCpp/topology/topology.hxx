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

#ifndef BIOCPP_TOPOLOGY_H
#define BIOCPP_TOPOLOGY_H

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>

namespace BioCpp{

template < typename atom_prop, typename bond_prop >
class topology{
  protected:
    Graph G;  
  public:
    typedef typename boost::adjacency_list<boost::listS, boost::vecS, boost::directedS, atom_prop, bond_prop > Graph;
    typedef typename boost::graph_traits<Graph>::vertex_descriptor vertex_t;
    typedef typename boost::graph_traits<Graph>::edge_descriptor edge_t;
    typedef typename boost::graph_traits<Graph>::vertex_iterator vertex_iterator;
    typedef typename boost::graph_traits<Graph>::edge_iterator edge_iterator;
  
    Graph& getGraph();
    
    template< typename atom_t, typename vertex_t, typename edge_t>
    friend class topology_constructor
};

template < typename atom_prop, typename bond_prop >
typename topology<atom_prop,bond_prop>::Graph& topology<atom_prop,bond_prop>::getGraph(){
  return G;
}

}

#endif

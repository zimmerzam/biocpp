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

#ifndef BIOCPP_PDB_STANDARD_RESIDUE_DICTIONARY_H
#define BIOCPP_PDB_STANDARD_RESIDUE_DICTIONARY_H

namespace BioCpp{
namespace standard{

template <typename atom_type>
struct vertex{
  typedef atom_type atom_t;
  atom_t* atom;
  void set( atom_t& atm ){
    atom = &atm;
  };
};

template <typename atom_t, typename bond_t>
struct edge{
  bond_t<atom_t> bond;
};

template <typename vertex_t, typename edge_t>
class topology_constructor : public BioCpp::topology_constructor<vertex_t::atom_t, vertex_t, edge_t>{
  
  BioCpp::topology_constructor<vertex_t::atom_t, vertex_t, edge_t>::topology_t 
  operator()( typename io::model<atom_t>::type& info, io::seqres_record& RseqRes, io::seqres_record& TseqRes, residue::dictionary_t& resdict ){
    topology<vertex_t, edge_t> topo;
    std::map< std::pair<int,int>, BioCpp::topology< vertex, edge >::vertex_t > added_vertex; // TODO add chain
    for(io::seqres_record::iterator seq = TseqRes.begin(); seq != TseqRes.end(); ++seq){
      for(std::string::iterator res = seq->second.begin(); res != seq->second.end(); ++res){
        int resid = resdict.string_to_id[*res];
        for(residue::dictionary_t::atom_list_t::iterator at = resdict.definition[resid].atomlist.begin(); at != resdict.definition[resid].atomlist.end(); ++at){
          BioCpp::topology< vertex, edge >::vertex_t u = boost::add_vertex(topo.G); // TODO add to added_vertex
        }
      }
    }
    
  }
  BioCpp::topology_constructor<vertex_t::atom_t, vertex_t, edge_t>::topology_t 
    operator()( typename io::model<atom_t>::type& info, io::seqres_record& RseqRes, residue::dictionary_t& resdict );
}

}
}

#endif

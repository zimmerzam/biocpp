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

#include <cctype>
#include <BioCpp/topology/topology_constructor.hxx>

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

template < typename bond_t >
struct edge{
  bond_t bond;
};

template <typename vertex_t, typename edge_t>
class topology_constructor : public BioCpp::base_topology_constructor<typename vertex_t::atom_t, vertex_t, edge_t >{
  
  public:
	  typename BioCpp::base_topology_constructor<typename vertex_t::atom_t, vertex_t, edge_t>::topology_t 
  	operator()( typename io::model<typename vertex_t::atom_t>::type& info, io::seqres_record& RseqRes, io::seqres_record& TseqRes, residue::dictionary_t& resdict ){
			typedef typename vertex_t::atom_t atom_t;
		
			std::map< std::pair<char,unsigned int>, std::set<int> > atom_list;
			std::map< std::pair< std::pair<char,unsigned int>, int>,  typename vertex_t::atom_t* > atom_ref;
			typename io::model<atom_t>::type::iterator atm = info.begin();
			
			char chainId = atm->chainId;
			unsigned int resSeq = 0;
			for(std::string::iterator res = RseqRes[chainId].begin(); res != RseqRes[chainId].end() and atm!=info.end(); ++resSeq){
				std::pair<char,unsigned int> resPos = std::make_pair( chainId,resSeq );
				atom_t base = *atm;
				if(*res == '-'){
					atom_list[ resPos ] = std::set<int>({});
					++res;
				}
				else if( not isalpha(*res) ){
					++res;
				}
				else{
					for(;atm!=info.end();++atm){
						if( atm->chainId!=base.chainId ){
							chainId = atm->chainId;
							res = RseqRes[chainId].begin();
							resSeq = -1;
							break;
						}
						else if( atm->resSeq!=base.resSeq or  atm->iCode!=base.iCode ){
							++res;
							break;
						}
						else{
							atom_list[ resPos ].insert( atm->id );
							atom_ref[ std::make_pair( resPos, atm->id ) ] = &(*atm);
						}
					}
				}
			}
			
			typename BioCpp::base_topology_constructor<typename vertex_t::atom_t, vertex_t, edge_t>::topology_t topo;
			for( io::seqres_record::iterator ch =  TseqRes.begin(); ch!=TseqRes.end(); ++ch){
				resSeq = 0;
				std::map< std::pair<unsigned int, int>, BioCpp::topology< vertex_t, edge_t >::vertex_t  > added_vertex;
				for(std::string::iterator res = ch->second.begin(); res != ch->second.end() ; ++res,++resSeq){
					residue::dictionary_t::definition_t::model_t res_model = ;// TODO  confronta la lista di atomi letti con i possibili modelli e restituisci il modello migliore
					std::pair<char,unsigned int> resPos = std::make_pair(*ch,resSeq);
					for( residue::dictionary_t::definition_t::model_t::atom_list_t at = res_model.atom_list.begin();  at != res_model.atom_list.end(); ++at){
						BioCpp::topology< vertex_t, edge_t >::vertex_t u = boost::add_vertex(topo.G);
          	added_vertex[ std::make_pair(resSeq,*at) ] = u;
          	std::pair< std::pair<char,unsigned int>, int > atmId = std::make_pair(resPos,*at);
          	std::map< std::pair< std::pair<char,unsigned int>, int>,  typename vertex_t::atom_t* >::iterator atm_ref = atom_ref.find( atmId );
          	if( atm_ref!=atom_ref.end() ){
	          	topo.getGraph()[u].set( *(atm_ref->second) );
	          }
	          else{
	          	// TODO aggiungi un atomo ad info...
	          }
					}
					for( residue::dictionary_t::definition_t::model_t::bond_list_t bn = res_model.bond_list.begin();  bn != res_model.bond_list.end(); ++bn){
						BioCpp::topology< vertex, edge >::edge_t e;
          	BioCpp::topology< vertex, edge >::vertex_t u = added_vertex[ std::make_pair(resSeq, bn->first) ];
          	BioCpp::topology< vertex, edge >::vertex_t v = added_vertex[ std::make_pair(resSeq, bn->second) ];
          	boost::tie(e,b) = boost::add_edge(u,v,G);
					}
					// TODO collega gli ultimi due residui...
				}
			
			return topo;
		}
		

  	typename BioCpp::base_topology_constructor<typename vertex_t::atom_t, vertex_t, edge_t>::topology_t 
    operator()( typename io::model<typename vertex_t::atom_t>::type& info, io::seqres_record& RseqRes, residue::dictionary_t& resdict ){
    	return (*this)( info, RseqRes, RseqRes, resdict );
    };
};

}
}

#endif

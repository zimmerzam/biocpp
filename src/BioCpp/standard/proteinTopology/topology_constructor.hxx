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
#include <float.h>
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
  private:
    residue::dictionary_t::definition_t::model_t getBestModel(residue::dictionary_t& resdict, int id, std::set<int> atom_list){ //TODO this is a stupid implementation...
      return *resdict.definition[id].model.begin();
    }
  public:
	  typename BioCpp::base_topology_constructor<typename vertex_t::atom_t, vertex_t, edge_t>::topology_t 
  	operator()( typename io::model<typename vertex_t::atom_t>::type& info, io::seqres_record& RseqRes, io::seqres_record& TseqRes, atom::dictionary_t& atmdict, residue::dictionary_t& resdict ){
			typedef typename vertex_t::atom_t atom_t;
		
			std::map< std::pair<char,unsigned int>, std::set<int> > atom_list;
			std::map< std::pair< std::pair<char,unsigned int>, int>,  typename vertex_t::atom_t* > atom_ref;
			typename io::model<atom_t>::type::iterator atm = info.begin();
			
			// create a list of all atoms contained in the structure
			char chainId = atm->chainId;
			unsigned int resSeq = 0;
			std::map<char,int> start_resSeq = { {chainId,atm->resSeq} };
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
							start_resSeq[chainId] = atm->resSeq;
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
			
//			for(std::map< std::pair<char,unsigned int>, std::set<int> >::iterator it = atom_list.begin(); it!=atom_list.end(); ++it ){
//			  std::cout << it->first.first << "  " << it->first.second << "  " <<resdict.id_to_string[ resdict.string_to_id[ std::string(1,RseqRes[it->first.first][it->first.second]) ] ] << "  [";
//			  for( std::set<int>::iterator at = it->second.begin(); at != it->second.end(); ++at ){
//			    std::cout << *at << ",";
//			  }
//			  std::cout <<"]"<< std::endl;
//			}
			
			// Build the topology using information from TseqRes.
			int new_serial = info.rbegin()->serial;
			typename BioCpp::base_topology_constructor<atom_t, vertex_t, edge_t>::topology_t topo;
			for( io::seqres_record::iterator ch =  TseqRes.begin(); ch!=TseqRes.end(); ++ch){
				unsigned int last_resSeq = -10;
				int last_end = 3;
				int atm_resSeq = start_resSeq[ch->first];
				resSeq = 0;
				std::map< std::pair<unsigned int, int>, typename BioCpp::topology< vertex_t, edge_t >::vertex_t  > added_vertex;
				for(std::string::iterator res = ch->second.begin(); res != ch->second.end() ; ++res,++resSeq,++atm_resSeq){
				  if(*res == '-' or *res == 'X'){
				    continue;
				  }
					residue::dictionary_t::definition_t::model_t res_model = getBestModel(resdict,resdict.string_to_id[std::string(1,*res)],atom_list[std::make_pair(ch->first,resSeq)]);
					std::pair<char,unsigned int> resPos = std::make_pair(ch->first,resSeq);
					for( residue::dictionary_t::definition_t::model_t::atom_list_t::iterator at = res_model.atom_list.begin();  at != res_model.atom_list.end(); ++at){
						typename BioCpp::topology< vertex_t, edge_t >::vertex_t u = boost::add_vertex(topo.getGraph());
          	added_vertex[ std::make_pair(resSeq,*at) ] = u;
          	std::pair< std::pair<char,unsigned int>, int > atmId = std::make_pair(resPos,*at);
          	typename std::map< std::pair< std::pair<char,unsigned int>, int>,  atom_t* >::iterator atm_ref = atom_ref.find( atmId );
            // if the atom exists, put a reference to it in the graph
          	if( atm_ref!=atom_ref.end() ){
	          	topo.getGraph()[u].set( *(atm_ref->second) );
	          }
	          // if the required atom does not exists, creates it, appends it to 'info' and creates a reference to it in the graph; 
	          else{
	          	atom_t tmp_atom;
	          	tmp_atom.serial = ++new_serial;
	          	tmp_atom.id = *at;
	          	tmp_atom.altLoc = ' ';
	          	tmp_atom.resName = resdict.string_to_id[std::string(1,*res)];
	          	tmp_atom.chainId = ch->first;
	          	tmp_atom.resSeq = atm_resSeq;
	          	tmp_atom.iCode = ' ';
	          	tmp_atom.coordinate = Eigen::Vector3d(999.,999.,999.);
	          	tmp_atom.occupancy = 0.0;
	          	tmp_atom.tempFactor = 0.0;
	          	tmp_atom.element = atmdict.definition[*at].element;
	          	tmp_atom.charge = 0.0;
	          	info.Append(tmp_atom.serial, tmp_atom);
	          	topo.getGraph()[u].set( info[tmp_atom.serial] );
	          }
					}
					// connects atoms in the same residue
					bool b = false;
					for( residue::dictionary_t::definition_t::model_t::bond_list_t::iterator bn = res_model.bond_list.begin();  bn != res_model.bond_list.end(); ++bn){
						typename BioCpp::topology< vertex_t, edge_t >::edge_t e;
          	typename BioCpp::topology< vertex_t, edge_t >::vertex_t u = added_vertex[ std::make_pair(resSeq, bn->first) ];
          	typename BioCpp::topology< vertex_t, edge_t >::vertex_t v = added_vertex[ std::make_pair(resSeq, bn->second) ];
          	boost::tie(e,b) = boost::add_edge(u,v,topo.getGraph());
					}
					// connects two consecutive residues
					if( last_resSeq == (resSeq-1) ){
  					typename BioCpp::topology< vertex_t, edge_t >::edge_t e;
            typename BioCpp::topology< vertex_t, edge_t >::vertex_t u = added_vertex[ std::make_pair(last_resSeq, last_end) ];
            typename BioCpp::topology< vertex_t, edge_t >::vertex_t v = added_vertex[ std::make_pair(resSeq, res_model.atom_start) ];
            boost::tie(e,b) = boost::add_edge(u,v,topo.getGraph());
          }
          last_resSeq = resSeq;
          last_end = res_model.atom_end;
				}
			}
			return topo;
		};

  	typename BioCpp::base_topology_constructor<typename vertex_t::atom_t, vertex_t, edge_t>::topology_t 
    operator()( typename io::model<typename vertex_t::atom_t>::type& info, io::seqres_record& RseqRes, atom::dictionary_t& atmdict, residue::dictionary_t& resdict ){
    	return (*this)( info, RseqRes, RseqRes, atmdict, resdict );
    };
};

}
}

#endif

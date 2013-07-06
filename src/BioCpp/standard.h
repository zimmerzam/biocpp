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

#ifndef BIOCPP_STANDARD_H
#define BIOCPP_STANDARD_H

#include <list>
#include "polimers/moiety_id.h"
#include "polimers/dpss/dpss_id.h"
#include "polimers/dpss/h_bridge_energy.h"
#include "pdb/pdb.h"
namespace BioCpp{
namespace standard{

/*! \brief A standard residue 

    Each atom is identified by an atom::id
*/
typedef BioCpp::base_container<BioCpp::atom::id, BioCpp::pdb_atom_info, BioCpp::amino_acid::id> residue;
/*! \brief A standard chain 

    Each residue is identified by an int (usually its seqRes)
*/
typedef BioCpp::base_container<int, residue, char> chain;
/*! \brief A standard complex

    Each chain is identified by a char (usually its chainId)
*/
typedef BioCpp::base_container<char, chain, std::string> complex;
/*! \brief A standard H_bridge_map */
typedef BioCpp::base_h_bridge_map<chain::iterator> h_bridge_map;

typedef std::pair< BioCpp::moiety::id, std::list<BioCpp::atom::id> > moiety_info;

}
}

/*!	\brief Simple constructor from a pdb_model
    
    \note You may want to overwrite this!
    \note This does not take the TseqRes into account
    
    \todo check!
    \todo use TseqRes
    \todo enhance with a more complete version
*/
namespace BioCpp{
template<>
template<>
base_container<char, standard::chain, std::string>::base_container(pdb_model& atom_list, pdb_seqres_record& TseqRes){
  this->Reserve(TseqRes.size());
  char prev_chain = ' ';
  int prev_residue = -9999;
  for(BioCpp::pdb_model::iterator it=atom_list.begin(); it!=atom_list.end(); ++it){
  	
  	if( it->chainId != prev_chain ){
  		standard::chain tmp_chain;
	    Append(it->chainId,tmp_chain);
	    (*this)[it->chainId].Reserve( TseqRes[it->chainId].size() );
	    (*this)[it->chainId].type = it->chainId;
	    prev_chain = it->chainId;
	  }
	  if( it->resSeq != prev_residue ){
	  	standard::residue tmp_res;
	  	tmp_res.type = it->resName;
	    (*this)[it->chainId].Append(it->resSeq,tmp_res);
  	  (*this)[it->chainId][it->resSeq].Reserve(30);
  	  prev_residue = it->resSeq;
  	}
    (*this)[prev_chain][prev_residue].Append(it->id, *it);
  }
}

/*!	\brief Constructor from a reference to a protein complex

    This function return a matrix whose entry (i,j) are True if there is an H bridge
    between i and j  
*/
template<>
template<>
base_h_bridge_map<standard::chain::iterator>::base_h_bridge_map(standard::complex& cmp){
	for(standard::complex::iterator ch1 = cmp.begin(); ch1 != cmp.end(); ++ch1){
		for(standard::chain::iterator res1 = ch1->begin(); res1 != ch1->end(); ++res1){
			Eigen::Vector3d c, o;
			if( not res1->exists(atom::C_) ){
//				std::cout << "missing backbone carbon" << std::endl;
				continue;
			}
			if(not res1->exists(atom::O_) ){
//				std::cout << "missing backbone oxygen" << std::endl;
				continue;
			}
			c = (*res1)[atom::C_].coordinate;
			o = (*res1)[atom::O_].coordinate;
			for(standard::complex::iterator ch2 = cmp.begin(); ch2 != cmp.end(); ++ch2){
				for(standard::chain::iterator res2 = ch2->begin(); res2 != ch2->end(); ++res2){
					if( res2!=ch2->begin() ){
						Eigen::Vector3d n, h;
						if(not res2->exists(atom::N_) ){
//							std::cout << "missing backbone nitrogen" << std::endl;
							continue;
						}
						if(not res2->exists(atom::H_) ){
//							std::cout << "missing backbone hydrogen" << std::endl;
							continue;
						}
						n = (*res2)[atom::N_].coordinate;
						h = (*res2)[atom::H_].coordinate;

						double en = dpss::h_bridge_energy(c, o, n, h);
						if(en<-1){
							(*this)[std::make_pair(res1, res2)] = true;
						}
					}
				}
			}
		}
	}
}

}
#endif

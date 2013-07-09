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

}
}

/*!	\brief Simple constructor from a pdb_model
    
    \note You may want to overwrite this!
*/
namespace BioCpp{
template<>
template<>
inline base_container<char, standard::chain, std::string>::base_container(pdb_model& atom_list, pdb_seqres_record& RseqRes, pdb_seqres_record& TseqRes){
  this->Reserve(TseqRes.size());
  BioCpp::pdb_seqres_record rseqres, tseqres;
  //replace gap with three virtual residues
  for( pdb_seqres_record::iterator ch = TseqRes.begin(); ch!=TseqRes.end(); ++ch ){
    tseqres[ch->first] = "";
    for(std::string::iterator s = ch->second.begin(); s != ch->second.end(); ++s){
      if((*s)!='-')
        tseqres[ch->first]+=(*s);
      else
        tseqres[ch->first]+="^^^";
    }
  }
  for( pdb_seqres_record::iterator ch = RseqRes.begin(); ch!=RseqRes.end(); ++ch ){
    rseqres[ch->first] = "";
    for(std::string::iterator s = ch->second.begin(); s != ch->second.end(); ++s){
      if((*s)!='-')
        rseqres[ch->first]+=(*s);
      else
        rseqres[ch->first]+="^^^";
    }
  }
  //for each chain...
	for( pdb_seqres_record::iterator ch = tseqres.begin(); ch!=tseqres.end(); ++ch ){
		standard::chain tmp_chain;
		tmp_chain.type = ch->first;
	  Append(ch->first, tmp_chain);
	  (*this)[ch->first].Reserve( ch->second.size() );
	  
	  BioCpp::pdb_model::iterator it=atom_list.begin(); //first atom in list
	  while( it->chainId != ch->first and it!=atom_list.end() ){ //skip atoms belonging to unwanted chains
	  	++it;
	  }
	  //now 'it' is the first atom of current chain
	  unsigned int res = 0;
	  int first_res = it->resSeq;
	  if( RseqRes[ch->first][0] == '^' ){ //count initial missing residue (in order to adjust first_res)
	  	while( rseqres[ch->first][res] == '^' ){
	  		--first_res;
	  		++res;
	  	}
	  }
	  res = 0;
		while( res < (ch->second.size()) ){
	  	if( rseqres[ch->first][res] == '^' ){
	  		standard::residue tmp_res;
        if( ch->second[res]=='^' ){
          tmp_res.type = BioCpp::amino_acid::UNK;
        }
        else{
  	  		std::stringstream sstype;
	    		sstype << ch->second[res];
	  	  	tmp_res.type = BioCpp::amino_acid::string_to_id[ sstype.str() ];
	  	  }
  	    (*this)[ch->first].Append(first_res+res,tmp_res);
		    ++res;
	  	}
	  	else{
	  		standard::residue tmp_res;
		  	std::stringstream sstype;
	  		sstype << ch->second[res];
		  	tmp_res.type = BioCpp::amino_acid::string_to_id[ sstype.str() ];
		    (*this)[ch->first].Append(first_res+res,tmp_res);
		    (*this)[ch->first][first_res+res].Reserve(40);
		    int prev_res = it->resSeq;
		    char prev_icode = it->iCode;
		    char prev_altloc = it->altLoc;
		  	while(it->resSeq == prev_res and it->iCode==prev_icode ){
		  	  if(it->altLoc!=prev_altloc){
		  	    ++it;
		  	    continue;
		  	  }
		  		(*this)[ch->first][first_res+res].Append(it->id, *it);
   				if( it != (atom_list.end()-1) )
  			  	++it;
  			  else
  			    break;
			  }
		    ++res;
	  	}
	  }
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

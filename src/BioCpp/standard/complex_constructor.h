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

#ifndef STANDARD_COMPLEX_CONSTRUCTOR_H
#define STANDARD_COMPLEX_CONSTRUCTOR_H

#include <BioCpp/polimers/proteins/base_complex_constructor.h>
#include <string>

#include "complex.h"

namespace BioCpp{

namespace standard{

template< typename atom_t >
class complex_constructor : public BioCpp::base_complex_constructor< atom_t, 
                                              typename residue<atom_t>::type::child_id_type, 
                                              typename residue<atom_t>::type::id_type, 
                                              typename chain<atom_t>::type::child_id_type, 
                                              typename chain<atom_t>::type::id_type,
                                              typename complex<atom_t>::type::child_id_type, 
                                              typename complex<atom_t>::type::id_type 
                                             >{
                                             
  public:
    typename standard::complex<atom_t>::type operator()( typename pdb::model<atom_t>::type& info, pdb::seqres_record& RseqRes, pdb::seqres_record& TseqRes ){
      typename standard::complex<atom_t>::type cmp;
    
      cmp.Reserve(TseqRes.size());
      BioCpp::pdb::seqres_record rseqres, tseqres;
      //replace gap with three virtual residues
      for( pdb::seqres_record::iterator ch = TseqRes.begin(); ch!=TseqRes.end(); ++ch ){
        tseqres[ch->first] = "";
        for(std::string::iterator s = ch->second.begin(); s != ch->second.end(); ++s){
          if((*s)!='-')
            tseqres[ch->first]+=(*s);
          else
            tseqres[ch->first]+="^^^";
        }
      }
      for( pdb::seqres_record::iterator ch = RseqRes.begin(); ch!=RseqRes.end(); ++ch ){
        rseqres[ch->first] = "";
        for(std::string::iterator s = ch->second.begin(); s != ch->second.end(); ++s){
          if((*s)!='-')
            rseqres[ch->first]+=(*s);
          else
            rseqres[ch->first]+="^^^";
        }
      }
      //for each chain...
	    for( pdb::seqres_record::iterator ch = tseqres.begin(); ch!=tseqres.end(); ++ch ){
	    	typename standard::chain<atom_t>::type tmp_chain;
	    	tmp_chain.type = ch->first;
	      cmp.Append(ch->first, tmp_chain);
	      cmp[ch->first].Reserve( ch->second.size() );
    
	      typename BioCpp::pdb::model<atom_t>::type::iterator it=info.begin(); //first atom in list
	      while( it->chainId != ch->first and it!=info.end() ){ //skip atoms belonging to unwanted chains
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
	      		typename standard::residue<atom_t>::type tmp_res;
            if( ch->second[res]=='^' ){
              tmp_res.type = BioCpp::amino_acid::UNK;
            }
            else{
      	  		std::stringstream sstype;
	        		sstype << ch->second[res];
	      	  	tmp_res.type = BioCpp::amino_acid::string_to_id[ sstype.str() ];
	      	  }
      	    cmp[ch->first].Append(first_res+res,tmp_res);
	    	    ++res;
	      	}
	      	else{
	      		typename standard::residue<atom_t>::type tmp_res;
	    	  	std::stringstream sstype;
	      		sstype << ch->second[res];
	    	  	tmp_res.type = BioCpp::amino_acid::string_to_id[ sstype.str() ];
	    	    cmp[ch->first].Append(first_res+res,tmp_res);
	    	    cmp[ch->first][first_res+res].Reserve(40);
	    	    int prev_res = it->resSeq;
	    	    char prev_icode = it->iCode;
	    	    char prev_altloc = it->altLoc;
	    	  	while(it->resSeq == prev_res and it->iCode==prev_icode ){
	    	  	  //if(it->altLoc!=prev_altloc){
	    	  	  //  ++it;
	    	  	  //  continue;
	    	  	  //}
	    	  	  if( not cmp[ch->first][first_res+res].exists(it->id) ){
	    		  		cmp[ch->first][first_res+res].Append(it->id, *it);
	    				}
       				if( it != (info.end()-1) )
      			  	++it;
      			  else
      			    break;
	    		  }
	    	    ++res;
	      	}
	      }
	    }
      return cmp;
    }

    typename standard::complex<atom_t>::type operator()( typename pdb::model<atom_t>::type& info, pdb::seqres_record& RseqRes){
      return (*this)(info, RseqRes);
    }

};
  
}
}
#endif

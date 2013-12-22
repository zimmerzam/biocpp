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

#include "parseSeqresFromFasta.hpp"

namespace BioCpp{
namespace pdb{

seqres_record parseSeqresFromFasta(std::string& fasta){
  fasta.erase(std::remove_if(fasta.begin(), fasta.end(), [](char ch){if( std::isalpha(ch) or ch=='>' or ch=='-'  ) return false; return true;}), fasta.end());
	seqres_record record;
	unsigned int i = 0;
	while(fasta[i]!='>' or not std::isalpha(fasta[i+1]) ){
		++i;
	}
	char cur_ch;
	for(; i< fasta.size(); ++i){
		if(fasta[i]=='>' and std::isalpha(fasta[i+1]) ){
			++i;
			cur_ch = fasta[i]; 
		}
		else{
			record[cur_ch]+=fasta[i];
		}
	}
	return record;
}

}
}

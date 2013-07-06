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

#ifndef PDB_SEQRES_FORMAT_H
#define PDB_SEQRES_FORMAT_H

#include <string>
#include <cstring>
#include <sstream>
#include <map>
#include <iomanip>
#include <stdexcept>
#include <algorithm>
#include <ctype.h>

#include "pdb_sections_and_records.h"
#include "../polimers/amino_acid_id.h"
namespace BioCpp{
/*! \brief Describes a SEQRES record of a pdb file
		
		@tparam char the chain id
		@tparam std::string the FASTA sequence corresponding to the given chain
		
		\see [FASTA sequence representation](http://en.wikipedia.org/wiki/FASTA_format) */
typedef std::map< char, std::string > pdb_seqres_record;

/*! \brief Print a SEQRES record
		
		Print a SEQRES record according to [pdb specification 3.3 for ATOM record](http://www.wwpdb.org/documentation/format33/sect3.html#SEQRES)	
		@param out the output stream (i.e. `std::cout`)
		@param chainId the chain identifier
		@param sequence the FASTA sequence */
std::ostream& print_pdb_seqres_record( std::ostream& out, char chainId, std::string sequence){
	int numRes = sequence.size();
	int line = numRes/13+1;
	for(int serNum = 0; serNum < line; ++serNum){
		out << "SEQRES ";
		out << std::setw(3) << serNum+1;
		out << " ";
		out << chainId;
		out << " ";
		out << std::setw(4) << numRes;
		out << "  ";
		std::string res;
		for(int i =0; i < std::min(13, numRes-13*serNum); ++i){
			try{
				res = sequence.substr(13*serNum + i,1);
			}
			catch(int e){
				break;
			}
			out << amino_acid::string_to_id[res];
			out << " ";
		}
		out << std::endl;
	}
	return out;
};

/*! \brief Print a SEQRES record 
		
		@param out the output stream (i.e. `std::cout`)
		@param seqres a pdb_seqres_record (i.e. BioCpp::pdb.TseqRes)
		\see print_pdb_seqres_record( std::ostream&, char, std::string), pdb
*/
std::ostream& print_pdb_seqres_record( std::ostream& out, const pdb_seqres_record& seqres){
	for(pdb_seqres_record::const_iterator seq=seqres.begin(); seq!=seqres.end(); ++seq){
		BioCpp::print_pdb_seqres_record(out, seq->first, seq->second);
	}
	return out;
}

/*! \brief Read a SEQRES record from a buffer
		
		\return a pdb_seqres_record object
		@param buffer a buffer containing a pdb SEQRES record
*/
pdb_seqres_record read_pdb_seqres_record( char buffer[] ){
	pdb_seqres_record record;
	char* c_line = strtok(buffer, "\n");
	while(c_line){
		std::string line(c_line, std::find(c_line, c_line + 80, '\0'));
		
		char chain_id;
		std::string s_res = "   ";
		if(get_pdb_record(line) == PDB_SEQRES){
			chain_id = isdigit( line.substr(11, 1).c_str()[0] ) ? char( 'A'+atoi(line.substr(11, 1).c_str()) ) : line.substr(11, 1).c_str()[0];
			chain_id = (chain_id==' ') ? 'A' : chain_id;
			for(int col = 19; col < 70; col+=4){
        try {
          s_res = line.substr(col,3);
        }
        catch (std::out_of_range& oor){
          break;
        }
        if( s_res != "   " ){
          amino_acid::id res = amino_acid::string_to_id[ s_res ];
          if( record.find( chain_id )!=record.end() )
            record[ chain_id ] += amino_acid::id_to_1_letter[res];
          else
            record.insert( std::make_pair(chain_id, amino_acid::id_to_1_letter[res] ) );
        }
      }
		}
		c_line = strtok(NULL, "\n");
	}
	return record;
}

/*! \brief Read a SEQRES record from a FASTA string 
		
		\return a pdb_seqres_record object
		@param fasta a string containing info about the sequences. Chain identifier has to be
		">chain_id_1", i.e. ">A"
*/
pdb_seqres_record read_fasta(std::string& fasta){
  fasta.erase(std::remove_if(fasta.begin(), fasta.end(), [](char ch){if( std::isalpha(ch) or ch=='>' or ch=='-'  ) return false; return true;}), fasta.end());
	pdb_seqres_record record;
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

std::string print_fasta(pdb_seqres_record& seqres){
  std::string fasta = "";
  for(pdb_seqres_record::iterator ch = seqres.begin(); ch!= seqres.end(); ++ch){
    fasta+=">";
    fasta+=ch->first;
    fasta+=ch->second;
  }
  return fasta;
}
} // end namespace
#endif

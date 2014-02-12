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

#include <fstream>
#include <float.h>
#include <BioCpp/base_atom/base_atom.hpp>
#include <BioCpp/io_files/model/model.hxx>
#include "file_reader.hpp"
#include "sections_and_records/sections_and_records.hpp"

namespace BioCpp{
namespace io{
namespace pdb{

file::file(const char* pdb_name, int init_flag ) : BioCpp::io::file( pdb_name, init_flag ){
  filename = pdb_name;
  n_models = 0;
  warning=WAR_NONE;
  error=ERR_NONE;
  std::ifstream file;
  file.open(filename, std::ios::binary);
  
  /* read whole file */
  file.seekg(0,std::ios::beg);
  std::streampos file_beg = file.tellg();
  file.seekg(0,std::ios::end);
  std::streampos file_end = file.tellg();
  file.seekg(0,std::ios::beg);
  int file_length = file_end-file_beg;
  buffer = new char [file_length];
  file.read(buffer,file_length);
  
  file.seekg(std::ios::beg);
  std::streampos prev_pos = file_beg;
  bool atom_section_found = false, seqres_section_found = false;
  std::string line, record;
  std::streamoff seqres_beg = file_beg, seqres_end = file_end;
  bool first_model=true;
  bool first_residue = true;
  Eigen::Vector3d prev_c(DBL_MAX,DBL_MAX,DBL_MAX), cur_n(0,0,0);
  int prev_resseq = -1;
  char prev_chain = 'Z';
  while (file.good()) {
    getline(file, line);
    /* SEQRES section */
    if(get_record(line) == SEQRES ){
      if(seqres_section_found==false){
        seqres_beg = prev_pos;
        seqres_section_found = true;
      }
    }
    /* MODEL section */
    else if(get_record(line) == MODEL){
      n_models++;
      model_beg_pos.insert( std::make_pair(n_models, prev_pos) );
      model_end_pos.insert( std::make_pair(n_models, file_end) );
      if( seqres_beg != file_beg and seqres_end == file_end)
         seqres_end = prev_pos;
    }
    /* end MODEL */
    else if(get_record(line) == ENDMDL){
      if( model_beg_pos.find(n_models)==model_beg_pos.end() ){
        model_beg_pos.insert( std::make_pair(n_models, file_beg) );
      }
      model_end_pos[n_models] = file.tellg();
      if(first_model)
        first_model=false;
      if(init_flag&BioCpp::io::INIT_FIRST_MODEL)
        break;
    }
    /* COORDINATE section */
    else if(get_record(line) == ATOM){ 
      if(atom_section_found==false)
        atom_section_found = true;
      if( seqres_beg != file_beg and seqres_end == file_end)
        seqres_end = prev_pos;
      /* read RseqRes if required */
      if( first_model and not (init_flag&BioCpp::io::INIT_FAST) ){
        BioCpp::base::atom atm=parseAtom<BioCpp::base::atom>(line);
        if(atm.chainId!=prev_chain){
          prev_chain=atm.chainId;
          prev_c=Eigen::Vector3d(DBL_MAX,DBL_MAX,DBL_MAX); 
          cur_n=Eigen::Vector3d(0,0,0);
          prev_resseq = -1;
          first_residue = true;
        }
        if(atm.id==atom::C and (prev_resseq!= -1) ){
          prev_c = atm.coordinate;
          first_residue = false;
        }
        if(atm.resSeq<=prev_resseq or atm.altLoc!=' ' )
          continue;
        else if(atm.id==atom::id::N){
          cur_n = atm.coordinate;
          double pept_bond_length = (cur_n-prev_c).norm();
          if( (pept_bond_length >= 1.5 or pept_bond_length <= 1.2) and not first_residue ){
            RseqRes[atm.chainId] += "-";
            warning = (BioCpp::warning)(warning|PDB_BACKBONE_HOLE);
          }
        }
        else{
          RseqRes[atm.chainId] += residue::dictionary.definition[ atm.resName ].one_letter_name;
          prev_resseq=atm.resSeq;
        }
      }
    }
    if( get_record(line) != SEQRES and seqres_section_found and seqres_end==file_end ){
      seqres_end = prev_pos;
    }
    prev_pos = file.tellg();
  }
  
  file.close();  
  
  /* adjust model counts */
  if(n_models==0 and atom_section_found==true){
    n_models=1;
    model_beg_pos.insert( std::make_pair(1, file_beg) );
    model_end_pos.insert( std::make_pair(1, file_end) );
  }
  else if(atom_section_found==false){
    error = (BioCpp::error)(error|PDB_COORDINATE_NOT_FOUND);
  }

  /* fill TseqRes */
  if( seqres_section_found==true ){
    int length = seqres_end-seqres_beg;
    char seqres[length];
    std::memcpy(seqres, buffer+seqres_beg-file_beg , length);
    TseqRes = parseSeqres( seqres );
  }
  else{
    warning = (BioCpp::warning)(warning|PDB_SEQRES_NOT_FOUND);
  }
}

file::~file(){}

} // end namespace
} // end namespace
} // end namespace

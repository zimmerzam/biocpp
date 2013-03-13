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

#ifndef PDB_SOURCE_H
#define PDB_SOURCE_H

#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <ctype.h>
#include <float.h>
#include <map>

#include "pdb_SEQRES_format.h"
#include "pdb_MODEL_format.h"
#include "pdb_ATOM_format.h"
#include "pdb_sections_and_records.h"
#include "../utils/errors_and_warnings.h"

namespace BioCpp{

/*! \brief These flags can be passed to pdb constructor to customize its behaviour while reading a file */
enum pdb_init_flags{
  PDB_INIT_COMPLETE    = 0,
  PDB_INIT_FAST        = (1<<1),  /*!< Do not read the structure sequence from the first model, but only read SeqRes section */
  PDB_INIT_FIRST_MODEL = (1<<2)   /*!< Read only the first model */
};

/*! \brief This is a `pdb file` object.

    Use this to read a pdb file and get a list of atom_info.
    \example pdb.cpp
*/
class pdb{
  private:
    const char* filename; /*!< The pdb filename. */
    char* buffer;    /*!< The pdb contents. */
    std::map<int, std::streampos> model_beg_pos; /*!< Positions (in the buffer stream) of the begin of each model. */
    std::map<int, std::streampos> model_end_pos; /*!< Positions (in the buffer stream) of the end of each model. */
  public:
    BioCpp::warning warning;
    BioCpp::error error;
    int n_models; /*!< Total number of models found in the file */
    pdb_seqres_record TseqRes; /*!< The sequence as read from SEQRES record */
    pdb_seqres_record RseqRes; /*!< The sequence as read from ATOM record */
    
    /*! \brief Read a pdb file 
        
        This read a pdb file line by line, saves the position of each models and eventually align/adjust the sequence of residues,
        according to the `init_flag` parameter.
        
        @param pdb_name the pdb filename
        @param init_flag can be any combination of pdb_init_flags
        \code
          #include <iostream>
          #include "BioCpp/BioCpp.h"

          int main( int argc, char* argv[] ){
            const char* filename = argc>1 ? argv[1] : "2RNM.pdb"; // if a pdb is passed, read that. else read an example pdb
  
            BioCpp::pdb PDB(filename, 0); // read all the models and get the sequence from ATOM lines of the first model.
            BioCpp::pdb PDB(filename, BioCpp::PDB_INIT_FIRST_MODEL); // read the first model only and get the sequence from ATOM lines.
            BioCpp::pdb PDB(filename, BioCpp::PDB_INIT_FAST); // read all the models and but don't get the sequence from ATOM lines of the first model. 
          }
        \endcode
    */
    pdb(const char* pdb_name, int init_flag);
    
    /*! \brief Destroy the pdb 
        
        \note The pdb file will not be deleted! 
    */
    ~pdb(){}
    
    /*! \brief Retrive informations about a specific model.
    
        \return a container of pdb_atom_info
        \note This can be useful for structure alignment and/or other operations that do not involve structural informations (residues, sequence, ..)
        but are atom-based.
    */
    pdb_model getModel(int mdl);
};

pdb::pdb(const char* pdb_name, int init_flag = (PDB_INIT_FAST|PDB_INIT_FIRST_MODEL) ){
  filename = pdb_name;
  n_models = 0;
  warning=WAR_NONE;
  error=ERR_NONE;
  std::ifstream file;
  file.open(filename, std::ios::binary);
  
  /* read whole file */
  std::streampos file_beg = file.tellg();
  file.seekg(0,std::ios::end);
  std::streampos file_end = file.tellg();
  file.seekg(0,std::ios::beg);
  int file_length = file_end-file_beg;
  buffer = new char [file_length];
  file.read(buffer,file_length);
  
  std::streampos prev_pos = file_beg;
  file.seekg(std::ios::beg);
  bool atom_section_found = false;
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
    if(get_pdb_record(line) == PDB_SEQRES ){
      if(seqres_beg == file_beg)
        seqres_beg = prev_pos;
    }
    /* MODEL section */
    else if(get_pdb_record(line) == PDB_MODEL){
      n_models++;
      model_beg_pos.insert( std::make_pair(n_models, prev_pos) );
      if( seqres_beg != file_beg and seqres_end == file_end)
         seqres_end = prev_pos;
    }
    /* end MODEL */
    else if(get_pdb_record(line) == PDB_ENDMDL){
      model_end_pos.insert( std::make_pair(n_models, prev_pos) );
      if(first_model)
        first_model=false;
      if(init_flag&PDB_INIT_FIRST_MODEL)
        break;
    }
    /* COORDINATE section */
    else if(get_pdb_record(line) == PDB_ATOM){ 
      if(atom_section_found==false)
        atom_section_found = true;
      if( seqres_beg != file_beg and seqres_end == file_end)
        seqres_end = prev_pos;
      /* read RseqRes if required */
      if( first_model and not (init_flag&PDB_INIT_FAST) ){
        pdb_atom_info atm=read_pdb_atom_line(line);
        if(atm.chainId!=prev_chain){
          prev_chain=atm.chainId;
          prev_c=Eigen::Vector3d(DBL_MAX,DBL_MAX,DBL_MAX); 
          cur_n=Eigen::Vector3d(0,0,0);
          prev_resseq = -1;
          first_residue = true;
        }
        if(atm.id==atom::C_ and (prev_resseq!= -1) ){
          prev_c = atm.coordinate;
          first_residue = false;
        }
        if(atm.resSeq<=prev_resseq)
          continue;
        else{
          RseqRes[atm.chainId] += amino_acid::id_to_1_letter[ atm.resName ];
          prev_resseq=atm.resSeq;
        }
        if(atm.id==atom::N_){
          cur_n = atm.coordinate;
          double pept_bond_length = (cur_n-prev_c).norm();
          if( (pept_bond_length >= 1.4 or pept_bond_length <= 1.27) and not first_residue ){
            RseqRes[atm.chainId] += "-";
            warning = (BioCpp::warning)(warning|PDB_BACKBONE_HOLE);
          }
        }
      }
    }
    else{ // line is different from MODEL, ENDMDL, SEQRES or ATOM
      if( seqres_beg != file_beg and seqres_end == file_end)
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
  if( seqres_beg!=file_beg or seqres_end!=file_end ){
    int length = seqres_end-seqres_beg;
    char seqres[length+1];
    std::memcpy(seqres, buffer + seqres_beg, length);
    TseqRes = read_pdb_seqres_record( seqres );
  }
  else{
    warning = (BioCpp::warning)(warning|PDB_SEQRES_NOT_FOUND);
  }
}

pdb_model pdb::getModel(int mdl){
  if(mdl>n_models){
    pdb_model record;
    return record;
  }
  int length = model_end_pos[mdl] - model_beg_pos[mdl];
  char mdl_buf[ length + 1 ];
  std::memcpy(mdl_buf, buffer + model_beg_pos[mdl], length);
  return read_pdb_model_record( mdl_buf );
}

} // end namespace

#endif

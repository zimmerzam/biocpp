#include <set>
#include <sstream>
#include <string>
#include <algorithm>
#include <getopt.h>
#include <BioCpp.h>

struct mutation{
  char chainId;
  int resSeq;
  BioCpp::amino_acid::id from;
  BioCpp::amino_acid::id to;
};

struct strict_residue_printer{
  BioCpp::pdb::print_atom_line printer;
  
  strict_residue_printer(std::ostream& dev): printer(dev) {};
  
  void operator()( BioCpp::standard::residue& res ){
    for(BioCpp::standard::residue::iterator at = res.begin(); at != res.end(); ++at){
      /* always print backbone atoms */
      if( at->id == BioCpp::atom::id::N_ or at->id == BioCpp::atom::id::CA or 
          at->id == BioCpp::atom::id::C_ or at->id == BioCpp::atom::id::OXT or
          at->id == BioCpp::atom::id::H_ or at->id == BioCpp::atom::id::O_ or
          at->id == BioCpp::atom::id::HA ){
        at->resName = res.type;
        printer(*at);
      }
      /* If residue is ALA */
      else if(res.type == BioCpp::amino_acid::id::ALA){
        if( at->id == BioCpp::atom::id::CB or at->id == BioCpp::atom::id::HB1 or
            at->id == BioCpp::atom::id::HB2 or at->id == BioCpp::atom::id::HB3 ) {
          at->resName = res.type;
          printer(*at);
        }
      }
      /* If residue is GLY */
      else if(res.type == BioCpp::amino_acid::id::ALA){
        if( at->id == BioCpp::atom::id::HA2 or at->id == BioCpp::atom::id::HA3 ) {
          at->resName = res.type;
          printer(*at);
        }
      }
      /* If residue is not ALA nor GLY and atom is a side-chain ones */
      else{
        at->resName = res.type;
        printer(*at);
      }
    }
  }
  
};

int main(int argc, char* argv[]){
  bool file_flag = false;
  bool mutation_flag = false;
  char* contactfile;
  std::set<std::string> mut;
  int c;
	while ((c = getopt (argc, argv, "f:m:")) != -1){
		switch (c){
			case 'f':
			  file_flag = true;
				contactfile = optarg;
				break;
      case 'm':
        mutation_flag = true;
        std::stringstream smut;
        smut << optarg;
        mut.insert( smut.str() );
        break;
		}
	}

  if( not file_flag or not mutation_flag){
    std::cout << "usage: .mutations -f 'structure.pdb' -m mutation1 -m mutation2" << std::endl
              << "Options: " << std::endl 
              << "\t-m:  mutation (ex: thrA237ala)" << std::endl;
    return 1;
  }
  bool exit_asap = false;
  std::vector<mutation> mutations;
  for(std::set<std::string>::iterator it = mut.begin(); it != mut.end(); ++it){
    if(it->size() < 8){
      exit_asap = true;
      std::cout << "Error: " << *it << " is not a valid mutation" << std::endl;
      continue;
    }
    std::string str = *it;
    std::transform(str.begin(), str.end(), str.begin(), ::toupper);
    mutation tmpmut;
    tmpmut.chainId = str.substr(3,1)[0];
    tmpmut.resSeq = atoi( str.substr(4,str.size()-7).c_str() );
    tmpmut.from = BioCpp::amino_acid::string_to_id[ str.substr(0,3) ];
    tmpmut.to = BioCpp::amino_acid::string_to_id[ str.substr(str.size()-3,3) ];
    if(tmpmut.to != BioCpp::amino_acid::id::ALA and tmpmut.to != BioCpp::amino_acid::id::GLY){
      exit_asap = true;
      std::cout << "Sorry! At the moment only mutations to ALA or GLY are supported" << std::endl;
      continue;
    }
    mutations.push_back(tmpmut);
  }
  if(mutations.size() == 0){
    std::cout << "Error: not even a valid mutation found. Bye!" << std::endl;
    return 1;
  }
  strict_residue_printer printer(std::cout);
  BioCpp::pdb::pdb PDB(contactfile, 0);
  for(int mdl = 1; mdl <= PDB.n_models; ++mdl){
    BioCpp::pdb::model all_info = PDB.getModel(mdl);
    BioCpp::standard::complex cmp( all_info, PDB.RseqRes, PDB.RseqRes );
    for( std::vector<mutation>::iterator mt = mutations.begin(); mt != mutations.end(); ++mt ){
      if( not cmp.exists(mt->chainId) ){
        std::cout << "Error: chain " << mt->chainId << " is not present in the target structure" << std::endl;
        exit_asap = true;
        continue;
      }
      if( not cmp[mt->chainId].exists(mt->resSeq) ){
        std::cout << "Error: residue " << mt->resSeq << " is not present in chain " << mt->chainId << std::endl;
        exit_asap = true;
        continue;
      }
      if( cmp[mt->chainId][mt->resSeq].type != mt->from ){
        exit_asap = true;
        std::cout << "Warning: residue " << mt->chainId << mt->resSeq << " is not an " << mt->from << std::endl;
      }
      cmp[mt->chainId][mt->resSeq].type = mt->to;
    }
    if(exit_asap){
      return 1;
    }
    std::cout << "MODEL     1" << std::endl;
    BioCpp::Iterate<BioCpp::standard::residue>(cmp, printer);    
    std::cout << "ENDMDL     " << std::endl;
  }

  return 0;
}

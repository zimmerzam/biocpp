/*!
    \file statistics_o_position.cpp
*/

#include <iostream>
#include <BioCpp/Core>
#include <BioCpp/Standard>

typedef BioCpp::base::atom atom;
typedef BioCpp::standard::base::residue residue;
typedef BioCpp::standard::base::chain chain;
typedef BioCpp::standard::base::complex complex;
typedef BioCpp::standard::base::complex_constructor complex_constructor;
typedef BioCpp::io::model<atom>::type model;

struct positions{
	Eigen::Vector3d ca;
	Eigen::Vector3d c;
	Eigen::Vector3d n;
	Eigen::Vector3d o;
};

int main( int argc, char* argv[] ){
	for( int file_num = 1; file_num!=argc; ++file_num ){
	  const char* filename = argv[file_num];
  	BioCpp::io::pdb::file PDB(filename, 0); // read the pdb file. 
  	bool is_good_structure = true;
  	
  	// gets the structure, align the sequence from ATOM lines to the one from SEQRES and builds the topology
		// if SEQRES has not been found and the sequence in ATOM has at least one hole raise an error. if the
		// sequence from ATOM lines has no holes no error will be raised and the topology will be built.
		model mdl = PDB.readModel<atom>(1); 
		BioCpp::error err = PDB.error;
  	BioCpp::warning war = PDB.warning;
		if(PDB.TseqRes.size() != 0){
  		for( BioCpp::io::seqres_record::iterator seqres = PDB.TseqRes.begin(); seqres!=PDB.TseqRes.end(); ++seqres ){
 	    	if( PDB.RseqRes.find(seqres->first)!=PDB.RseqRes.end() and seqres->second!=PDB.RseqRes[seqres->first]){
      		PDB.RseqRes[seqres->first].insert( PDB.RseqRes[seqres->first].begin(), '-' );
      		PDB.RseqRes[seqres->first].insert( PDB.RseqRes[seqres->first].end(), '-' );
        	double score = BioCpp::fasta::StrictNeedlemanWunsch(seqres->second, PDB.RseqRes[seqres->first], BioCpp::fasta::ZIMM1, err, war);
        	if(score < 0){
          	is_good_structure = false;
        	}
      	}
    	}
  	}
  	else{
   		for( BioCpp::io::seqres_record::iterator seqres = PDB.RseqRes.begin(); seqres!=PDB.RseqRes.end(); ++seqres ){
      	for( std::string::iterator str = seqres->second.begin(); str != seqres->second.end(); ++str ){
        	if((*str)=='-'){
          	is_good_structure = false;
        	}
      	}
    	}
  	}
  	if( not is_good_structure ){
  		break;
  	}
  	// Order atoms in 'all_info' by using informations stored in the pdb: missing residues are 
  	//detected and stored as empty residues 
  	BioCpp::standard::base::complex_constructor constr(BioCpp::residue::dictionary);
  	BioCpp::standard::base::complex cmp = constr(mdl, PDB.RseqRes, PDB.TseqRes);
  	
  	std::list<positions> position_list;
  	for( complex::iterator ch = cmp.begin(); ch!=cmp.end(); ++ch ){
  		chain::iterator last = ch->end(); --last;
  		for( chain::iterator res = ch->begin(); res!=last; ++res ){
  			chain::iterator sres = res; ++sres;
  			if( res->exists( BioCpp::atom::CA ) and res->exists( BioCpp::atom::C ) and res->exists( BioCpp::atom::O ) and sres->exists( BioCpp::atom::N ) ){
  				positions pos;
  				pos.ca = (*res)[BioCpp::atom::CA].coordinate;
  				pos.c =  (*res)[BioCpp::atom::C].coordinate;
  				pos.o =  (*res)[BioCpp::atom::O].coordinate;
  				pos.n = (*sres)[BioCpp::atom::N].coordinate;
  				position_list.push_back( pos );
  			} 
  		}
  	}
  	
  	for(std::list<positions>::iterator p = position_list.begin(); p != position_list.end(); ++p){
  		Eigen::Vector3d ca_c = p->c - p->ca;
  		Eigen::Vector3d c_n = p->n - p->c;
  		Eigen::Vector3d c_o = p->o - p->c;
  		
  		ca_c.normalize();
  		c_n.normalize();
  		Eigen::Vector x = ca_c+c_n;
  		Eigen::Vector y = ca_c-c_n;
  		Eigen::Vector z = c_n.cross(ca_c);
  		x.normalize();
  		y.normalize();
  		z.normalize();
  		
  	}
  	
  }
  
  return 0;
}

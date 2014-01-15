/*!
    \file topology_constructor.cpp
    \brief Read a pdb file, add missing atoms to the structure and print the structure back
    
    In this example we see how to build a topology from a pdb file. In order to see
    the effect of the building operation we print back the structure in pdb format. Eventual
    missing atoms should be added at the end of the stream.
*/

#include <iostream>
#include <BioCpp/base_atom/base_atom.hpp>
#include <BioCpp/standard/fileFormat/pdb/file_reader.hpp>
#include <BioCpp/topology/topology.hxx>
#include <BioCpp/standard/proteinTopology/topology_constructor.hxx>
#include <BioCpp/standard/fileFormat/pdb/printer/print_atom_t.hxx>
#include <BioCpp/base_container/Iterate_single.hxx>
#include <BioCpp/fasta/StrictNeedlemanWunsch.hpp>
#include <BioCpp/standard/alignmentMatrix/ZIMM1.hpp>
#include <boost/graph/breadth_first_search.hpp>

typedef BioCpp::base::atom atom;
typedef BioCpp::io::model<atom>::type model;
typedef BioCpp::element::dictionary_t element_dictionary;
typedef BioCpp::atom::dictionary_t atom_dictionary;
typedef BioCpp::residue::dictionary_t residue_dictionary;
typedef BioCpp::io::pdb::print_atom_t< element_dictionary, atom_dictionary, residue_dictionary > print_atom;

typedef BioCpp::standard::vertex<atom> vertex; 
typedef BioCpp::standard::edge< std::pair<atom,atom> > edge;
typedef BioCpp::topology<vertex, edge> topology;
typedef BioCpp::standard::topology_constructor<vertex,edge> topology_constructor;
typedef BioCpp::io::pdb::print_atom_t< element_dictionary, atom_dictionary, residue_dictionary > print_atom;

struct bfs_print_visitor : public boost::default_bfs_visitor {
	print_atom printer;
	bfs_print_visitor( element_dictionary& eleDict, atom_dictionary& atmDict, residue_dictionary& resDict ): printer(std::cout, eleDict, atmDict, resDict) {};
	
  void discover_vertex(topology::vertex_t u, const topology::Graph& G) {
  	printer( *(G[u].atom) );
  };
  
};  

int main( int argc, char* argv[] ){
  const char* pdb_filename = argv[1];
  const char* dic_filename = argv[2];
  
  // reading dictionaries
  libconfig::Config cfg;
  cfg.readFile(dic_filename);
  libconfig::Setting& root = cfg.getRoot();
  element_dictionary eleDict;
  eleDict.importSetting(root,{"elements"});
  atom_dictionary atmDict;
  atmDict.importSetting(root,{"atoms"});
  residue_dictionary resDict;
  resDict.importSetting(root,{"residues"}); 
  
  // reading pdb file
  BioCpp::io::pdb::file PDB(pdb_filename, 0); // read the pdb file. 
	
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
          std::cout << "Aligment failed!" << std::endl;
          return 1;
        }
      }
    }
  }
  else{
    for( BioCpp::io::seqres_record::iterator seqres = PDB.RseqRes.begin(); seqres!=PDB.RseqRes.end(); ++seqres ){
      for( std::string::iterator str = seqres->second.begin(); str != seqres->second.end(); ++str ){
        if((*str)=='-'){
          std::cout << "The structure has at least one hole and the SEQRES section is missing. A complete topology cannot be built." << std::endl;
          return 1;
        }
      }
    }
  } 
  topology_constructor topo_constr;
	topology topo = topo_constr( mdl, PDB.RseqRes, PDB.TseqRes, atmDict, resDict );
	
	topology::vertex_iterator v_beg, v_end, v;
  boost::tie(v_beg, v_end) = boost::vertices( topo.getGraph() );
	bfs_print_visitor printer(eleDict, atmDict, resDict);
	boost::breadth_first_search( topo.getGraph(), *v_beg, boost::visitor(printer) );
	
  return 0;
}

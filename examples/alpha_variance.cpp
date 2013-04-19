#include <map>
#include "../src/BioCpp_default.h"

typedef std::pair<char, int> residue;

struct alphaVariance{
	std::map<residue, std::vector< Eigen::Vector3d > > data;
	
	void operator()(BioCpp::standard::residue& res){
		if( not res.exists(BioCpp::atom::id::CA) ){
			return;
		}
		char chainId = res[ BioCpp::atom::id::CA ].chainId;
		int resSeq = res[ BioCpp::atom::id::CA ].resSeq;
		Eigen::Vector3d coordinate = res[ BioCpp::atom::id::CA ].coordinate;
		
		residue r = std::make_pair(chainId, resSeq);
		data[r].push_back( coordinate );
	}
	
	void compute_and_print(){
		for( std::map<residue, std::vector< Eigen::Vector3d > >::iterator r = data.begin(); r != data.end(); ++r ){
			Eigen::Vector3d mean = Eigen::Vector3d::Zero();
			for( std::vector< Eigen::Vector3d >::iterator c = r->second.begin(); c != r->second.end(); ++c ){
				mean += (*c);
			}
			mean /= r->second.size();
			double var = 0;
			for( std::vector< Eigen::Vector3d >::iterator c = r->second.begin(); c != r->second.end(); ++c ){
				var += ( (*c) - mean ).norm() * ( (*c) - mean ).norm()/( r->second.size()-1 );
			}
			std::cout << var << "  " << r->second.size() << std::endl;
		}
	}
};

int main(int argc, char* argv[]){
	bool file_flag = false;
  char* pdbfile;
  int c;
	while ((c = getopt (argc, argv, "f:")) != -1){
		switch (c){
			case 'f':
			  file_flag = true;
				pdbfile = optarg;
				break;
		}
	}

  if( not file_flag ){
    std::cout << "usage: .alpha_variance -f 'structure.pdb' " << std::endl;
    return 1;
  }

	alphaVariance variance;
  BioCpp::pdb PDB(pdbfile, 0);
  for(int mdl = 1; mdl <= PDB.n_models; ++mdl){
    BioCpp::pdb_model all_info = PDB.getModel(mdl);
    BioCpp::standard::complex cmp( all_info, PDB.RseqRes );
    BioCpp::Iterate< BioCpp::standard::residue >(cmp, variance);
  }
	variance.compute_and_print();
	
	return 0;
}

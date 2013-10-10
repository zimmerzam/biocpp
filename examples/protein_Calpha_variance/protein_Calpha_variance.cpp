#include <map>
#include <list>
#include <getopt.h>
#include <BioCpp.h>

typedef std::pair<char, int> residue;

struct alphaVariance{
	std::map<residue, std::list< Eigen::Vector3d > > data;
	
	void operator()(BioCpp::standard::base::residue& res){
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
		for( std::map<residue, std::list< Eigen::Vector3d > >::iterator r = data.begin(); r != data.end(); ++r ){
			Eigen::Vector3d mean = Eigen::Vector3d::Zero();
			for( std::list< Eigen::Vector3d >::iterator c = r->second.begin(); c != r->second.end(); ++c ){
				mean += (*c);
			}
			mean /= r->second.size();
			double var = 0;
			for( std::list< Eigen::Vector3d >::iterator c = r->second.begin(); c != r->second.end(); ++c ){
				var += ( (*c) - mean ).norm() * ( (*c) - mean ).norm()/( r->second.size()-1 );
			}
			std::cout << r->first.first << r->first.second << "  " << var << std::endl;
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
  BioCpp::pdb::pdb PDB(pdbfile, 0);
  for(int mdl = 1; mdl <= PDB.n_models; ++mdl){
    std::cout << "Now processing model " << mdl << "/" << PDB.n_models << std::endl;
    BioCpp::standard::base::model all_info = PDB.getModel<BioCpp::standard::base::atom>(mdl);
    BioCpp::standard::base::complex_constructor cmp_constr;
    BioCpp::standard::base::complex cmp = cmp_constr( all_info, PDB.RseqRes, PDB.RseqRes );
    BioCpp::Iterate< BioCpp::standard::base::residue >(cmp, variance);
  }
	variance.compute_and_print();
	
	return 0;
}

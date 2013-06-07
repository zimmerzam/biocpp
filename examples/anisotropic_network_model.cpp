#include <BioCpp_standard.h>
#include <geometry/Eigen/Dense>
#include <getopt.h>

struct anisotropic{
	typedef std::map< char, std::pair<unsigned int, unsigned int> > range_t;
	typedef std::map< char, int> offset_t;
	double p;
	double k1, k2;
	range_t range;
	offset_t offset;
	
	Eigen::MatrixXd matrix;
	
	anisotropic(unsigned int size, range_t ran):range(ran), matrix(Eigen::MatrixXd::Zero(size,size)){
		range_t::iterator it = ran.begin();
		offset[it->first] = 0;
		std::cout << it->first << "  " << offset[it->first] << "  " << it->second.first << "  " << it->second.second << std::endl;
		++it;
		for(range_t::iterator itm = ran.begin(); it != ran.end(); ++it, ++itm){
			offset[it->first] = offset[itm->first] + itm->second.second-itm->second.first;
			std::cout << it->first << "  " << offset[it->first] << "  " << it->second.first << "  " << it->second.second << std::endl;
		}
	}
	
	void operator()( BioCpp::standard::residue res1, BioCpp::standard::residue res2 ){
		char ch1  = res1[BioCpp::atom::id::CA].chainId;
		int nres1 = res1[BioCpp::atom::id::CA].resSeq;
		
		char ch2  = res2[BioCpp::atom::id::CA].chainId;
		int nres2 = res2[BioCpp::atom::id::CA].resSeq;
		
		int i = offset[ch1] + nres1 - range[ch1].first + 1;
		int j = offset[ch2] + nres2 - range[ch2].first + 1;
		
		Eigen::Vector3d c1 = res1[BioCpp::atom::id::CA].coordinate;
		Eigen::Vector3d c2 = res2[BioCpp::atom::id::CA].coordinate;
		
		double factor = 0;
		double dist_p = pow( (c1-c2).norm(),2.+ p );
		if( std::abs(j-i)==1 ){
			factor = -k1/dist_p;
		}
		else if( std::abs(j-i)>1 ){
			factor = -k2/dist_p;
		}
		
		matrix(3*i,3*j)     = factor*( c2(0)-c1(0) )*( c2(0)-c1(0) );
		matrix(3*i,3*j+1)   = factor*( c2(0)-c1(0) )*( c2(1)-c1(1) );
		matrix(3*i,3*j+2)   = factor*( c2(0)-c1(0) )*( c2(2)-c1(2) );
		matrix(3*i+1,3*j)   = factor*( c2(1)-c1(1) )*( c2(0)-c1(0) );
		matrix(3*i+1,3*j+1) = factor*( c2(1)-c1(1) )*( c2(1)-c1(1) );
		matrix(3*i+1,3*j+2) = factor*( c2(1)-c1(1) )*( c2(2)-c1(2) );
		matrix(3*i+2,3*j)   = factor*( c2(2)-c1(2) )*( c2(0)-c1(0) );
		matrix(3*i+2,3*j+1) = factor*( c2(2)-c1(2) )*( c2(1)-c1(1) );
		matrix(3*i+2,3*j+2) = factor*( c2(2)-c1(2) )*( c2(2)-c1(2) );
	}
};

int main(int argc, char* argv[]){
	char* pdbfilename;
	double p = 0;
	double k1 = 1., k2 = 1.;
  bool file_flag = false;
  bool usage_flag = false;
  bool eigenvalues_flag = true;
  bool eigenvectors_flag = false;
  bool mobility_flag = false;
  bool entropy_flag = false;
  
	int c;
	while ((c = getopt (argc, argv, "f:p:K:k:h")) != -1){
		switch (c){
			case 'f':
				file_flag = true;
				pdbfilename = optarg;
				break;
			case 'p':
				p = atof(optarg);
				break;
			case 'K':
				k1 = atof(optarg);
				break;
			case 'k':
				k2 = atof(optarg);
				break;
			case 'h':
				usage_flag = true;
				break;
		}
	}
	
	if( not file_flag or usage_flag){
    std::cout << "usage: .anisotropic_network_model -f file.pdb [options]" << std::endl
    					<< "\t-p: rescale distance exponent (default = 0)" << std::endl
    					<< "\t-K: spring constant for adjacent residues (default = 1.)" << std::endl
    					<< "\t-k: spring constant for other residues (default = 1.)" << std::endl
    					<< "\t-h: help" << std::endl;
    return 1;
  }
	
	BioCpp::pdb PDB(pdbfilename, BioCpp::PDB_INIT_COMPLETE);
	BioCpp::pdb_model all_info = PDB.getModel(1);
	unsigned int n_res = 0;
	for(BioCpp::pdb_seqres_record::iterator it = PDB.RseqRes.begin(); it!=PDB.RseqRes.end(); ++it){
		n_res+=it->second.size();
	}
	BioCpp::standard::complex cmp(all_info, PDB.RseqRes);
	
	std::map< char, std::pair<unsigned int, unsigned int> > range;
	for( BioCpp::standard::complex::iterator it = cmp.begin(); it != cmp.end(); ++it ){
		range[it->type] = std::make_pair( (*(it->begin()))[BioCpp::atom::id::CA].resSeq, (*(it->rbegin()))[BioCpp::atom::id::CA].resSeq );
	}
	
	anisotropic hessian( 3*n_res, range );
	hessian.p = p;
	hessian.k1 = k1;
	hessian.k2 = k2;
  
  BioCpp::Iterate<BioCpp::standard::residue, BioCpp::standard::residue>(cmp,cmp,hessian);

	int size = 3*n_res;

	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(hessian.matrix);
  if (eigensolver.info() != Eigen::Success) abort();
  // print eigenvalues
  if(eigenvalues_flag){
    std::cout << eigensolver.eigenvalues() << std::endl;
  }
  // print eigenvectors
  else if(eigenvectors_flag){
    std::cout << eigensolver.eigenvectors() << std::endl;
  }
  // print mobility
  else if(mobility_flag){
    Eigen::ArrayXd mobility = Eigen::VectorXd::Zero(size).array();
    for(int c = 1; c!= size; ++c){
      Eigen::ArrayXd mode = eigensolver.eigenvectors().col(c).array()*eigensolver.eigenvectors().col(c).array();
      mobility += mode/eigensolver.eigenvalues()[c];
    }
    std::cout << mobility << std::endl;
  }
  // print entropy
  else if(entropy_flag){
  	Eigen::ArrayXd entropy = Eigen::VectorXd::Zero(size).array();
    for(int c = 1; c!= size; ++c){
      Eigen::ArrayXd mode = eigensolver.eigenvectors().col(c).array()*eigensolver.eigenvectors().col(c).array();
      entropy -= mode*log( eigensolver.eigenvalues()[c] );
    }
    std::cout << entropy << std::endl;
  }
  
	return 0;
}

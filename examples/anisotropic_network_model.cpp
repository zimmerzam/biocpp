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
	Eigen::MatrixXd mass;
	
	std::map<BioCpp::amino_acid::id, double> residue_mass;
	
	anisotropic(unsigned int size, range_t ran):range(ran), matrix(Eigen::MatrixXd::Zero(size,size)), mass(Eigen::MatrixXd::Zero(size,size)){
		range_t::iterator it = ran.begin();
		offset[it->first] = 0;
//		std::cout << it->first << "  " << offset[it->first] << "  " << it->second.first << "  " << it->second.second << std::endl;
		++it;
		for(range_t::iterator itm = ran.begin(); it != ran.end(); ++it, ++itm){
			offset[it->first] = offset[itm->first] + itm->second.second-itm->second.first+1;
//			std::cout << it->first << "  " << offset[it->first] << "  " << it->second.first << "  " << it->second.second << std::endl;
		}
		residue_mass = std::map<BioCpp::amino_acid::id, double>
										 { { BioCpp::amino_acid::ALA, 71.079 },{ BioCpp::amino_acid::ARG, 156.188 },
										 { BioCpp::amino_acid::ASN, 114.104 },{ BioCpp::amino_acid::ASP, 115.089 },
         						 { BioCpp::amino_acid::CYS, 103.144 },{ BioCpp::amino_acid::GLN, 128.131 },
         						 { BioCpp::amino_acid::GLU, 129.116 },{ BioCpp::amino_acid::GLY, 57.052 },
        						 { BioCpp::amino_acid::HIS, 137.142 },{ BioCpp::amino_acid::ILE, 113.160 },
        						 { BioCpp::amino_acid::LEU, 113.160 },{ BioCpp::amino_acid::LYS, 128.174 },
         						 { BioCpp::amino_acid::MET, 131.198 },{ BioCpp::amino_acid::PHE, 147.177 },
         						 { BioCpp::amino_acid::PRO, 97.117 },{ BioCpp::amino_acid::SER, 87.078 },
        						 { BioCpp::amino_acid::THR, 101.105 },{ BioCpp::amino_acid::TRP, 186.213 },
        						 { BioCpp::amino_acid::TYR, 163.170 },{ BioCpp::amino_acid::VAL, 99.133 },
         					 };
	}
	
	void operator()( BioCpp::standard::residue res1, BioCpp::standard::residue res2 ){
		char ch1  = res1[BioCpp::atom::id::CA].chainId;
		int nres1 = res1[BioCpp::atom::id::CA].resSeq;
		
		char ch2  = res2[BioCpp::atom::id::CA].chainId;
		int nres2 = res2[BioCpp::atom::id::CA].resSeq;
		
		int i = offset[ch1] + nres1 - range[ch1].first;
		int j = offset[ch2] + nres2 - range[ch2].first;
		
		Eigen::Vector3d c1 = res1[BioCpp::atom::id::CA].coordinate;
		Eigen::Vector3d c2 = res2[BioCpp::atom::id::CA].coordinate;
		
		double factor = 0;
		double dist_p = pow( (c1-c2).norm(),2.+ p );
		if( j-i < 0 ){
			return;
		}
		else if(j-i==0){
			mass(3*i,3*i) = 1./sqrt( residue_mass[ res1[BioCpp::atom::id::CA].resName ] );
			mass(3*i+1,3*i+1) = mass(3*i,3*i);
			mass(3*i+2,3*i+2) = mass(3*i,3*i);
			return;
		}
		else if( std::abs(j-i)==1 ){ //TODO they have to belong to the same chain!!
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
		
		matrix(3*j,3*i)     = matrix(3*i,3*j);
		matrix(3*j,3*i+1)   = matrix(3*i,3*j+1);
		matrix(3*j,3*i+2)   = matrix(3*i,3*j+2);
		matrix(3*j+1,3*i)   = matrix(3*i+1,3*j);
		matrix(3*j+1,3*i+1) = matrix(3*i+1,3*j+1);
		matrix(3*j+1,3*i+2) = matrix(3*i+1,3*j+2);
		matrix(3*j+2,3*i)   = matrix(3*i+2,3*j);
		matrix(3*j+2,3*i+1) = matrix(3*i+2,3*j+1);
		matrix(3*j+2,3*i+2) = matrix(3*i+2,3*j+2);
	}
	
	void fill_diagonal(){
		int n_res = matrix.cols()/3;
	
		for( int i=0; i<n_res; ++i  ){
			for( int j=0; j<n_res; ++j  ){
				if(i==j){
					continue;
				}
				matrix(3*i,3*i)     -= matrix(3*i,3*j);
				matrix(3*i,3*i+1)   -= matrix(3*i,3*j+1);
				matrix(3*i,3*i+2)   -= matrix(3*i,3*j+2);
				matrix(3*i+1,3*i)   -= matrix(3*i+1,3*j);
				matrix(3*i+1,3*i+1) -= matrix(3*i+1,3*j+1);
				matrix(3*i+1,3*i+2) -= matrix(3*i+1,3*j+2);
				matrix(3*i+2,3*i)   -= matrix(3*i+2,3*j);
				matrix(3*i+2,3*i+1) -= matrix(3*i+2,3*j+1);
				matrix(3*i+2,3*i+2) -= matrix(3*i+2,3*j+2);
			}
		}
	}
	
	Eigen::VectorXd eigenvalues(){
		Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(mass*matrix*mass);
		return eigensolver.eigenvalues();
	}
	
	Eigen::MatrixXd eigenvectors(){
		Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(mass*matrix*mass);
		return mass*eigensolver.eigenvectors();
	}
};

int main(int argc, char* argv[]){
	char* pdbfilename;
	double p = 0;
	double k1 = 1., k2 = 1.;
  bool file_flag = false;
  bool usage_flag = false;
  bool eigenvalues_flag = false;
  bool eigenvectors_flag = false;
  bool mobility_flag = false;
  bool entropy_flag = false;
  bool other_flag = false;
  bool too_much_flags = false;
  
	int c;
	while ((c = getopt (argc, argv, "f:p:K:k:vVmeh")) != -1){
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
      case 'v':
        eigenvalues_flag=true;
        if(other_flag)
          too_much_flags=true;
        else
          other_flag=true;
        break;
      case 'V':
        eigenvectors_flag=true;
        if(other_flag)
          too_much_flags=true;
        else
          other_flag=true;
        break;
      case 'm':
        mobility_flag = true;
        if(other_flag)
          too_much_flags = true;
        else
          other_flag=true;
        break;
      case 'e':
        entropy_flag = true;
        if(other_flag)
          too_much_flags = true;
        else
          other_flag=true;
        break;
			case 'h':
				usage_flag = true;
				break;
		}
	}
	
	if( not file_flag or usage_flag or too_much_flags){
    std::cout << "usage: .anisotropic_network_model -f file.pdb [options]" << std::endl
    					<< "Options: " << std::endl
    					<< "\t-p: rescale distance exponent (default = 0)" << std::endl
    					<< "\t-K: spring constant for adjacent residues (default = 1.)" << std::endl
    					<< "\t-k: spring constant for other residues (default = 1.)" << std::endl
    					<< "\t-v: print eigenvalues" << std::endl
              << "\t-V: print eigenvectors" << std::endl
              << "\t-e: print entropy" << std::endl
              << "\t-m: print mobility" << std::endl
    					<< "\t-h: help" << std::endl
    					<< "No holes are allowed!! Residues in the same chain must have consecutive resSeq" << std::endl
    					<< "You can use only one option among 'v', 'V', 'm' and 'e'.";
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
	hessian.fill_diagonal();
	
	int size = 3*n_res;

  // print eigenvalues
  if(eigenvalues_flag){
    std::cout << hessian.eigenvalues() << std::endl;
  }
  // print eigenvectors
  else if(eigenvectors_flag){
    std::cout << hessian.eigenvectors() << std::endl;
  }
  // print mobility
  else if(mobility_flag){
    Eigen::ArrayXd mobility = Eigen::VectorXd::Zero(size).array();
    Eigen::VectorXd evalues = hessian.eigenvalues();
    Eigen::MatrixXd evector = hessian.eigenvectors();
    for(int c = 1; c!= size; ++c){
      Eigen::ArrayXd mode = evector.col(c).array()*evector.col(c).array();
      mobility += mode/evalues[c];
    }
    std::cout << mobility << std::endl;
  }
  // print entropy
  else if(entropy_flag){
  	Eigen::ArrayXd entropy = Eigen::VectorXd::Zero(size).array();
  	Eigen::VectorXd evalues = hessian.eigenvalues();
    Eigen::MatrixXd evector = hessian.eigenvectors();
    for(int c = 1; c!= size; ++c){
      Eigen::ArrayXd mode = evector.col(c).array()*evector.col(c).array();
      entropy -= mode*log( evalues[c] );
    }
    std::cout << entropy << std::endl;
  }
  
	return 0;
}

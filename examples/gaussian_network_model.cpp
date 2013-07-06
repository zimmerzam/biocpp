#include <string>
#include <sstream>
#include <istream>
#include <iterator>
#include <fstream>
#include <iostream>
#include <map>
#include <getopt.h>
#include <BioCpp/geometry/Eigen/Dense>

typedef std::map<char, std::pair<int, int> > contact_info;

struct tokens: std::ctype<char>{
  tokens(): std::ctype<char>(get_table()) {}

  static std::ctype_base::mask const* get_table(){
    typedef std::ctype<char> cctype;
    static const cctype::mask *const_rc= cctype::classic_table();

    static cctype::mask rc[cctype::table_size];
    std::memcpy(rc, const_rc, cctype::table_size * sizeof(cctype::mask));

    rc[','] = std::ctype_base::space; 
    return &rc[0];
  }
};

contact_info get_info(std::string line){
  contact_info info;
  std::stringstream sline(line);
  std::istream_iterator<std::string> begin(sline);
  std::istream_iterator<std::string> end;
  std::vector<std::string> stoken(begin, end);
  std::vector< std::vector<std::string> > all;
  for(std::vector<std::string>::iterator tk = stoken.begin(); tk != stoken.end(); ++tk){
    std::stringstream ss(*tk);
    ss.imbue(std::locale(std::locale(), new tokens()));
    std::istream_iterator<std::string> begin(ss);
    std::istream_iterator<std::string> end;
    std::vector<std::string> vstrings(begin, end);
    all.push_back(vstrings);
  }
  unsigned int nchains = all[0].size();
  for(unsigned int ch = 0; ch != nchains; ++ch ){
    std::pair<int, int> val = std::make_pair( atoi(all[1][ch].c_str()), atoi(all[2][ch].c_str()) );
    info[ all[0][ch][0] ] = val;
  }
  return info;
}

int main(int argc, char* argv[]){
  char* contactfile = 'x';
  double threshold = 0.;
  bool other_flag = false;
  bool too_much_flags = false;
  bool file_flag = false;
  bool eigenvalues_flag = false;
  bool eigenvectors_flag = false;
  bool mobility_flag = false;
  bool entropy_flag = false;
  bool fuzzy_flag = false;
  
  int c;
	while ((c = getopt (argc, argv, "f:t:vVme")) != -1){
		switch (c){
			case 'f':
				file_flag = true;
				contactfile = optarg;
				break;
			case 't':
				threshold = atof(optarg);
				break;
			case 'F':
				fuzzy_flag = true;
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
		}
	}

  if( not file_flag or not other_flag or too_much_flags ){
    std::cout << "usage: .gaussian_network_model -f <contact_file.dat> -t <threshold> - F [options (only one!)]" << std::endl
    					<< "\t-F: fuzzy; spring constant is proportional to entry in contact file (default=false)" << std::endl
    					<< "\t-t: threshold; consider only pair with entry in contact matrix greater than the threshold (default=0.)" << std::endl
              << "Options: " << std::endl 
              << "\t-v:  print eigenvalues" << std::endl
              << "\t-V:  print eigenvectors" << std::endl
              << "\t-e:  print entropy" << std::endl
              << "\t-m:  print mobility" << std::endl;
    return 1;
  }

  std::string line;
  std::ifstream infile(contactfile);
  /* read header */
  std::getline(infile, line); // title
  std::getline(infile, line); // x_info
  contact_info x_info_beg_end = get_info(line);
  std::getline(infile, line); // y_info
  contact_info y_info_beg_end = get_info(line);
  /* process header info */
  contact_info x_info_size_offset, y_info_size_offset;
  for(contact_info::iterator it = x_info_beg_end.begin(); it != x_info_beg_end.end(); ++it){
    if(it==x_info_beg_end.begin()){
       std::pair<int, int> val = std::make_pair( it->second.second - it->second.first + 1, 0 );
       x_info_size_offset[it->first] = val;
    }
    else{
      contact_info::iterator prev = it; --prev;
      int offset = x_info_size_offset[ prev->first ].first + x_info_size_offset[ prev->first ].second;
      x_info_size_offset[ it->first ] = std::make_pair( it->second.second - it->second.first + 1, offset );
    }
  }
  for(contact_info::iterator it = y_info_beg_end.begin(); it != y_info_beg_end.end(); ++it){
    if(it==y_info_beg_end.begin()){
       std::pair<int, int> val = std::make_pair( it->second.second - it->second.first + 1, 0 );
       y_info_size_offset[it->first] = val;
    }
    else{
      contact_info::iterator prev = it; --prev;
      int offset = y_info_size_offset[ prev->first ].first + y_info_size_offset[ prev->first ].second;
      y_info_size_offset[ it->first ] = std::make_pair( it->second.second - it->second.first + 1, offset );
    }
  }
  int x_size = x_info_size_offset.rbegin()->second.first + x_info_size_offset.rbegin()->second.second;
  int y_size = y_info_size_offset.rbegin()->second.first + y_info_size_offset.rbegin()->second.second;
  // fill the matrix
  Eigen::MatrixXd F(x_size, y_size);
  while( std::getline(infile, line) ){
    std::istringstream iss(line);
    char x_ch, y_ch;
    int x_am, y_am;
    double value;
    if (!(iss >> x_ch >> x_am >> y_ch >> y_am >> value ) ){
      break; 
    }
    if( x_info_beg_end.find(x_ch)!=x_info_beg_end.end() and y_info_beg_end.find(y_ch)!=y_info_beg_end.end() ){
      if( x_am >= x_info_beg_end[x_ch].first and x_am <= x_info_beg_end[x_ch].second and
          y_am >= y_info_beg_end[y_ch].first and y_am <= y_info_beg_end[y_ch].second ){
        int x = x_am - x_info_beg_end[x_ch].first + x_info_size_offset[x_ch].second;
        int y = y_am - y_info_beg_end[y_ch].first + y_info_size_offset[y_ch].second;
        if( not fuzzy_flag ){
        	if( value > threshold ){
        		value = 1.;
        	}
        	else{
        		value = 0.;
        	}
        }
        F(x,y) = -value;
      }
    }
  }
  for(int i = 0; i < x_size; ++i){
    F(i,i) = 0;
    for(int j = 0; j < y_size; ++j){
      if(j!=i){
        F(i,i) -= F(i,j);
      }
    }
  }
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(F);
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
    Eigen::ArrayXd mobility = Eigen::VectorXd::Zero(x_size).array();
    for(int c = 1; c!= y_size; ++c){
      Eigen::ArrayXd mode = eigensolver.eigenvectors().col(c).array()*eigensolver.eigenvectors().col(c).array();
      mobility += mode/eigensolver.eigenvalues()[c];
    }
    std::cout << mobility << std::endl;
  }
  else if(entropy_flag){
  	Eigen::ArrayXd entropy = Eigen::VectorXd::Zero(x_size).array();
    for(int c = 1; c!= y_size; ++c){
      Eigen::ArrayXd mode = eigensolver.eigenvectors().col(c).array()*eigensolver.eigenvectors().col(c).array();
      entropy -= mode*log( eigensolver.eigenvalues()[c] );
    }
    std::cout << entropy << std::endl;
  }
  return 0;
}

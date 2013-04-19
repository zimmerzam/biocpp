#include <map>
#include <iterator>
#include <sstream>
#include "../src/BioCpp_default.h"

typedef std::pair< char,int > residue;
typedef std::pair< residue, residue > contact;
typedef std::map< contact, double > contact_map_data;

struct contact_map{
  contact_map_data contacts;
  std::string title;
  std::set<char> xchainId, ychainId;
  std::map<char,int> xfirst, xlast;
  std::map<char,int> yfirst, ylast;
  bool use_hydrogens;
  double contact_distance;
  
  typedef std::map<double, int> contact_probability;
  
  void operator()( BioCpp::standard::residue& res1, BioCpp::standard::residue& res2 ){
    for( BioCpp::standard::residue::iterator r1 = res1.begin(); r1 !=res1.end(); ++r1 ){
      if( r1==res1.begin() ){
        xchainId.insert(r1->chainId);
        if( xfirst.find(r1->chainId)==xfirst.end() ){
          xfirst[r1->chainId] = 99999;
        }
        if( xlast.find(r1->chainId)==xlast.end() ){
          xlast[r1->chainId] = -99999;
        }
        if( r1->resSeq < xfirst[r1->chainId] ){
          xfirst[r1->chainId] = r1->resSeq;
        }
        else if( r1->resSeq > xlast[r1->chainId] ){
          xlast[r1->chainId] = r1->resSeq;
        }
      }
      for( BioCpp::standard::residue::iterator r2 = res2.begin(); r2 !=res2.end(); ++r2 ){
        if( r2==res2.begin() ){
          ychainId.insert(r2->chainId);
          if( yfirst.find(r2->chainId)==yfirst.end() ){
            yfirst[r2->chainId] = 99999;
          }
          if( ylast.find(r2->chainId)==ylast.end() ){
            ylast[r2->chainId] = -99999;
          }
          if( r2->resSeq < yfirst[r2->chainId] ){
            yfirst[r2->chainId] = r2->resSeq;
          }
          else if( r2->resSeq > ylast[r2->chainId] ){
            ylast[r2->chainId] = r2->resSeq;
          }
        }
        if( not use_hydrogens){
          if( r1->element==BioCpp::element::H or r2->element==BioCpp::element::H ){
            continue;
          }
        }
        residue resi1 = std::make_pair( r1->chainId, r1->resSeq );
        residue resi2 = std::make_pair( r2->chainId, r2->resSeq );
        contact cnt = std::make_pair( resi1, resi2 );
        if( contacts.find(cnt)==contacts.end() ){
          contacts[cnt]=0.;
        }
        if(&res1==&res2){
          return;
        }
        double distance = (r1->coordinate-r2->coordinate).norm();
        if( distance < contact_distance ){          
          contacts[cnt]+=1;
          return;
        }
      }
    }
  }
  void normalize(double fact){
    for( contact_map_data::iterator it = contacts.begin(); it != contacts.end(); ++it ){
      it->second/=fact;
    }
  }
  void chain_average(){
    bool failed = false;
    int first = xfirst.begin()->second;
    int last =  xlast.begin()->second;
    if( xchainId.size()==ychainId.size() ){
      for(std::set<char>::iterator i = xchainId.begin(), j=ychainId.begin(); i != xchainId.end(); ++i, ++j){
        if( (*i)!=(*j) ){
          failed=true;
        }
      }
    }
    else{
      failed=true;
    }
    for(std::map<char,int>::iterator i = xfirst.begin(); i != xfirst.end(); ++i){
      if(i->second!=first)
        failed=true;
    }
    for(std::map<char,int>::iterator i = yfirst.begin(); i != yfirst.end(); ++i){
      if(i->second!=first)
        failed=true;
    }
    for(std::map<char,int>::iterator i = xlast.begin(); i != xlast.end(); ++i){
      if(i->second!=last)
        failed=true;
    }
    for(std::map<char,int>::iterator i = ylast.begin(); i != ylast.end(); ++i){
      if(i->second!=last)
        failed=true;
    }
    if(failed){
      std::cout << "You can average each chain only if they have the same name, start and end" << std::endl
                << "The matrix has not been averaged" << std::endl;
      return;
    }
    
    for(int ri = first; ri <= last; ++ri){
      for(int rj = first; rj <= last; ++rj){
        for(std::set<char>::iterator j = xchainId.begin(); j != xchainId.end(); ++j){
          std::set<char>::iterator ti = xchainId.begin(), tj = j;
          int steps = std::distance( j, xchainId.end() );
          double val = 0.;
          for( ; tj != xchainId.end(); ++tj, ++ti ){
          residue resi = std::make_pair( *ti, ri );
              residue resj = std::make_pair( *tj, rj );
            val += contacts[ std::make_pair(resi, resj) ];
          }
          val /= double(steps);
          ti = xchainId.begin(); tj = j;
          for( ; tj != xchainId.end(); ++tj, ++ti ){
            residue resi = std::make_pair( *ti, ri );
            residue resj = std::make_pair( *tj, rj );
            contacts[ std::make_pair(resi, resj) ] = val;
            contacts[ std::make_pair(resj, resi) ] = val;
          }
        }
      }
    }
  }
  
  contact_probability probability_distribution(){
  	contact_probability prob;
  	for( contact_map_data::iterator it = contacts.begin(); it != contacts.end(); ++it ){
  		++prob[it->second];
  	}
  	return prob;
  }
  
};

std::ostream& operator << (std::ostream& out, contact_map map){
  out << map.title << std::endl;
  for( std::set<char>::iterator ch = map.xchainId.begin(); ch != map.xchainId.end(); ++ch ){
    if( ch!=map.xchainId.begin() ){
      out<< ",";
    }
    out << *ch;
  }
  out << "  ";
  for( std::set<char>::iterator ch = map.xchainId.begin(); ch != map.xchainId.end(); ++ch ){
    if( ch!=map.xchainId.begin() ){
      out<< ",";
    }
    out << map.xfirst[*ch];
  }
  out << "  ";
  for( std::set<char>::iterator ch = map.xchainId.begin(); ch != map.xchainId.end(); ++ch ){
    if( ch!=map.xchainId.begin() ){
      out<< ",";
    }
    out << map.xlast[*ch];
  }
  out << std::endl;
  
  for( std::set<char>::iterator ch = map.ychainId.begin(); ch != map.ychainId.end(); ++ch ){
    if( ch!=map.ychainId.begin() ){
      out<< ",";
    }
    out << *ch;
  }
  out << "  ";
  for( std::set<char>::iterator ch = map.ychainId.begin(); ch != map.ychainId.end(); ++ch ){
    if( ch!=map.ychainId.begin() ){
      out<< ",";
    }
    out << map.yfirst[*ch];
  }
  out << "  ";
  for( std::set<char>::iterator ch = map.ychainId.begin(); ch != map.ychainId.end(); ++ch ){
    if( ch!=map.ychainId.begin() ){
      out<< ",";
    }
    out << map.ylast[*ch];
  }
  out << std::endl;
  for( contact_map_data::iterator it = map.contacts.begin(); it != map.contacts.end(); ++it ){
    out << it->first.first.first << "  " << it->first.first.second << "  ";
    out << it->first.second.first << "  " << it->first.second.second << "  ";
    out << it->second;
    out << std::endl;
  }
  return out;
}

int main(int argc, char* argv[]){
  bool file_flag = false;
  bool average_flag = false;
  bool distribution_flag = false;
  bool map_flag = false;
  const char* contactfile;
  std::stringstream sstitle;
  
  contact_map map;
  map.use_hydrogens = false;
  map.contact_distance = 4.5;
  map.title = "Contact Probability";

  int c;
	while ((c = getopt (argc, argv, "f:t:d:hacp")) != -1){
		switch (c){
			case 'f':
			  file_flag = true;
				contactfile = optarg;
				break;
      case 't':
        sstitle << optarg;
        sstitle >> map.title;
        break;
      case 'd':
        map.contact_distance = atof(optarg);
        break;
      case 'h':
        map.use_hydrogens=true;
        break;
      case 'a':
        average_flag = true;
        break;
      case 'p':
        distribution_flag = true;
        break;
      case 'c':
      	map_flag = true;
      	break;
      
		}
	}

  if( not file_flag or not map_flag^distribution_flag ){
    std::cout << "usage: .contact_map -f 'structure.pdb' -d distance -h [-c -t title OR -p]" << std::endl
              << "Options: " << std::endl 
              << "\t-c:  print contact map (default false)" << std::endl
              << "\t-p:  print contact probability distribution (default false)" << std::endl
              << "\t-h:  use hydrogens (default false)" << std::endl
              << "\t-d:  set contact distance in Angstrom (default: 4.5)" << std::endl
              << "\t-t:  set map title (default: 'Contact Probability')" << std::endl;
    return 1;
  }

  BioCpp::pdb PDB(contactfile, 0);
  for(int mdl = 1; mdl <= PDB.n_models; ++mdl){
    BioCpp::pdb_model all_info = PDB.getModel(mdl);
    BioCpp::standard::complex cmp( all_info, PDB.RseqRes );
    BioCpp::Iterate<BioCpp::standard::residue, BioCpp::standard::residue>(cmp,cmp,map);
  }
  if(average_flag){
    map.chain_average();
  }
  map.normalize(PDB.n_models);
  if(map_flag){
	  std::cout << map;
	}
	else if(distribution_flag){
		contact_map::contact_probability prob = map.probability_distribution();
		for( contact_map::contact_probability::iterator it = prob.begin(); it != prob.end(); ++it ){
			std::cout << it->first << "  " << it->second << std::endl;
		}
	}
  return 0;
}

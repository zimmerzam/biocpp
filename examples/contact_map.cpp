#include <map>
#include <iterator>
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
  const char* filename = argv[1];
  
  contact_map map;
  map.use_hydrogens = bool(atoi(argv[2]));
  map.contact_distance = atof(argv[3]);
  map.title = "Contact Probability";
  
  BioCpp::pdb PDB(filename, 0);
  for(int mdl = 1; mdl <= PDB.n_models; ++mdl){
    BioCpp::pdb_model all_info = PDB.getModel(mdl);
    BioCpp::standard::complex cmp( all_info, PDB.RseqRes );
    BioCpp::Iterate<BioCpp::standard::residue, BioCpp::standard::residue>(cmp,cmp,map);
  }
  map.chain_average();
  map.normalize(PDB.n_models);
  std::cout << map;
  return 0;
}

#ifndef DICTIONARY_H
#define DICTIONARY_H

#include <map>
#include <list>
#include <string>
#include <libconfig.h++>
#include <iostream>

namespace BioCpp{

template <typename _definition_t>
class dictionary{
  public:
  	typedef _definition_t definition_t;
    std::map<int, std::string> id_to_string;
    std::map<std::string, int> string_to_id;
    std::map<int, definition_t> definition;
    
    dictionary( std::map<int, std::string> i2s = std::map<int, std::string>({}), std::map<std::string, int> s2i = std::map<std::string, int>({}), 
                      std::map<int, definition_t> des = std::map<int, definition_t>({}) );
    
    void importSetting(libconfig::Setting& setting);            
    void importSetting(libconfig::Setting& setting, std::list<std::string> import_lib);
    void writeSetting(std::string filename, std::string configname){};
};

template <typename _definition_t>
dictionary<_definition_t>::dictionary( std::map<int, std::string> i2s, 
                                      std::map<std::string, int> s2i,
                                      std::map<int, definition_t> des 
                                    ):
                                      id_to_string(i2s), 
                                      string_to_id(s2i), 
                                      definition(des) {};

template <typename _definition_t>
void dictionary<_definition_t>::importSetting(libconfig::Setting& setting){
  try{
    libconfig::Setting& import = setting["import"];
    int import_size = import.getLength();
    std::list<std::string> import_lib;
    for(int i = 0; i != import_size; ++i){
      import_lib.push_back( import[i] );
    }
    importSetting( setting, import_lib );
  }
  catch(...){
    std::cout << "Cannot find a list of libraries to be imported. Nothing has been done." << std::endl;
    return;
  }
};

template <typename definition_t>
void dictionary<definition_t>::importSetting(libconfig::Setting& setting, std::list<std::string> import_lib ){
  for( std::list<std::string>::iterator lib = import_lib.begin(); lib != import_lib.end(); ++lib ){
    try{
      libconfig::Setting& library = setting[*lib];
      int lib_size = library.getLength();
      for(int k = 0; k != lib_size; ++k){
        try{
          libconfig::Setting& item = library[k];
          int id = int(item["id"]);
          libconfig::Setting& strings = item["string"];
          int n_strings = strings.getLength();
          for(int i = 0; i != n_strings; ++i){
          	std::string tmp_string = strings[i];
            string_to_id[ tmp_string ] = id;
            if(i==0){
              id_to_string[id] = tmp_string;
            }
          }
          libconfig::Setting& definitions = item["definition"];
          if( definition.find(id)==definition.end() ){
            definition[id] = definition_t();
          }
          definition[id].importSetting(definitions);
        }
        catch(...){
          continue;
        };
      }
    }
    catch(...){
      continue;
    }
  }
};

}
#endif

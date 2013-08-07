/* ************************************************************************** */
/*                                                                            */
/*    Copyright 2013 Stefano Zamuner                                          */
/*                                                                            */
/*    This file is part of BioCpp.                                            */
/*                                                                            */
/*    BioCpp is free software: you can redistribute it and/or modify          */
/*    it under the terms of the GNU General Public License as published by    */
/*    the Free Software Foundation, either version 3 of the License, or       */
/*    (at your option) any later version.                                     */
/*                                                                            */
/*    BioCpp is distributed in the hope that it will be useful,               */
/*    but WITHOUT ANY WARRANTY; without even the implied warranty of          */
/*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           */
/*    GNU General Public License for more details.                            */
/*                                                                            */
/*    You should have received a copy of the GNU General Public License       */
/*    along with BioCpp.  If not, see <http://www.gnu.org/licenses/>.         */
/*                                                                            */
/* ************************************************************************** */

#ifndef BIOCPP_BASE_DICTIONARY_H
#define BIOCPP_BASE_DICTIONARY_H

#include <map>
#include <list>
#include <string>
#include <boost/algorithm/string.hpp>

namespace BioCpp{

class base_dictionary{
  protected:
    typedef std::map< int, std::string > 2str;
    std::map< std::string, int > data_name_to_id;
    std::map< int, 2str > data_to_string;
    std::map< std::string, int > data_to_id;
    
    // TODO check if word.size==length
    void readStringsLine( std::string& line, int id ){
      if( line.substr(0,7)!="STRINGS" ){
        return;
      }
      int length = atoi( line[8] );
      if( data_to_string[length].find(id)!=data_to_string[length].end() ){
        return;
      }
      std::list<std::string> words;
      boost::split(words, line.substr(10), boost::is_any_of(" \t"), boost::token_compress_on);
      data_to_string[length][id] = words.front().substr(1,length-2);
      for( std::list<std::string>::iterator w = words.begin(); w!=words.end(); ++w; ){
        data_to_id[w->substr(1,length-2) ] = id;
      }
    }
    
  public:
    virtual void readFromFile( std::string& filename ) = 0;
    virtual void read FromFile( std::list<std::string> filenames ){
      for( std::list<std::string>::iterator f = filenames.begin(); f!=filenames.end(); ++f ){
        readFromFile( *f );
      }
    }
  
    int to_id( std::string str ){
      return data_to_id.find(str)!=data_to_id.end() ? data_to_id[str] : -1;
    };
    std::string to_string( int length, int id ){
      if( data_to_string.find(length)==data_to_string.end() ){
        return std::string(length, ' ');
      }
      if( data_to_string[length].find(id)==data_to_string[length].end() ){
        return std::string(length, ' ');
      }
      return data_to_string[length][id];
    };
    int get_id_from_name(std::string name){
      return data_name_to_id.find(name)!=data_name_to_id.end() ? data_name_to_id[name] : -1;
    }
};


}
#endif

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

#ifndef BIOCPP_ATOM_DICTIONARY_H
#define BIOCPP_ATOM_DICTIONARY_H

namespace BioCpp{
namespace residue{

class dictionary : public base_dictionary {
  protected:
    std::map< int, std::list< std::vector<int> > > children;
    std::map< int, std::list< std::pair<int,int> > > connections;
  public:
    void readFromFile( std::string filename ){
    
    };
    
    std::pair<int, std::map< int, int > > get_children( std::vector<int>& child );
};

}
}
#endif

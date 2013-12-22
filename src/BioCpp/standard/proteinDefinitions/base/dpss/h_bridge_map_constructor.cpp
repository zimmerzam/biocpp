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

#include "h_bridge_map_constructor.hpp"
#include <BioCpp/standard/chemicalComponentDictionary/atom_dictionary_standard.hpp>
#include <BioCpp/proteins/dpss/h_bridge_energy.hpp>

namespace BioCpp{
namespace dpss{

template<>
class base_h_bridge_map_constructor<BioCpp::standard::base::chain::iterator,BioCpp::standard::base::complex>{
  public:
    BioCpp::standard::base::h_bridge_map operator()( BioCpp::standard::base::complex& cmp ){
      
      BioCpp::standard::base::h_bridge_map data;
      
      for(BioCpp::standard::base::complex::iterator ch1 = cmp.begin(); ch1 != cmp.end(); ++ch1){
        for(BioCpp::standard::base::chain::iterator res1 = ch1->begin(); res1 != ch1->end(); ++res1){
          
          Eigen::Vector3d c, o;
          if( not res1->exists(atom::C) ){
            continue;
          }
          if(not res1->exists(atom::O) ){
            continue;
          }
          c = (*res1)[atom::C].coordinate;
          o = (*res1)[atom::O].coordinate;
          for(BioCpp::standard::base::complex::iterator ch2 = cmp.begin(); ch2 != cmp.end(); ++ch2){
            for(BioCpp::standard::base::chain::iterator res2 = ch2->begin(); res2 != ch2->end(); ++res2){
              if( res2!=ch2->begin() ){
                Eigen::Vector3d n, h;
                if(not res2->exists(atom::N) ){
                  continue;
                }
                if(not res2->exists(atom::H) ){
                  continue;
                }
                n = (*res2)[atom::N].coordinate;
                h = (*res2)[atom::H].coordinate;
                double en = h_bridge_energy(c, o, n, h);
                if(en<-1){
                  data[std::make_pair(res1, res2)] = true;
                }
              }
            }
          }
        }
      }
      return data;
    }                        
};

}
}

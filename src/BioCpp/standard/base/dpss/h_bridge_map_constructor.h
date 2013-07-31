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

#ifndef BIOCPP_STANDARD_BASE_DPSS_H_BRIDGE_MAP
#define BIOCPP_STANDARD_BASE_DPSS_H_BRIDGE_MAP

#include <BioCpp/standard/base/chain.h>
#include <BioCpp/polimers/proteins/dpss/base_h_bridge_map_constructor.h>
#include <BioCpp/polimers/proteins/dpss/h_bridge_energy.h>

namespace BioCpp{
namespace standard{
namespace base{
namespace dpss{

typedef BioCpp::dpss::base_h_bridge_map_constructor<BioCpp::standard::base::chain::iterator, 
                                                    BioCpp::standard::base::complex>
  h_bridge_map_constructor;

}
}
}

namespace dpss{

template<>
class base_h_bridge_map_constructor<BioCpp::standard::base::chain::iterator,BioCpp::standard::base::complex>{
  public:
    BioCpp::standard::base::dpss::h_bridge_map operator()( BioCpp::standard::base::complex& cmp ){
      
      BioCpp::standard::base::dpss::h_bridge_map data;
      
      for(BioCpp::standard::base::complex::iterator ch1 = cmp.begin(); ch1 != cmp.end(); ++ch1){
        for(BioCpp::standard::base::chain::iterator res1 = ch1->begin(); res1 != ch1->end(); ++res1){
          
          Eigen::Vector3d c, o;
          if( not res1->exists(atom::C_) ){
            continue;
          }
          if(not res1->exists(atom::O_) ){
            continue;
          }
          c = (*res1)[atom::C_].coordinate;
          o = (*res1)[atom::O_].coordinate;
          for(BioCpp::standard::base::complex::iterator ch2 = cmp.begin(); ch2 != cmp.end(); ++ch2){
            for(BioCpp::standard::base::chain::iterator res2 = ch2->begin(); res2 != ch2->end(); ++res2){
              if( res2!=ch2->begin() ){
                Eigen::Vector3d n, h;
                if(not res2->exists(atom::N_) ){
                  continue;
                }
                if(not res2->exists(atom::H_) ){
                  continue;
                }
                n = (*res2)[atom::N_].coordinate;
                h = (*res2)[atom::H_].coordinate;
                double en = h_bridge_energy(c, o, n, h);
                if(en<-1){
                  data[std::make_pair(res1, res2)] = true;
//                  std::cout << res1->begin()->chainId << res1->begin()->resSeq << std::endl;
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
#endif

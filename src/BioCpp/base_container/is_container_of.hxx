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

#ifndef IS_CONTAINER_OF_H
#define IS_CONTAINER_OF_H

#include "base_container.hxx"

namespace BioCpp{

/*! \brief template struct to be used to determine if a type is contained in
    a container
    
    This is similar to the standard library std::is_base and std::is_same
    
    \code
      #include <iostream>
      #include "BioCpp/BioCpp.h"
            
      int main(){
        std::cout << "char in double " 
                  << BioCpp::is_container_of<double, char>::value 
                  << std::endl; // output: 0
        std::cout << "char in int " 
                  << BioCpp::is_container_of<int, char>::value 
                  << std::endl; // output: 0
        std::cout << "int in double " 
                  << BioCpp::is_container_of<int, double>::value 
                  << std::endl; // output: 0
        return 0;
      }
    \endcode
*/
template <typename, typename>
struct is_container_of{
	static bool const value = false; /*!< by default its value is false */
};

/*! \brief template struct to be used to determine if a type is contained in
    a container
    
    This is a specialization of 
    template <typename, typename> struct is_container_of for the case in which 
    a container is checked.
    Its value is true if the container or the containers in it contains the
    given type.
    
    \code
      #include <iostream>
      #include "BioCpp/BioCpp.h"
      
      typedef BioCpp::base::container< int, char, int > word;
      typedef BioCpp::base::container< int, word, int > sentence;
      
      int main(){
        std::cout << "char in word " 
                  << BioCpp::is_container_of<word, char>::value 
                  << std::endl; // output: 1
        std::cout << "char in sentence " 
                  << BioCpp::is_container_of<sentence, char>::value 
                  << std::endl; // output: 1
        std::cout << "double in word " 
                  << BioCpp::is_container_of<sentence, double>::value 
                  << std::endl; // output: 0
        return 0;
      }
    \endcode
    
*/
template < template <typename, typename, typename> class container, typename item, typename child_id, typename obj, typename id>
struct is_container_of< container<child_id, obj, id>, item >{
	static bool const value = std::is_same<obj, item>::value or is_container_of<obj,item>::value;
};


}
#endif

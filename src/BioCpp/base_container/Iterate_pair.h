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

#ifndef CONTAINER_ITERATE_PAIR_H
#define CONTAINER_ITERATE_PAIR_H

#include "base_container.h"
#include "is_container_of.h"

namespace BioCpp{

/*! \brief Iterate over a pair of container and apply a functor to each possible pair of (item1,item2) 
    
    This function efficiently iterate pair of items of the required type.
    
    \code{.cpp}
      #include <iostream>
      #include <utility>

      #include "../src/polimers/Iterate_single.h"
      #include "../src/polimers/Iterate_pair.h"

      // print functor. It can accept both char or iterator to char
      struct print{
        void operator()(char s){
          std::cout << s << std::endl;
        }
        void operator()(char s1, char s2){
          std::cout << s1 << "  " << s2 << std::endl;
        }
        template <typename T,
                  typename = typename std::enable_if< std::is_same<typename std::iterator_traits<T>::value_type, char>::value, char>::type >
        void operator()(T s){
          std::cout << *s << "  " << "iter" << std::endl;
        }
        template <typename T,
                  typename = typename std::enable_if< std::is_same<typename std::iterator_traits<T>::value_type, char>::value, char>::type >
        void operator()(T s1, T s2){
          std::cout << *s1 << "  " << *s2 << "  " << "iter" << std::endl;
        }
      };

      typedef BioCpp::base_container< int, char, int > word; // define a word
      typedef BioCpp::base_container< int, word, int > sentence; // define a sentence

      int main(){
        word s1; // build a word containing 'A', 'B' and 'C'
        s1.Reserve(3);
        s1.Append(1, 'A');
        s1.Append(2, 'B');
        s1.Append(3, 'C');
        word s2; // build a word containing 'F', 'G' and 'H'
        s2.Reserve(3);
        s2.Append(11, 'F');
        s2.Append(12, 'G');
        s2.Append(13, 'H');

        sentence b1; // build a sentence containing the two words defined above
        b1.Append(1,s1);
        b1.Append(2,s2);
        
        std::cout << "b1 size: " << b1.size() << std::endl; // output: 2
        std::cout << "s1 size: " << s1.size() << std::endl; // output: 3
        std::cout << "s2 size: " << s2.size() << std::endl; // output: 3
        
        std::cout << "s1[2] is " << s1[2] << std::endl; // output: 'B'
        std::cout << "s2[13] is " << s2[13] << std::endl; // output: 'H'
        
        std::cout << "s1 contains 3: " << s1.exists(1) << std::endl; // output: 1
        std::cout << "s1 contains 11: " << s1.exists(11) << std::endl; // output: 0
        
        std::cout << "******************" << std::endl;
        
        print p; // print functor
        BioCpp::Iterate<char,char>(b1, s1, p);  // print each possible pair of chars
       
        std::cout << "******************" << std::endl;
       
        BioCpp::Iterate<char,char>(b1, b1, p); // print each possible pair of chars
        
        return 0;
      }
    \endcode
*/
template< typename item1, typename item2,
          typename Function,
          typename container1, typename container2,
          typename is_cont1 = typename std::enable_if< is_container_of< container1, item1 >::value, item1 >::type,
          typename is_cont2 = typename std::enable_if< is_container_of< container2, item2 >::value, item2 >::type,
          typename is_parent1 = typename std::enable_if< std::is_same< typename container1::child_type, item1 >::value, item1 >::type,
          typename is_parent2 = typename std::enable_if< std::is_same< typename container2::child_type, item2 >::value, item2 >::type,
          typename not_parent1 = void,
          typename not_parent2 = void >
void Iterate( container1& cont1, container2& cont2, Function& todo );

template< typename item1, typename item2,
          typename Function,
          typename container1, typename container2,
          typename is_cont1 = typename std::enable_if< is_container_of< container1, item1 >::value, item1 >::type,
          typename is_cont2 = typename std::enable_if< is_container_of< container2, item2 >::value, item2 >::type,
          typename is_parent1 = typename std::enable_if< std::is_same< typename container1::child_type, item1 >::value, item1 >::type,
          typename is_parent2 = void,
          typename not_parent1 = void,
          typename not_parent2 = typename std::enable_if< not std::is_same< typename container2::child_type, item2 >::value, item2 >::type >
void Iterate( container1& cont1, container2& cont2, Function& todo, bool a = true );

template< typename item1, typename item2,
          typename Function,
          typename container1, typename container2,
          typename is_cont1 = typename std::enable_if< is_container_of< container1, item1 >::value, item1 >::type,
          typename is_cont2 = typename std::enable_if< is_container_of< container2, item2 >::value, item2 >::type,
          typename is_parent1 = void,
          typename is_parent2 = typename std::enable_if< std::is_same< typename container2::child_type, item2 >::value, item2 >::type,
          typename not_parent1 = typename std::enable_if< not std::is_same< typename container1::child_type, item1 >::value, item1 >::type,
          typename not_parent2 = void >
void Iterate( container1& cont1, container2& cont2, Function& todo, bool a = true, bool b = true );

template< typename item1, typename item2,
          typename Function,
          typename container1, typename container2,
          typename is_cont1 = typename std::enable_if< is_container_of< container1, item1 >::value, item1 >::type,
          typename is_cont2 = typename std::enable_if< is_container_of< container2, item2 >::value, item2 >::type,
          typename is_parent1 = void,
          typename is_parent2 = void,
          typename not_parent1 = typename std::enable_if< not std::is_same< typename container1::child_type, item1 >::value, item1 >::type,
          typename not_parent2 = typename std::enable_if< not std::is_same< typename container2::child_type, item2 >::value, item2 >::type >
void Iterate( container1& cont1, container2& cont2, Function& todo, bool a = true, bool b = true, bool c = true );


/*! \brief Iterate over a pair of container and apply a functor to each possible pair of (iterator to item1, iterator to item2) 
    
    This function efficiently iterate pair of items of the required type. The function is applied to iterators such that the value_type of the iterator
    is of the required type
    
    \code{.cpp}
      #include <iostream>
      #include <utility>

      #include "../src/polimers/Iterate_single.h"
      #include "../src/polimers/Iterate_pair.h"

      // print functor. It can accept both char or iterator to char
      struct print{
        void operator()(char s){
          std::cout << s << std::endl;
        }
        void operator()(char s1, char s2){
          std::cout << s1 << "  " << s2 << std::endl;
        }
        template <typename T,
                  typename = typename std::enable_if< std::is_same<typename std::iterator_traits<T>::value_type, char>::value, char>::type >
        void operator()(T s){
          std::cout << *s << "  " << "iter" << std::endl;
        }
        template <typename T,
                  typename = typename std::enable_if< std::is_same<typename std::iterator_traits<T>::value_type, char>::value, char>::type >
        void operator()(T s1, T s2){
          std::cout << *s1 << "  " << *s2 << "  " << "iter" << std::endl;
        }
      };

      typedef BioCpp::base_container< int, char, int > word; // define a word
      typedef BioCpp::base_container< int, word, int > sentence; // define a sentence

      int main(){
        word s1; // build a word containing 'A', 'B' and 'C'
        s1.Reserve(3);
        s1.Append(1, 'A');
        s1.Append(2, 'B');
        s1.Append(3, 'C');
        word s2; // build a word containing 'F', 'G' and 'H'
        s2.Reserve(3);
        s2.Append(11, 'F');
        s2.Append(12, 'G');
        s2.Append(13, 'H');

        sentence b1; // build a sentence containing the two words defined above
        b1.Append(1,s1);
        b1.Append(2,s2);
        
        std::cout << "b1 size: " << b1.size() << std::endl; // output: 2
        std::cout << "s1 size: " << s1.size() << std::endl; // output: 3
        std::cout << "s2 size: " << s2.size() << std::endl; // output: 3
        
        std::cout << "s1[2] is " << s1[2] << std::endl; // output: 'B'
        std::cout << "s2[13] is " << s2[13] << std::endl; // output: 'H'
        
        std::cout << "s1 contains 3: " << s1.exists(1) << std::endl; // output: 1
        std::cout << "s1 contains 11: " << s1.exists(11) << std::endl; // output: 0
        
        std::cout << "******************" << std::endl;
        
        print p; // print functor
        BioCpp::Iterate_iter<char,char>(b1, s1, p);  // print each possible pair of chars
       
        std::cout << "******************" << std::endl;
       
        BioCpp::Iterate_iter<char,char>(b1, b1, p); // print each possible pair of chars
        
        return 0;
      }
    \endcode
*/
template< typename item1, typename item2,
          typename Function,
          typename container1, typename container2,
          typename is_cont1 = typename std::enable_if< is_container_of< container1, item1 >::value, item1 >::type,
          typename is_cont2 = typename std::enable_if< is_container_of< container2, item2 >::value, item2 >::type,
          typename is_parent1 = typename std::enable_if< std::is_same< typename container1::child_type, item1 >::value, item1 >::type,
          typename is_parent2 = typename std::enable_if< std::is_same< typename container2::child_type, item2 >::value, item2 >::type,
          typename not_parent1 = void,
          typename not_parent2 = void >
void Iterate_iter( container1& cont1, container2& cont2, Function& todo );

template< typename item1, typename item2,
          typename Function,
          typename container1, typename container2,
          typename is_cont1 = typename std::enable_if< is_container_of< container1, item1 >::value, item1 >::type,
          typename is_cont2 = typename std::enable_if< is_container_of< container2, item2 >::value, item2 >::type,
          typename is_parent1 = typename std::enable_if< std::is_same< typename container1::child_type, item1 >::value, item1 >::type,
          typename is_parent2 = void,
          typename not_parent1 = void,
          typename not_parent2 = typename std::enable_if< not std::is_same< typename container2::child_type, item2 >::value, item2 >::type >
void Iterate_iter( container1& cont1, container2& cont2, Function& todo, bool a = true );

template< typename item1, typename item2,
          typename Function,
          typename container1, typename container2,
          typename is_cont1 = typename std::enable_if< is_container_of< container1, item1 >::value, item1 >::type,
          typename is_cont2 = typename std::enable_if< is_container_of< container2, item2 >::value, item2 >::type,
          typename is_parent1 = void,
          typename is_parent2 = typename std::enable_if< std::is_same< typename container2::child_type, item2 >::value, item2 >::type,
          typename not_parent1 = typename std::enable_if< not std::is_same< typename container1::child_type, item1 >::value, item1 >::type,
          typename not_parent2 = void >
void Iterate_iter( container1& cont1, container2& cont2, Function& todo, bool a = true, bool b = true );

template< typename item1, typename item2,
          typename Function,
          typename container1, typename container2,
          typename is_cont1 = typename std::enable_if< is_container_of< container1, item1 >::value, item1 >::type,
          typename is_cont2 = typename std::enable_if< is_container_of< container2, item2 >::value, item2 >::type,
          typename is_parent1 = void,
          typename is_parent2 = void,
          typename not_parent1 = typename std::enable_if< not std::is_same< typename container1::child_type, item1 >::value, item1 >::type,
          typename not_parent2 = typename std::enable_if< not std::is_same< typename container2::child_type, item2 >::value, item2 >::type >
void Iterate_iter( container1& cont1, container2& cont2, Function& todo, bool a = true, bool b = true, bool c = true );

}// end namespace


template< typename item1, typename item2,
          typename Function,
          typename container1, typename container2,
          typename is_cont1,
          typename is_cont2,
          typename is_parent1,
          typename is_parent2,
          typename not_parent1,
          typename not_parent2 >
void BioCpp::Iterate( container1& cont1, container2& cont2, Function& todo ){
  for( typename container1::iterator it1 = cont1.begin(); it1!=cont1.end(); ++it1 ){
    for( typename container2::iterator it2 = cont2.begin(); it2!=cont2.end(); ++it2 ){
      todo( *it1, *it2 );
    }
  }
}

template< typename item1, typename item2,
          typename Function,
          typename container1, typename container2,
          typename is_cont1,
          typename is_cont2,
          typename is_parent1,
          typename is_parent2,
          typename not_parent1,
          typename not_parent2 >
void BioCpp::Iterate( container1& cont1, container2& cont2, Function& todo, bool a = true ){
  for( typename container2::iterator it2 = cont2.begin(); it2!=cont2.end(); ++it2 ){
    BioCpp::Iterate<item1, item2>( cont1, *it2, todo );
  }
}

template< typename item1, typename item2,
          typename Function,
          typename container1, typename container2,
          typename is_cont1,
          typename is_cont2,
          typename is_parent1,
          typename is_parent2,
          typename not_parent1,
          typename not_parent2 >
void BioCpp::Iterate( container1& cont1, container2& cont2, Function& todo, bool a = true, bool b = true ){
  for( typename container1::iterator it1 = cont1.begin(); it1!=cont1.end(); ++it1 ){
    BioCpp::Iterate<item1, item2>( *it1, cont2, todo );
  }
}

template< typename item1, typename item2,
          typename Function,
          typename container1, typename container2,
          typename is_cont1,
          typename is_cont2,
          typename is_parent1,
          typename is_parent2,
          typename not_parent1,
          typename not_parent2 >
void BioCpp::Iterate( container1& cont1, container2& cont2, Function& todo, bool a = true, bool b = true, bool c = true ){
  for( typename container1::iterator it1 = cont1.begin(); it1!=cont1.end(); ++it1 ){
    for( typename container2::iterator it2 = cont2.begin(); it2!=cont2.end(); ++it2 ){
      BioCpp::Iterate<item1, item2>( *it1, *it2, todo );
    }
  }
}


template< typename item1, typename item2,
          typename Function,
          typename container1, typename container2,
          typename is_cont1,
          typename is_cont2,
          typename is_parent1,
          typename is_parent2,
          typename not_parent1,
          typename not_parent2 >
void BioCpp::Iterate_iter( container1& cont1, container2& cont2, Function& todo ){
  for( typename container1::iterator it1 = cont1.begin(); it1!=cont1.end(); ++it1 ){
    for( typename container2::iterator it2 = cont2.begin(); it2!=cont2.end(); ++it2 ){
      todo( it1, it2 );
    }
  }
}

template< typename item1, typename item2,
          typename Function,
          typename container1, typename container2,
          typename is_cont1,
          typename is_cont2,
          typename is_parent1,
          typename is_parent2,
          typename not_parent1,
          typename not_parent2 >
void BioCpp::Iterate_iter( container1& cont1, container2& cont2, Function& todo, bool a = true ){
  for( typename container2::iterator it2 = cont2.begin(); it2!=cont2.end(); ++it2 ){
    BioCpp::Iterate_iter<item1, item2>( cont1, *it2, todo );
  }
}

template< typename item1, typename item2,
          typename Function,
          typename container1, typename container2,
          typename is_cont1,
          typename is_cont2,
          typename is_parent1,
          typename is_parent2,
          typename not_parent1,
          typename not_parent2 >
void BioCpp::Iterate_iter( container1& cont1, container2& cont2, Function& todo, bool a = true, bool b = true ){
  for( typename container1::iterator it1 = cont1.begin(); it1!=cont1.end(); ++it1 ){
    BioCpp::Iterate_iter<item1, item2>( *it1, cont2, todo );
  }
}

template< typename item1, typename item2,
          typename Function,
          typename container1, typename container2,
          typename is_cont1,
          typename is_cont2,
          typename is_parent1,
          typename is_parent2,
          typename not_parent1,
          typename not_parent2 >
void BioCpp::Iterate_iter( container1& cont1, container2& cont2, Function& todo, bool a = true, bool b = true, bool c = true ){
  for( typename container1::iterator it1 = cont1.begin(); it1!=cont1.end(); ++it1 ){
    for( typename container2::iterator it2 = cont2.begin(); it2!=cont2.end(); ++it2 ){
      BioCpp::Iterate_iter<item1, item2>( *it1, *it2, todo );
    }
  }
}

#endif

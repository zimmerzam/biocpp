/*!
    \file iterate.cpp
    \brief Example usage of iteration methods over base_container
    
    In this example we define two new types:
      - word (a container of chars)
      - sentence (a container of words)
    The example shows how to build the containers and how to appropriately 
    iterate over them.
*/

#include <iostream>
#include <utility>

#include <BioCpp/base_container/Iterate_single.h>
#include <BioCpp/base_container/Iterate_pair.h>

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
  BioCpp::Iterate<char>(b1, p);  // print each char in sentence b1
 
  std::cout << "******************" << std::endl;
 
  BioCpp::Iterate_iter<char,char>(b1, b1, p); // print each possible pair of chars (iterator overload)
  
  return 0;
}

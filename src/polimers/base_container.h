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

#ifndef BASE_CONTAINER_H
#define BASE_CONTAINER_H

#include <map>
#include <vector>

namespace BioCpp{

/*! \brief This class represent a generic container (residue, chain, complex, ...)
    
    This is a very flexible class that allows the creation of custom containers.
    @tparam childID the identifier type of its children
    @tparam childType the type of its children
    @tparam ID a type describing it.
    
    \example container.cpp
    Example usage of base_container methods
*/
template <typename childID, typename childType, typename ID=int>
class base_container {
  private:
    typedef typename std::vector<childType> child_list;
  protected:
    /*! \brief a list of pointers to child objects */
    std::map< childID, childType* > child; 
    /*! \brief child object list */
    child_list children;
  public:
    typedef childID child_id_type;
    typedef childType child_type;
    typedef ID id_type;
    
    /*! \brief Describes the container 'class' */
    ID type;
    
    /*! \brief A pointer to the container that contains this one */
    base_container* parent;

    /*! void constructor */
    base_container(){};

    /*! \brief constructor from a vector of children
        
        @param children a vector of `std::pair` containing child identifier and child object.
        \note This method is safer than base_container::Append, but is also slower. Use this if
        the number of children is large, otherwise use base_container::Append. Internally uses
        base_container::Reserve and base_container::Append
        \todo check if it works
    */
    base_container(std::vector< std::pair<childID, childType> >& children);
    
    /*! \brief General purpose constructor, to be specialized
        \todo write an example code
        \todo write a version with two info arguments
    */
    template <typename list, typename info>
    base_container(list& lst, info& i);
    
    /*! \brief General purpose constructor, to be specialized
        \todo write an example code
        \todo write a version with two info arguments
    */
    template <typename list, typename info>
    base_container(list& lst, info& i, info& j);
    
    typedef typename child_list::iterator iterator; /*!< iterator over children */
    typedef typename child_list::const_iterator const_iterator; /*!< const iterator over children */
    typedef typename child_list::reverse_iterator reverse_iterator; /*!< reverse iterator over children */
    
    iterator begin(){return children.begin();}; /*!< iterator to the first item */
    iterator end(){return children.end();}; /*!< this is not the last item */
    reverse_iterator rbegin(){return children.rbegin();}; /*!< iterator to the last item */
    reverse_iterator rend(){return children.rend();}; /*!< this is not the first item */

    /*! \ retrieve the size of the container
        \return the number of children 
    */
    unsigned int size(){return children.size();};

    /*! \brief Reserve the memory for at least 'n_child' children
        
        \param n_child the final number of children
        \note Always use this if possible.
        
        \code{.cpp}
          #include <iostream>
          #include "BioCpp/BioCpp.h"
          
          typedef BioCpp::base_container< int, char, int > word;
          
          int main(){
            word s1; // build a word containing 'A', 'B' and 'C'
            s1.Reserve(3); // You can now safely append three items to this container
            return 0;
          }
        \endcode
    */
    void Reserve(int n_child);
       
    /*! \brief Append a child 
        
        \param id child identifier
        \param ch child object
        
        \code{.cpp}
          #include <iostream>
          #include "BioCpp/BioCpp.h"
          
          typedef BioCpp::base_container< int, char, int > word;
          
          int main(){
            word s1; // build a word containing 'A', 'B' and 'C'
            s1.Reserve(3);
            s1.Append(1, 'A');
            s1.Append(2, 'B');
            s1.Append(3, 'C');
            s1.Append(4, 'D'); // Attention!! You have reserved memory for three items only!!
            return 0;
          }
        \endcode
    */
    void Append(childID id, childType ch);
    
    /*! \brief Check if a child with a given id exists
        
        \return false if there is no element with the specified id, true otherwise
        \param id child identifier

        \code{.cpp}
          #include <iostream>
          #include "BioCpp/BioCpp.h"
          
          typedef BioCpp::base_container< int, char, int > word;
          
          int main(){
            word s1; // build a word containing 'A', 'B' and 'C'
            s1.Reserve(3);
            s1.Append(1, 'A');
            s1.Append(2, 'B');
            s1.Append(3, 'C');
            std::cout << s1.exists(1) << std::endl; //output: 1
            std::cout << s1.exists(4) << std::endl; //output: 0
            return 0;
          }
        \endcode
    */
    bool exists(childID id);

    /*! \brief Get the child with a given identifier
        
        \return a reference to the child with the given `id`
        \param id the child identifier
        \note if the desired child has been added with base_container::Append
        the returned reference can be broken. If the total number of children is known in advance, 
        use base_container::Reserve method before appending any child, 
        or use base_container(std::vector< std::pair<childID, childType> >&)
        
        \code{.cpp}
          #include <iostream>
          #include "BioCpp/BioCpp.h"
          
          typedef BioCpp::base_container< int, char, int > word;
          
          int main(){
            word s1; // build a word containing 'A', 'B' and 'C'
            s1.Reserve(3);
            s1.Append(1, 'A');
            s1.Append(2, 'B');
            s1.Append(3, 'C');
            std::cout << s1[1] << "  " << s1[2] << "  " << s1[3] << std::endl; //output: 'A'  'B'  'C'
            return 0;
          }
        \endcode
    */
    childType& operator[](childID id);
};

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
template <typename childID, typename childType, typename ID>
inline base_container<childID, childType, ID>::base_container(std::vector< std::pair<childID, childType> >& children){
  Reserve(children.size());
  for( typename std::vector< std::pair<childID, childType> >::iterator ch = children.begin(); ch < children.end(); ++ch ){
    Append(ch->first, ch->second);
  }
}

/* operators */
template <typename childID, typename childType, typename ID>
inline void base_container<childID, childType, ID>::Reserve(int size){
  children.reserve(size);
}

template <typename childID, typename childType, typename ID>
inline void base_container<childID, childType, ID>::Append(childID id, childType ch){
  if(child.find(id)==child.end()){
    children.push_back( ch );
    child.insert(std::make_pair(id, &(*(rbegin())) ));
  }
}

template <typename childID, typename childType, typename ID>
inline bool base_container<childID, childType, ID>::exists(childID id){
  if(child.find(id)==child.end())
    return false;
  return true;
};

template <typename childID, typename childType, typename ID>
inline childType& base_container<childID, childType, ID>::operator[](childID id){
  return *(child[id]);
};

} // end namespace

#endif

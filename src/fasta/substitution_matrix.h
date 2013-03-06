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

#ifndef ALIGNMENT_MATRIX_H
#define ALIGNMENT_MATRIX_H

#include <map>

namespace BioCpp{
namespace fasta{

class substitution_matrix{
  typedef std::map< char, std::map< char , double > > align_map;
  private:
    align_map map;
  public:
    substitution_matrix( align_map map1 ){map=map1;};
    double operator()(char s1, char s2);
};

inline double substitution_matrix::operator()(char s1, char s2){
  return map[s1][s2];
}

}// end fasta
}// end biocpp

#endif

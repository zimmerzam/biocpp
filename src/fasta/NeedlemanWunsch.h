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

#ifndef NEEDLEMAN_WUNSCH_H
#define NEEDLEMAN_WUNSCH_H

#include <string>
#include "substitution_matrix.h"
#include "../geometry/Eigen/Core"
#include "../utils/errors_and_warnings.h"

namespace BioCpp{
namespace fasta{

double NeedlemanWunsch(std::string& S1, std::string& S2, substitution_matrix& matrix, BioCpp::error& errors, BioCpp::warning& warnings){

  // input string have to be uppercase
  std::transform(S1.begin(), S1.end(), S1.begin(), ::toupper);
  std::transform(S2.begin(), S2.end(), S2.begin(), ::toupper);

  // control if input strings are valid
  std::string s1 = S1, s2 = S2;
  s1.erase(std::remove_if(s1.begin(), s1.end(), [](char ch){if(ch=='-' or ch=='^') return true; return false;}), s1.end());
  s2.erase(std::remove_if(s2.begin(), s2.end(), [](char ch){if(ch=='-' or ch=='^') return true; return false;}), s2.end());
  std::replace_if( s1.begin(), s1.end(), 
       [&warnings](char ch){if(ch=='O' or ch=='U'){ warnings=(BioCpp::warning)(ALIGN_NOT_A_VALID_SEQUENCE|warnings); return true;} return false;}, 'X');
  std::replace_if( s2.begin(), s2.end(), 
       [&warnings](char ch){if(ch=='O' or ch=='U'){ warnings=(BioCpp::warning)(ALIGN_NOT_A_VALID_SEQUENCE|warnings); return true;} return false;}, 'X');
  
  unsigned int size1 = s1.size(), size2 = s2.size();
  
  // fill alignment matrix F
  Eigen::MatrixXd F = Eigen::MatrixXd::Zero(size1+1, size2+1);
  for(unsigned int t = 0; t < size1+1; ++t){
    F(t,0) = t*matrix(s1[t], 'X');
  }
  for(unsigned int r = 0; r < size2+1; ++r){
    F(0,r) = r*matrix('X', s2[r]);
  }
  for(unsigned int t = 0; t < size1; ++t){
    for(unsigned int r = 0; r < size2; ++r){
      double match = F(t,r) + matrix(s1[t], s2[r]);
      double deletion = F(t,r+1) + matrix(s1[t], 'X');
      double insertion = F(t+1,r) + matrix('X', s2[r]);
      F(t+1,r+1) = std::max( match, std::max(deletion, insertion) );
    }
  }
  // save score
  double score = F(size1, size2);
  std::cout << score << std::endl; 
  if(score < 0){
    errors=(BioCpp::error)(ALIGN_FAILED|errors);
    return score;
  }
  // backtracking
  std::string align1 = "", align2 = "";
  while( size1>0 or size2>0 ){
    if( size1>0 and size2>0 and F(size1,size2) == F(size1-1,size2-1) + matrix(s1[size1-1], s2[size2-1]) ){
      align1 = s1[size1-1] + align1;
      align2 = s2[size2-1] + align2;
      --size1;
      --size2;
    }
    else if (size1 > 0 and F(size1,size2) == F(size1-1,size2) + matrix(s1[size1-1], 'X') ){
      align1 = s1[size1-1] + align1;
      align2 = '^' + align2;
      --size1;
    }
    else if (size2 > 0 and F(size1,size2) == F(size1,size2-1) + matrix('X', s2[size2-1]) ){
      align1 = '^' + align1;
      align2 = s2[size2-1] + align2;
      --size2;
    }
  } 

  //update input sequences
  S1 = align1;
  S2 = align2;
  return score;
}

}// end biocpp
}// end fasta

#endif

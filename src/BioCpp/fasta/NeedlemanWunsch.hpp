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

#ifndef NEEDLEMAN_WUNSCH_HPP
#define NEEDLEMAN_WUNSCH_HPP

#include <string>
#include "substitution_matrix.hpp"
#include <Eigen/Core>
#include "../utils/errors_and_warnings/errors_and_warnings.hpp"

namespace BioCpp{
namespace fasta{

/*!
  This is an implementation of the Needlemann-Wunsch alignment algorithm.
  
  \param S1 the first sequence
  \param S2 the second sequence
  \param matrix the substitution matrix to be used during the alignment
  \param errors a BioCpp::error object. This is going to be modified
  \param warnings a BioCpp::warning object. This is going to be modified
  
  \note input sequences are transformed to uppercase!
*/
double NeedlemanWunsch(std::string& S1, std::string& S2, substitution_matrix& matrix, BioCpp::error& errors, BioCpp::warning& warnings);

}// end biocpp
}// end fasta

#endif

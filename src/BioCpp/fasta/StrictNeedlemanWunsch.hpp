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

#ifndef STRICT_NEEDLEMAN_WUNSCH_ALIGNMENT_HPP
#define STRICT_NEEDLEMAN_WUNSCH_ALIGNMENT_HPP

#include <string>
#include <algorithm>

#include "substitution_matrix.hpp"
#include "../utils/errors_and_warnings/errors_and_warnings.hpp"

namespace BioCpp{
namespace fasta{

/*!
    This is a strict version of the Needlemann-Wunsch alignment algorithm: here insertion and deletion
    can happen only in those positions where the input sequences exhibit a gap character -. In this case 
    the score associated is that of the current amino-acid with X.
*/
double StrictNeedlemanWunsch(std::string& S1, std::string& S2, substitution_matrix& matrix, BioCpp::error& errors, BioCpp::warning& warnings);

} // end fasta
} // end biocpp
#endif

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

#ifndef DPSS_ID_H
#define DPSS_ID_H

#include <map>
#include "base_h_bridge_map.h"

namespace BioCpp{
namespace dpss{

/*! \brief A list of DPSS secondary structure classification.

		For more informations read: 'Dictionary
		of Protein Secondary Structure: Pattern Recognition of Hydrogen-Bonded and
		Geometrical Features' by Wolfgang Kabsch and Christian Sander
    \note Only parallel bridge, antiparallel bridge and alpha helix is supported at the moment. */
enum id{ 
	PARA_BRIDGE, /*!< \brief a parallel bridge */
	ANTI_BRIDGE, /*!< \brief an antiparallel bridge */
	FOUR_HELIX,  /*!< \brief an alpha-helix */
	NOT_A_SECONDARY_STRUCTURE /*!< \brief not a recognized secondary structure */
	};

std::map< id, std::string > id_to_string = {
    {NOT_A_SECONDARY_STRUCTURE, "NOT A SECONDARY STRUCTURE"},
    {PARA_BRIDGE, "PARALLEL BRIDGE"},
    {ANTI_BRIDGE, "ANTIPARALLEL BRIDGE"},
    {FOUR_HELIX, "FOUR HELIX"}
  };

/*! \brief To be used if a classification of protein secondary structure is needed.

		\return the dpss class from a 'contact binary representation'
		@param contact describes the interaction matrix between two generic objects. 
		It is an 8 bit number and each bit is described by the following
		table:
		\htmlonly
		  <table>
		  	<tr>
		  		<td>(i-1,j)</td><td>(i-1,j+1)</td><td>(i,j)</td><td>(i,j+1)</td>
		  		<td>(j-1,i)</td><td>(j-1,i+1)</td><td>(j,i+1)</td><td>(j,i)</td> 
		  	</tr>
	  	</table>
		  \endhtmlonly
		\latexonly
		  \newline
		  \begin{tabular}{ |c|c|c|c|c|c|c|c| }
		    \hline
		    bit 7   & bit 6     & bit 5 & bit 4   & bit 3   & bit 2     & bit 1   & bit 0 \\
		    \hline
        (i-1,j) & (i-1,j+1) & (i,j) & (i,j+1) & (j-1,i) & (j-1,i+1) & (j,i+1) & (j,i) \\
        \hline
      \end{tabular}
    \vskip 1mm
		\endlatexonly
		@param distance is the distance of the two objects along the chain
		\note According to DPSS definition:
		\htmlonly
		  <table>
			  <tr><td>Dpss_id (def. number)</td><td>binary representation for `contact`</td> <td>decimal equivalent</td></tr>
			  <tr><td>PARA_BRIDGE (1) </td><td>10000010</td> <td>130</td></tr>
			  <tr><td>PARA_BRIDGE (2) </td><td>00011000</td> <td>24 </td></tr>
			  <tr><td>ANTI_BRIDGE (1) </td><td>00100001</td> <td>33 </td></tr>
			  <tr><td>ANTI_BRIDGE (2) </td><td>01000100</td> <td>68 </td></tr>
			  <tr><td>FOUR_HELIX      </td><td>10010000</td> <td>144</td></tr>
		  </table>
		\endhtmlonly
		\latexonly
		  \newline
		  \vskip 1mm
		  \begin{tabular}{ |c|c|c| }
		    \hline
		    Dpss id (def number) & binary representation for contact & decimal equivalent \\
        \hline
        PARA BRIDGE ($1$)  & $10000010$ & $130$ \\
			  PARA BRIDGE ($2$)  & $00011000$ & $24$  \\
			  ANTI BRIDGE ($1$)  & $00100001$ & $33$  \\
			  ANTI BRIDGE ($2$)  & $01000100$ & $68$  \\
			  FOUR HELIX         & $10010000$ & $144$ \\
			  \hline
      \end{tabular}
    \vskip 1mm
		\endlatexonly
		*/
inline id getSecondaryStructure(short contact, int distance){
	if(distance < 3)
		return NOT_A_SECONDARY_STRUCTURE;
	else if( (contact&130)==130 or (contact&24)==24 )
		return PARA_BRIDGE;
	else if( (contact&33)==33 or (contact&68)==68 )
		return ANTI_BRIDGE;
	else if( distance==3 and (contact&144)==144 )
		return FOUR_HELIX;
	return NOT_A_SECONDARY_STRUCTURE;
}

/*! \brief To be used if a classification of protein secondary structure is needed.

		\return the [dpss_id](@ref dpss_id) classification for the pair of objects i and j
		@tparam T has to support the operator `T+int`. Moreover `T-T` has to return int. Typically 
		`int` or `iterator` are used.
		@param map the H-bridge map of the system.
		@param  i   the first object of the pair. It may be one of the keys of [map](@ref h_bridge_map) or not.
		@param  j   the second object of the pair. It may be one of the keys of [map](@ref h_bridge_map) or not.
		\see getDpssClass */
template <typename T, typename bridge_map >
inline id getSecondaryStructure( bridge_map& map, T& i, T& j ){
	short contact = (map[std::make_pair(i-1,j)]<<7)|(map[std::make_pair(i-1,j+1)]<<6)|(map[std::make_pair(i,j)]<<5)|
									(map[std::make_pair(i,j+1)]<<4)|(map[std::make_pair(j-1,i)]<<3)|(map[std::make_pair(j-1,i+1)]<<2)|
									(map[std::make_pair(j,i+1)]<<1)|(map[std::make_pair(j,i)]<<0);
	return getSecondaryStructure(contact, abs(j-i) );
}

} //end namespace
} //end namespace

/*! \brief Print a string describing the id
*/
std::ostream& operator << ( std::ostream& out, BioCpp::dpss::id dpss_id){
  out << BioCpp::dpss::id_to_string[dpss_id];
	return out;
}
#endif

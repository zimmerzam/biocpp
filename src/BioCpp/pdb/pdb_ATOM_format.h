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

#ifndef PDB_ATOM_FORMAT_H
#define PDB_ATOM_FORMAT_H

#include <string>
#include <sstream>
#include <map>
#include <stdexcept>

#include <Eigen/Core>
#include "../polimers/amino_acid_id.h"
#include "../polimers/atom_id.h"
#include "../polimers/element_id.h"
#include "pdb_sections_and_records.h"

namespace BioCpp{
namespace pdb{

/*! \brief ATOM line in pdb file

		According to [pdb specification 3.3 for ATOM record ](http://www.wwpdb.org/documentation/format33/sect9.html#ATOM ) 
*/
struct atom_info{
  int serial;                     /*!< serial (progressive) number */
 	atom::id id;                     /*!< atom identifier (CA, CB, ...) */
 	char altLoc;                    /*!< alternative location */
 	amino_acid::id resName;          /*!< residue name */
 	char chainId;                   /*!< chain identifier */
 	int resSeq;                     /*!< residue number */
 	char iCode;                     /*!< code for insertion of residues */
 	Eigen::Vector3d coordinate;     /*!< coordinate */
 	double occupancy;               /*!< occupancy */
 	double tempFactor;              /*!< temperature factor */
 	element::id element;             /*!< element type (carbon, oxygen, ...) */
 	double charge;                  /*!< charge */
};

/*! \brief Get atom line info from a pdb ATOM line 
		
		\return a atom_info containing all informations from the ATOM line
		@param line the ATOM line
*/
atom_info read_atom_line(std::string& line){
	atom_info info;
	if(get_record(line) == ATOM){
	  info.serial = atoi(line.substr(6, 5).c_str()); 
  	info.id = atom::string_to_id[line.substr(12, 4)];
  	info.altLoc = line.substr(16, 1).c_str()[0]=='A' ? ' ' : line.substr(16, 1).c_str()[0];
  	info.resName = amino_acid::string_to_id.find( line.substr(17, 3) )!=amino_acid::string_to_id.end() ? amino_acid::string_to_id[line.substr(17, 3)] : amino_acid::UNK;
  	info.chainId = isdigit( line.substr(21, 1).c_str()[0] ) ? char( 'A'+atoi(line.substr(21, 1).c_str()) ) : line.substr(21, 1).c_str()[0];
  	if(info.chainId == ' ') info.chainId = 'A';
  	info.resSeq = atoi( line.substr(22, 4).c_str() );
  	info.iCode = line.substr(26, 1).c_str()[0];
  	info.coordinate = Eigen::Vector3d( atof(line.substr(30, 8).c_str()), atof(line.substr(38, 8).c_str()),
  	                           atof(line.substr(46, 8).c_str()));
	
  	try{
  	  info.occupancy = atof(line.substr(54, 6).c_str()) ? atof(line.substr(54, 6).c_str()) : 0.;
  	}
  	catch(std::out_of_range& oor){
  	  info.occupancy = 0.0;
  	}
	
  	try{
  	  info.tempFactor = atof(line.substr(60, 6).c_str()) ? atof(line.substr(60, 6).c_str()) : 0.;
  	}
  	catch (std::out_of_range& oor){
  	  info.tempFactor = 0.0;
  	}
  	info.element = element::string_to_id[line.substr(12, 4)]; //damn! I cannot always perform element::string_to_id[line.substr(76, 2)]
  	info.charge = 0;
  	try {
  	  if( isdigit( line.substr(78, 2).c_str()[0] ) and ( line.substr(78, 2).c_str()[1]=='-' or line.substr(78, 2).c_str()[1]=='+') ){
  	    info.charge = line.substr(78, 2).c_str()[1]=='-' ? -atof(line.substr(78, 1).c_str()) : atof(line.substr(78, 1).c_str());
  	  }
  	}
  	catch (std::out_of_range& oor){
  	  info.charge = 0;
  	}
  }
  return info;
};

/*! \brief Print an ATOM line from a atom_info
		
		@param out the output stream (i.e. `std::cout`)
		@param info the atom_info to be printed 
*/
struct print_atom_line{
  std::ostream& out;

  print_atom_line(std::ostream& dev): out(dev) {};

  /*! \brief Print an ATOM line
		@param out the output stream (i.e. `std::cout`)
		@param serial atom serial number
		@param id atom identifier
		@param altLoc alternative location
		@param resName residue name
		@param chainId chain identifier
		@param resSeq residue number
		@param iCode code for insertion of residues
		@param coordinate atom coordinate
		@param occupancy occupancy
		@param tempFactor temperature factor
		@param element element type
		@param charge atom charge 
  */
  std::ostream& operator()( int& serial, atom::id& id, char& altLoc, 
														amino_acid::id& resName, char& chainId, int& resSeq, char& iCode, 
														Eigen::Vector3d& coordinate, double& occupancy, 
														double& tempFactor, element::id& element, double& charge ){
	  out << "ATOM  ";
	  out << std::setw(5) << serial;
	  out << " ";
	  out << BioCpp::atom::id_to_string[id];
	  out << altLoc;
	  out << BioCpp::amino_acid::id_to_string[resName];
	  out << " ";
	  out << chainId;
	  out << std::setw(4) << resSeq;
	  out << iCode;
	  out << "   ";
	  out << std::fixed << std::setprecision(3);
	  out << std::setw(8) << coordinate(0);
	  out << std::setw(8) << coordinate(1);
	  out << std::setw(8) << coordinate(2);
	  out << std::fixed << std::setprecision(2);
	  out << std::setw(6) << occupancy;
	  out << std::setw(6) << tempFactor;
	  out << "          "; // blank
	  out << BioCpp::element::id_to_string[element];
	  if( charge == 0. )
		  out << "  ";
	  else{
		  out << std::setprecision(0) << fabs(charge);
		  if( charge > 0 )
			  out << '+';
		  else if(charge < 0){
			  out << '-';
		  }
	  }
	  out << std::endl;
	  return out;
  };

  std::ostream& operator()( atom_info& info){
	  return (*this)( info.serial, info.id, info.altLoc, 
									  info.resName, info.chainId, info.resSeq, 
									  info.iCode, info.coordinate, info.occupancy, 
									  info.tempFactor, info.element, info.charge );
  };

};

} // end namespace
} // end namespace

std::ostream& operator << (std::ostream& out, BioCpp::pdb::atom_info& info ){
  BioCpp::pdb::print_atom_line print(out);
	return print( info);
}

#endif

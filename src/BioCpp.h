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

/*! \mainpage BioCpp Documentation

    BioCpp is a simple C++ template library for basic operation on protein structures.

    \section get_it How to get it
    You can download BioCpp from https://github.com/zimmer/biocpp .
    \section Installation
    BioCpp is an header-only library: just copy the downloaded folder into your 
    project directory and include "/path/to/your/project/BioCpp/BioCpp.h" into
    your project file
*/

#ifndef BIOCPP_H
#define BIOCPP_H

#include "polimers/element_id.h"
#include "polimers/atom_id.h"
#include "polimers/amino_acid_id.h"
#include "polimers/base_container.h"
#include "polimers/is_container_of.h"
#include "polimers/Iterate_single.h"
#include "polimers/Iterate_pair.h"

#include "polimers/dpss/dpss_id.h"
#include "polimers/dpss/h_bridge_energy.h"

#include "polimers/morphology/surface_area_lcpo.h"

#include "polimers/reconstruction/bb_hydrogen.h"

#include "pdb/pdb.h"

#include "fasta/PAM30.h"
#include "fasta/PAM70.h"
#include "fasta/BLOSUM45.h"
#include "fasta/BLOSUM62.h"
#include "fasta/BLOSUM80.h"
#include "fasta/ZIMM1.h"
#include "fasta/NeedlemanWunsch.h"
#include "fasta/StrictNeedlemanWunsch.h"

/*! \example dpss.cpp */
/*! \example iterate.cpp */
/*! \example pdb.cpp */
/*! \example fasta_align.cpp */
/*! \example container.cpp */
/*! \example ids.cpp */

#endif

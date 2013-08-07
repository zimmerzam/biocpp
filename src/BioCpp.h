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
    You can download BioCpp from https://github.com/zimmerzam/biocpp .
    \section Installation
    BioCpp is an header-only library: just copy the downloaded folder into your 
    project directory and include "/path/to/your/project/BioCpp/BioCpp.h" into
    your project file
*/

#ifndef BIOCPP_H
#define BIOCPP_H

#include "BioCpp/base_container/base_container.h"
#include "BioCpp/base_container/is_container_of.h"
#include "BioCpp/base_container/Iterate_single.h"
#include "BioCpp/base_container/Iterate_pair.h"
#include "BioCpp/base_atom/base_atom.h"

/* include atom,element,... identifiers */
#if defined BIOCPP_INCLUDE_ID || defined BIOCPP_INCLUDE_ALL
  #include "BioCpp/standard/ids/element_id.h"
  #include "BioCpp/standard/ids/atom_id.h"
  #include "BioCpp/standard/ids/amino_acid_id.h"
  #include "BioCpp/standard/ids/moiety_id.h"
#endif

/* Include dssp definitions */
#if defined BIOCPP_INCLUDE_DPSS || defined BIOCPP_INCLUDE_ALL
  #include "BioCpp/polimers/proteins/dpss/dpss_id.h"
  #include "BioCpp/polimers/proteins/dpss/h_bridge_energy.h"
#endif

/* Include header files for protein morphology */
#if defined BIOCPP_INCLUDE_MORPHOLOGY || defined BIOCPP_INCLUDE_ALL
  #include "BioCpp/morphology/surface_area_lcpo.h"
  #if not defined BIOCPP_DONT_INCLUDE_CGAL
    #define BIOCPP_INCLUDE_CGAL
  #endif
#endif
#if defined BIOCPP_INCLUDE_CGAL
  #include "BioCpp/morphology/cgal_extensible_kernel.h"
  #include "BioCpp/morphology/cgal_alpha_shape_3/alpha_vertices.h"
#endif

/* Include header files for reconstructing atom coordinates */
#if defined BIOCPP_INCLUDE_RECONSTRUCTION || defined BIOCPP_INCLUDE_ALL
  #include "BioCpp/polimers/proteins/reconstruction/bb_hydrogen.h"
#endif

/* Include header files for reading/writing files in pdb format */
#if defined BIOCPP_INCLUDE_PDB || defined BIOCPP_INCLUDE_ALL
  #include "BioCpp/pdb/pdb.h"
#endif

/* Include header files for sequence alignment */
#if defined BIOCPP_INCLUDE_FASTA_SUBST_MATRIX_PAM30 || defined BIOCPP_INCLUDE_FASTA_SUBST_MATRIX_PAM70
  #define BIOCPP_INCLUDE_FASTA_SUBST_MATRIX
#endif
#if defined BIOCPP_INCLUDE_FASTA_SUBST_MATRIX_BLOSUM45 || defined BIOCPP_INCLUDE_FASTA_SUBST_MATRIX_BLOSUM62 || defined BIOCPP_INCLUDE_FASTA_SUBST_MATRIX_BLOSUM80
  #define BIOCPP_INCLUDE_FASTA_SUBST_MATRIX
#endif
#if defined BIOCPP_INCLUDE_FASTA_SUBST_MATRIX_ZIMM1 || defined BIOCPP_INCLUDE_FASTA_SUBST_MATRIX_ALL
  #define BIOCPP_INCLUDE_FASTA_SUBST_MATRIX
#endif
#if defined BIOCPP_INCLUDE_FASTA_SUBST_MATRIX || defined BIOCPP_INCLUDE_FASTA_ALIGNMENT
  #define BIOCPP_INCLUDE_FASTA
#endif
#if defined BIOCPP_INCLUDE_FASTA || defined BIOCPP_INCLUDE_ALL
  #if defined BIOCPP_INCLUDE_FASTA_SUBST_MATRIX || defined BIOCPP_INCLUDE_ALL
    #if defined BIOCPP_INCLUDE_FASTA_SUBST_MATRIX_PAM30 || defined BIOCPP_INCLUDE_FASTA_SUBST_MATRIX_ALL || defined BIOCPP_INCLUDE_ALL
      #include "BioCpp/fasta/PAM30.h"
    #endif
    #if defined BIOCPP_INCLUDE_FASTA_SUBST_MATRIX_PAM70 || defined BIOCPP_INCLUDE_FASTA_SUBST_MATRIX_ALL || defined BIOCPP_INCLUDE_ALL
    #include "BioCpp/fasta/PAM70.h"
    #endif
    #if defined BIOCPP_INCLUDE_FASTA_SUBST_MATRIX_BLOSUM45 || defined BIOCPP_INCLUDE_FASTA_SUBST_MATRIX_ALL || defined BIOCPP_INCLUDE_ALL
    #include "BioCpp/fasta/BLOSUM45.h"
    #endif
    #if defined BIOCPP_INCLUDE_FASTA_SUBST_MATRIX_BLOSUM62 || defined BIOCPP_INCLUDE_FASTA_SUBST_MATRIX_ALL || defined BIOCPP_INCLUDE_ALL
    #include "BioCpp/fasta/BLOSUM62.h"
    #endif
    #if defined BIOCPP_INCLUDE_FASTA_SUBST_MATRIX_BLOSUM80 || defined BIOCPP_INCLUDE_FASTA_SUBST_MATRIX_ALL || defined BIOCPP_INCLUDE_ALL
    #include "BioCpp/fasta/BLOSUM80.h"
    #endif    
    #if defined BIOCPP_INCLUDE_FASTA_SUBST_MATRIX_ZIMM1 || defined BIOCPP_INCLUDE_FASTA_SUBST_MATRIX_ALL || defined BIOCPP_INCLUDE_ALL
    #include "BioCpp/fasta/ZIMM1.h"
    #endif
  #endif
  #if defined BIOCPP_INCLUDE_FASTA_ALIGNMENT || defined BIOCPP_INCLUDE_ALL
    #include "BioCpp/fasta/NeedlemanWunsch.h"
    #include "BioCpp/fasta/StrictNeedlemanWunsch.h"
  #endif
#endif

#include "BioCpp/utils/sgn.h"

/* Include standard definitions of residue, chain, ... */
#if defined BIOCPP_INCLUDE_STANDARD || defined BIOCPP_INCLUDE_ALL
  #include "BioCpp/standard.h"
#endif

/* Use topology */
#if defined BIOCPP_INCLUDE_TOPOLOGY || defined BIOCPP_INCLUDE_ALL
  #include "BioCpp/topology/topology.h"
#endif

/* Use Denavit-Hartenberg notation */
#if defined BIOCPP_INCLUDE_DENAVIT_HARTENBERG || defined BIOCPP_INCLUDE_ALL
  #include "BioCpp/topology/denavit_hartenberg/dh_matrix.h"
  #include "BioCpp/topology/denavit_hartenberg/dh_edge.h"
  #include "BioCpp/topology/denavit_hartenberg/dh_chain.h"
#endif

/*! \example dpss.cpp */
/*! \example iterate.cpp */
/*! \example pdb.cpp */
/*! \example fasta_align.cpp */
/*! \example container.cpp */
/*! \example ids.cpp */

#endif

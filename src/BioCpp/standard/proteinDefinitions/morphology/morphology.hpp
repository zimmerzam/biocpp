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

#ifndef BIOCPP_STANDARD_MORPHOLOGY_DEFINITION_HPP
#define BIOCPP_STANDARD_MORPHOLOGY_DEFINITION_HPP

#include <BioCpp/base_atom/base_atom.hpp>
#include <BioCpp/io_files/model/model.hxx>
#include <BioCpp/standard/proteinDefinitions/templated.hxx>
#include <BioCpp/standard/chemicalComponentDictionary/residue_dictionary_standard.hpp>

#include <BioCpp/morphology/cgal_extensible_kernel.hxx>
#include <CGAL/Filtered_kernel.h>
#include <CGAL/Alpha_shape_3.h>
#include <CGAL/Delaunay_triangulation_3.h>

namespace BioCpp{
namespace standard{
namespace morphology{

typedef BioCpp::morphology::cgal_extensible_kernel::point<BioCpp::base::atom> atom;
typedef BioCpp::standard::model<atom>::type model;
typedef BioCpp::standard::residue<atom>::type residue;
typedef BioCpp::standard::chain<atom>::type chain;
typedef BioCpp::standard::complex<atom>::type complex;
typedef BioCpp::standard::complex_constructor<atom, BioCpp::residue::dictionary_t> complex_constructor;

typedef BioCpp::morphology::cgal_extensible_kernel::kernel<BioCpp::base_atom, double>  extKernel;
typedef CGAL::Filtered_kernel_adaptor<extKernel>                                       kernel;
typedef CGAL::Alpha_shape_cell_base_3<kernel>                                          cell_base;
typedef CGAL::Alpha_shape_vertex_base_3<kernel>                                        vertex_base;
typedef CGAL::Triangulation_data_structure_3<vertex_base,cell_base>                    triangulation_data_structure;
typedef CGAL::Delaunay_triangulation_3<kernel, triangulation_data_structure>           triangulation_3;
typedef CGAL::Alpha_shape_3<triangulation_3>                                           alpha_shape_3;

}
}
}

#endif

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

#ifndef BIOCPP_STANDARD_MORPHOLOGY_KERNEL_H
#define BIOCPP_STANDARD_MORPHOLOGY_KERNEL_H

#include <BioCpp/base_atom/base_atom.h>
#include <BioCpp/morphology/cgal_extensible_kernel.h>
#include <CGAL/Filtered_kernel.h>

namespace BioCpp{
namespace standard{
namespace morphology{

typedef BioCpp::morphology::cgal_extensible_kernel::kernel<BioCpp::base_atom, double>  extKernel;
typedef CGAL::Filtered_kernel_adaptor<extKernel>                                       kernel;

}
}
}

#endif

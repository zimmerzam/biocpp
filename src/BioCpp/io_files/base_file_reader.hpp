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

#ifndef BIOCPP_BASE_FILE_READER_HPP
#define BIOCPP_BASE_FILE_READER_HPP

#include <map>
#include <BioCpp/utils/errors_and_warnings/errors_and_warnings.hpp>

#include "model/model.hxx"
#include "seqres/seqres_record.hpp"
#include "file_reader_init_flags.hpp"

namespace BioCpp{
namespace io{

/*! \brief This is a `structure file` object.

    Use this to read a pdb file and get a list of atom_info.
*/
class file{
  protected:
    const char* filename;                        /*!< The pdb filename. */
    char* buffer;                                /*!< The pdb contents. */
    std::map<int, std::streampos> model_beg_pos; /*!< Positions (in the buffer stream) of the begin of each model. */
    std::map<int, std::streampos> model_end_pos; /*!< Positions (in the buffer stream) of the end of each model.   */
  public:
    BioCpp::warning warning;
    BioCpp::error error;
    int n_models;                                /*!< Total number of models found in the file */
    seqres_record TseqRes;                       /*!< The sequence as read from SEQRES record  */
    seqres_record RseqRes;                       /*!< The sequence as read from ATOM record    */
    
    /* \brief Default constructor is void */
    file(const char* spp_name, int init_flag){};
    
    /*! \brief Destroy the file 
        \note The file will not be deleted! 
    */
    ~file(){};
    
    template<typename atom_t>
    typename model<atom_t>::type readModel(int mdl);
    
    template<typename atom_t, typename eleDict, typename atmDict, typename resDict>
    typename model<atom_t>::type readModel(int mdl, 
                                           eleDict& ele_dict, atmDict& atm_dict,
                                           resDict& res_dict);
};

} // end namespace
} // end namespace

#endif

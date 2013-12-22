#include <map>
#include <string>
#include "errors_and_warnings.hpp"

namespace BioCpp{

std::map< BioCpp::warning, std::string > warning_to_string = {
    {WAR_NONE, "No warnings "},
    {PDB_SEQRES_NOT_FOUND, "No SEQRES section found"},
    {PDB_BACKBONE_HOLE, "Chain break detected: one ore more residues may be missing"},
    {ALIGN_SEQUENCE_NEQ_FASTA, "Protein sequence in the model structure is different from the primary sequence"},
    {ALIGN_NOT_A_VALID_SEQUENCE, "Not a valid sequence"}
  };

std::map< BioCpp::error, std::string > error_to_string = {
    {ERR_NONE, "Success"},
    {PDB_COORDINATE_NOT_FOUND, "No ATOM section found"},
    {ALIGN_MISSING_PRIMARY_SEQUENCE, "Alignment is not possible due to missing target sequence"},
    {ALIGN_FAILED, "Alignment failed"},
    {UNSPECIFIED, "Unspecified problem. Sorry"}
  };

} //end namespace

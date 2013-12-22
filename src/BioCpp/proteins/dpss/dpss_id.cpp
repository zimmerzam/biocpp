#include "dpss_id.hpp"

namespace BioCpp{
namespace dpss{

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

template <typename T, typename cnt_map >
inline id getSecondaryStructure( cnt_map& map, T& i, T& j ){
	short contact = (map[std::make_pair(i-1,j)]<<7)|(map[std::make_pair(i-1,j+1)]<<6)|(map[std::make_pair(i,j)]<<5)|
									(map[std::make_pair(i,j+1)]<<4)|(map[std::make_pair(j-1,i)]<<3)|(map[std::make_pair(j-1,i+1)]<<2)|
									(map[std::make_pair(j,i+1)]<<1)|(map[std::make_pair(j,i)]<<0);
	return getSecondaryStructure(contact, abs(j-i) );
}

} //end namespace
} //end namespace

std::ostream& operator << ( std::ostream& out, BioCpp::dpss::id dpss_id){
  out << BioCpp::dpss::id_to_string[dpss_id];
	return out;
}

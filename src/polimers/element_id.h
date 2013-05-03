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

#ifndef ELEMENT_ID
#define ELEMENT_ID

#include <iostream>
#include <cstring>
#include <map>

#include "../utils/list_of_type.h"

namespace BioCpp{
namespace element{

/*! \brief id list of the standard atom elements

*/
enum id { Unk, Ac, Ag, Al, Am, Ar, As, At, Au, B, Ba, Be, Bh, Bi, Bk, Br,
          C, Ca, Cd, Ce, Cf, Cl, Cm, Cn, Co, Cr, Cs, Cu, Db, Ds, Dy,
          Er, Es, Eu, F, Fe, Fm, Fr, Ga, Gd, Ge, H, He, Hf, Hg, Ho,
          Hs, I, In, Ir, K, Kr, La, Li, Lr, Lu, Md, Mg, Mn, Mo, Mt,
          N, Na, Nb, Nd, Ne, Ni, No, Np, O, Os, P, Pa, Pb, Pd, Pm,
          Po, Pr, Pt, Pu, Ra, Rb, Re, Rf, Rg, Rh, Rn, Ru, S, Sb, Sc,
          Se, Sg, Si, Sm, Sn, Sr, Ta, Tb, Tc, Te, Th, Ti, Tl, Tm, U,
          Uuh, Uuo, Uup, Uuq, Uus, Uut, V, W, Xe, Y, Yb, Zn, Zr
};

/*! \brief Map a string to an element id

    All strings that do not correspond to a valid element nomenclature are 
    mapped onto X_ by default

    \code
      #include <iostream>
      #include "../src/BioCpp.h"

      int main(){

        std::cout << BioCpp::element::string_to_id["HE22"] << std::endl; // output: " H"
        std::cout << BioCpp::element::string_to_id["Rb"] << std::endl; // output: " H"

        return 0;
      }
    \endcode

*/
std::map< std::string, id > string_to_id = 
        map_list_of_type< std::string, id >
    /* these map PDB element entries (chars 76 and 77) to element_id */
    ("Unk", Unk)
    ("Ac", Ac)("Ag", Ag)("Al", Al)("Am", Am)("Ar", Ar)("As", As)("At", At)
    ("Au", Au)(" B",  B)("Ba", Ba)("Be", Be)("Bh", Bh)("Bi", Bi)("Bk", Bk)
    ("Br", Br)(" C",  C)("Ca", Ca)("Cd", Cd)("Ce", Ce)("Cf", Cf)("Cl", Cl)
    ("Cm", Cm)("Cn", Cn)("Co", Co)("Cr", Cr)("Cs", Cs)("Cu", Cu)("Db", Db)
    ("Ds", Ds)("Dy", Dy)("Er", Er)("Es", Es)("Eu", Eu)(" F",  F)("Fe", Fe)
    ("Fm", Fm)("Fr", Fr)("Ga", Ga)("Gd", Gd)("Ge", Ge)(" H",  H)("He", He)
    ("Hf", Hf)("Hg", Hg)("Ho", Ho)("Hs", Hs)(" I",  I)("In", In)("Ir", Ir)
    (" K",  K)("Kr", Kr)("La", La)("Li", Li)("Lr", Lr)("Lu", Lu)("Md", Md)
    ("Mg", Mg)("Mn", Mn)("Mo", Mo)("Mt", Mt)(" N",  N)("Na", Na)("Nb", Nb)
    ("Nd", Nd)("Ne", Ne)("Ni", Ni)("No", No)("Np", Np)(" O",  O)("Os", Os)
    (" P",  P)("Pa", Pa)("Pb", Pb)("Pd", Pd)("Pm", Pm)("Po", Po)("Pr", Pr)
    ("Pt", Pt)("Pu", Pu)("Ra", Ra)("Rb", Rb)("Re", Re)("Rf", Rf)("Rg", Rg)
    ("Rh", Rh)("Rn", Rn)("Ru", Ru)(" S",  S)("Sb", Sb)("Sc", Sc)("Se", Se)
    ("Sg", Sg)("Si", Si)("Sm", Sm)("Sn", Sn)("Sr", Sr)("Ta", Ta)("Tb", Tb)
    ("Tc", Tc)("Te", Te)("Th", Th)("Ti", Ti)("Tl", Tl)("Tm", Tm)(" U",  U)
    ("uh", Uuh)("uo", Uuo)("up", Uup)("uq", Uuq)("us", Uus)("ut", Uut)
    (" V",  V)(" W",  W)("Xe", Xe)(" Y",  Y)("Yb", Yb)("Zn", Zn)("Zr", Zr)
    /* these map PDB atom id entries (chars 12, 13, 14 and 15) to element_id */
    ( " C  ", C)( " CA ", C)( " CB ", C)( " CD ", C)
    ( " CD1", C)( " CD2", C)( " CD3", C)( " CE ", C)
    ( " CE1", C)( " CE2", C)( " CE3", C)( " CG ", C)
    ( " CG1", C)( " CG2", C)( " CG3", C)( " CZ ", C)
    ( " CZ1", C)( " CZ2", C)( " CZ3", C)( " CH1", C)
    ( " CH2", C)( " H  ", H)
    ( " HA ", H)( " HA1", H)( "1HA ", H)( " HA2", H)( "2HA ", H)( " HA3", H)( "3HA ", H)
    ( " HB ", H)( " HB1", H)( "1HB ", H)( " HB2", H)( "2HB ", H)( " HB3", H)( "3HB ", H)
    ( " HD ", H)( " HD1", H)( "1HD ", H)( " HD2", H)( "2HD ", H)( " HD3", H)( "3HD ", H)
    ( "HD11", H)( "1HD1", H)( "HD12", H)( "2HD1", H)( "HD13", H)( "3HD1", H)
    ( "HD21", H)( "1HD2", H)( "HD22", H)( "2HD2", H)( "HD23", H)( "3HD2", H)
    ( " HE ", H)( " HE1", H)( "1HE ", H)( " HE2", H)( "2HE ", H)( " HE3", H)( "3HE ", H)
    ( "HE11", H)( "1HE1", H)( "HE12", H)( "2HE1", H)( "HE13", H)( "3HE1", H)
    ( "HE21", H)( "1HE2", H)( "HE22", H)( "2HE2", H)( "HE23", H)( "3HE2", H)
    ( " HG ", H)( " HG1", H)( "1HG ", H)( " HG2", H)( "2HG ", H)( " HG3", H)( "3HG ", H)
    ( "HG11", H)( "1HG1", H)( "HG12", H)( "2HG1", H)( "HG13", H)( "3HG1", H)
    ( "HG21", H)( "1HG2", H)( "HG22", H)( "2HG2", H)( "HG23", H)( "3HG2", H)
    ( " HH ", H)( " HH1", H)( "1HH ", H)( " HH2", H)( "2HH ", H)( " HH3", H)( "3HH ", H)
    ( "HH11", H)( "1HH1", H)( "HH12", H)( "2HH1", H)( "HH13", H)( "3HH1", H)
    ( "HH21", H)( "1HH2", H)( "HH22", H)( "2HH2", H)( "HH23", H)( "3HH2", H)
    ( " HN ", H)( " HN1", H)( "1HN ", H)( " HN2", H)( "2HN ", H)( " HN3", H)( "3HN", H)
    ( "HN11", H)( "1HN1", H)( "HN12", H)( "2HN1", H)( "HN13", H)( "3HN1", H)
    ( "HN21", H)( "1HN2", H)( "HN22", H)( "2HN2", H)( "HN23", H)( "3HN2", H)
    ( " HZ ", H)( " HZ1", H)( "1HZ ", H)( " HZ2", H)( "2HZ ", H)( " HZ3", H)( "3HZ ", H)
    ( " H1 ", H)( " H2 ", H)( " H3 ", H)
    ( " N  ", N)( " ND ", N)
    ( " ND1", N)( " ND2", N)( " NE ", N)( " NE1", N)
    ( " NE2", N)( " NH ", N)( " NH1", N)( " NH2", N)
    ( " NZ ", N)( " O  ", O)( " OD ", O)( " OD1", O)
    ( " OD2", O)( " OE ", O)( " OE1", O)( " OE2", O)
    ( " OH ", O)( " OH1", O)( " OH2", O)( " OXT", O)
    ( " OG ", O)( " OG1", O)( " OG2", O)
    ( " SD ", S)( " SG ", S);

/*! \brief Map element id to a two-letter string

    \code
      #include <iostream>
      #include "../src/BioCpp.h"

      int main(){

        std::cout << BioCpp::element::id_to_string[BioCpp::element::N] << std::endl; // output: " N"
  
        return 0;
      }
    \endcode

*/
std::map< id, std::string > id_to_string = 
        map_list_of_type< id, std::string >
    (Unk, "Unk")
    ( Ac, "Ac")( Ag, "Ag")( Al, "Al")( Am, "Am")( Ar, "Ar")( As, "As")( At, "At")( Au, "Au")
    ( B, " B")( Ba, "Ba")( Be, "Be")( Bh, "Bh")( Bi, "Bi")( Bk, "Bk")( Br, "Br")( C, " C")
    ( Ca, "Ca")( Cd, "Cd")( Ce, "Ce")( Cf, "Cf")( Cl, "Cl")( Cm, "Cm")( Cn, "Cn")( Co, "Co")
    ( Cr, "Cr")( Cs, "Cs")( Cu, "cu")( Db, "Db")( Ds, "Ds")( Dy, "Dy")( Er, "Er")( Es, "Es")
    ( Eu, "Eu")( F, " F")( Fe, "Fe")( Fm, "Fm")( Fr, "Fr")( Ga, "Ga")( Gd, "Gd")( Ge, "Ge")
    ( H, " H")( He, "He")( Hf, "Hf")( Hg, "Hg")( Ho, "Ho")( Hs, "Hs")( I, " I")( In, "In")
    ( Ir, "Ir")( K, " K")( Kr, "Kr")( La, "La")( Li, "Li")( Lr, "Lr")( Lu, "Lu")( Md, "Md")
    ( Mg, "Mg")( Mn, "Mn")( Mo, "Mo")( Mt, "Mt")( N, " N")( Na, "Na")( Nb, "Nb")( Nd, "Nd")
    ( Ne, "Ne")( Ni, "Ni")( No, "No")( Np, "Np")( O, " O")( Os, "Os")( P, " P")( Pa, "Pa")
    ( Pb, "Pb")( Pd, "Pd")( Pm, "Pm")( Po, "Po")( Pr, "Pr")( Pt, "Pt")( Pu, "Pu")( Ra, "Ra")
    ( Rb, "Rb")( Re, "Re")( Rf, "Rf")( Rg, "rg")( Rh, "Rh")( Rn, "Rn")( Ru, "Ru")( S, " S")
    ( Sb, "Sb")( Sc, "Sc")( Se, "Se")( Sg, "Sg")( Si, "Si")( Sm, "Sm")( Sn, "Sn")( Sr, "Sr")
    ( Ta, "Ta")( Tb, "Tb")( Tc, "Tc")( Te, "Te")( Th, "Th")( Ti, "Ti")( Tl, "Tl")( Tm, "Tm")
    ( U, " U")( Uuh, "uh")( Uuo, "uo")( Uup, "up")( Uuq, "uq")( Uus, "us")( Uut, "ut")( V, " V")
    ( W, " W")( Xe, "Xe")( Y, " Y")( Yb, "Yb")( Zn, "Zn")( Zr, "Zr");

} // end element
} // end BioCpp

/*! \brief Print an id as a two-letter string
    
    \code
      #include <iostream>
      #include "../src/BioCpp.h"

      int main(){

        std::cout << BioCpp::element::C << std::endl; // output: " C"
        std::cout << BioCpp::element::string_to_id["HE22"] << std::endl; // output: " H"
        std::cout << BioCpp::element::id_to_string[BioCpp::element::N] << std::endl; // output: " N"
  
        return 0;
      }
    \endcode
*/
std::ostream& operator << (std::ostream& out, BioCpp::element::id element_id ){
  out << BioCpp::element::id_to_string[element_id];
	return out;
}

#endif

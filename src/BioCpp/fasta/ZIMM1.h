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

#ifndef ZIMM1_H
#define ZIMM1_H

#include "substitution_matrix.h"
namespace BioCpp{
namespace fasta{

substitution_matrix ZIMM1( {
    {'A', {{'A',4.0},{'R',-10000.0},{'N',-10000.0},{'D',-10000.0},{'C',-10000.0},{'Q',-10000.0},{'E',-10000.0},{'G',-10000.0},{'H',-10000.0},{'I',-10000.0},{'L',-10000.0},{'K',-10000.0},{'M',-10000.0},{'F',-10000.0},{'P',-10000.0},{'S',-10000.0},{'T',-10000.0},{'W',-10000.0},{'Y',-10000.0},{'V',-10000.0},{'B',-10000.0},{'J',-10000.0},{'Z',-10000.0},{'X',1.0},{'*',-10000.0},{'-',-10000.0}}},
    {'R', {{'A',-10000.0},{'R',5.0},{'N',-10000.0},{'D',-10000.0},{'C',-10000.0},{'Q',-10000.0},{'E',-10000.0},{'G',-10000.0},{'H',-10000.0},{'I',-10000.0},{'L',-10000.0},{'K',-10000.0},{'M',-10000.0},{'F',-10000.0},{'P',-10000.0},{'S',-10000.0},{'T',-10000.0},{'W',-10000.0},{'Y',-10000.0},{'V',-10000.0},{'B',-10000.0},{'J',-10000.0},{'Z',-10000.0},{'X',1.0},{'*',-10000.0},{'-',-10000.0}}},
    {'N', {{'A',-10000.0},{'R',-10000.0},{'N',6.0},{'D',-10000.0},{'C',-10000.0},{'Q',-10000.0},{'E',-10000.0},{'G',-10000.0},{'H',-10000.0},{'I',-10000.0},{'L',-10000.0},{'K',-10000.0},{'M',-10000.0},{'F',-10000.0},{'P',-10000.0},{'S',-10000.0},{'T',-10000.0},{'W',-10000.0},{'Y',-10000.0},{'V',-10000.0},{'B',4.0},{'J',-10000.0},{'Z',-10000.0},{'X',1.0},{'*',-10000.0},{'-',-10000.0}}},
    {'D', {{'A',-10000.0},{'R',-10000.0},{'N',-10000.0},{'D',6.0},{'C',-10000.0},{'Q',-10000.0},{'E',-10000.0},{'G',-10000.0},{'H',-10000.0},{'I',-10000.0},{'L',-10000.0},{'K',-10000.0},{'M',-10000.0},{'F',-10000.0},{'P',-10000.0},{'S',-10000.0},{'T',-10000.0},{'W',-10000.0},{'Y',-10000.0},{'V',-10000.0},{'B',4.0},{'J',-10000.0},{'Z',-10000.0},{'X',1.0},{'*',-10000.0},{'-',-10000.0}}},
    {'C', {{'A',-10000.0},{'R',-10000.0},{'N',-10000.0},{'D',-10000.0},{'C',9.0},{'Q',-10000.0},{'E',-10000.0},{'G',-10000.0},{'H',-10000.0},{'I',-10000.0},{'L',-10000.0},{'K',-10000.0},{'M',-10000.0},{'F',-10000.0},{'P',-10000.0},{'S',-10000.0},{'T',-10000.0},{'W',-10000.0},{'Y',-10000.0},{'V',-10000.0},{'B',-10000.0},{'J',-10000.0},{'Z',-10000.0},{'X',1.0},{'*',-10000.0},{'-',-10000.0}}},
    {'Q', {{'A',-10000.0},{'R',-10000.0},{'N',-10000.0},{'D',-10000.0},{'C',-10000.0},{'Q',5.0},{'E',-10000.0},{'G',-10000.0},{'H',-10000.0},{'I',-10000.0},{'L',-10000.0},{'K',-10000.0},{'M',-10000.0},{'F',-10000.0},{'P',-10000.0},{'S',-10000.0},{'T',-10000.0},{'W',-10000.0},{'Y',-10000.0},{'V',-10000.0},{'B',-10000.0},{'J',-10000.0},{'Z',4.0},{'X',1.0},{'*',-10000.0},{'-',-10000.0}}},
    {'E', {{'A',-10000.0},{'R',-10000.0},{'N',-10000.0},{'D',-10000.0},{'C',-10000.0},{'Q',-10000.0},{'E',5.0},{'G',-10000.0},{'H',-10000.0},{'I',-10000.0},{'L',-10000.0},{'K',-10000.0},{'M',-10000.0},{'F',-10000.0},{'P',-10000.0},{'S',-10000.0},{'T',-10000.0},{'W',-10000.0},{'Y',-10000.0},{'V',-10000.0},{'B',-10000.0},{'J',-10000.0},{'Z',4.0},{'X',1.0},{'*',-10000.0},{'-',-10000.0}}},
    {'G', {{'A',-10000.0},{'R',-10000.0},{'N',-10000.0},{'D',-10000.0},{'C',-10000.0},{'Q',-10000.0},{'E',-10000.0},{'G',6.0},{'H',-10000.0},{'I',-10000.0},{'L',-10000.0},{'K',-10000.0},{'M',-10000.0},{'F',-10000.0},{'P',-10000.0},{'S',-10000.0},{'T',-10000.0},{'W',-10000.0},{'Y',-10000.0},{'V',-10000.0},{'B',-10000.0},{'J',-10000.0},{'Z',-10000.0},{'X',1.0},{'*',-10000.0},{'-',-10000.0}}},
    {'H', {{'A',-10000.0},{'R',-10000.0},{'N',-10000.0},{'D',-10000.0},{'C',-10000.0},{'Q',-10000.0},{'E',-10000.0},{'G',-10000.0},{'H',8.0},{'I',-10000.0},{'L',-10000.0},{'K',-10000.0},{'M',-10000.0},{'F',-10000.0},{'P',-10000.0},{'S',-10000.0},{'T',-10000.0},{'W',-10000.0},{'Y',-10000.0},{'V',-10000.0},{'B',-10000.0},{'J',-10000.0},{'Z',-10000.0},{'X',1.0},{'*',-10000.0},{'-',-10000.0}}},
    {'I', {{'A',-10000.0},{'R',-10000.0},{'N',-10000.0},{'D',-10000.0},{'C',-10000.0},{'Q',-10000.0},{'E',-10000.0},{'G',-10000.0},{'H',-10000.0},{'I',4.0},{'L',-10000.0},{'K',-10000.0},{'M',-10000.0},{'F',-10000.0},{'P',-10000.0},{'S',-10000.0},{'T',-10000.0},{'W',-10000.0},{'Y',-10000.0},{'V',-10000.0},{'B',-10000.0},{'J',3.0},{'Z',-10000.0},{'X',1.0},{'*',-10000.0},{'-',-10000.0}}},
    {'L', {{'A',-10000.0},{'R',-10000.0},{'N',-10000.0},{'D',-10000.0},{'C',-10000.0},{'Q',-10000.0},{'E',-10000.0},{'G',-10000.0},{'H',-10000.0},{'I',-10000.0},{'L',4.0},{'K',-10000.0},{'M',-10000.0},{'F',-10000.0},{'P',-10000.0},{'S',-10000.0},{'T',-10000.0},{'W',-10000.0},{'Y',-10000.0},{'V',-10000.0},{'B',-10000.0},{'J',3.0},{'Z',-10000.0},{'X',1.0},{'*',-10000.0},{'-',-10000.0}}},
    {'K', {{'A',-10000.0},{'R',-10000.0},{'N',-10000.0},{'D',-10000.0},{'C',-10000.0},{'Q',-10000.0},{'E',-10000.0},{'G',-10000.0},{'H',-10000.0},{'I',-10000.0},{'L',-10000.0},{'K',5.0},{'M',-10000.0},{'F',-10000.0},{'P',-10000.0},{'S',-10000.0},{'T',-10000.0},{'W',-10000.0},{'Y',-10000.0},{'V',-10000.0},{'B',-10000.0},{'J',-10000.0},{'Z',-10000.0},{'X',1.0},{'*',-10000.0},{'-',-10000.0}}},
    {'M', {{'A',-10000.0},{'R',-10000.0},{'N',-10000.0},{'D',-10000.0},{'C',-10000.0},{'Q',-10000.0},{'E',-10000.0},{'G',-10000.0},{'H',-10000.0},{'I',-10000.0},{'L',-10000.0},{'K',-10000.0},{'M',5.0},{'F',-10000.0},{'P',-10000.0},{'S',-10000.0},{'T',-10000.0},{'W',-10000.0},{'Y',-10000.0},{'V',-10000.0},{'B',-10000.0},{'J',2.0},{'Z',-10000.0},{'X',1.0},{'*',-10000.0},{'-',-10000.0}}},
    {'F', {{'A',-10000.0},{'R',-10000.0},{'N',-10000.0},{'D',-10000.0},{'C',-10000.0},{'Q',-10000.0},{'E',-10000.0},{'G',-10000.0},{'H',-10000.0},{'I',-10000.0},{'L',-10000.0},{'K',-10000.0},{'M',-10000.0},{'F',6.0},{'P',-10000.0},{'S',-10000.0},{'T',-10000.0},{'W',-10000.0},{'Y',-10000.0},{'V',-10000.0},{'B',-10000.0},{'J',-10000.0},{'Z',-10000.0},{'X',1.0},{'*',-10000.0},{'-',-10000.0}}},
    {'P', {{'A',-10000.0},{'R',-10000.0},{'N',-10000.0},{'D',-10000.0},{'C',-10000.0},{'Q',-10000.0},{'E',-10000.0},{'G',-10000.0},{'H',-10000.0},{'I',-10000.0},{'L',-10000.0},{'K',-10000.0},{'M',-10000.0},{'F',-10000.0},{'P',7.0},{'S',-10000.0},{'T',-10000.0},{'W',-10000.0},{'Y',-10000.0},{'V',-10000.0},{'B',-10000.0},{'J',-10000.0},{'Z',-10000.0},{'X',1.0},{'*',-10000.0},{'-',-10000.0}}},
    {'S', {{'A',-10000.0},{'R',-10000.0},{'N',-10000.0},{'D',-10000.0},{'C',-10000.0},{'Q',-10000.0},{'E',-10000.0},{'G',-10000.0},{'H',-10000.0},{'I',-10000.0},{'L',-10000.0},{'K',-10000.0},{'M',-10000.0},{'F',-10000.0},{'P',-10000.0},{'S',4.0},{'T',-10000.0},{'W',-10000.0},{'Y',-10000.0},{'V',-10000.0},{'B',-10000.0},{'J',-10000.0},{'Z',-10000.0},{'X',1.0},{'*',-10000.0},{'-',-10000.0}}},
    {'T', {{'A',-10000.0},{'R',-10000.0},{'N',-10000.0},{'D',-10000.0},{'C',-10000.0},{'Q',-10000.0},{'E',-10000.0},{'G',-10000.0},{'H',-10000.0},{'I',-10000.0},{'L',-10000.0},{'K',-10000.0},{'M',-10000.0},{'F',-10000.0},{'P',-10000.0},{'S',-10000.0},{'T',5.0},{'W',-10000.0},{'Y',-10000.0},{'V',-10000.0},{'B',-10000.0},{'J',-10000.0},{'Z',-10000.0},{'X',1.0},{'*',-10000.0},{'-',-10000.0}}},
    {'W', {{'A',-10000.0},{'R',-10000.0},{'N',-10000.0},{'D',-10000.0},{'C',-10000.0},{'Q',-10000.0},{'E',-10000.0},{'G',-10000.0},{'H',-10000.0},{'I',-10000.0},{'L',-10000.0},{'K',-10000.0},{'M',-10000.0},{'F',-10000.0},{'P',-10000.0},{'S',-10000.0},{'T',-10000.0},{'W',11.0},{'Y',-10000.0},{'V',-10000.0},{'B',-10000.0},{'J',-10000.0},{'Z',-10000.0},{'X',1.0},{'*',-10000.0},{'-',-10000.0}}},
    {'Y', {{'A',-10000.0},{'R',-10000.0},{'N',-10000.0},{'D',-10000.0},{'C',-10000.0},{'Q',-10000.0},{'E',-10000.0},{'G',-10000.0},{'H',-10000.0},{'I',-10000.0},{'L',-10000.0},{'K',-10000.0},{'M',-10000.0},{'F',-10000.0},{'P',-10000.0},{'S',-10000.0},{'T',-10000.0},{'W',-10000.0},{'Y',7.0},{'V',-10000.0},{'B',-10000.0},{'J',-10000.0},{'Z',-10000.0},{'X',1.0},{'*',-10000.0},{'-',-10000.0}}},
    {'V', {{'A',-10000.0},{'R',-10000.0},{'N',-10000.0},{'D',-10000.0},{'C',-10000.0},{'Q',-10000.0},{'E',-10000.0},{'G',-10000.0},{'H',-10000.0},{'I',-10000.0},{'L',-10000.0},{'K',-10000.0},{'M',-10000.0},{'F',-10000.0},{'P',-10000.0},{'S',-10000.0},{'T',-10000.0},{'W',-10000.0},{'Y',-10000.0},{'V',4.0},{'B',-10000.0},{'J',2.0},{'Z',-10000.0},{'X',1.0},{'*',-10000.0},{'-',-10000.0}}},
    {'B', {{'A',-10000.0},{'R',-10000.0},{'N',4.0},{'D',4.0},{'C',-10000.0},{'Q',-10000.0},{'E',-10000.0},{'G',-10000.0},{'H',-10000.0},{'I',-10000.0},{'L',-10000.0},{'K',-10000.0},{'M',-10000.0},{'F',-10000.0},{'P',-10000.0},{'S',-10000.0},{'T',-10000.0},{'W',-10000.0},{'Y',-10000.0},{'V',-10000.0},{'B',4.0},{'J',-10000.0},{'Z',-10000.0},{'X',1.0},{'*',-10000.0},{'-',-10000.0}}},
    {'J', {{'A',-10000.0},{'R',-10000.0},{'N',-10000.0},{'D',-10000.0},{'C',-10000.0},{'Q',-10000.0},{'E',-10000.0},{'G',-10000.0},{'H',-10000.0},{'I',3.0},{'L',3.0},{'K',-10000.0},{'M',2.0},{'F',-10000.0},{'P',-10000.0},{'S',-10000.0},{'T',-10000.0},{'W',-10000.0},{'Y',-10000.0},{'V',2.0},{'B',-10000.0},{'J',3.0},{'Z',-10000.0},{'X',1.0},{'*',-10000.0},{'-',-10000.0}}},
    {'Z', {{'A',-10000.0},{'R',-10000.0},{'N',-10000.0},{'D',-10000.0},{'C',-10000.0},{'Q',4.0},{'E',4.0},{'G',-10000.0},{'H',-10000.0},{'I',-10000.0},{'L',-10000.0},{'K',-10000.0},{'M',-10000.0},{'F',-10000.0},{'P',-10000.0},{'S',-10000.0},{'T',-10000.0},{'W',-10000.0},{'Y',-10000.0},{'V',-10000.0},{'B',-10000.0},{'J',-10000.0},{'Z',4.0},{'X',1.0},{'*',-10000.0},{'-',-10000.0}}},
    {'X', {{'A',1.0},{'R',1.0},{'N',1.0},{'D',1.0},{'C',1.0},{'Q',1.0},{'E',1.0},{'G',1.0},{'H',1.0},{'I',1.0},{'L',1.0},{'K',1.0},{'M',1.0},{'F',1.0},{'P',1.0},{'S',1.0},{'T',1.0},{'W',1.0},{'Y',1.0},{'V',1.0},{'B',1.0},{'J',1.0},{'Z',1.0},{'X',1.0},{'*',-10000.0},{'-',-10000.0}}},
    {'*', {{'A',-10000.0},{'R',-10000.0},{'N',-10000.0},{'D',-10000.0},{'C',-10000.0},{'Q',-10000.0},{'E',-10000.0},{'G',-10000.0},{'H',-10000.0},{'I',-10000.0},{'L',-10000.0},{'K',-10000.0},{'M',-10000.0},{'F',-10000.0},{'P',-10000.0},{'S',-10000.0},{'T',-10000.0},{'W',-10000.0},{'Y',-10000.0},{'V',-10000.0},{'B',-10000.0},{'J',-10000.0},{'Z',-10000.0},{'X',-10000.0},{'*',1.0},{'-',-10000.0}}},
    {'-', {{'A',-10000.0},{'R',-10000.0},{'N',-10000.0},{'D',-10000.0},{'C',-10000.0},{'Q',-10000.0},{'E',-10000.0},{'G',-10000.0},{'H',-10000.0},{'I',-10000.0},{'L',-10000.0},{'K',-10000.0},{'M',-10000.0},{'F',-10000.0},{'P',-10000.0},{'S',-10000.0},{'T',-10000.0},{'W',-10000.0},{'Y',-10000.0},{'V',-10000.0},{'B',-10000.0},{'J',-10000.0},{'Z',-10000.0},{'X',-10000.0},{'*',-10000.0},{'-',1.0}}}
  } );
  
}//end fasta
}//end biocpp
#endif

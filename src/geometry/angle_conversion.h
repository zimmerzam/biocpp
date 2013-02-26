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

#ifndef ANGLE_CONVERSION_H
#define ANGLE_CONVERSION_H

#include <math.h>

#define PI 3.14159265358979323846

/*! \brief Convert degrees in radians

		\return the value in radians
		@param deg the angle in degrees */
inline double DegToRad(double deg){
	return deg*PI/180.;
}

/*! \brief Convert radians in degrees 
		\return the value in degrees
		@param rad the angle in radians */
inline double RadToDeg(double rad){
	return rad*180./PI;
}

#endif

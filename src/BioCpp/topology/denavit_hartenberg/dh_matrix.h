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

#ifndef BIOCPP_DENAVIT_HARTENBERG_MATRIX
#define BIOCPP_DENAVIT_HARTENBERG_MATRIX

#include <Eigen/Core>
#include <float.h>

namespace BioCpp{
namespace DenavitHartenberg{

Eigen::Matrix4d Matrix( double d, double r, double alpha, double theta ){
	Eigen::Matrix4d result;
	double sa = sin(alpha);
	double ca = cos(alpha);
	double st = sin(theta);
	double ct = cos(theta);
	result << ct , -st*ca , st*sa  , r*ct ,
	          st , ct*sa  , -ct*sa , r*st ,
	          0  , sa     , ca     , d    ,
	          0. , 0.     , 0.     , 1.   ;
	
	return result;
}

/*******************************************************/
/* Denavit-Hartenberg matrix associated with the joint */
/*******************************************************/
Eigen::Matrix4d MatrixZ( double theta, double d ){
	Eigen::Matrix4d result;
	double ct = cos(theta);
	double st = sin(theta);
	result << ct , -st , 0 , 0 ,
	          st ,  ct , 0 , 0 ,
	          0  ,  0  , 1 , d ,
	          0  ,  0  , 0 , 1 ;
	return result;
}

/* first derivative with respect to theta */
Eigen::Matrix4d MatrixZt( double theta, double d ){
	Eigen::Matrix4d result;
	double ct = cos(theta);
	double st = sin(theta);
	result << -st , -ct , 0 , 0 ,
	          ct  , -st , 0 , 0 ,
	          0   ,  0  , 0 , 0 ,
	          0   ,  0  , 0 , 0 ;
	return result;
}

/* second derivative with respect to theta */
Eigen::Matrix4d MatrixZtt( double theta, double d ){
	Eigen::Matrix4d result;
	double ct = cos(theta);
	double st = sin(theta);
	result << -ct ,  st , 0 , 0 ,
	          -st , -ct , 0 , 0 ,
	          0   ,  0  , 0 , 0 ,
	          0   ,  0  , 0 , 0 ;
	return result;
}

/******************************************************/
/* Denavit-Hartenberg matrix associated with the link */
/******************************************************/
Eigen::Matrix4d MatrixX( double alpha, double r ){
	Eigen::Matrix4d result;
	double ca = cos(alpha);
	double sa = sin(alpha);
	result << 1 , 0  , 0  , r ,
	          0 , ca , -sa, 0 ,
	          0 , sa , ca , 0 ,
	          0 , 0  , 0  , 1 ;
	return result; 
}

}
}
#endif

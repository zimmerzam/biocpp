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

#ifndef BIOCPP_DENAVIT_HARTENBERG_CHAIN
#define BIOCPP_DENAVIT_HARTENBERG_CHAIN

#include <Eigen/Core>
#include <Eigen/Dense>
#include <float.h>
#include <vector>
#include <iostream>

namespace BioCpp{
namespace DenavitHartenberg{

template <typename atom_t>
class chain{
  protected:
    int n_edges;
    std::vector< edge<atom_t> > edges;
    
    std::vector<double> v_alpha;
		std::vector<double> v_theta;
		std::vector<double> v_r;
		std::vector<double> v_d;
		
	public:
	  chain(){};
		chain(std::vector< edge<atom_t> >& links);
		
		double theta( typename std::vector< edge<atom_t> >::iterator& link);
		double alpha( typename std::vector< edge<atom_t> >::iterator& link);
		double d( typename std::vector< edge<atom_t> >::iterator& link);
		double r( typename std::vector< edge<atom_t> >::iterator& link);		
		
		int get_n_edges();
		double getAlpha(int i);
		double getTheta(int i);
		double getR(int i);
		double getD(int i);
		
		void setAlpha(int i, double el);
		void setTheta(int i, double el);
		void setR(int i, double el);
		void setD(int i, double el);
		
		Eigen::Matrix4d M( );
		template <class V> Eigen::Matrix4d M( V& theta );
		template <class V> Eigen::Matrix4d Mp(int i, V& theta);
};

template <typename atom_t>
chain<atom_t>::chain(std::vector< edge<atom_t> >& links){
	/* Copy the links */
	for( typename std::vector< edge<atom_t> >::iterator lk = links.begin(); lk != links.end(); ++lk){
		edges.push_back(*lk);
	}
	/* Store the number of links */
	n_edges=links.size();
	/* Set the base frame */
	edges.begin()->origin = edges.begin()->base;
	edges.begin()->xdir = Eigen::Vector3d(0,1,0).cross( edges.begin()->zdir );
	edges.begin()->xdir /= edges.begin()->xdir.norm();

	/* Set a coordinate system for other links */
	for( typename std::vector< edge<atom_t> >::iterator lk = edges.begin()+1; lk != edges.end(); ++lk){
		typename std::vector< edge<atom_t> >::iterator plk = lk-1;
		lk->computeFrame(*plk);
	}
	
	/* Compute DH parameters */
	for( typename std::vector< edge<atom_t> >::iterator lk = edges.begin()+1; lk < edges.end(); ++lk){
		v_alpha.push_back( alpha(lk) );
		v_theta.push_back( theta(lk) );
		v_d.push_back( d(lk) );
		v_r.push_back( r(lk) );
	}
}

template <typename atom_t>
int chain<atom_t>::get_n_edges(){
	return n_edges;
}

template <typename atom_t>
double chain<atom_t>::theta( typename std::vector< edge<atom_t> >::iterator& link){
	typename std::vector< edge<atom_t> >::iterator pr_link = link-1;
	Eigen::Vector3d cross = (pr_link->xdir).cross(link->xdir);
	double theta = fabs( cross.norm() ) < 1. ? asin( cross.norm() ) : asin(1.);
	if( pr_link->xdir.dot( link->xdir ) < 0 ) theta = M_PI-theta;
	if(cross.dot( pr_link->zdir ) > 0) return theta;
	else return -theta;
}

template <typename atom_t>
double chain<atom_t>::alpha( typename std::vector< edge<atom_t> >::iterator& link){
	typename std::vector< edge<atom_t> >::iterator pr_link = link-1;
	Eigen::Vector3d crossv = ( pr_link->zdir ).cross( link->zdir );
	double alpha = fabs( crossv.norm() ) < 1. ? asin( crossv.norm() ) : asin(1.);
	if( pr_link->zdir.dot( link->zdir ) < 0 ) alpha = M_PI-alpha;
	if( crossv.dot( link->xdir ) > 0 ) return alpha;
	else return -alpha;
}

template <typename atom_t>
double chain<atom_t>::d( typename std::vector< edge<atom_t> >::iterator& link){
	typename std::vector< edge<atom_t> >::iterator pr_link = link-1;
	return (link->origin - pr_link->origin).dot( pr_link->zdir );
}

template <typename atom_t>
double chain<atom_t>::r( typename std::vector< edge<atom_t> >::iterator& link){
	typename std::vector< edge<atom_t> >::iterator pr_link = link-1;
	double dot = (link->origin - pr_link->origin).dot( link->xdir );
//  edge<atom_t> tmp( link->origin, link->origin+link->xdir);
	return dot;
}

template <typename atom_t>
double chain<atom_t>::getAlpha(int i){
	if( i >=n_edges or i < 0){
		std::cout << "alpha is not defined for i=" << i << "." << std::endl;
		exit(1);
	}
	std::vector<double>::iterator it = v_alpha.begin() + i;
	return *it;
}

template <typename atom_t>
double chain<atom_t>::getTheta(int i){
	if( i >=n_edges or i < 0){
		std::cout << "theta is not defined for i=" << i << "." << std::endl;
		exit(1);
	}
	std::vector<double>::iterator it = v_theta.begin() + i;
	return *it;
}

template <typename atom_t>
double chain<atom_t>::getR(int i){
	if( i >=n_edges or i < 0){
		std::cout << "r is not defined for i=" << i << "." << std::endl;
		exit(1);
	}
	std::vector<double>::iterator it = v_r.begin() + i;
	return *it;
}

template <typename atom_t>
double chain<atom_t>::getD(int i){
	if( i >=n_edges or i < 0){
		std::cout << "d is not defined for i=" << i << "." << std::endl;
		exit(1);
	}
	if( i==(n_edges-1) ){
		Eigen::Vector3d vec = edges[i].target - edges[i].base;
		return vec.norm();
	}
	std::vector<double>::iterator it = v_d.begin() + i;
	return *it;
}

template <typename atom_t>		
void chain<atom_t>::setAlpha(int i, double el){
	if( i >=n_edges or i < 0){
		std::cout << "alpha is not defined for i=" << i << "." << std::endl;
		exit(1);
	}
	std::vector<double>::iterator it = v_alpha.begin() + i;
	*it = el;
}

template <typename atom_t>
void chain<atom_t>::setTheta(int i, double el){
	if( i >=n_edges or i < 0){
		std::cout << "theta is not defined for i=" << i << "." << std::endl;
		exit(1);
	}
	std::vector<double>::iterator it = v_theta.begin() + i;
	*it = el;
}

template <typename atom_t>
void chain<atom_t>::setR(int i, double el){
	if( i >=n_edges or i < 0){
		std::cout << "r is not defined for i=" << i << "." << std::endl;
		exit(1);
	}
	std::vector<double>::iterator it = v_r.begin() + i;
	*it = el;
}

template <typename atom_t>
void chain<atom_t>::setD(int i, double el){
	if( i >=n_edges or i < 0){
		std::cout << "d is not defined for i=" << i << "." << std::endl;
		exit(1);
	}
	std::vector<double>::iterator it = v_d.begin() + i;
	*it = el;
}


template <typename atom_t>
Eigen::Matrix4d chain<atom_t>::M( ){
	return M( v_theta );
}

template <typename atom_t>
template <class V>
Eigen::Matrix4d chain<atom_t>::M( V& theta ){
	int i;
	Eigen::Matrix4d result, Z, X;
	result.setIdentity();
	
	for(i = 0; i < n_edges -1; i++){
		Z = MatrixZ( theta[i], getD(i) );
		X = MatrixX( getAlpha(i), getR(i) );
		result*=Z;
		result*=X;
	}
	return result;
}

template <typename atom_t>
template <class V>
Eigen::Matrix4d chain<atom_t>::Mp( int i, V& theta ){
	if(i < 0 or i > n_edges - 1){
		std::cout << "cannot derive by theta_" << i << "." << std::endl;
		exit(1);
	}
	
	int k;
	Eigen::Matrix4d dhmtx, Z, X;
	dhmtx.setIdentity();
	for(k = 0; k < n_edges - 1 ; k++){
		if( k==i ){
			Z = MatrixZt( theta(k), getD(k) );
			X = MatrixX( getAlpha(k), getR(k) );		
		}
		else{
			Z = MatrixZ( theta(k), getD(k) );
			X = MatrixX( getAlpha(k), getR(k) );
		}
		dhmtx*=Z;
		dhmtx*=X;
	}

	return dhmtx;
}

}
}
#endif

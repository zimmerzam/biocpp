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

#ifndef BIOCPP_DENAVIT_HARTENBERG_EDGE
#define BIOCPP_DENAVIT_HARTENBERG_EDGE

#include <Eigen/Core>
#include <float.h>

namespace BioCpp{
namespace DenavitHartenberg{

template <typename atom_t>
class chain;

template <typename atom_t>
class edge{
  private:
    atom_t* atom_base;
    atom_t* atom_target;
  protected:
    Eigen::Vector3d base; 
		Eigen::Vector3d target; 
		Eigen::Vector3d origin;
		Eigen::Vector3d zdir;  
		Eigen::Vector3d xdir;  
  public:
    edge(){};
    void set( atom_t& base_atm, atom_t& target_atm );
    void set();
    
    Eigen::Vector3d computeFrame(edge<atom_t>& link);         /* compute new frame origin and xdir. */
	
		bool operator || (edge<atom_t>& link);        /* return true if link is parallel to *this */
		Eigen::Vector3d operator |= (edge<atom_t>& link);    /* return the intersection point or (MAX_DBL, MAX_DBL, MAX_DBL) if there is no intersection */
		
		friend class chain<atom_t>;
    
};

template <typename atom_t>
void edge<atom_t>::set( atom_t& base_atm, atom_t& target_atm ){   /* Constructor from two points */
	atom_base = &base_atm;
	atom_target = &target_atm;
  set();
}

template <typename atom_t>
void edge<atom_t>::set(  ){   /* Constructor from two points */
	this->base = atom_base->coordinate;
	this->target = atom_target->coordinate;
	this->zdir= target-base;
	this->zdir/=this->zdir.norm();
}

template <typename atom_t>
Eigen::Vector3d edge<atom_t>::computeFrame(edge<atom_t>& link){          /* Origin and x direction */
	if( (*this)|| link ){                         /* case: parallel links */
		this->origin = this->base;
		Eigen::Vector3d diff = this->base - link.base;
		Eigen::Vector3d par = link.zdir*( diff.dot( link.zdir ) );
		this->xdir = diff - par;
	}
	else if( ( (*this)|= link ).norm() < DBL_MAX/2. ){ /* case: intersecting links */
		this->origin = (*this)|= link;
		this->xdir = link.zdir.cross( this->zdir );
	}
	else{	                                         /* otherwise..	 */
		Eigen::Vector3d diff = link.base-this->base;
		double vec1_diff = diff.dot( link.zdir );
		double vec2_diff = diff.dot( this->zdir );
		double vec1_vec2 = link.zdir.dot( this->zdir );
		double mu1 = ( vec2_diff * vec1_vec2 - vec1_diff ) / 
  	             ( 1. - vec1_vec2 * vec1_vec2 );
		double mu2 = ( vec2_diff + vec1_vec2 * mu1);
		this->origin = this->base + this->zdir*mu2;
		Eigen::Vector3d pr_intersect = link.base + link.zdir*mu1;
		this->xdir = this->origin - pr_intersect;
	}
	this->xdir/=this->xdir.norm();
	return this->xdir;
}

template <typename atom_t>
bool edge<atom_t>::operator || (edge<atom_t>& link){         /* Parallel */
	if( link.zdir == this->zdir or link.zdir == -this->zdir )
		return true;
	else return false;
}

template <typename atom_t>
Eigen::Vector3d edge<atom_t>::operator |= (edge<atom_t>& link){     /* Intersection */
	if(this->base == link.target)
		return this->base;

	double EPS = DBL_EPSILON;

	Eigen::Vector3d crossv = this->zdir.cross( link.zdir );
	crossv /= crossv.norm();
	double distance = crossv.dot( this->base - link.base );

	if( fabs( distance ) < EPS ){
		Eigen::Vector3d R = link.base - this->base;
		double Rb = R.dot( link.zdir );
		double Ra = R.dot( this->zdir );
		double ab = link.zdir.dot( this->zdir );
		double mu = ( Rb - Ra*ab )/( ab*ab - 1. );
		return link.base + link.zdir*mu; 
	}
	else return Eigen::Vector3d( DBL_MAX, DBL_MAX, DBL_MAX );
}

}
}

#endif

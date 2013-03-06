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

#include "../geometry/Eigen/Core"

template<typename T>
class has_coordinate{
    struct Fallback { Eigen::Vector3d coordinate; };
    struct Derived : T, Fallback { };
 
    template<typename U, U> struct Check;
 
    typedef char ArrayOfOne[1];
    typedef char ArrayOfTwo[2];
 
    template<typename U> 
    static ArrayOfOne & func(Check<Eigen::Vector3d Fallback::*, &U::coordinate> *);
 
    template<typename U> 
    static ArrayOfTwo & func(...);
 
  public:
    typedef has_coordinate type;
    enum { value = sizeof(func<Derived>(0)) == 2 };
};

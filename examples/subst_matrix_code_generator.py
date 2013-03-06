# ########################################################################## #
#                                                                            #
#    Copyright 2013 Stefano Zamuner                                          #
#                                                                            #
#    This file is part of BioCpp.                                            #
#                                                                            #
#    BioCpp is free software: you can redistribute it and/or modify          #
#    it under the terms of the GNU General Public License as published by    #
#    the Free Software Foundation, either version 3 of the License, or       #
#    (at your option) any later version.                                     #
#                                                                            #
#    BioCpp is distributed in the hope that it will be useful,               #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of          #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           #
#    GNU General Public License for more details.                            #
#                                                                            #
#    You should have received a copy of the GNU General Public License       #
#    along with BioCpp.  If not, see <http://www.gnu.org/licenses/>.         #
#                                                                            #
# ########################################################################## #

import sys

name = sys.argv[1]
filename = sys.argv[2]

file = open(filename,'r')
chars = []
print "substitution_matrix "+name+"( map_list_of_type< char, std::map< char, double > >",
for idx,line in enumerate(file):
  if(idx==0):
    chars = line.split(";")
  else:
    token = line.split(";")
    first = token[0][1:2]
    print "('" + first + "', map_list_of_type< char, double >",
    for index, value in enumerate(token[1:]):
      second = chars[index+1][1:2]
      print "('" + second + "'," + str(float(value)) + ")",
    print " )"
print ");"

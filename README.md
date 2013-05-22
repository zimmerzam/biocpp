# BioCpp #

This is a simple template library for basic operations on proteins. 

# System requirements #

+ Operating system: GNU/Linux
+ C++ compiler: gcc4.6 (or newer) or equivalent

**BioCpp** is C++0x standard.

## Installation and setup ##

BioCpp uses Eigen3 for atom coordinates, geometrical transformations and this 
kind of stuffs. So you need to download itfrom [here](http://eigen.tuxfamily.org/). 
Eigen3 needs no istallation: simply extract the sources in 
*path-to-BioCpp/src/geometry/Eigen/* and you're done.  

In order to use BioCpp in your project simply

+ add `export CPLUS_INCLUDE_PATH=$CPLUS_INCLUDE_PATH:/path-to/BioCpp/src/` 
to your `.bashrc` and/or `.bash_profile` file
+ import BioCpp header file

```c++
#include <BioCpp.h>
```
or
```c++
#include <BioCpp_default.h>
```

## Documentation ##

Documentation is available in both [html](http://biocpp.zimlotech.com/html/) and 
[pdf](http://biocpp.zimlotech.com/pdf/refman.pdf) format.


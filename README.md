# BioCpp #

This is a simple template library for basic operations on proteins. 

# System requirements #

+ Operating system: GNU/Linux
+ C++ compiler: gcc4.6 (or newer) or equivalent

**BioCpp** is C++0x standard.

## Installation and setup ##

BioCpp uses Eigen3 for atom coordinates, geometrical transformations and this 
kind of stuffs. So you need to download it from [here](http://eigen.tuxfamily.org/).  

Eigen3 needs no installation: simply extract the sources in 
*path-to-BioCpp/src/geometry/Eigen/* and you're done.  

In order to use BioCpp in your project simply

+ add the path containing the source code to your CPLUS_INCLUDE_PATH. In order 
to do this you can modify the `C_INC_PATH` flag to point to your 
`/path-to/BioCpp/src/`. Alternatively you can export the path at login by adding

     
```bash
$ export CPLUS_INCLUDE_PATH=$CPLUS_INCLUDE_PATH:/path-to/BioCpp/src/
``` 
to your `.bashrc` and/or `.bash_profile` file if you are using *bash* or *sh* 
shell or, alternatively

```bash
$ if ( $?CPLUS_INCLUDE_PATH ) then
$   setenv CPLUS_INCLUDE_PATH {$CPLUS_INCLUDE_PATH}:/path-to/BioCpp/src/
$ else
$   setenv CPLUS_INCLUDE_PATH /path-to/BioCpp/src/
$ endif
```
to your `.cshrc` file if you are using *tcsh* or *csh*.

+ import BioCpp header file
```c++
#include <BioCpp.h>
```
or
```c++
#include <BioCpp_standard.h>
```

## Documentation ##

Documentation is available in both [html](http://biocpp.zimlotech.com/html/) and 
[pdf](http://biocpp.zimlotech.com/pdf/refman.pdf) format.

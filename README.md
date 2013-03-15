#BioCpp#

This is a simple template library for basic operations on proteins. 

##Where to start##
BioCpp uses Eigen3 for atom coordinates, geometrical transformations and this kind of stuffs. So you need to download it
from [here](http://eigen.tuxfamily.org/). Eigen3 needs no istallation: simply extract the sources in BioCpp/src/geometry/Eigen/
and you're done.

In order to use BioCpp in your project simply include BioCpp.h

```c++
#include "/path-to-BioCpp-folder/src/BioCpp.h"
```
or
```c++
#include "/path-to-BioCpp-folder/src/BioCpp_default.h"
```

You need gcc4.6 or newer in order to be able to compile it!

##Documentation##
Looking for detailed documentation or examples? Please visit the documentation page at [biocpp-doc](http://biocpp.zimlotech.com/).

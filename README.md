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
export CPLUS_INCLUDE_PATH=$CPLUS_INCLUDE_PATH:/path-to/BioCpp/src/
``` 
to your `.bashrc` and/or `.bash_profile` file if you are using *bash* or *sh* 
shell or, alternatively,
```bash
if ( $?CPLUS_INCLUDE_PATH ) then
  setenv CPLUS_INCLUDE_PATH {$CPLUS_INCLUDE_PATH}:/path-to/BioCpp/src/
else
  setenv CPLUS_INCLUDE_PATH /path-to/BioCpp/src/
endif
```
to your `.cshrc` file if you are using *tcsh* or *csh*.

+ import BioCpp header file
```c++
#include <BioCpp.h>
```

### Preprocessor directives ###

Including the whole library in your project is usually not necessary. This is why
*BioCpp.h* provides only few basic features. If your project requires more features, 
add them by using one (or more) of the following *directives*:  

| Directive                                       | Description                                                    |
| :---------------------------------------------: | :------------------------------------------------------------: |
| BIOCPP\_INCLUDE\_ID                             | Add standard identifiers for amino acids, atoms, elements, ... |
| BIOCPP\_INCLUDE\_DPSS                           | Add definitions of secondary stuctures according to dssp       |
| BIOCPP\_INCLUDE\_MORPHOLOGY                     | Compute solvent-accessible surface area                        |
| BIOCPP\_INCLUDE\_RECONSTRUCTION                 | Compute position of missing atoms                              |
| BIOCPP\_INCLUDE\_PDB                            | Read and write *pdb* files                                     |
| BIOCPP\_INCLUDE\_FASTA_ALIGN                    | Sequence alignment using *NeedlemanWunsch* algorithm           |
| BIOCPP\_INCLUDE\_FASTA\_SUBST\_MATRIX\_PAM30    | Use PAM30 substitution matrix                                  |
| BIOCPP\_INCLUDE\_FASTA\_SUBST\_MATRIX\_PAM70    | Use PAM70 substitution matrix                                  |
| BIOCPP\_INCLUDE\_FASTA\_SUBST\_MATRIX\_BLOSUM45 | Use BLOSUM45 substitution matrix                               |
| BIOCPP\_INCLUDE\_FASTA\_SUBST\_MATRIX\_BLOSUM62 | Use BLOSUM62 substitution matrix                               |
| BIOCPP\_INCLUDE\_FASTA\_SUBST\_MATRIX\_BLOSUM80 | Use BLOSUM80 substitution matrix                               |
| BIOCPP\_INCLUDE\_FASTA\_SUBST\_MATRIX\_ZIMM1    | Use a strict version of BLOSUM62                               |
| BIOCPP\_INCLUDE\_FASTA\_SUBST\_MATRIX\_MATRIX   | Use all the substitution matrices                              |
| BIOCPP\_INCLUDE\_STANDARD                       | Add standard definitions of residue, chain, complex and moiety |
| BIOCPP\_INCLUDE\_SUBST\_MATRIX\_ALL             | Add all the features provided by *BioCpp*                      |

By default, if no directives are defined, only *base_container* and *iterators* 
are included (see [Documentation][] and examples).  

In order to add *directives* you can both  
```c++
#define DIRECTIVE_NAME
```  
**before** including *BioCpp.h*, **or** pass the directive name to your compiler.
If you use *g++*, simply  
```bash
g++ -D DIRECTIVE_NAME ...
```  
(see Makefile provided in the *examples* folder).

## Warning ##

This is a **developing code**. Newer versions of BioCpp may not be compatible with
older ones.  
If you distribute code based on BioCpp, please also provide the version of 
BioCpp on which it is based.  

## Documentation ##

Documentation is available in both [html](http://biocpp.zimlotech.com/html/) and 
[pdf](http://biocpp.zimlotech.com/pdf/refman.pdf) format.

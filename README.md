## The ANETO Library

Project homepage: https://github.com/DSantosO/anetolib

Copyright (c) 2017 Daniel Santos-Oliván and Carlos F. Sopuerta

__ANETO__ is under GNU General Public License ("GPL").

The __ANETO__ (Arbitrary precisioN solvEr with pseudo-specTral methOds) library's
main purpose is to provide a tool to perform simple one-dimensional problems
or evolution ones with the Method of lines with __Arbitrary Precision__ using
__Pseudo-Spectral Collocation methods__ (PSC).

The exponential convergence of spectral methods makes them the best option to go beyond
the standard double precision (64-bits).
For this precision, maximum accuracy is usually reached with a very low number of
discretization points so it is not that computational expensive to go further.

The classes implemented here can be used with arbitrary data types allowing us to
control the accuracy of our numerical computations until the level we need / can
computationally afford.

A paper with a basic review of PSC methods and showing some of the possibilities
of the library will be published soon.

For comments, questions or suggestions about new functionalities, the user can
contact anetolib@gmail.com

The authors thank [Lluís Gesa](https://github.com/Hlod-Wig) and Víctor Martín 
for suggestions and improvements in the code of the library.



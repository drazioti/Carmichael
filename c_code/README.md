#  C/C++ implementation of Product Subset Problem and application for finding Carmichael numbers

![GPLv2][]

[GPLv2]: https://img.shields.io/badge/license-GPLv2-lightgrey.svg

Initial Author : E. Tiganourias (etiganou@csd.auth.gr)

We apply this problem to find Carmichael numbers by usinng Erdos algorithm.

Applications of MSPP to cryptography see https://github.com/drazioti/NSK-birthday-attack

## Requirements
Linux

GMP

## Compile
```
g++ carmi1.cpp Combinations3.cpp -lgmpxx -lgmp
```
## Contribute
First fork this repository. Make the changes you want (e.g. update some tables, correct a bug to the code etc)

Contribute by using pull request to this repo. 

## TODO
- Add hash-trick

- Add c/c++ code

- Add c/c++/openmp  code for the parallel attack
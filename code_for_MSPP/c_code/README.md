#  C/C++ implementation of Product Subset Problem and application for finding Carmichael numbers with many primes

![GPLv2][]

[GPLv2]: https://img.shields.io/badge/license-GPLv2-lightgrey.svg

Initial Author 						  : E. Tiganourias (etiganou@csd.auth.gr)

Refactoring (openmp and exponents)    : V. Martidis    (vamartid@yandex.com)

We apply this problem to find Carmichael numbers by usinng Erdos algorithm.

Applications of MSPP to cryptography see https://github.com/drazioti/NSK-birthday-attack

## Requirements
Linux/GMP/OpenSSL

## Compile
```
g++ carmi1.cpp Combinations3.cpp -lgmpxx -lgmp -lcrypto
```

## Contribute
First fork this repository. Make the changes you want (e.g. update some tables, correct a bug to the code etc)

Contribute by using pull request to this repo. 

## TODO
- Add hash

- Add c/c++/openmp  code for the parallel attack

## Issues

- unordered map?
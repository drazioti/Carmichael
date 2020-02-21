#  C/C++ parallel implementation of Product Subset Problem and application for finding Carmichael numbers with many primes

![GPLv2][]

[GPLv2]: https://img.shields.io/badge/license-GPLv2-lightgrey.svg

Initial Author 						  : E. Tiganourias (etiganou@csd.auth.gr)

Refactoring (openmp,exponents,args)    : V. Martidis    (vamartid@yandex.com)

We apply this problem to find Carmichael numbers by usinng Erdos algorithm.

Applications of MSPP to cryptography see https://github.com/drazioti/NSK-birthday-attack

## Requirements
Linux/gmp/OpenSSL/Openmp/args

For the library args :
```
$git clone https://github.com/Taywee/args.git
$cd args
$sudo make install
```

## Compile
Makefile
```
make all #results carmichael.out (serial) and carmicael_par.out(parallel)
make serial # results carmichael.out (serial) and compiled files
make parallel # carmicael_par.out(parallel) and compiled files
...
```

## Run
```
./a.out L --ham 15 -b 32 -f 1
```

Where      
          

	   -b is the cardinality of the sets I_1 and I_2
	   
	   -f is the fragmentation parameter
	   
	   -q is the length of hash stored

	  --ham is the local hamming weight (we choose h1,h2 such that h1+h2=ham)

	    L is written in the form a1 a2 ... ar where aj are the expoents for construction Lambda

E.g.
```
./a.out 20 5 4 1 1 --ham 15 -b 32 -f 1 -q 10
```
Note that you have to hard-code the parameter Q in the carmi.cpp, main().

## Contribute
First fork this repository. Make the changes you want (e.g. update some tables, correct a bug to the code etc)

Contribute by using pull request to this repo. 

## TODO
- add comments
- Find a Carmichael number with many prime factors!


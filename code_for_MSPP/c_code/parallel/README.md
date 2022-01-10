#  C/C++ parallel implementation of Product Subset Problem and application for finding Carmichael numbers with many primes

![GPLv2][]

[GPLv2]: https://img.shields.io/badge/license-GPLv2-lightgrey.svg

Initial Author 						  : E. Tiganourias (etiganou@csd.auth.gr)

Refactoring (openmp,exponents,args)    : V. Martidis    (vamartid@yandex.com)

We apply this problem to find Carmichael numbers by using Erdos algorithm.

Applications of MSPP to cryptography see https://github.com/drazioti/NSK-birthday-attack

## Requirements
Linux/gcc/g++/make
and gmp/OpenSSL/Openmp/args

In Debian based systems, for gcc
```
$sudo apt-get install build-essential
$sudo apt-get install manpages-dev
```

Except openssl you will need
```
$sudo apt-get install libssl-dev
```

For gmp in Debian systems you can use
```
$sudo apt-get install libgmp3-dev
```

For the library args :
```
$git clone https://github.com/Taywee/args.git
$cd args
$sudo make install
```

## Compile
```
$make all
```

### run single core
```
$./carmichael.out 20 5 4 1 1 --ham 15 -b 32 -f 1 -q 10
```

### run multi core
```
$./carmichael_par.out L --ham 15 -b 32 -f 1 -q 10
```
Where,      
```
	   L is written in the form a1 a2 ... ar where aj are the exponents for construction of Lambda (a1>=a2>=...>=ar)

	   --ham is the local hamming weight (we choose h1,h2 such that h1+h2=ham)

	   -b is the cardinality of the sets I_1 and I_2 (|I_1|=|I_2|=b)
	   
	   -f is the fragmentation parameter
	   
	   -q is the length of hash stored
```

## Contribute
First fork this repository. Make the changes you want (e.g. update some tables, correct a bug to the code etc)

Contribute by using pull request to this repo. 

## TODO
- add comments
- Find a Carmichael number with many prime factors!


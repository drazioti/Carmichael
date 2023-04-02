#  Modular Subset Product Problem (MSPP)

![GPLv2][]

[GPLv2]: https://img.shields.io/badge/license-GPLv2-lightgrey.svg

In this repository we provide some code for solving Modular Subset Product Problem.<br>
We apply this problem to find Carmichael numbers by using Erdos algorithm (which uses MSPP), see section 3 of [paper](https://arxiv.org/abs/2002.07095). <br>
Applications of MSPP to cryptography see [here](https://github.com/drazioti/NSK-birthday-attack)<br>
Published as:<br>
*K. A. Draziotis, V. Martidis and S. Tiganourias, Product Subset Problem : Applications to number theory and cryptography, Book Chapter, Analysis, Cryptography and Information Science, Chapter 5, Vol. 10, World Scientific, 2023.* [https://www.worldscientific.com/worl...](https://www.worldscientific.com/worldscibooks/10.1142/13296#t=aboutBook)

## The code
We provide Sagemath code, which we used to build all the (small instances) tables.

We also provide c++ single and multi core variants. Using the single core case we managed
to produce a Carmichael number with 19589 prime factors in a I5/16Gb Linux PC in three hours.


## Tables for Carmichael Numbers
We provide tables for Carmichael numbers in directory Tables.

For instance,

*5-80.txt* contains Carmichael numbers with 5 to 80 prime factors.

*2542.txt* contains a Carmichael number with 2542 prime factors.

## Get the tables
You can download the source code with the tables by using

```sh
git clone https://github.com/drazioti/Carmichael.git
```

## Contribute
First fork this repository. Make the changes you want (e.g. update some tables, correct a bug to the code etc)

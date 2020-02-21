// $g++ --std=c++11 carmi.cpp Combinations.cpp -lgmpxx -lgmp -fopenmp -lcrypto
// $nohup ./a.out > script.out 2>&1 &
#ifndef carmi_H
#define carmi_H

#include <fcntl.h> // open
// #include <unistd.h>
// #include <stdio.h>
// #include <stdlib.h>
#include <iostream>
#include <new>
#include <math.h>
#include <time.h>
#include <signal.h>
#include <list>
#include <algorithm>
#include <ctime>
// #include <fstream>
#include <args.hxx> // argument lib args by Taywee
#include <gmpxx.h> // gmp
#include "subset_product.h"
using namespace std;

//BASIC TEMPORARY UTILITY

void from_list_to_array(std::list<unsigned char*> &source, unsigned char** &dest, int r);

//FUNCTION (2): CALCULATING the modulus LAMBDA

void Lambda(int* Q,int* H,int r, mpz_class &Lambda);

//FUNCTION (3): CHOOSING FIRST r PRIMES FOR Q SET

void set_Q(int r,int* &Q);

//FUNCTION (4): FUNCTIONS FOR MAKING THE P SET

void dup_array(int* A,int* B,int size, int** &out);
void multiply_list(int* H,int r, mpz_class &result);
void divisors(int* P, int* H, int r,mpz_class &size, unsigned char** &divrep);
void make_P_set(int* Q, int* H,int r,mpz_class &Lambda, std::list<unsigned char*> &P);

//THIS CODE ABOVE IS USED FOR THE OPTIMIZED VERSION OF STORING THE
//THE P SET

//FUNCTION(6): MAKE T_SET

int T_set(int* Q, int r,unsigned char** &P, mpz_class &sizeP, mpz_class &Lambda, mpz_class** &I, mpz_class &sizeI,int local_hamming_weight, double total_time, int frag, char Q_bytes);

//FUNCTION(7): PRODUCING THE CARMICHAEL NUMBER FROM THE SOL1,SOL2 SETS
//	       THIS FUNCTION IMPLEMENTS THE S SET OF ERDOS ALGORITHM	

void extract_number(mpz_class* &P, mpz_class** &I,mpz_class &sizeI, std::list<int*> &sol1, std::list<int*> &sol2, int h1, int h2, mpz_class*  &numbers);
void euler_totient(int* Q, int* H, int size, mpz_class &result);
void density(int* Q, int* H, int size, int Psize, double &result1, double &density1 );

#endif
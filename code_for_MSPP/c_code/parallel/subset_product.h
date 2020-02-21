
#ifndef subset_product_H
#define subset_product_H

// #include <stdio.h>
// #include <stdlib.h>
#include <iostream> // cout
#include <string.h>
#include <algorithm>
#include <vector>
#include <list>
#include <fstream> // file access
#include "unordered_map" // hash table
#include <openssl/md5.h> // open ssl md5
#include <sys/resource.h> // definition for XSI resource operations
#include <sys/time.h> // time types
#include <sys/types.h> // data types
#include <unistd.h> // read, close
// #include <inttypes.h>
#include <signal.h> // raise(SIGINT)
#include <gmpxx.h> // gmp
#include <omp.h> // openmp
#include "Combinations.h"
using namespace std;


//############################BASIC FUNCTIONS###########################

//FUNCTIONS (1): IMPLAMENTATION OF SIMPLE USEFULL FUNCTONS FOR MPZ OBJECTS
unsigned long long mpz_2_ull(mpz_class z);
int chrcmp(const char chr1, const char chr2);
void arrcpy(mpz_class* &dest, mpz_class* &source, mpz_class &size);
mpz_class fact(mpz_class &n);
void comb_size(mpz_class &n, mpz_class k, mpz_class &size);

//FUNCTION (1.2): MEMORY CALCULATION FOR UPPER BOUND IN U1

int getmem();

//FUNCTION (1.3): MAKING BOUND TO DIVIDE HASHMAP
//NEEDS DEVELOPMENT to calcuate on Runtime optimal bound

void U_bound(mpz_class size, mpz_class &bound,int frag);

//FUNCTION (2): MD5 HASH FUNCTION FOR MPZ_CLASS OBJECTS

string to_md5_f6_str(mpz_class number,char Q_bytes);

//FUCNTION (3): IS_CARMICHAEL FUNCTION TO DEAL WITH UNWANTED COLLISIONS

bool is_intersection(mpz_class &Lambda, int* Q_s, int r, unsigned char** &P, mpz_class &u2, int* &E1, int h1, mpz_class* &I, mpz_class &sizeI, mpz_class &c);
int is_carmichael(mpz_class &n, mpz_class* &factors, mpz_class &sizef);

//FUCNTION (4): CODE FOR "func1"

mpz_class func1(mpz_class* &P, int* E, mpz_class &Lambda, mpz_class &c, int flag, int size);

//FUNCTION (5): CODE FOR "func2"

void func2(mpz_class* &P, int** E, mpz_class &Lambda, mpz_class &c, int flag, unordered_multimap <string,int*> &Map, mpz_class &sizeI, int h1, mpz_class &begin, mpz_class &end, char hash_flag, char Q_bytes);

//FUNCTION (6): FUNCTIONS FOR U1

mpz_class get_P_element(int* &Q, unsigned char* &H, int r);
void U1(int* &Q_s, int r, unsigned char** &P, mpz_class* &I, int** E, mpz_class &Lambda,mpz_class &sizeI,int h1, unordered_multimap<string, int*> &Map, mpz_class &begin, mpz_class &end, char hash_flag, char Q_bytes);

//FUNCTION (7): FUNCTION FOR SAVING SUBSET AFTER FINDING INTERSECTION

void sol(std::list<int*> &sol1, std::list<int*> &sol2, int* E1, int* E2, int h1, int h2);

//FUNCTION (8): INTERSECTION FUNCTIONS aka HASHMAP SEARCH

int intersection(int* &Q_s, int r, unsigned char** &P,mpz_class &sizeP, mpz_class** &I, int* E, mpz_class &Lambda, mpz_class &c, mpz_class &sizeE, mpz_class &sizeI, unordered_multimap <string, int*> &Map, mpz_class &count, std::list<int*> &sol1, std::list<int*> &sol2,int h1, int h2, double total_time, char hash_flag, char Q_bytes);

//FUNCTION (8): GENERATE I SET FUNCTION

void randomize_I(gmp_randclass &rr, unsigned int seed);
void gen_I(mpz_class &n, mpz_class &b, int flag, mpz_class** &I, gmp_randclass &rr, unsigned int seed);

//FUNCTION(9): PRODUCT SUBSET ATTACK PHASE 1

int product_attack_1(int* &Q, int r, unsigned char** &P,mpz_class &sizeP, mpz_class &Lambda, mpz_class &c, mpz_class** &I,int local_hamming_weight,mpz_class &sizeI, std::list<int*> &sol1, std::list<int*> &sol2, mpz_class &count, double total_time, int frag, char Q_bytes);

#endif
// $g++ --std=c++11 carmi.cpp Combinations3.cpp -lgmpxx -lgmp -fopenmp -lcrypto
// $nohup ./a.out > script.out 2>&1 &

#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <gmpxx.h>
#include "iostream"	
#include <new>
#include <math.h>
#include <time.h>
#include <signal.h>
#include "subset_product.cpp"
#include <list>
#include <algorithm>
#include <ctime>
#include <fstream>
using namespace std;

//BASIC TEMPORARY UTILITY

void from_list_to_array(std::list<unsigned int*> &source, unsigned int** &dest, int r){
	unsigned long j=0;
	for(std::list<unsigned int*>::iterator it=source.begin();it != source.end(); ++it){
		for (unsigned long i=0;i<r;i++) 
			dest[j][i] = *(*it+i);
		j++;
        }
}

//FUNCTION (2): CALCULATING the modulus LAMBDA

void Lambda(int* Q,int* H,int r, mpz_class &Lambda){
	
	mpz_class exp;
	for (int i=0;i<r;i++){
		exp = pow(Q[i], H[i]); 
		Lambda *= exp;
	}
	return;
}

//FUNCTION (3): CHOOSING FIRST r PRIMES FOR Q SET

void set_Q(int r,int* &Q){
	if (r<3){
		cout << "At least 3 primes factors" << endl;
		return ;	
	}
	int q=3;
	Q[0] = 2;
	//cout << "DEBUGGIN THE Q SET FUNCTION" << endl;
	mpz_class temp;
	for (int i=1;i<r;i++){
		temp = q;
                while (mpz_probab_prime_p(temp.get_mpz_t(), 5) != 1 && mpz_probab_prime_p(temp.get_mpz_t(), 5) != 2)
                {
			q+=2;
			temp = q;
                }	       	
                Q[i] = q;
                q +=2;
        }
}

//FUNCTION (4): FUNCTIONS FOR MAKING THE P SET

void dup_array(int* A,int* B,int size, int** &out){

        for (int i=0;i<size;i++){
               // cout << "P["<<i<<"] is: " << A[i] << endl;
               // cout << "H["<<i<<"] is: " << B[i] << endl;
                out[0][i] = A[i];
                out[1][i] = B[i];
        }
}

void multiply_list(int* H,int r, mpz_class &result){
      
        for (int i=0;i<r;i++){
                result *= (H[i]+1);
        }
	return;
}

void divisors(int* P, int* H, int r,mpz_class &size, unsigned int** &divrep){
		// size to ull
		unsigned long long size_ull=mpz_2_ull(size);

        //getting the 2d array needed
        int** dup;
        dup = new int*[2];
        dup[0] = new int[size_ull];
        dup[1] = new int[size_ull];
        dup_array(P,H,r,dup);

        //-------------------------
		mpz_class* divs = new mpz_class[size_ull];
        divs[0] = 1;
        
		// create divisor array with "size" rows and "r" columns
		// unsigned int** divrep; 			// init pointer
		divrep = new unsigned int*[r];	// allocate columns
		for(int i=0;i<r;i++){
			// for each column allocate "size" rows
			divrep[i] = new unsigned int[size_ull];
		}
		// in divrep array has "size" elements
		// each element is represented by a row of "r" values
		// each value of row is the exponent of equivalent P element
		// so if P is [2,3,5] and the first row of divreps is [0,2,1]
		// the element represented is pow(2,0)*pow(3,2)*pow(5,1) = 45

		unsigned long long out_count =1;
        for (int i=0;i<r;i++){
                mpz_class* prev;
                prev = new mpz_class[size_ull];
                mpz_class pn = 1;
                arrcpy(prev, divs, size);
                for (int j=0;j<dup[1][i];j++){
                        pn *= dup[0][i];
                        unsigned long long k = 0;
                        while (prev[k] != 0 && k<size_ull){
                                divs[out_count] = prev[k] * pn;
                                // create a copy of the last multiplied element
								for(int y=0;y<r;y++){
									divrep[y][out_count]=divrep[y][k];
								}
								// set the new power according to what we multiplied 
								divrep[i][out_count]=1+j;
								out_count++;
                                k++;
                        }
                }
		delete[] prev;
        }

	delete[] divs;
	delete[] dup[0];
	delete[] dup[1];
	delete[] dup;
}

void make_P_set(int* Q, int* H,int r,mpz_class &Lambda, std::list<unsigned int*> &P){
    mpz_class size_d=1; 
	multiply_list(H,r,size_d);
        unsigned int** divs;			// init pointer (divs is divrep now)
        divisors(Q, H, r, size_d, divs);
		unsigned long long size_d_ull = mpz_2_ull(size_d);
		// cout << " SIZE OF P SET " << size_d_ull << " " << size << endl;
        for (unsigned long long i=0;i<size_d_ull;i++){
			mpz_class num = 1;
			// create the number
			for(int j=0;j<r;j++){
				num = num * pow(Q[j], divs[j][i]) ;
			}
			num+=1; // increase the number to be odd 
			if (mpz_probab_prime_p(num.get_mpz_t(), 5) !=0 && Lambda%num!=0 && num<Lambda){
				unsigned int* numrep = new unsigned int[r];
				for(int k=0;k<r;k++){
					numrep[k]=divs[k][i];
				}
				P.insert(P.end(), numrep);
			}
        }
	delete[] divs;
	return;
}

//THIS CODE ABOVE IS USED FOR THE OPTIMIZED VERSION OF STORING THE
//THE P SET

int is_carmichael(mpz_class &n, mpz_class* &factors, mpz_class &sizef);

//FUNCTION(6): MAKE T_SET

int T_set(int* Q, int r,unsigned int** &P, mpz_class &sizeP, mpz_class &Lambda, mpz_class** &I, mpz_class &sizeI,int local_hamming_weight, double total_time){
	mpz_class b=1;
	for (unsigned long long i=0;i<sizeP.get_ui();i++)
	{
		b*= get_P_element(Q, P[i], r);
		mpz_mod(b.get_mpz_t(), b.get_mpz_t(), Lambda.get_mpz_t());
	}
	mpz_class count=0;
	std::list<int*> sol1;
	std::list<int*> sol2;
	product_attack_1(Q, r, P,sizeP,Lambda, b, I, local_hamming_weight,sizeI, sol1, sol2, count, total_time);
	
	cout << "T set product attack finished " << endl;	
	if(sol1.begin() != sol1.end() && sol2.begin() != sol2.end())
	{
		cout << "Found " << count << " carmichael numbers with " << sizeP - local_hamming_weight << " factors" << endl;
		sol1.clear();
        	sol2.clear();
		return 1;
	}
	sol1.clear();
	sol2.clear();
	return 0;
}	

//FUNCTION(7): PRODUCING THE CARMICHAEL NUMBER FROM THE SOL1,SOL2 SETS
//	       THIS FUNCTION IMPLEMENTS THE S SET OF ERDOS ALGORITHM	

void extract_number(mpz_class* &P, mpz_class** &I,mpz_class &sizeI, std::list<int*> &sol1, std::list<int*> &sol2, int h1, int h2, mpz_class*  &numbers){
	mpz_class* subset1;
	subset1 = new mpz_class[mpz_get_ui(sizeI.get_mpz_t())];
	mpz_class* subset2;
	subset2 = new mpz_class[mpz_get_ui(sizeI.get_mpz_t())];	 //SUBSETS OF P BASED ON I1,I2
	
	for (mpz_class i=0;i<sizeI;++i)
	{
		mpz_class index = I[0][mpz_get_ui(i.get_mpz_t())];
		subset1[mpz_get_ui(i.get_mpz_t())] = P[mpz_get_ui(index.get_mpz_t())];
		index = I[1][mpz_get_ui(i.get_mpz_t())];
		subset2[mpz_get_ui(i.get_mpz_t())] = P[mpz_get_ui(index.get_mpz_t())];
	}	
	std::list<int*>::iterator it1 = sol1.begin();
	std::list<int*>::iterator it2 = sol2.begin();
	cout << endl << "--------------"<< endl;
	//cout << "---DEBUG: SUBSETS DONE --- " << endl;
	mpz_class counter=0;

	for(;it1 != sol1.end() && it2!=sol2.end() && counter<sol1.size(); ++it1, ++it2, ++counter)
	{

//Prime FACTORS FOR the CARMICHAEL FUNCTION

		mpz_class* factors;
		mpz_class fsize = h1 + h2;
		factors = new mpz_class[mpz_get_ui(fsize.get_mpz_t())];
		int fcount =0;
		mpz_class number=1;

		for (int i =0;i<h1;++i)
		{
			int index = *(*it1 + i);
			//cout << index << ' ';
			number *= subset1[index];
			factors[fcount] = subset1[index];
			fcount++;
		}
		for (int i=0;i<h2;++i)
		{
			int index = *(*it2 + i);
			//cout << index << ' ';
			number *= subset2[index];
			factors[fcount] = subset1[index];
			fcount++;
		}
			if (is_carmichael(number,factors,fsize)==1)
			{
			numbers[mpz_get_ui(counter.get_mpz_t())] = number;
			cout << number << " is a carmichael number"<<endl;
		}
	}	
	cout << "SET S IS COMPLETE" << endl;
}

void euler_totient(int* Q, int* H, int size, mpz_class &result)
// we compute euler function but we know the factorization.
{    
	mpz_class exp1;
	mpz_class exp2;
	for(int i=0;i<size;i++){
		exp1 =  pow(Q[i], H[i]);
		exp2 =  pow(Q[i], H[i]-1);
		result *= (exp1-exp2);
	}
	return;	
}


void density(int* Q, int* H, int size, int Psize, double &result1, double &density1 )
// we compute euler function but we know the factorization.
// We also compute the density of the Product subset problem.
{
	mpz_class exp1;
	mpz_class exp2;
	mpz_class diff;
	double difflogsum ;
	double difflog = 0;
	for(int i=0;i<size;i++){
		exp1 =  pow(Q[i], H[i]);
		exp2 =  pow(Q[i], H[i]-1);
		diff = exp1 - exp2;
		difflog = log2( mpz_get_d(diff.get_mpz_t()));
		difflogsum +=difflog;
	}
	result1 = difflogsum;
    density1 = Psize/result1;
	return;
}


//--------------------------------------------------------------//

int main(){
	clock_t startP, endP;
	mpz_class L=1;
	
//-----CHANGE THESE PARAMETERS TO RUN a new instance-------//	
// Also there is the parameter frag in subset_product.cpp //
	
	int r=5;		    //number of first primes
	int hamming =15;	//the hamming weight	
	mpz_class b =32;    //the bound
	int H[r];
	
	//INITIALIZING exponents H TO ONES. Which is the default value.
	
	for (int i=0;i<r;i++){
		H[i] =1;
	}

	H[0]=19;
	H[1]=5;
	H[2]=4;
	H[3]=1;
	H[4]=1;
	H[5]=1;
	H[6]=1;
	//H[7]=5;

//-----------------------------------------//

//START DEBUGGING

	int* Q;
    	Q = new int[r];

	set_Q(r, Q);
	Lambda(Q, H, r, L);
	cout << "Lambda is : " << L << endl; 
    
    //Here we print some basic things we need to know

	std::list<unsigned int*> P;
	make_P_set(Q,H,r,L,P);

	int Psize = P.size();
	cout << "hamming   : " << hamming << endl;
	cout << "P size is : " << Psize << endl;
	cout << "looking for Carmichaels with " << Psize - hamming << " factors" << endl;
	printf("\n");
//	mpz_class R=1;
//	euler_totient(Q,H,r,R);
//	cout << "euler_phi of Lambda is     : " << R << endl;

    	double euler_phi_log;
    	double density_of_the_problem;
    	density(Q,H,r,Psize,euler_phi_log,density_of_the_problem);
//    cout << "log_2 (euler_phi(Lambda) ) : " << fixed << euler_phi_log<< endl;
    	cout << "density of the problem     : " << fixed << density_of_the_problem << endl;
    	printf("\n");
//WHOLE TESTING
	
	unsigned long list_size = P.size();
	unsigned int** P2;
	P2 = new unsigned int*[list_size];
	for (unsigned long i=0;i<list_size;i++)
		P2[i] = new unsigned int[r];
	from_list_to_array(P, P2,r);

//GENERATING I set
	
	mpz_class n=list_size;
	int found = 0;				
	double begin = omp_get_wtime();
	
//START THE TEST
	while(found==0){
//	for(int ite=0;ite<100;ite++){ 		
		mpz_class** I;			
		I = new mpz_class*[2];		//RESULTS
		I[0] = new mpz_class[mpz_get_ui(b.get_mpz_t())];
		I[1] = new mpz_class[mpz_get_ui(b.get_mpz_t())];
		gen_I(n,b,1, I);
		mpz_class count =0;		//THIS IS A COUNTER FOR HOW MANY INTERSECTIONS WE GET

		//AT THIS STAGE WE HAVE THE COMBINATIONS IN SOL1, SOL2 THAT 
		//PRODUCE THE CARMICHAEL NUMBERS
		//SO WE NEED TO EXTRACT THESE NUMBERS AND CHECK IF THEY 
		//ARE INDEED CARMICHAEL
		
		double total_time = omp_get_wtime();		//TOTAL TIMER

		found = T_set(Q,r,P2, n, L, I, b, hamming, total_time);
		delete[] I[0];
		delete[] I[1];
		delete[] I;
		if (found==1)
			break;
}
	double end = omp_get_wtime();

	cout << endl;
	if(found==0)
		cout << "The algorithm did ALL The ITERATIONS WITHOUT any SUCCESS" << endl;	
	delete[] Q;
	delete[] P2;	
	return 0;	
}


#include <iostream>
#include <inttypes.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <gmpxx.h>
#include "iostream"	
#include <math.h>
#include <time.h>
#include <fstream>
#include <list>
using namespace std;

void Lambda(int* Q,int* H,int r, mpz_class &Lambda){
	
	mpz_class exp;
	for (int i=0;i<r;i++){
		exp = pow(Q[i], H[i]); 
		Lambda *= exp;
	}
	return;
}

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


void arrcpy(mpz_class* &dest, mpz_class* &source, mpz_class &size){
	for (mpz_class i=0;i<size;i++)
	{	
		dest[mpz_get_ui(i.get_mpz_t())] = source[mpz_get_ui(i.get_mpz_t())];
	}

}

void dup_array(int* A,int* B,int size, int** &out){

        for (int i=0;i<size;i++){
               // cout << "P["<<i<<"] is: " << A[i] << endl;
               // cout << "H["<<i<<"] is: " << B[i] << endl;
                out[0][i] = A[i];
                out[1][i] = B[i];
        }
}


void divisors(int* P, int* H, int r,mpz_class &size, mpz_class* &divs){

        //getting the 2d array needed
        int** dup;
        dup = new int*[2];
        dup[0] = new int[mpz_get_ui(size.get_mpz_t())];
        dup[1] = new int[mpz_get_ui(size.get_mpz_t())];
        dup_array(P,H,r,dup);

        //-------------------------

        divs[0] = 1;
        mpz_class out_count =1;
        for (int i=0;i<r;i++){
                mpz_class* prev;
                prev = new mpz_class[mpz_get_ui(size.get_mpz_t())];
                mpz_class pn = 1;
                arrcpy(prev, divs, size);
                for (int j=0;j<dup[1][i];j++){
                        pn *= dup[0][i];
                        mpz_class k = 0;
                        while (prev[mpz_get_ui(k.get_mpz_t())] != 0 && k<size){
                                divs[mpz_get_ui(out_count.get_mpz_t())] = prev[mpz_get_ui(k.get_mpz_t())] * pn;
                                out_count++;
                                k++;
                        }
                }
		delete[] prev;
        }
	delete[] dup[0];
	delete[] dup[1];
	delete[] dup;
}


void multiply_list(int* H,int r, mpz_class &result){
      
        for (int i=0;i<r;i++){
                result *= (H[i]+1);
        }
	return;
}

void make_P_set(int* Q, int* H, int r,mpz_class &Lambda, std::list<mpz_class> &P){
    mpz_class size_d=1; 
	multiply_list(H,r,size_d);
        mpz_class* divs;
        divs = new mpz_class[mpz_get_ui(size_d.get_mpz_t())];
        divisors(Q, H, r, size_d, divs);
        for (mpz_class i=0;i<size_d;i++){
                mpz_class num = divs[mpz_get_ui(i.get_mpz_t())] + 1;
                if (mpz_probab_prime_p(num.get_mpz_t(), 5) !=0 && Lambda%num!=0 && num<Lambda){
			P.insert(P.end(), num);
                }
        }
	delete[] divs;
	return;
}


void euler_totient(int* Q, int* H, int size, mpz_class &result)
// L <-- Lambda
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

int main(){
	mpz_class L=1;
	int r=3;		    //number of first primes
	int hamming =9;	//the hamming weight	
	mpz_class b =40;    //the bound
	int H[r];

	//INITIALIZING H TO ONES
	
	for (int i=0;i<r;i++){
		H[i] =1;
	}

	H[0]=8;
	H[1]=5;
	H[2]=3;
	H[3]=1;
	H[4]=1;
	
	int* Q;
    Q = new int[r];
	//cout << "1) STARTED RUNNING" << endl; 		    //DEBUG POINT (1)
	set_Q(r, Q);
	//cout << "Q set done" << endl;					//DEBUG POINT(2)
	Lambda(Q, H, r, L);
	std::list<mpz_class> P;
	make_P_set(Q,H,r,L,P);
	int Psize = P.size();
	cout << "P size is : " << Psize << endl;
	cout << "Lambda is : " << L << endl; 
	mpz_class R=1;
	euler_totient(Q,H,r,R);
	cout << "euler_phi of Lambda is     : " << R << endl;

	double LOG2 = log2( mpz_get_ui(R.get_mpz_t()) );
	cout << "log_2 (euler_phi(Lambda) ) : " << LOG2 << endl;

	double density = Psize/LOG2;
	cout << "density of the problem : " << density << endl;


	//cout << "density of Lambda          : " << density << endl; 

}
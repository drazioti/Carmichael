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

void from_list_to_array(std::list<mpz_class> &source, mpz_class* &dest){
	unsigned long j=0;
	for(std::list<mpz_class>::iterator it=source.begin();it != source.end(); ++it){
		dest[j] = *it ;
		j++;
        }
}


//FUCNTION (1): FUNCTIONS FOR SIZE OF P

//NOT USED

//FUNCTION (2): CALCULATING LAMBDA

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
	cout << "DEBUGGIN THE Q SET FUNCTION" << endl;
	mpz_class temp;
	for (int i=1;i<r;i++){
		temp = q;
		//cout << "First loop entered" << endl;
		//cout << "Temp value is : " << temp << endl;
		//raise(SIGINT);
                while (mpz_probab_prime_p(temp.get_mpz_t(), 5) != 1 && mpz_probab_prime_p(temp.get_mpz_t(), 5) != 2)
                {
			//cout  << "Probab prime value: " << mpz_probab_prime_p(temp.get_mpz_t(), 5) <<  endl;
                        //cout  << "And temp value : " << temp <<endl;
			q+=2;
			temp = q;
                }
	       	
                Q[i] = q;
                q +=2;
        }
}


//FUNCTION (4): FUNCTIONS FOR MAKING THE P SET


//void arrcpy(mpz_class* &dest, mpz_class* source,mpz_class size){
//        for (mpz_class i=0;i<size;i++){
//                dest[mpz_get_ui(i.get_mpz_t())] = source[mpz_get_ui(i.get_mpz_t())];
//        }
//}



void dup_array(int* A,int* B,int size, int** &out){

        for (int i=0;i<size;i++){
                cout << "P["<<i<<"] is: " << A[i] << endl;
                cout << "H["<<i<<"] is: " << B[i] << endl;
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
                        //printf("Round P[%d]=%d, current power: %d \n",i,P[i],j);
                        pn *= dup[0][i];
                        //cout << "Pn this round got:" << pn << endl;
                        mpz_class k = 0;
                        while (prev[mpz_get_ui(k.get_mpz_t())] != 0 && k<size){
                                //cout << "Got into the last while in round" << i << " and current prev is:" << prev[mpz_get_ui(k.get_mpz_t())] << endl;
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


void make_P_set(int* Q, int* H,int r,mpz_class &Lambda, std::list<mpz_class> &P){

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

//FUNCTION(5): FUNCTION FOR GETTING MPZ_CLASS NUMBER P FROM
//	       LIST OF EXPONENTS

mpz_class get_P_element(int* Q, int* H, int r){
	mpz_class p=1;
	for (int i=0;i<r;i++){
		int exp =  pow(Q[i], H[i]);
		p *= exp;
	}
	return p;
}

//THIS CODE ABOVE IS USED FOR THE OPTIMIZED VERSION OF STORING THE
//THE P SET
int is_carmichael(mpz_class &n, mpz_class* &factors, mpz_class &sizef);

//FUNCTION(6): MAKE T_SET

int T_set(mpz_class* &P, mpz_class &sizeP, mpz_class &Lambda, mpz_class** &I, mpz_class &sizeI,int local_hamming_weight){
	mpz_class b=1;
	for (mpz_class i=0;i<sizeP;i++)
	{
		b*= P[mpz_get_ui(i.get_mpz_t())];
		mpz_mod(b.get_mpz_t(), b.get_mpz_t(), Lambda.get_mpz_t());
	}
	mpz_class hamming =8;
	mpz_class count=0;
	std::list<int*> sol1;
	std::list<int*> sol2;
	product_attack_1(P,sizeP,Lambda, b, I, local_hamming_weight,sizeI, sol1, sol2, count);
	
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
	subset2 = new mpz_class[mpz_get_ui(sizeI.get_mpz_t())];					//SUBSETS OF P BASED ON I1,I2
	
	//cout << "---DEBUG: INSIDE EXTRACT NUMBER----" <<endl;
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
//FACTORS FOR IS CARMICHAEL FUNCTION
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

//--------------------------------------------------------------
int main(){
	clock_t startP, endP;
	mpz_class L=1;
	
//-----CHANGE THESE PARAMETERS TO RUN-------//	
	int r=10;		//number of first primes
	int hamming =21;		
	mpz_class b =29;
	
	int H[r];
	//INITIALIZING H TO ONES FOR SIMPLICITY
	for (int i=0;i<r;i++){
		H[i] =1;
	}
	H[0]=8;
	H[1]=3;
	H[2]=3;
	H[3]=3;
	H[4]=2;
	//H[5]=7;
	//H[6]=5;
	//H[7]=5;

//-----------------------------------------//

//START DEBUGGING
	int* Q;
    Q = new int[r];

	cout << "1) STARTED RUNNING" << endl; 				//DEBUG POINT (1)
	set_Q(r, Q);
	cout << "Q set done" << endl;					//DEBUG POINT(2)
	Lambda(Q, H, r, L);
	cout << "Lambda done" << endl;
	cout << "Lambda is : " << L << endl; 
	startP = clock();
  //DEBUG POINT(3)
	std::list<mpz_class> P;
	make_P_set(Q,H,r,L,P);

	endP = clock();
	printf("\n\n\n");

	printf("\n");

	//print out P set
	//for(std::list<mpz_class>::iterator it=P.begin();it != P.end(); ++it){	
	//	cout << "P element: "<< *it << endl;
	//}	
	//print Q set
	//printf("\n\n");
	for (int i=0;i<r;i++){
		cout << "Q["<<i<<"] is: "<< Q[i] << endl; 
	}
	printf ("\nTime for P set is : %f seconds", (double) (endP-startP)/1000000);
//WHOLE TESTING
	
	unsigned long list_size = P.size();
	mpz_class* P2;
	P2 = new mpz_class[list_size];
	from_list_to_array(P, P2);

//GENERATING I set
	
	mpz_class n=list_size;
	int found = 0;				//FOR TESTING WE WILL CHOOSE
	clock_t begin = clock();
	
//START THE TEST

	while(found==0){
//	for(int ite=0;ite<100;ite++){ 
		
		
		mpz_class** I;			
		I = new mpz_class*[2];		//RESULTS
		I[0] = new mpz_class[mpz_get_ui(b.get_mpz_t())];
		I[1] = new mpz_class[mpz_get_ui(b.get_mpz_t())];
	//b=2;
		gen_I(n,b,1, I);
		//cout << "I size is : " << n << endl;
		//for(int i=0;i<b;i++){
              	//	cout << "I[" << i<< "] are: "<< I[0][i] << " and " << I[1][i] << endl;
		//}
		//for(int i=0;i<n;i++){
		//	cout << "Array P[" << i <<"] element: " << P2[i] << endl;
		//}
		mpz_class c=1;
		//std::list<int*> sol1;
		//std::list<int*> sol2;
		mpz_class count =0;						//THIS IS A COUNTER FOR HOW MANY INTERSECTIONS WE GET
	//product_attack_1(P2, L, c, I,8,b, sol1, sol2, count);    

//AT THIS STAGE WE HAVE THE COMBINATIONS IN SOL1, SOL2 THAT 
//PRODUCE THE CARMICHAEL NUMBERS
//SO WE NEED TO EXTRACT THESE NUMBERS AND CHECK IF THEY 
//ARE INDEED CARMICHAEL

//!!!!!BAD IMPLEMENTATION NEED CHANGE FOR H1,H2!!!!!!!!
		int h1;
        int h2;
		int local_hamming_weight = 8;
    	if(local_hamming_weight%2==1){
            	h1 = local_hamming_weight/2;
            	h2 = h1+1;
    	}
    	else{
            	h1 = local_hamming_weight/2;
            	h2=h1;
    	}

		//mpz_class* numbers;
		//cout << "ALLOCATING " << count << " BOXES FOR NUMBERS" << endl;
		//numbers = new mpz_class[mpz_get_ui(count.get_mpz_t())];
		
		//extract_number(P2, I, b, sol1, sol2, h1, h2, numbers);

		found = T_set(P2, n, L, I,b,hamming);	
		delete[] I[0];
		delete[] I[1];
		delete[] I;
		cout << "P size is : " << P.size() << endl;
		if (found==1)
			break;
//	}
}
	clock_t end = clock();
	cout << "Time elapsed: " << double(end - begin)/CLOCKS_PER_SEC << endl;
	if(found==0)
		cout << "DID ALL ITERATIONS WITHOUT SUCCESS " << endl;
	//sol1.merge(sol2);
//TEST SET
	delete[] Q;
	delete[] P2;	
	return 0;	
}

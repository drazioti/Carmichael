#include <iostream>
#include <string.h>
#include <gmpxx.h>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <vector>
#include "Combinations.h"
#include <list>
#include "unordered_map"
using namespace std;


//############################BASIC FUNCTIONS###########################

//FUNCTIONS (1): IMPLAMENTATION OF SIMPLE USEFULL FUNCTONS FOR MPZ OBJECTS


int chrcmp(const char chr1, const char chr2) {
	size_t length1, length2;

	char s1[2] = {chr1, '\0'};
	char s2[2] = {chr2, '\0'};

	length1 = strlen(s1);
	length2 = strlen(s2);

	if (length1 == 1 && length2== 1){
		if (strcmp(s1,s2) == 0){
			return 1;
		}
		else{
			return 0;
		}
	}
	else{
		return 0;
	}
}

void arrcpy(mpz_class* &dest, mpz_class* &source, mpz_class &size){
	for (mpz_class i=0;i<size;i++)
	{	
		dest[mpz_get_ui(i.get_mpz_t())] = source[mpz_get_ui(i.get_mpz_t())];
	}

}


mpz_class fact(mpz_class &n){
	mpz_class f=1;
	for(mpz_class i=n;i>0;i--){
		f*=i;
	}
	return f;	
}

void comb_size(mpz_class &n, mpz_class k, mpz_class &size)
{
	mpz_class tempN,tempK,temp;
	tempN = fact(n);
	tempK = fact(k);
	temp = n-k;
	temp = tempK*fact(temp);
       	size = tempN/temp;	
}



//FUCNTION (4): CODE FOR "func1"

mpz_class func1(mpz_class* &P, int* E, mpz_class &Lambda, mpz_class &c, int flag, int size){
		
	if (c<=0)
	{
		cout << "c must be positive" << endl;
		return 0;
	}
	mpz_class prod =1;
	
	if (flag==1){		//element of set U1
		//cout << "Func1 size : " << size << endl;
		for (mpz_class i=0;i<size;i++){
			//cout << "Prod element: " << P[mpz_get_ui(i.get_mpz_t())] << endl;
			int index = E[mpz_get_ui(i.get_mpz_t())];
			prod *= P[index]; 					//modulo after every prod
			mpz_mod(prod.get_mpz_t(), prod.get_mpz_t(), Lambda.get_mpz_t());
			
		}
	}
	else if (flag==2){		//element of set U2
		
		prod *= c;
		for (mpz_class i=0;i<size;i++) 		//for loop to calculate product
		{
				int index = E[mpz_get_ui(i.get_mpz_t())];
				mpz_class temp;
				if(mpz_invert(temp.get_mpz_t(), P[index].get_mpz_t(), Lambda.get_mpz_t())!=0)
				{
					//cout << "Inverse modulo of " << P[mpz_get_ui(i.get_mpz_t())] << " is: " << temp << endl;
					prod *= temp;	
					mpz_mod(prod.get_mpz_t(), prod.get_mpz_t(), Lambda.get_mpz_t());
				}
				else{cout << "Number "<< P[mpz_get_ui(i.get_mpz_t())].get_mpz_t() << " doesnt have inverse modulo "<< Lambda << endl;
				}	
	}
	}
	return prod;

}

//FUNCTION (5): CODE FOR "func2"

void func2(mpz_class* &P, int** E, mpz_class &Lambda, mpz_class &c, int flag, mpz_class &sizeE, unordered_map <string,int*> &Map, mpz_class &sizeI, int h1){
	cout <<"Entered func2\n";
	//cout << "Number of rounds or sizeE : " << sizeE << endl;
	for (mpz_class i=0;i<sizeE;i++){
		mpz_class temp=func1(P,E[mpz_get_ui(i.get_mpz_t())],Lambda,c,flag,h1);
		Map[temp.get_str()] = E[mpz_get_ui(i.get_mpz_t())]; 	
	}


}

//FUNCTION (6): FUNCTIONS FOR U1

void U1(mpz_class* &P, mpz_class* &I, int** E, mpz_class &Lambda,mpz_class &sizeE,mpz_class &sizeI,int h1, unordered_map<string, int*> &Map){

	mpz_class* Q;
	Q = new mpz_class[mpz_get_ui(sizeI.get_mpz_t())];
	//cout << "Size of I is : " << mpz_get_ui(sizeI.get_mpz_t()) << endl;
	//cout << "Q1 countains: ";
	for (mpz_class i=0;i<sizeI;i++)
	{
		mpz_class index = I[mpz_get_ui(i.get_mpz_t())];
		Q[mpz_get_ui(i.get_mpz_t())] = P[mpz_get_ui(index.get_mpz_t())];
	//	cout << Q[mpz_get_ui(i.get_mpz_t())] << " " ;
	}
	//cout <<"."<<endl;
	//cout << endl;
	mpz_class c=1;
	cout << "ENTERED U1\n";	
	func2(Q,E,Lambda,c,1,sizeE, Map, sizeI,h1);
	delete[] Q;
}
//FUNCTION (7): FUNCTION FOR SAVING SUBSET AFTER FINDING INTERSECTION

void sol(std::list<int*> &sol1, std::list<int*> &sol2, int* E1, int* E2)
{
	//cout << "Entered SOL FUCNTION" << endl;
	sol1.push_back(E1);
	sol2.push_back(E2);
}

//FUNCTION (8): INTERSECTION FUNCTIONS aka HASHMAP SEARCH

void intersection(mpz_class* &P, mpz_class* &I, int ** E, mpz_class &Lambda, mpz_class &c, mpz_class &sizeE, mpz_class &sizeI, unordered_map <string, int*> &Map, mpz_class &count, std::list<int*> &sol1, std::list<int*> &sol2, int h2){
	mpz_class* Q;
	Q = new mpz_class[mpz_get_ui(sizeI.get_mpz_t())];
	//cout << "Q2 contains: ";
	for (mpz_class i=0;i<sizeI;i++)
	{
		mpz_class index = I[mpz_get_ui(i.get_mpz_t())];
		Q[mpz_get_ui(i.get_mpz_t())] = P[mpz_get_ui(index.get_mpz_t())];
	//	cout << Q[mpz_get_ui(i.get_mpz_t())] << ", ";
	}
	cout << endl;
	
	//ON THE FLY CALCULATION OF U2 ELEMENTS
	for (mpz_class i=0;i<sizeE;i++){
                mpz_class temp=func1(Q,E[mpz_get_ui(i.get_mpz_t())],Lambda,c,2,h2);	
	
		//EACH ELEMENT OF U2 IS STORED IN A TEMPORARY
		//VARIABLE TEMP AND TESTED FOR INTERSECTION

		unordered_map<string, int*>::const_iterator got = Map.find(temp.get_str());
		if (got != Map.end()){
			cout << "!!!!!!!!FOUND INTERSECTION!!!!!!!!" <<  endl;
			//cout << "Product : " << got->first << " with combination ";
			//for (int k=0;k<h2;k++)
				//cout << ' '  <<got->second[k];
			//cout << endl;
			count++;
			sol(sol1,sol2,got->second,E[mpz_get_ui(i.get_mpz_t())]); 
			//cout << "RETURNED BEFORE FINISHING SIZE_E"<<endl;
			delete[] Q;
			return ;
		}
	}
	printf("\n\n");
	cout << "!!!!!!!!!!!FOUND " << count << " CARMICHAEL NUMBERS!!!!!!!!!! " << endl;
	printf ("\n\n\n");
	cout << "END OF INTERSECTION" << endl; 
	delete[] Q;
//	return 0;

}

//FUNCTION (8): GENERATE I SET FUNCTION

void gen_I(mpz_class &n, mpz_class &b, int flag, mpz_class** &I){
        if (b>=n){
                cout << "bound must be smaller than I size" << endl;
                return;
        }
        mpz_class bound;
        if (flag==0){                   //bound is set to n/2
                cout << "0-option sets the bound to n/2" << endl;
                bound= n/2;
        }
        else if (flag==1){
                cout << "1-option sets the bound to b" << endl;
                bound = b;
        }
        else{
                cout << "flag can be only 1 or 0" << endl;
                return;
        }

        mpz_class temp;
        mpz_sqrt(temp.get_mpz_t(), n.get_mpz_t());
        if (2*b>temp){
                mpz_class* G;
                G = new mpz_class[mpz_get_ui(n.get_mpz_t())];
                for(mpz_class i=0;i<n;i++)
                        G[mpz_get_ui(i.get_mpz_t())] = i;
                random_shuffle(&G[0], &G[mpz_get_ui(n.get_mpz_t())]);
                for (mpz_class i=0;i<b;i++){
                        mpz_class index = i + b;
                        I[0][mpz_get_ui(i.get_mpz_t())] = G[mpz_get_ui(i.get_mpz_t())];
                        I[1][mpz_get_ui(i.get_mpz_t())] = G[mpz_get_ui(index.get_mpz_t())];
                }
		delete[] G;
        }
        else{
		unordered_map <unsigned int, unsigned int> S;
                gmp_randstate_t state;
                gmp_randinit_default(state);
                mpz_class i=0;
                while(i<2*b)
                {
                        mpz_class choice;
                        mpz_urandomm(choice.get_mpz_t(), state, n.get_mpz_t());
                        unordered_map<unsigned int, unsigned int>::const_iterator got = S.find(mpz_get_ui(choice.get_mpz_t()));
                        if (got == S.end())
                        {
                                S[mpz_get_ui(i.get_mpz_t())] = mpz_get_ui(choice.get_mpz_t());
                                i++;
                        }
                }
                for (mpz_class i=0;i<b;i++){
                        mpz_class index = i+b;
                        I[0][mpz_get_ui(i.get_mpz_t())] = S[mpz_get_ui(i.get_mpz_t())];
                        I[1][mpz_get_ui(i.get_mpz_t())] = S[mpz_get_ui(index.get_mpz_t())];
                }
		S.clear();
	}
	
}



//FUNCTION(9): PRODUCT SUBSET ATTACK PHASE 1

int product_attack_1(mpz_class* &P, mpz_class &Lambda, mpz_class &c, mpz_class** &I,int local_hamming_weight,mpz_class &sizeI, std::list<int*> &sol1, std::list<int*> &sol2, mpz_class &count)
{
	int h1;
	int h2;
	if(local_hamming_weight%2==1){
		h1 = local_hamming_weight/2;
		h2 = h1+1;
	}
	else{
		h1 = local_hamming_weight/2;
		h2=h1;
	}

        mpz_class sizeE1;
	mpz_class sizeE2;
	comb_size(sizeI, h1, sizeE1);
	comb_size(sizeI, h2, sizeE2);
	int** E1;
        E1 = new int*[mpz_get_ui(sizeE1.get_mpz_t())];
        for (mpz_class i=0;i<sizeE1;i++)
                E1[mpz_get_ui(i.get_mpz_t())] = new int[h1]; //MAKE ARRAY OF INDEXES OF COLUMN SIZE = comb_size
	int** E2;				             //AND ROW SIZE = H1 OR H2
        E2 = new int*[mpz_get_ui(sizeE2.get_mpz_t())];
        for (mpz_class i=0;i<sizeE2;i++)
                E2[mpz_get_ui(i.get_mpz_t())] = new int[h2]; //
	cout << "DEBUG BEFORE PERMS" << endl;
        //-----------COMBINATIONS-------------------//
	//for E1

	Combinations* c_obj1;
	c_obj1	=new Combinations(mpz_get_si(sizeI.get_mpz_t()),h1);
	
	cout << "SizeI is : " << sizeI << endl;
	cout << "h1 is : " << h1 << endl;
	
	for (unsigned long i=0;i<mpz_get_ui(sizeE1.get_mpz_t());i++){
		std::vector<short> cmb;
		cmb = c_obj1->next_combination();
		int j=0;
		for (std::vector<short>::iterator it = cmb.begin();it != cmb.end(); ++it)
			{
				E1[i][j] = *it;
				j++;
			}
	}
	//for E2
		
	Combinations* c_obj2;
        c_obj2  =new Combinations(mpz_get_si(sizeI.get_mpz_t()),h2);

	for (unsigned long i=0;i<mpz_get_ui(sizeE2.get_mpz_t());i++){
                std::vector<short> cmb;
                cmb = c_obj2->next_combination();
		int j=0;
                for (std::vector<short>::iterator it = cmb.begin();it != cmb.end(); ++it)
                        {
                                E2[i][j] = *it;
                                j++;
                        }
        }

	
	cout << "Perms done ! " << endl;
	//
	unordered_map <string, int*> U;
	
	cout <<endl;
	cout <<endl;
	U1(P, I[0], E1, Lambda, sizeE1, sizeI,h1, U);
       	cout << "U1 done" << endl;
//TESTING THE MAP
	cout << endl;
//SOLUTION SETS
//cout  << "SOLUTION SET1 " << endl;	
	
//TESTING THE INSTERSECTION
	intersection(P, I[1], E2, Lambda, c, sizeE2, sizeI, U, count,sol1,sol2,h2);	
	cout << endl;
	for (mpz_class i=0;i<sizeE1;i++)
                delete[] E1[mpz_get_ui(i.get_mpz_t())];
	for (mpz_class i=0;i<sizeE2;i++)
                delete[] E2[mpz_get_ui(i.get_mpz_t())];
        delete[] E1;
	delete[] E2;
	cout <<"Deleted E1,E2 try to delete c_obj"<<endl;
	delete c_obj1;
	delete c_obj2;
	U.clear();
	cout << endl;
	return 0;	
		
}	








//##################################TESTING###########################
//int main(){
//	mpz_class* L ;
//	L = new mpz_class[5];
//	L[0]=1;
//	L[1]=0;
//	L[2]=1;
//	L[3]=0;
//	L[4]=1;
	//string str = "abc";
	//mpz_class* ones;
	//ones = new mpz_class[5];
	//ones = positions_of_1(L,5);
	// n = str.size();
	//perms(str,0,n-1);
	//for (int i=0;i<5;i++)
	//{
	//	cout << ones[i] << endl;
	//}


//PERMS TESTING : START
	
//	mpz_class count=0;
//	char** E;
//	E = new char*[10];
//	for (int i=0;i<10;i++)
//		E[i] = new char[5];
//	for (int i=0;i<10;i++){
//                for (int j=0;j<5;j++){
//                        E[i][j] = '0';
//                }
//        }

	//perms(5,3,count, E);
	//for (int i=0;i<10;i++){
	//	cout << "Bitstring " <<i<< " is: ";
	//	for (int j=0;j<5;j++){
	//		cout << E[i][j];
	//	}
	//	cout << endl;
	//}
	


//PERMS TESTING : END


//GEN_I TESTING: START
	
//	mpz_class n=100;			// n is the size of |P|
//	mpz_class b=15;			// b must be smaller or equal to n
//	mpz_class** I;
//	I = new mpz_class*[2];
//	I[0] = new  mpz_class[mpz_get_ui(b.get_mpz_t())];
//	I[1] = new mpz_class[mpz_get_ui(b.get_mpz_t())];
//	gen_I(n,b,1,I);

//GET_I TESTING: END


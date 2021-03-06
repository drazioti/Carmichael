// $g++ --std=c++11 carmi.cpp Combinations.cpp -lgmpxx -lgmp -lcrypto

// $nohup ./a.out > script.out 2>&1 &
#include <fcntl.h>
#include <unistd.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <gmpxx.h>
#include "iostream"	
#include <new>
#include <math.h>
#include <sys/time.h>
#include <sys/types.h>
#include <signal.h>
#include "subset_product.cpp"
#include <list>
#include <algorithm>
#include <ctime>
#include <fstream>
#include <args.hxx>
using namespace std;

//BASIC TEMPORARY UTILITY

void from_list_to_array(std::list<mpz_class> &source, mpz_class* &dest){
	unsigned long j=0;
	for(std::list<mpz_class>::iterator it=source.begin();it != source.end(); ++it){
		dest[j] = *it ;
		j++;
        }
}

//FUNCTION (2): CALCULATING the modulus LAMBDA

void Lambda(int* Q,int* H,int r, mpz_class &Lambda){
	
	mpz_class exp;
	for (int i=0;i<r;i++){
		
		mpz_ui_pow_ui(exp.get_mpz_t(),Q[i], H[i]); 
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
		mpz_class exp;
		mpz_ui_pow_ui(exp.get_mpz_t(),Q[i], H[i]);
		p *= exp;
	}
	return p;
}

//THIS CODE ABOVE IS USED FOR THE OPTIMIZED VERSION OF STORING THE
//THE P SET

int is_carmichael(mpz_class &n, mpz_class* &factors, mpz_class &sizef);

//FUNCTION(6): MAKE T_SET

  int T_set(mpz_class* &P, mpz_class &sizeP, mpz_class &Lambda, mpz_class** &I, mpz_class &sizeI,int local_hamming_weight,double total_time, int frag, char Q_bytes){
	mpz_class b=1;
	for (mpz_class i=0;i<sizeP;i++)
	{
		b*= P[mpz_get_ui(i.get_mpz_t())];
		mpz_mod(b.get_mpz_t(), b.get_mpz_t(), Lambda.get_mpz_t());
	}
	mpz_class count=0;
	std::list<int*> sol1;
	std::list<int*> sol2;
    
	product_attack_1(P,sizeP,Lambda, b, I, local_hamming_weight,sizeI, sol1, sol2, count,frag, total_time, Q_bytes);
	
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

int main(int argc, char **argv){
        //double total_time = CLOCKTIME();
        cout<<"Setting up arguments"<< endl;
        mpz_class L=1;

//-----CHANGE THESE PARAMETERS TO RUN a new instance-------//   

        //declare Argument Parser
    args::ArgumentParser parser("This is an implementation of a propability attack\n to the Subsect Product Problem for the generation of Carmichael numbers.", "Make sure to place arguments right.");
    args::Group exponents_group(parser, "Exponents", args::Group::Validators::All);
        args::Group attack_group(parser, "Attack Values", args::Group::Validators::All);
    args::Group help_group(parser, "Help:", args::Group::Validators::DontCare);

//-----------------------------------------//
// args::ValueFlag<int> frag(exponents_group,"integer > 0", "Number of first primes",{'r',"firstprimes"});
args::PositionalList<int> Harray(exponents_group,"integers","List of 'r' exponets");//{'e',"exponents"}
args::ValueFlag<int> hammingWeight(attack_group,"integer", "Hamming weight",{"ham","hamming"});
args::ValueFlag<int> bound(attack_group,"integer", "Bound",{'b',"bound"});
args::ValueFlag<int> hashlenght(attack_group,"integer", "Q",{'q',"lenght of hash"});
args::ValueFlag<int> fragmentationSubsets(attack_group,"integer > 0", "Subsets of the set stored in memory",{'f',"fragmentation"});
args::HelpFlag help(help_group, "Help", "Display help menu", {'h', "help"});
// Manage arguments
    try
    {
        parser.ParseCLI(argc, argv);
    }
    catch (const args::Completion& e)
    {
        std::cout << e.what();
        return 0;
    }
    catch (const args::Help&)
    {
        std::cout << parser;
        return 0;
    }
    catch (const args::ParseError& e)
    {
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
        return 1;
    }
    catch (args::ValidationError e)
    {
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
        return 1;
    }
// SETTING VARS FROM ARGS
        int r = args::get(Harray).size(); // exponent array size example 5
        int H[r]; // exponent array example 20 5 4 1 1
        //cout << " EXPONENTS: ";
        for (int i=0;i<r;i++){
                H[i] =args::get(Harray)[i];
                //cout << H[i] << ", ";
        }
        // cout << endl;

        mpz_class b = args::get(bound); // bound example 32
        int fragmentation = args::get(fragmentationSubsets); // frag vallue
        int hamming = args::get(hammingWeight); // hamming example 15

        cout << "H values : ";
        for(int i=0;i<r;i++){
                cout << H[i] << " ";
                if ( i<r-1 && H[i]<H[i+1] ){
                        cout << endl << "H values were not valid.\nMake sure they meet requirements"<< endl << endl;
                        return 0;
                }
        }cout<<endl;

        // set Q hash lenght
        char Q_bytes = args::get(hashlenght);

        cout<<"Bound : "<<b<<endl;
        cout<<"Parameter Q : "<< (int)Q_bytes <<endl;
        cout<<"Fragmentation : "<<fragmentation<<endl;
//START DEBUGGING
    int* Q;
    Q = new int[r];


        set_Q(r, Q);
        Lambda(Q, H, r, L);
        cout << "Lambda is : " << L << endl;

    //Here we print some basic things we need to know

        std::list<mpz_class> P;
        make_P_set(Q,H,r,L,P);

        int Psize = P.size();
        cout << "Hamming   : " << hamming << endl;
        cout << "P size is : " << Psize << endl;
        cout << endl <<"Looking for Carmichael numbers with " << Psize - hamming << " factors" << endl;
        //printf("\n");
//      mpz_class R=1;
//      euler_totient(Q,H,r,R);
//      cout << "euler_phi of Lambda is     : " << R << endl;

        double euler_phi_log;
        double density_of_the_problem;
        density(Q,H,r,Psize,euler_phi_log,density_of_the_problem);
//    cout << "log_2 (euler_phi(Lambda) ) : " << fixed << euler_phi_log<< endl;
        cout << "Density of the problem : " << fixed << density_of_the_problem << endl;
        //printf("\n");
//WHOLE TESTING
        unsigned long list_size = P.size();
        mpz_class* P2;
        P2 = new mpz_class[list_size];
        from_list_to_array(P, P2);

//GENERATING I set

        mpz_class n=list_size;
        int found = 0;

        double total_time = omp_get_wtime();

        int ite =0;
//START THE TEST
        int randfile = open("/dev/random", O_RDONLY);
        unsigned int seed = 12345;  //fixed for no randomized option
        char ran=0;
	if (randfile>0)
       	{
		if (ran==1)
               		size_t check = read(randfile, &seed, sizeof(seed));
                //cout << "seed is: " << seed << endl;
                //cout << endl;
        }
        else
        {       cout << "Failed to open /dev/random" << endl;
                	exit(0);
       	}
        gmp_randclass rr(gmp_randinit_mt);
        randomize_I(rr,seed);                   //RANDOMIZE ONCE
        close(randfile);
        while(found==0){
//      for(int ite=0;ite<100;ite++){

//CHOOSING RANDOM I
                mpz_class** I;
                I = new mpz_class*[2];          //RESULTS

                I[0] = new mpz_class[mpz_get_ui(b.get_mpz_t())];
                I[1] = new mpz_class[mpz_get_ui(b.get_mpz_t())];

                gen_I(n,b,1, I, rr,seed);
                mpz_class count =0;             //THIS IS A COUNTER FOR HOW MANY INTERSECTIONS WE GET
//I CHECKS

                //AT THIS STAGE WE HAVE THE COMBINATIONS IN SOL1, SOL2 THAT
                //PRODUCE THE CARMICHAEL NUMBERS
                //SO WE NEED TO EXTRACT THESE NUMBERS AND CHECK IF THEY
                //ARE INDEED CARMICHAEL

                found = T_set(P2, n, L, I, b, hamming,total_time, fragmentation, Q_bytes);
                delete[] I[0];
                delete[] I[1];
                delete[] I;
                if (found==1)
                        break;
}

        // double end = CLOCKTIME();
        // cout << begin - end<<endl;


        cout << endl;
        if(found==0)
  cout << "The algorithm did ALL The ITERATIONS WITHOUT any SUCCESS" << endl;
        delete[] Q;
        delete[] P2;
        return 0;
}




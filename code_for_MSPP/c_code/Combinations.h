// Combinations3.h
#ifndef Combinations3_H
#define Combinations3_H


#include <vector>
#include <gmpxx.h>//gmp c++ interface

class Combinations{
	public:
		/**
		 *	init n,r,vector throw exc if input is 0
		 */
		Combinations(int n,int r);
		/**
		 * set mpz_class result as the number of different compinations
		 */
		void get_size(mpz_class& result);
		/** 
		 *	set vector as the new combination
		 */
		std::vector<int> next_combination();
	private:
		int n,r;// n,r
		std::vector<int> cmb_v;//vec with last result
		/**
		 *  set mpz_class t as the factorial of n
		 */
		void fact(mpz_class& t,int n);
};
#endif

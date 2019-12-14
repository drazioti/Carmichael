// Combinations.cpp
#include "Combinations.h"

Combinations::Combinations(int n,int r){
    //if inputs are 0
    if (n==0){
        throw std::exception();
    }
    this->n=n;
    this->r=r;
    if(this->r==0){
        cmb_v.push_back(-1);
    }
    
    // for (unsigned i=0; i<cmb_v.size(); i++)
    //     cmb_v.at(i)=i;

    //push to the vector 0,1,2,...,b
    for (int i = 0; i < this->r; i++){
        cmb_v.push_back(i);
    }
    // std::cout << "c -> n : "<< this->n << std::endl;
    // std::cout << "c -> r : "<< this->r << std::endl;
    // std::cout << "size : "<< cmb_v.size() << std::endl;
};

/**
 * set result as the number of different compinations
 */
void Combinations::get_size(mpz_class &result){
        //declare
    mpz_class nf;
    mpz_class rf;
    mpz_class nrf;
    //init
    fact(nf,n);//n!
    fact(rf,r);//r!
    fact(nrf,n-r);//(n-r)!
    //print n! r! (n-r)!
    // std::cout<<nf<<std::endl;
    // std::cout<<rf<<std::endl;
    // std::cout<<nrf<<std::endl;
    //calc
    nf = nf/rf; //nf=nf/rf
    result = nf/nrf; //result=nf/nrf
};

/** 
 *	set vector as the new combination
 */
std::vector<int> Combinations::next_combination(){
    if(r==0){
        return cmb_v;
    }
    auto prev = cmb_v;
    
    int i = cmb_v.size()-1;

    bool changed = false;
    do{
        int x = cmb_v[i];
        if (x == n-1){
            i--;
        }else if (x + (cmb_v.size() - i) < n){  
            cmb_v[i]++;
            x = cmb_v[i++] + 1;
            for (; i < cmb_v.size(); i++){
                cmb_v[i] = x++;
            }
            changed = true;
            break;
        }else{
            i--;
        }
    }while (i >= 0);

    if(!changed){
        // reset?
        for (int i = 0; i < cmb_v.size(); i++)
            cmb_v[i] = i;
    }
    return prev;
};

/**
 *  set mpz_t t as the factorial of n
 */
void Combinations::fact(mpz_class& t,int n){
    mpz_class t_tmp=1;
    for (int i = 1; i < n+1; ++i){
        t_tmp=t_tmp*i;
    }
    t=t_tmp;
};

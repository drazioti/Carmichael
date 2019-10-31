'''
Sagemath (ver.8.1).
Code for Modular Subset Produc Problem
Initial author : K.A.Draziotis (2019), drazioti@gmail.com
Licence : GPL v.2


### An example for product subset problem
iterations = 60
hamming = [8..11]
N = int(1e6)
print "N:",N
print "bits of N:",bits(N)
P = [2..N] # P is the set consisting from the first N natural numbers
Lambda = next_prime(N+5) # the modulus
print "modulus=",Lambda
c = 190238      # target number
print "Target number:",c
SET = product_subset_attack_second_phase(P,Lambda,c,hamming,iterations) # birthday attack
if SET!='attack failed':
    prod(SET)%Lambda==c # verification


'''

def bits(n):
    return floor(log(n,2))+1
             
             
##################  Product Subset Problem - Birthday Attack


def positions_of_1(L): # L is a list
    A = []
    for i in range(len(L)):
        if L[i]==1:
            A.append(i)
    return A
 
        
def bin2int(a):
    '''
    convert a binary string to integer
    sage: bin2int('0101')
          5
    '''
    
    return int("".join(str(x) for x in a),2)
      

def func1(P, E, Lambda, c, flag):
    ''' 
    P is a list say P=[2,3,5,7]
    E is the list of exponents, from the function : perms_list(n,hamming,flag)
    Lambda is the modulus and c=1 or c>1
    then it returns the product P[i]^(E[i]) if c = 1
    or product c*P[i]^(-E[i]) if c>1
    
    '''
    n = len(P)
    if c<=0:
        return " c must be positive. "
    if flag == 1: # we construct the set U1 (does not have 'c')
        A = str(  int( mod(prod(P[i]^(E[i]) for i in range(n)),Lambda) ) )
    if flag == 2:  # for the set U2, it needs 'c'
        A = str( int( mod(c*prod(P[i]^(-E[i]) for i in range(n)),Lambda) ) )      
    return A
 


def func2(P,E,Lambda,c,flag): 
    # E is a list of lists (list of binary exponents)
    # P is the set of primes
    # Lambda is the modulus
    # the function returns a list of the form [ ['a',integer1],['b',integer2], etc ]
    
    return [ [func1(P,E[i],Lambda,c,flag),bin2int(E[i]) ] for i in range(len(E))]
        

        
def U1(I,P,E,Lambda):
    U = []
    Q = [P[i] for i in I]
    L = func2(Q,E,Lambda,1,1)
    for x in L:
        U.append(int(x[0]))
    return U
        
             
def U2(I,P,E,Lambda,c):
    U = []
    Q = [P[i] for i in I]
    L = func2(Q,E,Lambda,c,2)
    for x in L:
        U.append(int(x[0]))
    return U

   
def perms(b,hamming,flag): 
    ''' 
    Input : b, hamming, flag, hamming < b, flag = 0 or \not =0
    Output: Returns all the permutations of b-bit words with specific Hamming weight (variable hamming). 
    if flag!=0, then it returns all the m-bit words (the algorithm ignores the variable hamming).
    Finally, it returns a set of strings not lists.
    For instance, 
    sage : list(perms(4,2,0))
    
    [['1100', '1010', '1001', '0110', '0101', '0011']]
    
    sage : list(perms(3,2,1)) # all 2^3=8, 3-bits binary words, it ignores hamming = 1
    
    ['000', '001', '010', '011', '100', '101', '110', '111']
    
    '''
    
    import itertools
    if not b:
        return
    if flag == 0:
        result = []
        for bits in itertools.combinations(range(b), hamming):
            s = ['0'] * b
            for bit in bits:
                s[bit] = '1'
            result.append(''.join(s))
        yield result
        
    if flag!=0:
        for i in xrange(2**b):
            s = bin(i)[2:]
            s = "0" * (m-len(s)) + s
            yield s

  
def perms_list(n,hamming,flag):
    # the same as the previous function, but it returns the elements as lists (not strings)
    s = list(perms(n,hamming,flag))
    if flag!=0:
        return [[int(s[j][i]) for i in range(n)] for j in range(len(s))]
    if flag==0:
        s = s[0]
        return [[int(s[j][i]) for i in range(n)] for j in range(len(s))]

def gen_I(n,b,flag):
    ''' this function return two disjoint subsets of {0,1,2,...,n} with b elements'''
    
    import hashlib,random
    a1 = 0
    b1 = n
    L = []
    I2 = []
    I1 = []
    bound = floor(b)
    if b >= b1 - a1:
        print "bound must be smaller than ",b1-a1
        return
    if flag not in [0,1]:
        print "flag is either 0 or 1"
        return         
    if flag==0: # when you choose flag=0 then bound = n/2. The function ignores the parameter b.
        #hash = hashlib.md5(os.urandom(2*n)).digest()
        #random.seed(hash)
        L = random.sample(range(a1,b1), n)
        I1 = L[0:bound]
        I2 = L[bound:len(L)]        
    if flag==1:
        #hash = hashlib.md5(os.urandom(2*n)).digest()
        #random.seed(hash)
        L = random.sample(range(a1,b1), 2*bound)
        I1 = L[0:bound]
        I2 = L[bound:2*bound]
    return I1,I2

def product_subset_attack_first_phase(P,Lambda,c,I1,I2,local_hamming_weight):
    '''
    Input :
    =======
    
    P      : the original set set
    Lambda : modulus
    c      : target number
    I1, I2 disjoint subsets of {0,1,2,...,n}
    
    Output :
    ========
    
    Null or the subset of P say L = [a1,a2,..,ar] (r = 2 * local_hamming_weight) such that
    a1*a2*...*ar = c (mod Lambda)

    '''
    def extract_set(P,U_1,U_2,I1,I2,E1,E2,sol):
        index_1 = U_1.index(sol[0])
        index_2 = U_2.index(sol[0])
        Q1 = [P[i] for i in I1]
        Q2 = [P[i] for i in I2]
        R1 = [Q1[i] for i in positions_of_1(E1[index_1])]
        R2 = [Q2[i] for i in positions_of_1(E2[index_2])]
        M  = R1+R2
        return M

    h1 = floor(local_hamming_weight/2.)
    h2 = ceil(local_hamming_weight/2.)
    #print "h1,h2:",h1,h2
    E1 = perms_list(len(I1),h1,0)     # the set of binary exponents having hamming weight h1
    E2 = perms_list(len(I2),h2,0)     # the set of binary exponents having hamming weight h2
    U_1 = U1(I1,P,E1,Lambda)
    U_1 = U1(I1,P,E1,Lambda)         # the first set U1
    U_2 = U2(I2,P,E2,Lambda,c)       # the second set U2
    sol = list((Set(U_1)).intersection((Set(U_2)))) # calculate the intersection
    if len(sol)>=1:
        return extract_set(P,U_1,U_2,I1,I2,E1,E2,sol)
    else:
        return "No solution found"



def product_subset_attack_second_phase(P,Lambda,c,hamming,iterations):
    '''
    Input :
    =======
    
    P      : the original set set
    Lambda : modulus
    c      : target number
    hamming: is a list of the possible hamming weights ton check
    
    Output :
    ========
    
    Null or the subset of P say L = [a1,a2,..,ar] (r = 2 * local_hamming_weight) such that
    a1*a2*...*ar = c (mod Lambda)

    '''
    for local_hamming_weight in hamming:
        print "H(c):",local_hamming_weight
        for i in range(iterations):  # TODO :  how to choose the number of iterations?
            I1,I2 =  gen_I(len(P),local_hamming_weight,1)
            #print "I1,I2:",I1,I2
            #print local_hamming_weight
            Out = product_subset_attack_first_phase(P,Lambda,c,I1,I2,local_hamming_weight)           
            if Out!='No solution found':
                print Out
                return Out
    print "attack failed" 
    return "attack failed"

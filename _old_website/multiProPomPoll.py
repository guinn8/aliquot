import multiprocessing as mp
from math import expm1, exp, log, pi, floor, factorial
from sympy import divisors, divisor_count
from functools import reduce, partial

def driver(y):
    ks = [0,1,2,3,4,5,6,7,8,9,10]
    for k in ks:
        
        print("Delta" , k , " = " , multiPro(y, k))


def multiPro(y, k):
    
    z = floor(y/9)
    divisions = [[1 + x *z , z + x *z ] for x in range(0, 9)]
   
    if __name__ == '__main__':
        with mp.Pool(10) as p:
            parts =  sum(p.map(partial(preDeltaK, k), divisions)) * 1/log(y)
            p.close()
            return parts
    
    

def preDeltaK(k,argList):
  
    acc = 0
   
    for a in range(argList[0], argList[1]+1):
        #if(a % 2 == 0 ): acc += (a**k/s(a)**k) * (1/( a * factorial(k))) * exp(-a/s(a)) 
        if(a % 2 == 0 ): acc += 1.0561531762023264**k * (1/( a * factorial(k))) * exp(-a/s(a))
        #if(a % 2 == 0 ): acc += (1/( a * factorial(k))) * exp(-a/s(a))

    return acc

def s(n):
    return reduce(lambda a, b : a + b, divisors(n)) - n 

def averageSum(y):
    
    for k in range(1, 10):
        acc = 0
        evens = 0 
        for a in range(2,y):
            
            if a % 2 == 0:
                evens+=1
                acc += (a/s(a))**k
        print("where k = ", k, " average (a/s(a))^k = ", acc/evens  )
                

    return acc/evens

if __name__ == '__main__':
    bound = input("Enter integer upper bound:")
    #driver(int(bound))
    print(averageSum(int(bound)))




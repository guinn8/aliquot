import numpy as np
import sys
from sympy.ntheory.factor_ import totient 
from sympy import divisors, divisor_count
from functools import reduce
from math import expm1, exp, log, pi, floor
import matplotlib.pyplot as plt 
import csv
import multiprocessing as mp

def data(yMax):
    stops = []
    acc = 50000
    while acc <= yMax: 
        stops.append(acc)
        acc += 50000
    
    
    sheet = []

    with open('C:\\Users\\GavinPC\\OneDrive\\Desktop\\guinn8.github.io\\aliqData.csv', mode='w') as dat:
        writer = csv.writer(dat, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        for y in stops:
            print(y)
            row = []
            A_y = np.uint64(np.lcm.reduce(list(range(1, y+1))))
            div_Ay = divisors(int(A_y))
            
            row.append(y)
            row.append(conjDelta_0(y))
            row.append(conjDelta_1(y))
            row.append(conjDelta_2(y))
            row.append(conjDelta_3(y))
            row.append(conjDelta_4(y))
            sheet.append(row)
            
        for row in sheet:
            writer.writerow(row)

           
def test(y):
    accEven = 0
    accOdd = 0
    denom = 0
    for a in range(2, y+1): 
        denom +=1
        if(a % 2 == 0):
            accEven += a/s(a)
        else: 
            accOdd += a/s(a)  
    print("even average: ", accEven/denom, " odd average ", accOdd/denom)

def driver(y):
    A_y = np.uint64(np.lcm.reduce(list(range(1, y+1))))
    div_Ay = divisors(int(A_y))
    print("Conj delta 0 :", conjDelta_0(y))
    print("Conj delta 1 :", conjDelta_1(y))
   
    print("Conj delta 2 :", conjDelta_2(y))
    print("Conj delta 3 :", conjDelta_3(y))
    print("Conj delta 4 :", conjDelta_4(y))
    print("Conj delta 5 :", conjDelta_5(y))

    print("Density of 0 parent numbers:", delta_0(A_y, div_Ay))
    print("Density of 1 parent numbers:", delta_1(A_y, div_Ay))
    print("Density of 2 parent numbers:", delta_2(A_y, div_Ay))
    print("Density of 3 parent numbers:", delta_3(A_y, div_Ay))

def testDriver(y):
    print("Conj delta 1 :", conjDelta_1(y))
    print("Conj delta 1 test :", conjDelta_1_test(y))
    print("Conj delta 2 :", conjDelta_2(y))
    print("Conj delta 2 test :", conjDelta_2_test(y))

#Implementation of (3.4)
def conjDelta_0(y):
    acc = 0
    for a in range(2, y+1):
        if(a % 2 == 0):
            acc += (1/a) * exp(-a/s(a))
    return (1 / log(y)) * acc



def conjDelta_1(y):
    acc = 0
    for a in range(2, y+1):
            if(a % 2 ==0 ): acc +=  (1/s(a)) * exp(-a/s(a))
    return (1 / log(y)) * acc
    
   
def conjDelta_2(y):
    acc = 0
    for a in range(2, y+1):
        if(a % 2 ==0 ): acc += (a/(2*(s(a)**2))) * exp(-a/s(a))
    return (1 / log(y)) * acc

def conjDelta_3(y):
    acc = 0
    for a in range(2, y+1):
        
        if(a % 2 ==0 ): acc += ((a**2)/(6*(s(a)**3))) * exp(-a/s(a))
    return (1 / log(y)) * acc

def conjDelta_4(y):
    acc = 0
    for a in range(2, y+1):
        
        if(a % 2 ==0 ): acc += ((a**3)/(24*(s(a)**4))) * exp(-a/s(a))
    return (1 / log(y)) * acc

def conjDelta_5(y):
    acc = 0
    for a in range(2, y+1):
        
        if(a % 2 ==0 ): acc += ((a**4)/(120*(s(a)**5))) * exp(-a/s(a))
    return (1 / log(y)) * acc

def conjDelta_1_test(list):
    acc = 0
    
    for a in range(list[0], list[1]+1):
        
        if(a % 2 == 0 ): acc +=  (1/a) * exp(-a/s(a))
    return (1 / log(list[1])) * acc

def conjDelta_1_test_multi(list):
    acc = 0
    
    for a in range(list[0], list[1]+1):
        if(a % 2 == 0 ): acc +=  (1/(a)) * exp(-a/s(a))
    return acc

def conjDelta_2_test_multi(list):
    acc = 0
    for a in range(list[0], list[1]+1):
        if(a % 2 == 0 ): acc +=  1.11555*(1/(2*a)) * exp(-a/s(a))
    return acc

def multiPro(y):
    
    z = floor(y/9)
    divisions = [[1 + x *z , z + x *z ] for x in range(0, 9)]
    if __name__ == '__main__':
        with mp.Pool(10) as p:
            parts = (p.map(conjDelta_2_test_multi, divisions))
            print (sum(parts)* (1/(log(y))))



#Implementation of (3.1)
def delta_0(A_y, div_Ay):
    acc = 0
    for a in div_Ay:
        if(a % 2 ==0 ): acc +=  (1/a) * exp(-a/s(a))
    return  (totient (A_y) / A_y) * acc

def delta_1(A_y, div_Ay):
    acc = 0
    for a in div_Ay:
        if(a % 2 ==0  ): acc +=  (1/s(a)) * exp(-a/s(a))
    return  (totient (A_y) / A_y) * acc

def delta_2(A_y, div_Ay):
    acc = 0
    for a in div_Ay:
        if(a % 2 ==0  ): acc +=  (a/(2*(s(a)**2))) * exp(-a/s(a))
    return  (totient (A_y) / A_y) * acc

def delta_3(A_y, div_Ay):
    acc = 0
    for a in div_Ay:
        if(a % 2 ==0  ): acc += ((a**2)/(6*(s(a)**3))) * exp(-a/s(a))
    return  (totient (A_y) / A_y) * acc

def s(n):
    return reduce(lambda a, b : a + b, divisors(n)) - n 

#test(int(sys.argv[1]))   
#data(int(sys.argv[1]))
driver(int(sys.argv[1]))
#plot(int(sys.argv[1]))
#testDriver(int(sys.argv[1]))
#print(conjDelta_1_test([1, int(sys.argv[1])]))
#multiPro(int(sys.argv[1]))
#print(conjDelta_2(int(sys.argv[1])))
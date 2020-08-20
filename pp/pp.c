#include <stdio.h>
#include <stdlib.h>
#include <assert.h> 
#include <flint/arith.h>

#define UPPERPARENTS 16



int numChunks = 10;
int chunk_size;
int buffer_size;

unsigned long s(unsigned long n);
unsigned long factorial(unsigned long n); 
double accumulator( unsigned long start_a, char k,  unsigned long * propSumDiv, long double * acc);

//This program implements an generaliztion of Conj. 1.4 of Pollack/Pomerance "Some problems of Erdos on the Sum of Divisors Function"
//Instead of estimating the natural density of only aliqout orphans this program also estimates the density of k-parent aliquot numbers
//n is a k-parent aliqout number iff there are k distinct natural numbers m st s(m) = n
//An aliquot orphan is a 0-parent aliquot number
//let delta-k be the estimated density of k-parent aliquot numbers and s(n) be the sum-of-proper-divisors function
//delta-k = 1/log(max_bound) * sum(forall a <= max_bound)( (a^(k-1) * e^(-a/s(a)) / k! * s(a)^k) )
int main(int argc, char *argv[]){
    unsigned long max_bound;

    long double acc[UPPERPARENTS] = {0};

    if(argc < 1){
        printf("./[max_bound]");
        exit(0);
    }

    max_bound = atol(argv[1]);
    chunk_size = max_bound/numChunks;
    buffer_size = chunk_size /2;
    
    //loop through the chunks for multi-threading
    for (int i = 0; i < numChunks; i++){

        unsigned long m = (chunk_size * i) + 2; //first even number in chunk 
        unsigned long * sigma =  malloc(buffer_size * sizeof(unsigned long));
        
        for(int j = 0; j < buffer_size; j++){
           sigma[j] = s(m + 2*j);
        }
       
        for(int g = 0; g < UPPERPARENTS; g++){
            accumulator(m, g, sigma, acc);
        }
    }
    
    for(int i = 0; i < UPPERPARENTS; i++){
        acc[i] *= 1/(log(max_bound) * factorial(i));
        printf("delta %d = %Lf\n", i, acc[i]);
    }
}

unsigned long s(unsigned long n){
    fmpz_t res, num;
    fmpz_init(res);
    fmpz_init_set_ui(num,n);

    arith_divisor_sigma(res, num, 1);

    fmpz_sub(res, res, num);
    return fmpz_get_ui(res);
}

unsigned long factorial(unsigned long n) { 
    if (n == 0) return 1; 
    else return n * factorial(n - 1); 
} 


double accumulator( unsigned long start_a, char k,  unsigned long * propSumDiv, long double * acc){
    double numer;
    double denom;
    double a;
    unsigned long s_a;
   
    for(unsigned long i = 0; i < buffer_size; i++){
       
        a = (2 * i) + start_a;
        s_a = propSumDiv[i];

        //This needs to be treated casewise because computers dont understand neg. exp. correctly
        if(k == 0){
            numer =exp(-a / s_a);
            denom = a;
        }else{
            numer = pow(a, k-1) * exp(-a / s_a);
            denom = pow(s_a, k);
        }
        
        acc[k] += numer / denom;
    }
}
  
#include "sieve.h"
#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <assert.h> 

//void factorial(mpz_t result, unsigned long input);
unsigned long max_bound;
int numChunks = 10;
int chunk_size;
long factorial(int n);


const int initWheel [3] = {2, 3, 5 };
const int inc [8] = {4, 2, 4, 2, 4, 6, 2, 6};

unsigned long wheelDivSum(unsigned long n){

    if(n == 1 || n == 0) return 0;
    else if (n < 0){
        printf("\nValue must be non-neg \n");
        exit(EXIT_FAILURE);
    }

    unsigned long savedN = n;
    unsigned long curr_sum = 1;
    unsigned long curr_term = 1;
    unsigned long res = 1;
  
    for (int i = 0; i < 3; i++){
        if(n % initWheel[i] == 0){
            do{
                n = n/initWheel[i];
                curr_term *= initWheel[i];
                curr_sum += curr_term;
            }while(n % initWheel[i] == 0);

            res *= curr_sum;
            curr_sum = curr_term = 1;
        }
    }

    unsigned long k = 7;
    int i = 0;
    while(k <= n){
        if(n % k == 0){
            do{
                n = n / k;
                curr_term *= k;
                curr_sum += curr_term;
            }while(n % k == 0);
            res *= curr_sum;
            curr_sum =  curr_term = 1;
        } else {
            k += inc[i];
            if(i < 7) i++;
            else i = 0;
        }
    }
    res *= curr_sum;
    return res - savedN;
}

double kParentEst( unsigned long start_a, char k,  unsigned long * sigma, long double * acc){
    unsigned long numer;
    unsigned long denom;
    unsigned long a;
    unsigned long s_a;
    for(unsigned long i = 0; i < chunk_size; i++){
       
        a = i + start_a;
        s_a = sigma[i] -a;
      
       // printf("s(%lu) = %lu\n", a, s_a );
       // assert(a % 2 == 0);
       // assert(s_a % 2 == 0);


        //numer =exp(-a / s_a);
        //numer = pow(a, k-1) * exp(-a / s_a);
        //denom = pow(s_a, k);
        
        //acc[k] = numer / a;
    }
}

int main(int argc, char *argv[]){
    long double acc[140] = {0};

    if(argc < 1){
        printf("./[max_bound]");
        exit(0);
    }

    max_bound = atol(argv[1]);
    chunk_size = max_bound/numChunks;

    //These buffers are nessicary for the prime seive
    //Most of this is pulled straight from Anton's implementation
    //so I really dont have a good idea how it work
    // const unsigned long max_prime = (unsigned long) sqrt((double) (max_bound << 1));
	// const unsigned prime_bound = (unsigned long) (1.25506 * (max_prime + 1) / log(max_prime + 1));
	// unsigned int * primes = (unsigned int *) malloc(sizeof(unsigned int) * (prime_bound + 1));
	// prime_sieve(max_prime, primes);

    for (int i = 0; i < numChunks; i++){
        unsigned long m = (chunk_size * i)+1;
        unsigned long * sigma =  malloc(chunk_size * sizeof(unsigned long));
        
        printf("m = %lu\n", m);
        for(int j = 1; j < chunk_size; j+=2){
            wheelDivSum(m+j)
        }
       // sum_of_divisors(chunk_size , m, sigma, primes);

        kParentEst(m, 0, sigma, acc);
    }
    
    acc[0] *= 1 / log(max_bound);
    printf("\n acc[0] = %Lf\n", acc[0]);

}

long factorial(int n)
{
   int c;
   long result = 1;
 
   for( c = 1 ; c <= n ; c++ )
         result = result*c;
 
   return ( result );
}


// void factorial(mpz_t result, unsigned long input) {
//     mpz_set_ui(result, 1);
//     while (input > 1) {
//         mpz_mul_ui(result, result, input--);
//     }
// }
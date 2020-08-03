#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

unsigned long max_bound;

unsigned long wheelDivSum(unsigned long n);
void recParents(char characArr[max_bound/2]);

int main(int argc, char *argv[]){
    max_bound = atol(argv[1]);
    char * localCharac = malloc(max_bound / 2);
    unsigned long s_n ;

    for(unsigned long i = 2; i <= max_bound*2; i+=2){
        s_n = wheelDivSum(i);
        printf("i: %lu out of: %lu\n", i, max_bound*2);
        if(s_n <= max_bound && s_n % 2 == 0 ) {
            localCharac[(s_n / 2)-1]++;
        }
       
    }

    for(unsigned long i = 1; i <= max_bound; i+= 2){
        s_n = wheelDivSum(i*i);
        printf("i: %lu out of: %lu\n", i, max_bound);
        if(s_n <= max_bound && s_n % 2 == 0 ) {
            localCharac[(s_n / 2)-1]++;
        }
       
    }

    recParents(localCharac);
    
 
}

const int initWheel [3] = {2, 3, 5 };
const int inc [8] = {4, 2, 4, 2, 4, 6, 2, 6};
unsigned long wheelDivSum(unsigned long n){

    if(n == 1) return 0;
    else if (n <= 0){
        printf("\nValue must be non-positive \n");
        exit(EXIT_FAILURE);
    }

    unsigned long savedN = n;
    unsigned long curr_sum = 1;
    unsigned long curr_term = 1;
    unsigned long res = 1 ;


    for (int i = 0; i < 3; i++){
        if(n % initWheel[i] == 0){
            
            do{
                n = n/initWheel[i];
                curr_term *= initWheel[i];
                curr_sum += curr_term;
            }while(n % initWheel[i] == 0);

            res *= curr_sum;
            curr_sum = 1;
            curr_term = 1;
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
            curr_sum = 1;
            curr_term = 1;
        } else {
            k += inc[i];
            if(i < 7) i++;
            else i = 0;
        }
    }

    res *= curr_sum;
    return res -savedN;
}

void recParents(char characArr[max_bound/2]){

    unsigned long accPreimages[256] = {0};

    int count = 0;
    accPreimages[0]++; //accounts for 5 which is the
    for(unsigned long i = 0;  i < max_bound/2; i++){
        accPreimages[characArr[i]]++;
    }

    printf("\nBound = %lu\n", max_bound);
    for(int j = 0;  j < 8; j++){
       printf("\n%d count: %lu, density = %f", j, accPreimages[j], (float)accPreimages[j]/max_bound);
    }
}
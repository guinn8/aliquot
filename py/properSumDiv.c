#include "properSumDiv.h"

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
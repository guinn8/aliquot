#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>

const unsigned long max_bound = 1000000000;
const int countMax = 256;

unsigned long wheelDivSum(unsigned long n);
void recParents(char * characArr, unsigned long max_bound);

int main(int argc, char *argv[]){
    remove("preimagesParents.csv");

   
    unsigned long current_bound = 1000;
    unsigned long last_bound = 0;

    FILE * fp = fopen("preimagesParents.csv","a");
    fprintf(fp, "bound, 0");
    for(int i = 1; i < countMax; i++){
       fprintf(fp, ",%d", i);
    }
    fprintf(fp, "\n");
    fclose(fp);

    char * localCharac = malloc(max_bound / 2);
    for(int i = 0; i < max_bound / 2; i++){
            localCharac[i] = 0;
    }

   
    while(current_bound <= max_bound){
        

        #pragma omp parallel
        {
            #pragma omp section
            {
                #pragma omp for 
                for(unsigned long i = 2; i <= current_bound*2; i+=2){
                    unsigned long s_n = wheelDivSum(i);
                    //printf("i: %lu out of: %lu\n", i, max_bound*2);
                    if(s_n > last_bound && s_n <= current_bound && s_n % 2 == 0 ) {
                        #pragma omp atomic
                        localCharac[(s_n / 2)-1]++;
                    }   
                }
            }
        
       
            #pragma omp section
            {
                for(unsigned long i = 1; i <= current_bound; i+= 2){
                    unsigned long s_n = wheelDivSum(i*i);
                    //printf("i: %lu out of: %lu\n", i, max_bound);
                    if(s_n > last_bound && s_n <= current_bound && s_n % 2 == 0 ) {
                        #pragma omp atomic
                        localCharac[(s_n / 2)-1]++;
                    }
                }    
            }
        }
       

        recParents(localCharac, current_bound);
        last_bound = current_bound;
        current_bound += 10000;
    }
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

void recParents(char * characArr, unsigned long max_bound){

    FILE * fp = fopen("preimagesParents.csv","a");

    unsigned long accPreimages[256] = {0};

    int count = 0;
    accPreimages[0]++; //accounts for 5 which is the
    for(unsigned long i = 0;  i < max_bound/2; i++){
        accPreimages[characArr[i]]++;
    }

    fprintf(fp,"%lu", max_bound);
    for(int j = 0;  j < countMax; j++){
       fprintf(fp, ",%lu",  accPreimages[j] );
    }
    fprintf(fp, "\n");

    fclose(fp);
}
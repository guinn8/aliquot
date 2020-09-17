#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include "properSumDiv.h"

#define FILENAME "counts.csv"

const unsigned long max_bound = 1000000000;
const int countMax = 256;


void recParents(char * characArr, unsigned long max_bound);

int main(int argc, char *argv[]){
    remove(FILENAME);

    unsigned long current_bound = 1000;
    unsigned long last_bound = 0;

    FILE * fp = fopen(FILENAME, "a");
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


void recParents(char * characArr, unsigned long max_bound){

    FILE * fp = fopen(FILENAME,"a");

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
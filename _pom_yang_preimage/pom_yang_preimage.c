// Author: Gavin Guinn
// Contact: gavinguinn1@gmail.com

//#define FLINT
//#define ASSERT //for any production runs ensure assert is not defined
#define TIMING

#define FILENAME "counts.csv"

#include "sieve.h"
#include <omp.h>
#include <assert.h> 
#include <string.h>
#include <stdint.h>
//#include "/home/gavin.guinn/flint-2.6.2/arith.h"

#ifdef FLINT
#include <flint/arith.h>
uint64_t s(uint64_t n);
uint64_t s_sqInput(uint64_t n);
#else
#include "properSumDiv.h"
#endif

//uints64 are used ensure consistant behaviour about over unsigned wrap
//We are safe as long as all values are below 18,446,744,073,709,551,615
//The largest number that can possibly be produced is: max_bound ^ 4/3
//This means that approx 281474976700000 (2.81x10^14) is the largest safe value for max_bound

void writeBuffer(uint64_t * imageBuffer, uint32_t bufferInd, uint8_t * characFunc);
void tabStats(uint8_t * characArr);
uint32_t writeImage(uint64_t image, uint64_t * imageBuffer, uint32_t bufferInd, uint8_t * characFunc);
uint32_t chooseNumChunks(uint64_t max_bound);

uint32_t buffer_size; //size of buffer alloced to each chunk
uint64_t max_bound;

//These 2 values are used to estimate the total percentage we are through the computation
#ifdef TIMING
uint64_t numPreimages;
uint64_t writtenPreimages = 0;
#endif

//Purpose: Computed the counts of EVEN k-parent aliquot numbers  [2, max_bound]
//Parameters: 
//  -[max_bound] number to compute upto (inclusive)
//  -(optional)[buffer_size] controls how big of a buffer each thread gets before writing to central characterstic function
//Preconditions: Valid input
//Postconditions: the counts of every k-parent EVEN value are appended to counts.csv
void main(uint32_t argc, uint8_t *argv[]){

    double startTime = omp_get_wtime();

    if(argc == 0){
        printf("What tf are you trying to input\n");
        printf("./[max_bound] OR ./[max_bound][buffer_size]");
        exit(0);
    }else if(argc == 2){
        max_bound = atol(argv[1]);
        buffer_size = 10000; 
    }else if(argc == 3){
        max_bound = atol(argv[1]);
        buffer_size = atol(argv[2]);
    }else{
        printf("What tf are you trying to input\n");
        printf("./[max_bound] OR ./[max_bound][buffer_size]");
        exit(0);
    }

    #ifdef TIMING
    numPreimages = max_bound * .65;
    #endif

    assert(max_bound % 2 == 0);
    assert(max_bound < 281474976700000); //max safe value for max_bound
    assert(buffer_size > 0);

    uint32_t numChunks = chooseNumChunks(max_bound); //number of pieces we break the problem into for threading
    uint64_t chunk_size = max_bound/numChunks; //amount of numbers over which one thread will run the PY algorithm 

    //byte array that accumates parent information
    //characFunc[i] holds the number of parents for (i + 1) * 2
    uint8_t * characFunc = malloc(max_bound/2);
    memset(characFunc, 0, max_bound/2 );

    //These buffers are nessicary for the prime seive
    //Most of this is pulled straight from Anton's implementation
    //so I really dont have a good idea how it work
    const uint64_t max_prime = (uint64_t) sqrt(2 * max_bound);
	const uint64_t prime_bound = (uint64_t) (1.25506 * (max_prime + 1) / log(max_prime + 1));
	uint32_t * primes = ( uint32_t *) malloc(sizeof(uint32_t) * (prime_bound + 1));
	prime_sieve(max_prime, primes);

    //These vars handle the odd composite squares
    //these numbers are pushed into a buffer to be processed in the second part of the algor
    const uint64_t compositeBound = (uint64_t)pow(max_bound, .6666666666666666666);
    uint64_t * compSquares =  malloc(sizeof(uint64_t) * (compositeBound));
    uint64_t compSquareCounter = 0;

    #pragma omp parallel
    {
        //sigma holds the sum-of-divisors returned by antons sum_of_divisors_odd
        uint64_t * sigma =  malloc ((chunk_size /2)*sizeof(uint64_t));

        //this buffer holds preimages before they are written to disk in a batch
        uint64_t * imageBuffer = malloc (sizeof(uint64_t)* buffer_size);
        uint32_t bufferInd; //counter for the imageChunk buffer
       
        //splits the problem into chunks for threading
        #pragma omp for schedule(auto) 
        for(uint32_t i = 0; i < numChunks; i++){

            #ifdef TIMING
            double threadStart = omp_get_wtime();
            #endif

            bufferInd = 0;

            uint64_t m = (i * chunk_size);//the number currently being processed by algo

            sum_of_divisors_odd(chunk_size, m, sigma, primes);//Antons again, very fast sum of divisors
      
            //runs through the thread specific chunk
            for(uint64_t j = 0; j < chunk_size / 2; j++){
              
                //sigma(m) = sigma[j] 
                m = (i * chunk_size + 1) + (2 * j); //cannot remember where tf this offset comes from

                #ifdef ASSERT
                assert(s(m)+m == sigma[j]);
                #endif

                //catches evens values of sigma[j] and runs them through recurrance
                if(!(sigma[j] & 1)){ 

                    uint64_t  t = 3 * sigma[j] - (2 * m);
                   
                    #ifdef ASSERT
                    uint32_t pow = 2;
                    #endif

                    while (t <= max_bound){

                        #ifdef ASSERT
                        //printf("s(%lu) = %lu and t = %lu\n", pow*m, s(pow * m), t);
                        assert(s(pow * m) == t);
                        pow *= 2;
                        #endif

                        bufferInd = writeImage(t, imageBuffer, bufferInd, characFunc);
                        t = (2 * t) + sigma[j];
                    }
                }

                //catches primes are records them appropriatly 
                if(sigma[j] == m+1){

                    #ifdef ASSERT
                    if(s_sqInput(m)  != sigma[j]) printf("   m: %lu s(%lu) = %lu != sigma[j] = %lu\n", m,m*m,s_sqInput(m), sigma[j] );
                    //assert(s(m*m) == sigma[j]);
                    #endif

                    bufferInd = writeImage(m+1, imageBuffer, bufferInd, characFunc);
                }

                //catches odd composite values of m which are pushed to
                //buffer for later processing
                else if(m <= compositeBound && m != 1){ // We are counting the number of preimages for even numbers [1,max_bound] since s(1) = 0 it falls outside that range and should be excluded
                    #pragma omp critical (writeComp)
                    {
                        compSquares[compSquareCounter] = m;
                        compSquareCounter++;
                    }
                }
            }

            //finish up any leftover info
           if(bufferInd > 0) {
               writeBuffer(imageBuffer, bufferInd, characFunc);
               bufferInd = 0;
           }
        }

        free(sigma);

        //Im not sure this nowait is safe but it seems to work
        #pragma omp for nowait
        for(uint64_t  i = 0; i < compSquareCounter; i++){

	        #ifdef FLINT
            uint64_t s_mSq = s(compSquares[i]*compSquares[i]);
	        #else
            uint64_t s_mSq = wheelDivSum(compSquares[i]*compSquares[i]);
	        #endif

            #ifdef ASSERT
            assert(compSquares[i] % 2 != 0);
            #endif

            if(s_mSq <= max_bound) {
                bufferInd = writeImage(s_mSq, imageBuffer, bufferInd, characFunc);
            }
        }

        if(bufferInd > 0){
            writeBuffer(imageBuffer, bufferInd, characFunc);
            bufferInd = 0;
        } 
    }

    tabStats(characFunc);

    #ifdef TIMING
    printf("\n\nFinished in %f seconds after writing %lu images\n", omp_get_wtime()-startTime, writtenPreimages);
    #else
    printf("\n\nFinished in %f seconds \n", omp_get_wtime()-startTime);
    #endif
}

//Purpose: writes a image of s(n) to the chunks buffer and if the buffer is full it is written to the charac function
//Parameters:
//  -image: s(n) = image found by the PY algorithm
//  -* imageBuffer = buffer of images being held for a later write to the charac func
//  -bufferInd = index currently being filled by preimages
//  -* characFunc = charac function where statistics for each number is accumulated
//Preconditions: valid input
//Postconditions: the current chunkcount is returned
uint32_t writeImage(uint64_t image, uint64_t * imageBuffer, uint32_t bufferInd, uint8_t * characFunc){

    imageBuffer[bufferInd] = image;
    bufferInd++;
    
    if(bufferInd == buffer_size){
        writeBuffer(imageBuffer, bufferInd, characFunc);
        bufferInd = 0;
    }

    return bufferInd;
}

//Purpose: This function writes a buffer of preimages to the characFile
//Parameters:
//  -* imageBuffer: buffer full of images to be acc'ed in characFunc
//  -bufferInd = index currently being filled by preimages
//  -* characFunc = charac function where statistics for each number is accumulated
//Preconditions: Valid input
//Postconditions: buffer is acc'ed in characFunc
void writeBuffer( uint64_t * imageBuffer, uint32_t bufferInd, uint8_t * characFunc){

    #ifdef TIMING
    double chunkTime  = omp_get_wtime();
    #endif

    for(uint32_t i = 0; i < bufferInd; i++){
        #ifdef ASSERT
        assert(imageChunk[i] <= max_bound);
        assert (imageChunk[i] > 0);
        assert(imageChunk[i] % 2 == 0);
        #endif

        #pragma omp atomic
        characFunc[(imageBuffer[i]/2)-1]++;
    }
    
    #ifdef TIMING
    #pragma omp atomic
    writtenPreimages += bufferInd;
    printf("%%%-.2f || %d preimages written in %f seconds\n", (float)100 * writtenPreimages/numPreimages ,bufferInd, omp_get_wtime()-chunkTime);
    #endif
}

//Purpose: This functions processes the characArray and tabulates statisics about k-parent numbers characArr must be max_bound/2 in size
//Parameters: 
//  -characArr: char array where the number of preimages for each evenb number is written
//Preconditions: processed characArr
//Postconditions: None
void tabStats(uint8_t * characArr){

    FILE * fp = fopen(FILENAME,"a");

    uint64_t accPreimages[256] = {0};

    for(uint64_t i = 0;  i < max_bound/2; i++){
        accPreimages[characArr[i]]++;
    }

    fprintf(fp,"%lu", max_bound);
    for(uint32_t j = 0;  j < 255; j++){
       fprintf(fp, ",%lu",  accPreimages[j] );
    }
    fprintf(fp, "\n");

    fclose(fp);
}

#ifdef FLINT

//Purpose: compute the sum of proper divisors for input
//Parameters:
//  -n: natural number input for s(n) 
//Precondtions: None
//Postconditions: s(n) is returned
uint64_t s(uint64_t n){
    fmpz_t res, num;
    fmpz_init(res);
    fmpz_init_set_ui(num,n);

    arith_divisor_sigma(res, num, 1);

    fmpz_sub(res, res, num);
    return fmpz_get_ui(res);
}

//Purpose: compute the sum of proper divisors for input squared
//Parameters:
//  -n: natural number input for s(n*n) 
//Precondtions: None
//Postconditions: s(n*n) is returned
uint64_t s_sqInput(uint64_t n){
    fmpz_t res, num;
    
    fmpz_init(res);
    fmpz_init_set_ui(num,n);
    fmpz_mul(num, num, num);
    arith_divisor_sigma(res, num, 1);

    fmpz_sub(res, res, num);
    return fmpz_get_ui(res);
}
#endif

//Purpose: When threading in this manner it is nessicary to choose a number of chunks to 
//  split the problem into such that each thread can work through the problem a chunk at a time
//  chooses the first divisor greater than the sqrt of max bound to be the number of chunks
//Parameters:
//  -max_bound: as described
//Preconditions: Valid input
//Postconditions: Returns a valid number of chunks to break the problem into
uint32_t chooseNumChunks(uint64_t max_bound){
    uint32_t chunks;
    uint32_t root = sqrt(max_bound);
    for (uint32_t i=1; i <= max_bound; i++){
        if (max_bound % i==0){
            
            if(i >= root && (max_bound / i) % 2 == 0){
                chunks = i;
                return chunks;
            }
        } 
    }
    printf("Could not find a suitable divisor i of maxbound st maxbound/i is even");
    exit(EXIT_FAILURE);
}
#include "mmapArray.h"
#include "properSumDiv.h"
#include "sieve.h"
#include <omp.h>
#include <assert.h> 

#define TIMING

unsigned long max_bound;
unsigned long chunk_size;

const int numChunks = 1000;

const int subBufferSize = 1000000;
const int subBuffers = 100;
int subBufRange;

unsigned long totalWrites;

void writePreimage(ulong preimage, ulong * imageChunk, int * chunkCount, char * characFunc, omp_lock_t subRangeLocks[subBuffers]);
void writeBuffer( ulong * imageChunk, int  chunkCount, int subBuffer, char * characFunc, omp_lock_t subRangeLocks[subBuffers]);
void tabStats(unsigned char * characArr);
int cmpfunc (const void * a, const void * b);







int main(int argc, char *argv[]){
    

   // #ifdef TIMING
    double startTime = omp_get_wtime();
  //  #endif
    
   

    max_bound = atol(argv[1]);
    chunk_size = max_bound/numChunks;

    assert(argc > 1);
    assert(max_bound % 10 == 0 );
    assert(max_bound % subBuffers == 0);

    subBufRange = max_bound / subBuffers;

    //byte array that accumates parent information
    //characFunc[i] holds the number of parents for (i + 1) * 2
    unsigned char * characFunc = createByteArray("charac", max_bound/2);

    //These buffers are nessicary for the prime seive
    //Most of this is pulled straight from Anton's implementation
    //so I really dont have a good idea how it work
    const unsigned long max_prime = (unsigned long) sqrt((double) (max_bound << 1));
	const unsigned prime_bound = (unsigned long) (1.25506 * (max_prime + 1) / log(max_prime + 1));
	unsigned int * primes = (unsigned int *) malloc(sizeof(unsigned int) * (prime_bound + 1));
	prime_sieve(max_prime, primes);

    //These vars handle the odd composite squares
    //these numbers are pushed into a buffer to be processed in the second part of the algor
    unsigned long compositeBound = (unsigned long)pow(max_bound, .6666666666666);
    unsigned long * compSquares =  malloc(sizeof(unsigned long) * (compositeBound));
    unsigned long compSquareCounter = 0;

    omp_lock_t subRangeLocks[subBuffers];
    
    for(int i = 0; i < subBuffers; i++){
        omp_init_lock(&(subRangeLocks[i]));
    }

    #pragma omp parallel
    {
        //sigma holds the sum-of-divisors returned by antons sum_of_divisors_odd
        unsigned long * sigma =  malloc ((chunk_size /2)*sizeof(unsigned long));

        //this buffer holds preimages before they are written to disk in a batch
        //We are also splitting the preimages into subbuffers that cover a sorted portion of the charac function 
        unsigned long * imageChunk = malloc (sizeof(unsigned long)* subBufferSize * subBuffers);
        int chunkCount[subBuffers]; //counter for each subbuffer 
       
        //splits the problem into chunks for threading
        #pragma omp for schedule(auto) 
        for(int i = 0; i < numChunks; i++){

            #ifdef TIMING
            double threadStart = omp_get_wtime();
            #endif

            //zero out each subbuffer counter
            for(int i = 0; i < subBuffers; i++){
                chunkCount[i] = 0;
            }

            unsigned long m = (i * chunk_size );//the number currently being processed by algo

            sum_of_divisors_odd(chunk_size, m, sigma, primes);//Antons again, very fast sum of divisors

            //runs through the thread specific chunk
            for(int j = 0; j < chunk_size / 2; j++){
              
                //sigma(m) = sigma[j] ? 
                m = (i * chunk_size + 1) + (j << 1); //cannot remember where tf this offset comes from

                //catches evens values of sigma[j] and runs them through recurrance
                if(!(sigma[j] & 1)){ 

                    unsigned long t = 3 * sigma[j] - (m << 1);

                    while (t <= max_bound){
                        writePreimage(t, imageChunk, chunkCount, characFunc, subRangeLocks);
                        t = (t << 1) + sigma[j];
                    }
                }

                //catches primes are records them appropriatly 
                if(sigma[j] == m+1){
                    writePreimage(m+1, imageChunk, chunkCount, characFunc, subRangeLocks);
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
            for(int i = 0; i < subBuffers; i++){
                if(chunkCount[i] > 0) writeBuffer(imageChunk, chunkCount[i], i, characFunc, subRangeLocks);
                chunkCount[i] = 0;
            }
        }

        free(sigma);

        //Im not sure this nowait is safe but it seems to work
        #pragma omp for nowait
        for(int i = 0; i < compSquareCounter; i++){
            unsigned long s_mSq = wheelDivSum(compSquares[i]*compSquares[i]);

            if(s_mSq <= max_bound) {
                writePreimage(s_mSq, imageChunk, chunkCount, characFunc, subRangeLocks);
            }
        }

        //finish up any leftover info
        for(int i = 0; i < subBuffers; i++){
            if(chunkCount[i] > 0) writeBuffer(imageChunk, chunkCount[i], i, characFunc, subRangeLocks);
            chunkCount[i] = 0;
        }
    }
    
    tabStats(characFunc);
    closeByteArray(characFunc, max_bound/2);

     for(int i = 0; i < subBuffers; i++){
        omp_destroy_lock(&subRangeLocks[i]);
    }


   //#ifdef TIMING
    printf("\n\nFinished in %f seconds at %lu images per millisecond\n", omp_get_wtime()-startTime, (unsigned long)(totalWrites/ 10*(omp_get_wtime()-startTime)));
   // #endif
}
 
//writes a preimage to the chunks buffer
//if the buffer is full it is written to the charac function
void writePreimage(ulong preimage, ulong * imageChunk, int * chunkCount, char * characFunc, omp_lock_t subRangeLocks[subBuffers]){    
    
    for(int i = 0; i < subBuffers; i++){
        if(chunkCount[i] == subBufferSize){
                writeBuffer(imageChunk, chunkCount[i], i, characFunc, subRangeLocks);
                chunkCount[i] = 0;
        }
        if(preimage <= subBufRange * (i + 1) && preimage > subBufRange * i){
           // printf("wrote: %lu to subbuffer: %d\n", preimage, i);
            // assert(preimage > 0);
            // assert(chunkCount[i] < subBufferSize);
            // assert((i * subBufferSize) + chunkCount[i] < subBufferSize * subBuffers);
            // assert((i * subBufferSize) + chunkCount[i] >= 0);
            imageChunk[(i * subBufferSize) + chunkCount[i]] = preimage;
           
            chunkCount[i]++;
            break;
        } 
    }
}

//This function writes a buffer of preimages to the characFile
void writeBuffer( ulong * imageChunk, int chunkCount, int subBuffer, char * characFunc, omp_lock_t subRangeLocks[subBuffers]){

    #ifdef TIMING
    double chunkTime  = omp_get_wtime();
    double waitTime;
    #endif
    
    qsort(&imageChunk[subBuffer * subBufferSize], chunkCount, sizeof(*imageChunk), cmpfunc);

    omp_set_lock(&subRangeLocks[subBuffer]);

    #ifdef TIMING
        waitTime  = omp_get_wtime();
    #endif
    for(int i = 0; i < chunkCount; i++){
        // printf("subBuf: %d chunkCount: %lu imagechunk [%d] = %lu\n", subBuffer,chunkCount,i, imageChunk[(subBuffer * subBufRange) + i]);
        // assert(imageChunk[(subBuffer * subBufferSize) + i] <= max_bound);
        // assert (imageChunk[(subBuffer * subBufferSize) + i] > 0);
        // assert(imageChunk[(subBuffer * subBufferSize) + i] % 2 == 0);
        
        characFunc[(imageChunk[(subBuffer * subBufferSize) + i]/2)-1]++;
    }

    omp_unset_lock(&subRangeLocks[subBuffer]);
    
    #ifdef TIMING
    printf("preimages:%-10d  rate(pre's/second):%-10d    wait-time:%f\n", chunkCount, (int)( chunkCount/(10 *( omp_get_wtime() - chunkTime))), omp_get_wtime() - waitTime  );
    #endif
}

//This functions processes the characArray and tabulates statisics about k-parent numbers
//characArr must be max_bound/2 in size
void tabStats(unsigned char * characArr){
    totalWrites = 0;
    unsigned long accPreimages[256] = {0};
    
    accPreimages[0]++; //accounts for 5 which is the only odd untouchable
    for(ulong i = 0;  i < max_bound/2; i++){
        accPreimages[characArr[i]]++;
        totalWrites += characArr[i];
    }

    printf("\nBound = %lu\n", max_bound);
    for(int j = 0;  j < 8; j++){
       printf("\n%d count: %lu, density = %f", j, accPreimages[j], (float)accPreimages[j]/max_bound);
    }
}

int cmpfunc (const void * a, const void * b) {
   return ( *(int*)a - *(int*)b );
}

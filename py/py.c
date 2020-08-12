//#include "mmapArray.h"
#include "properSumDiv.h"
#include "sieve.h"
#include <omp.h>
#include <assert.h> 
#include "/home/gavin.guinn/flint-2.6.2/arith.h"

#define FLINT
//#include <flint/arith.h>
//#define FLINT
//#define TIMING
//#define ASSERT

unsigned long s(unsigned long n);
int writePreimage(ulong preimage, ulong * imageChunk, int chunkCount, char * characFunc);
void writeBuffer( ulong * imageChunk, int chunkCount, char * characFunc);
void tabStats(unsigned char * characArr);
//int cmpfunc (const void * a, const void * b) ;

unsigned long max_bound;
unsigned long chunk_size;

const int buffer_size = 100000;
const int numChunks = 256;

int main(int argc, char *argv[]){

//    #ifdef TIMING
    double startTime = omp_get_wtime();
//    #endif

    if(argc < 1){
        printf("./[max_bound]");
        exit(0);
    }

    max_bound = atol(argv[1]);
    chunk_size = max_bound/numChunks;

    if(max_bound % numChunks != 0 ){
        printf("max_bound must be divisable by numChunks");
        exit(0);
    }

    //byte array that accumates parent information
    //characFunc[i] holds the number of parents for (i + 1) * 2
    unsigned char * characFunc = malloc(max_bound/2);//creasteByteArray("charac", max_bound/2);
  

    //These buffers are nessicary for the prime seive
    //Most of this is pulled straight from Anton's implementation
    //so I really dont have a good idea how it work
    const unsigned long max_prime = (unsigned long) sqrt((double) (2 * max_bound));
	const unsigned prime_bound = (unsigned long) (1.25506 * (max_prime + 1) / log(max_prime + 1));
	unsigned int * primes = (unsigned int *) malloc(sizeof(unsigned int) * (prime_bound + 1));
	prime_sieve(max_prime, primes);

    //These vars handle the odd composite squares
    //these numbers are pushed into a buffer to be processed in the second part of the algor
    unsigned long compositeBound = (unsigned long)pow(max_bound, .6666666666666);
    unsigned long * compSquares =  malloc(sizeof(unsigned long) * (compositeBound));
    unsigned long compSquareCounter = 0;

    #pragma omp parallel
    {
        //sigma holds the sum-of-divisors returned by antons sum_of_divisors_odd
        unsigned long * sigma =  malloc ((chunk_size /2)*sizeof(unsigned long));

        //this buffer holds preimages before they are written to disk in a batch
        unsigned long * imageChunk = malloc (sizeof(unsigned long)* buffer_size);
        int chunkCount; //counter for the imageChunk buffer
       
        //splits the problem into chunks for threading
        #pragma omp for schedule(auto) 
        for(int i = 0; i < numChunks; i++){

            #ifdef TIMING
            double threadStart = omp_get_wtime();
            #endif

            chunkCount = 0;

            unsigned long  m = (i * chunk_size );//the number currently being processed by algo

            sum_of_divisors_odd(chunk_size, m, sigma, primes);//Antons again, very fast sum of divisors
      

            //runs through the thread specific chunk
            for(unsigned long j = 0; j < chunk_size / 2; j++){
              
                //sigma(m) = sigma[j] 
                m = (i * chunk_size + 1) + (2 * j); //cannot remember where tf this offset comes from
                

                #ifdef ASSERT
                assert(s(m)+m == sigma[j]);
                #endif

                //catches evens values of sigma[j] and runs them through recurrance
                if(!(sigma[j] & 1)){ 

                    unsigned long  t = 3 * sigma[j] - (2 * m);
                   
                    #ifdef ASSERT
                    int pow = 2;
                    #endif

                    while (t <= max_bound){

                        #ifdef ASSERT
                        //printf("s(%lu) = %lu and t = %lu\n", pow*m, s(pow * m), t);
                        assert(s(pow * m) == t);
                        pow *= 2;
                        #endif

                        chunkCount = writePreimage(t, imageChunk, chunkCount, characFunc);
                        t = (2 * t) + sigma[j];
                    }
                }

                //catches primes are records them appropriatly 
                if(sigma[j] == m+1){

                    #ifdef ASSERT
                    assert(s(m*m) == sigma[j]);
                    #endif

                    chunkCount = writePreimage(m+1, imageChunk, chunkCount, characFunc);
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
           if(chunkCount > 0) {
               writeBuffer(imageChunk, chunkCount, characFunc);
               chunkCount = 0;
           }
        }

        free(sigma);

        //Im not sure this nowait is safe but it seems to work
        #pragma omp for nowait
        for(unsigned long  i = 0; i < compSquareCounter; i++){

	        #ifdef FLINT
            unsigned long s_mSq = s(compSquares[i]*compSquares[i]);
	        #else
            unsigned long s_mSq = wheelDivSum(compSquares[i]*compSquares[i]);
	        #endif

            #ifdef ASSERT
            assert(compSquares[i] % 2 != 0);
            #endif

            if(s_mSq <= max_bound) {
                chunkCount = writePreimage(s_mSq, imageChunk, chunkCount, characFunc);
            }
        }
        if(chunkCount > 0){
            writeBuffer(imageChunk, chunkCount, characFunc);
            chunkCount = 0;
        } 
    }

    tabStats(characFunc);
    //closeByteArray(characFunc, max_bound/2);

    //#ifdef TIMING
    printf("\n\nFinished in %f seconds\n", omp_get_wtime()-startTime);
   // #endif
}

//writes a preimage to the chunks buffer
//if the buffer is full it is written to the charac function
int writePreimage(ulong preimage, ulong * imageChunk, int chunkCount, char * characFunc){

    imageChunk[chunkCount] = preimage;
    chunkCount++;
    
    if(chunkCount == buffer_size){
        writeBuffer(imageChunk, chunkCount, characFunc);
        chunkCount = 0;
    }

    return chunkCount;
}

//This function writes a buffer of preimages to the characFile
void writeBuffer( ulong * imageChunk, int chunkCount, char * characFunc){

    #ifdef TIMING
    double chunkTime  = omp_get_wtime();
    #endif

  //  qsort(imageChunk, chunkCount, sizeof(*imageChunk), cmpfunc);

    for(int i = 0; i < chunkCount; i++){

        #ifdef ASSERT
        assert(imageChunk[i] <= max_bound);
        assert (imageChunk[i] > 0);
        assert(imageChunk[i] % 2 == 0);
        #endif

        //printf("imagechunk [%d] = %lu\n", i, imageChunk[i]);
        if(characFunc[(imageChunk[i]/2)-1] > 20) printf("%lu has %d pre's\n", imageChunk[i], characFunc[(imageChunk[i]/2)-1]);
        #pragma omp atomic
        characFunc[(imageChunk[i]/2)-1]++;
    }
    
    #ifdef TIMING
    printf("%d preimages written in %f seconds\n", chunkCount,omp_get_wtime() - chunkTime );
    #endif


}

//This functions processes the characArray and tabulates statisics about k-parent numbers
//characArr must be max_bound/2 in size
void tabStats(unsigned char * characArr){

    FILE * fp = fopen("counts.csv","a");

    unsigned long accPreimages[256] = {0};

    accPreimages[0]++; //accounts for 5 which is the
    for(unsigned long i = 0;  i < max_bound/2; i++){
       // if(characArr[i] > 10)printf("%lu has %d preimages\n", (i+1)*2, characArr[i]);
        accPreimages[characArr[i]]++;
    }

    fprintf(fp,"%lu", max_bound);
    for(int j = 0;  j < 255; j++){
       fprintf(fp, ",%lu",  accPreimages[j] );
    }
    fprintf(fp, "\n");

    fclose(fp);
}


unsigned long s(unsigned long n){
    fmpz_t res, num;
    fmpz_init(res);
    fmpz_init_set_si(num,n);

    arith_divisor_sigma(res, num, 1);
    ulong result  = fmpz_get_si(res);
    return result -n;
}

#include <stdio.h>
#include <omp.h>
#include <stdlib.h>
#include "sieve.h"
#include <sys/mman.h>
#include <unistd.h>
#include <fcntl.h>
#include <string.h>
#include <sys/stat.h>
#include <dirent.h>
#include <ctype.h>

//#define printFile

unsigned long max_bound;
unsigned long chunk_size;
uint * primes;

const int buffer_size = 10000000;
const int numChunks = 200;

unsigned long fileTicker = 0;




const int initWheel [3] = {2, 3, 5 };
const int inc [8] = {4, 2, 4, 2, 4, 6, 2, 6};
unsigned long wheelDivSum(unsigned long n);

void updateCharac(int fileNum);
void createCharac();
long int findSize(const char *file_name);
void recParents(unsigned char *characArr);
void u_read( char * fileName, size_t fileSize, ulong * buffer);
void readCharac();

void openCharac(char * localCharac);
int writePreimage(ulong preimage, ulong * imageChunk, int chunkCount);
void u_write( ulong * buffer, size_t bufferSize);

int main(int argc, char *argv[]){
    

    double startTime =  omp_get_wtime();
    createCharac();

    if(argc < 1){
        printf("./[max_bound]");
        exit(0);
    }

    max_bound = atol(argv[1]);
    chunk_size = max_bound/numChunks;

    if(max_bound % 10 != 0 ){
        printf("max_bound must be divisable by 10");
        exit(0);
    }

    
    unsigned long compositeBound = (unsigned long)pow(max_bound, .6666666666666);

    const unsigned long max_prime = (ulong) sqrt((double) (max_bound << 1));
	const unsigned prime_bound = (ulong) (1.25506 * (max_prime + 1) / log(max_prime + 1));

	primes = (uint *) malloc(sizeof(uint) * (prime_bound + 1));
	prime_sieve(max_prime, primes);

    unsigned long * compSquares =  malloc(sizeof(unsigned long) * (compositeBound));
    unsigned long compSquareCounter = 0;

    

    #pragma omp parallel
    {

        unsigned long * sigma =  malloc ((chunk_size /2)*sizeof(ulong));
        ulong * imageChunk = malloc (sizeof(ulong)* buffer_size);
        int chunkCount;
       
        #pragma omp for schedule(dynamic) 
        for(int i = 0; i < numChunks; i++){

            double threadStart = omp_get_wtime();

            chunkCount = 0;

            unsigned long m = (i * chunk_size );
            sum_of_divisors_odd(chunk_size, m, sigma, primes);

            for(int j = 0; j < chunk_size/2; j++){

                unsigned long t;
                m = ( i * chunk_size + 1)+(j << 1);

                if(!(sigma[j] & 1)){ //true if even

                    t = 3 * sigma[j] - (m << 1);

                    while (t <= max_bound){
                        chunkCount = writePreimage(t, imageChunk, chunkCount);
                        t = (t << 1) + sigma[j];
                    }
                }

                if(sigma[j] == m+1){
                    chunkCount = writePreimage(m+1, imageChunk, chunkCount);
                }

                else if(m <= compositeBound && m != 1){ //s(1) = 0 is a special case which I need to think more about
                    #pragma omp critical (writeComp)
                    {
                        compSquares[compSquareCounter] = m;
                        compSquareCounter++;
                    }
                }
            }

            //finish up any leftover info
            if(chunkCount > 0) u_write( imageChunk, chunkCount);
              
        }

        free(sigma);

        double part2Start = omp_get_wtime();
        chunkCount = 0;

        //Im not sure this nowait is safe but it seems to work
        #pragma omp for nowait
        for(int i = 0; i < compSquareCounter; i++){

            unsigned long s_mSq = wheelDivSum(compSquares[i]*compSquares[i]);
            if(s_mSq <= max_bound) {
                 //if(s_mSq==0)printf("found, %lu", compSquares[i]*compSquares[i]);
                chunkCount = writePreimage(s_mSq, imageChunk, chunkCount);
            }
        }

        if(chunkCount > 0) u_write( imageChunk, chunkCount);
    }
    
    readCharac(); 
    printf("\n\nFinished in %f seconds\n", omp_get_wtime()-startTime);
}

long int findSize(const char *file_name){
    struct stat st; 
    if(stat(file_name,&st)==0) return st.st_size;
    else return -1;
}

unsigned long wheelDivSum(unsigned long n){

    if(n == 1 || n == 0) return 0;
    else if (n<0){
        printf("\nValue must be non-neg \n");
        exit(EXIT_FAILURE);
    }

    unsigned long savedN = n;
    unsigned long curr_sum, curr_term, res = 1;
  
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

void recParents(unsigned char * characArr){

    unsigned long accPreimages[256] = {0};
   // memset(accPreimages, 0, 256);

    
    accPreimages[0]++; //accounts for 5 which is the
    for(ulong i = 0;  i < max_bound/2; i++){
        
        accPreimages[characArr[i]]++;
    }

    printf("\nBound = %lu\n", max_bound);
    for(int j = 0;  j < 8; j++){
       printf("\n%d count: %lu, density = %f", j, accPreimages[j], (float)accPreimages[j]/max_bound);
    }
}

void u_read( char * fileName, size_t fileSize, ulong * buffer){
    size_t bufferSize = fileSize/sizeof(ulong);
    int fd;
    ulong *map;  /* mmapped array of int's */

    fd = open(fileName, O_RDONLY);
    if (fd == -1) {
	    perror("Error opening file %d for writing");
        printf("%s", fileName);
	    exit(EXIT_FAILURE);
    }

    #ifdef printFile
    printf("\n opening file %s with size %lu and %lu images\n", fileName, fileSize, bufferSize);
    #endif

    map = mmap(NULL, fileSize , PROT_READ , MAP_SHARED, fd, 0);
    close(fd);

    if (map == MAP_FAILED) {
	    perror("Error mmapping the file");
	    exit(EXIT_FAILURE);
    }

    for (ulong i = 0; i < bufferSize; i++){
        buffer[i] = map[i];
    }

    if (munmap(map, fileSize) == -1){
        perror("Error un-mmapping the file");
        exit(EXIT_FAILURE);
    } 
}


void u_write( ulong * buffer, size_t bufferSize){

    ulong grabbedTicker;
    char fileName[200];

    #pragma omp critical
    {
        
        grabbedTicker = fileTicker;
        fileTicker ++;
       
    }

    sprintf(fileName, "data/w%lu", grabbedTicker);

	size_t fileSize = (bufferSize)*sizeof(ulong)  ;

    int fd;
    int result;
    ulong *map;  /* mmapped array of int's */

    fd = open(fileName, O_RDWR | O_CREAT | O_TRUNC, (mode_t)0600);
    if (fd == -1) {
	    perror("Error opening file for writing");
	    exit(EXIT_FAILURE);
    }

    result = lseek(fd, (fileSize-1), SEEK_SET);
    if (result == -1) {
        close(fd);
	    perror("Error calling lseek() to 'stretch' the file");
	    exit(EXIT_FAILURE);
    }

    #ifdef printFile
    printf("\n writing file %s with size %lu and %lu images\n", fileName, fileSize, bufferSize);
    #endif

    result = write(fd, "", 1);
    if (result != 1) {
	    close(fd);
	    perror("Error writing last byte of the file");
	    exit(EXIT_FAILURE);
    }

    map = mmap(NULL, fileSize , PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
    close(fd);

    if (map == MAP_FAILED) {
	    perror("Error mmapping the file");
	    exit(EXIT_FAILURE);
    }

    for (int i = 0; i < bufferSize; i++) {
        map[i] = buffer[i];
       
    } 
       
    if (munmap(map, bufferSize) == -1) {
	    perror("Error un-mmapping the file");
        exit(EXIT_FAILURE);
    }
    char newName[200];

    sprintf(newName, "data/%lu", grabbedTicker);
    rename(fileName, newName);
   
}

void createCharac(){

    size_t bufferSize = (max_bound / 2);
	size_t fileSize = buffer_size*sizeof(char);

    int fd;
    int result;
    char *map;  /* mmapped array of int's */

    fd = open("data/charac", O_RDWR | O_CREAT | O_TRUNC, (mode_t)0600);
    if (fd == -1) {
	    perror("Error opening file for writing");
	    exit(EXIT_FAILURE);
    }

    result = lseek(fd, (fileSize-1), SEEK_SET);
    if (result == -1) {
        close(fd);
	    perror("Error calling lseek() to 'stretch' the file");
	    exit(EXIT_FAILURE);
    }

    #ifdef printFile
    printf("\n writing file %s with size %lu and %lu images\n", fileName, fileSize, bufferSize);
    #endif

    result = write(fd, "", 1);
    if (result != 1) {
	    close(fd);
	    perror("Error writing last byte of the file");
	    exit(EXIT_FAILURE);
    }

    map = mmap(NULL, fileSize , PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
    close(fd);

    if (map == MAP_FAILED) {
	    perror("Error mmapping the file");
	    exit(EXIT_FAILURE);
    }

    for (int i = 0; i < bufferSize; i++) {
        map[i] = 0;
       
    } 
       
    if (munmap(map, fileSize) == -1) {
	    perror("Error un-mmapping the file");
        exit(EXIT_FAILURE);
    } 
   
}

void updateCharac(int fileNum){

    
    
    size_t bufferSize = (max_bound / 2);
	size_t fileSize = buffer_size*sizeof(char);

    int fd;
    unsigned char *map;  /* mmapped array of int's */ 

    fd = open("data/charac", O_RDWR);
    if (fd == -1) {
	    perror("Error opening file %d for writing");
        printf("%s", "data/charac");
	    exit(EXIT_FAILURE);
    }

    #ifdef printFile
    printf("\n opening file %s with size %lu and %lu images\n", fileName, fileSize, bufferSize);
    #endif

    map = mmap(NULL, fileSize , PROT_READ | PROT_WRITE,  MAP_SHARED, fd, 0);
    //printf("%d\n", map[10]);
    close(fd);

    //recParents(map);

    if (map == MAP_FAILED) {
	    perror("Error mmapping the file");
	    exit(EXIT_FAILURE);
    }

    for (ulong i = 0; i < inputsize; i++){
       
        map[(input[i]/2)-1]++;
    }
   
    msync(map, fileSize, MS_SYNC);

    if (munmap(map, fileSize) == -1){
        perror("Error un-mmapping the file");
        exit(EXIT_FAILURE);
    } 
     
}

void readCharac(){
    size_t bufferSize = (max_bound / 2);
	size_t fileSize = buffer_size*sizeof(char);

    int fd;
    unsigned char *map;  /* mmapped array of int's */

    fd = open("data/charac", O_RDWR);
    if (fd == -1) {
	    perror("Error opening file %d for writing");
        printf("%s", "data/charac");
	    exit(EXIT_FAILURE);
    }

    #ifdef printFile
    printf("\n opening file %s with size %lu and %lu images\n", fileName, fileSize, bufferSize);
    #endif

    map = mmap(NULL, fileSize , PROT_READ | PROT_WRITE,  MAP_SHARED, fd, 0);
    //printf("%d\n", map[10]);
    close(fd);

    //recParents(map);

    if (map == MAP_FAILED) {
	    perror("Error mmapping the file");
	    exit(EXIT_FAILURE);
    }
     
    recParents(map);

    if (munmap(map, fileSize) == -1){
        perror("Error un-mmapping the file");
        exit(EXIT_FAILURE);
    } 
     
}

int writePreimage(ulong preimage, ulong * imageChunk, int chunkCount){

    imageChunk[chunkCount] = preimage;
    chunkCount++;
    
    if(chunkCount == buffer_size){
        u_write( imageChunk, chunkCount);
        chunkCount = 0;
    }

    

    return chunkCount;
}

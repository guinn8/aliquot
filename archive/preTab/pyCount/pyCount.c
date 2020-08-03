#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <fcntl.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/mman.h>
#include <math.h>
#include "../config/config.h"
#include "sieve.h"
#include <pari/pari.h>

//#define ONLY_STATS
#ifndef ulong
#define ulong unsigned long
#endif

#ifndef uint
#define uint unsigned int
#endif

#define TEMPSIZE 32768 

ulong max_solutions;
ulong blocksize;
ulong max_bound;


void writeToCharac(unsigned char * characFunc, ulong * buffer, ulong * preimage, size_t bufferSize);
ulong divSum(ulong num);
void u_write(const char * filename, ulong * buffer, ulong * preimage, size_t total, ulong min_value, omp_lock_t u_lock);
void runCharac(unsigned char * characFunc);

/*
***********************************************************
Huge big note:
When this program encounters a prime p it will record s(p) = p+1
THIS IS NOT TRUE!!!!!!
It is the case that s(p^2) = p+1
I have left this bug in because we would need 64bit ints to account for it
***********************************************************

This is the implementation of Pomerance-Yang Algorithm.
One provides a collection of [files] files containing sigma(m) for odd m.
The total number of odd m's in a single file is [blocksize].
The output is [files] files with even values of s(n) < 2*[limit].

I Pulled this from a different implmentation

This is the implementation of the Pomerance-Yang Algorithm.
See the paper of Pomerance and Yang, " Variant of a theorem of Erdős on
the sum-of-proper-divisors function", in Mathematics of Computation 83,
pp. 1903-1913, 2014.

Parameters:
1: [blocksize] - number of n's in the file;
2: [files_start] - start processing from this file;
3: [files_end] - stop processing at this file (excluded);
4: [sigma_files] - total number of files;
5: [max_solutions] - max number of n such that s(m) = n is recorded
6: [max_bound] - the first bound to record;

./py3Simple 10  32  1000000 (Successfully compiles)
./py3      [1]    [2][3] [4] [5] [6] 
*/
int main(int argc, char * argv[]){

	if(system("exec rm -r ../data/*") == -1) printf("error");

	if (argc < 3){
		printf("./py3Simple  [sigma_files] [max_solutions] [max_bound] ...\n");
		exit(1);
	}

	const ulong sigma_files = atol(argv[1]);
	max_solutions = atol(argv[2]);
	 max_bound = atol(argv[3]);
	blocksize = max_bound/sigma_files;
	
	const char * u_folder = "../data";


	

	if ((blocksize & 3) != 0){
	 	perror("The quantity [max_bound] / [sigma_files] has to be divisible by 4.\n");
	 	exit(1);
	}

	#ifdef ONLY_STATS
	//since we are only doing even numbers accArray[i] = i*2 preimages
	unsigned char accArray[max_bound/2];
	for(int i = 0; i < max_bound/2; i++){
		accArray[max_bound/2] = 0;
	}
	#endif



	// threads
	const int n = omp_get_max_threads();
	printf("Total threads: %d.\n", n);
	printf("Max preimages: %lu.\n", max_solutions);

	//write to config file
	const ulong arg_vec[6]={blocksize, 00, 00, sigma_files, max_solutions, max_bound};
	writeConfig(arg_vec);

	// enumeration of untouchable numbers
	omp_lock_t u_lock[sigma_files];
	ulong u_min[sigma_files];

	char filename[240];

	for (ulong f = 0; f < sigma_files; f++){

		sprintf(filename, "%s/u%lu", u_folder, f);

		if (access(filename, F_OK) == -1){
			close(creat(filename, 0755));
		}

		int truncRet = truncate(filename, max_solutions * blocksize * sizeof(ulong)); //I believe this function sets the size of the filename
		if(truncRet == -1) printf("Something went wrong");

		omp_init_lock(&u_lock[f]);
		u_min[f] = ((f * blocksize) << 1) + 2;
	}

	
	const ulong max_prime = (ulong) sqrt((double) (2 * blocksize * sigma_files));
	const ulong prime_bound = (ulong) (1.25506 * (max_prime + 1) / log(max_prime + 1));

	uint * primes = (uint *) malloc(sizeof(uint) * (prime_bound + 1));
	prime_sieve(max_prime, primes);

	void * reallocRet = realloc(primes, sizeof(uint) * (primes[0] + 1));
	if(reallocRet == NULL) printf("Something went wrong");

	struct timeval begin_loop, end_loop;
	double exec_time_loop;
	gettimeofday(&begin_loop, NULL);
	
	#pragma omp parallel  //This block of code is executed in parallel
	{	
		ulong sigma[TEMPSIZE];
		ulong t, f0, i, j, n_cur, m, u_cur_file;
		ulong u_count[sigma_files];

		// pointer to array of size sigma_files (unalloced) array of unsigned longs
		ulong * imageBuffer[sigma_files]; 
		ulong * preimageBuffer[sigma_files];
	
		for (f0 = 0; f0 < sigma_files; f0++){
			//u_buffer is being filled up with pointers to arrays? of unsigned longs
			//these arrays have size TEMPSIZE * unsigned long, temp size is a macro for 32768
			imageBuffer[f0] = (ulong *) malloc(TEMPSIZE * sizeof(ulong));
			preimageBuffer[f0] = (ulong *) malloc(TEMPSIZE * sizeof(ulong));
		}

		char filename[240];

		

		#pragma omp for schedule(dynamic) //dynamic opens up a reuntime queue to pass around the loop seems slightly faster

		for (int ind = 0; ind < sigma_files/2; ind++){ //This loop needed to cover all m's in bound

			// Fills the array at u_count with sigma_files * sizeof(ulong) number of zeros
			// sizeof(ulong) = 8 (bytes)
			memset(u_count, 0, sigma_files * sizeof(ulong)); 

			for (f0 = 0; f0 < sigma_files; f0++) memset(imageBuffer[f0], 0, TEMPSIZE * sizeof(ulong));
				
			m = ((ind * blocksize) << 1) + 1;
			
			for (i = 0; i <= blocksize / TEMPSIZE; i++){
				
				//This function writes the sum of diviors every odd number to array sigma
				//sigma is TEMPSIZE large
				sum_of_divisors_odd(TEMPSIZE << 1, m - 1, sigma, primes);
				
				for (j = 0; j < TEMPSIZE; j++, m += 2){

					
					// m is either a) even and neither a square nor twice a square; OR b) m is an odd square
					if ((sigma[j] & 1) == 0){ // checks if sum_div(m) is odd
					
						// m is prime
						if (sigma[j] == m + 1 && (m + 1) <= max_bound){
						
							

							
							t = m + 1;

							u_cur_file = ((t >> 1) - 1) / blocksize;
							
							if (u_count[u_cur_file] == TEMPSIZE){		
							

								#ifdef ONLY_STATS
								writeToCharac(accArray, imageBuffer[u_cur_file], preimageBuffer[u_cur_file], TEMPSIZE);	
								#else						
								sprintf(filename, "%s/u%lu", u_folder, u_cur_file); //formats the file name and writes to filename
								u_write(filename, imageBuffer[u_cur_file], preimageBuffer[u_cur_file], TEMPSIZE, u_min[u_cur_file], u_lock[u_cur_file]);			
								#endif
								u_count[u_cur_file] = 0;
							}

							imageBuffer[u_cur_file][u_count[u_cur_file]] = t;

							//********************************************
							//This should really be:  u_preimage[u_cur_file][u_count[u_cur_file]++] = m*m;
							//I have left it as is, It would require 64bit ints 
							//********************************************
							preimageBuffer[u_cur_file][u_count[u_cur_file]++] = m; 
						}

						
						n_cur = (m << 1);
					
						t = 3L * sigma[j] - (m << 1);
					
						while (t <= max_bound){ //This runs through that odd recurrance
							//printf("t: %lu m: %lu\n", t, m);
							
							u_cur_file = ((t >> 1) - 1) / blocksize;

							if (u_cur_file < (ind >> 1)){

								#pragma omp critical
								{
									printf("Error: the value of u_cur_file is too small.\nn=%lu, s(n)=%lu, u_cur_file=%lu\n", n_cur, t, u_cur_file);
									fflush(stdout);
								}
							}

							if (u_count[u_cur_file] == TEMPSIZE){
								#ifdef ONLY_STATS
								writeToCharac(accArray, imageBuffer[u_cur_file], preimageBuffer[u_cur_file], TEMPSIZE);	
								#else
								sprintf(filename, "%s/u%lu", u_folder, u_cur_file);
								u_write(filename, imageBuffer[u_cur_file], preimageBuffer[u_cur_file], TEMPSIZE, u_min[u_cur_file], u_lock[u_cur_file]);	
								#endif
								u_count[u_cur_file] = 0;
							}

							imageBuffer[u_cur_file][u_count[u_cur_file]] = t;
							preimageBuffer[u_cur_file][u_count[u_cur_file]++] = n_cur;
							
							n_cur <<= 1;
							t = (t << 1) + sigma[j];
						}
					}
				}
			}
		
			//I believe this simply cleans up any leftover files that were not maxed in the recurrance
			for (int i = 0; i < sigma_files; i++){ // Sigma_files is the number of files written to as output
			
				// u_count is an array initalized to be sigma_files large
				// what does u_count contain
				if (u_count[i] > 0){	
					#ifdef ONLY_STATS
					writeToCharac(accArray, imageBuffer[i], preimageBuffer[i], u_count[i]);	
					#else
					sprintf(filename, "%s/u%d", u_folder, i);
					u_write(filename,  imageBuffer[i], preimageBuffer[i], u_count[i], u_min[i], u_lock[i]);	
					#endif				
				}
			}
		}
	}
	gettimeofday(&end_loop, NULL);
	exec_time_loop = (end_loop.tv_sec*1e6 + end_loop.tv_usec) - (begin_loop.tv_sec*1e6 + begin_loop.tv_usec);
	printf("\nPart 1 of Algorithm computed in %.3f sec.\n\n",exec_time_loop / 1000000.0);


	//timing
	struct timeval begin_total, end_total;
	double exec_time_total;
	gettimeofday(&begin_total, NULL);

	//This block runs through s(m^2) where m <= max_bound^(2/3)
	//Where m is odd and composite, See the P-Y algor. for more deets
	//My goal is to buffer the I/O into blocks
	int bufferBlock = 10000;
	int bound = pow(max_bound, .666666666666);
	ulong preimageBuffer[sigma_files][bufferBlock];//contains preimages s(preimage[i][j]) = imageBuffer[i][j]
	ulong imageBuffer[sigma_files][bufferBlock];
	int sizeCount[sigma_files];//counts the amount written to each array
	pari_init(max_bound, 0);
	for(int i = 0; i < sigma_files; i++) sizeCount[i] = 0;

	//run the odd m's upto max_bound^(2/3)
	//every odd comp. number is ensured have an even image under s	
	for(ulong m = 3; m < bound; m+= 2){
		if(isprime(stoi(m))) continue;
		ulong preimage = m*m;
		ulong image = gtolong(sumdiv(stoi(preimage))) - preimage;
		
		if(image > max_bound) continue; //check if within range or if prime 
		
		ulong fileNum = ((image >> 1) - 1) / blocksize;
		//printf("m: %lu s(%lu) = %lu blocksize: %lu fileNum: %lu\n",m, preimage,image, blocksize, fileNum);
		
		preimageBuffer[fileNum][sizeCount[fileNum]] = preimage;
		imageBuffer[fileNum][sizeCount[fileNum]] = image;
		sizeCount[fileNum]++;

		if(sizeCount[fileNum] == bufferBlock){

			#ifdef ONLY_STATS
			writeToCharac(accArray, imageBuffer[fileNum], preimageBuffer[fileNum], bufferBlock);	
			#else
			sprintf(filename, "%s/u%lu", u_folder, fileNum);
			u_write(filename, imageBuffer[fileNum], preimageBuffer[fileNum], bufferBlock, u_min[fileNum], u_lock[fileNum]);
			#endif

			for(int i = 0; i < bufferBlock; i++){
				preimageBuffer[fileNum][i] = 0;
				imageBuffer[fileNum][i] = 0;
			} 
			sizeCount[fileNum] = 0;
		}
	}
	pari_close();

	for(ulong i = 0 ; i <sigma_files && sizeCount[i] > 0; i++){
		#ifdef ONLY_STATS
			writeToCharac(accArray, imageBuffer[i], preimageBuffer[i], sizeCount[i]);
		#else
			sprintf(filename, "%s/u%lu", u_folder, i);
			u_write(filename,  imageBuffer[i], preimageBuffer[i],  sizeCount[i], u_min[i], u_lock[i]);
		#endif	
	}
	
	

	free(primes);


	gettimeofday(&end_total, NULL);
	exec_time_total = (end_total.tv_sec*1e6 + end_total.tv_usec) - (begin_total.tv_sec*1e6 + begin_total.tv_usec);
	printf("\nPart 2 of Algorithm computed in %.3f sec.\n\n",exec_time_total / 1000000.0);

	#ifdef ONLY_STATS
	runCharac(accArray);
	#endif
}

//This function actually writes the data to the output file
//const char * filename: non-const pointer to a const char (possibly array?)
//size_t size: unsigned integral data type returned by sizeof()
//ulong * buffer: Pointer to unsigned long (array?)
//ulong * preimage: Pointer to unsigned long (array?)
//ulong max_solutions: Direct pass of unsigned long
//size_t total: unsigned integral data type returned by sizeof()
//ulong min_value: Direct pass of unsigned long

//s(preimage[i]) = buffer[i]?
void u_write(const char * filename, ulong * buffer, ulong * preimage, size_t total, ulong min_value, omp_lock_t u_lock){

	omp_set_lock(&u_lock);
	int fd = open(filename, O_RDWR);
	int size = max_solutions * blocksize * sizeof(ulong);

	if (fd == -1){
		char err[360];
		sprintf(err, "%d: unable to open the file %s for writing.\n", omp_get_thread_num(), filename);
		perror(err);
		exit(1);
	}

	unsigned long * file = (unsigned long *) mmap(0, size, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);


	if (file == MAP_FAILED){
		char err[360];
		sprintf(err, "%d: unable to map the file %s for writing.\n", omp_get_thread_num(), filename);
		perror(err);
		exit(1);
	}

	ulong j, index;
	

	for (ulong i = 0; i < total; i++){
		index = max_solutions * ((buffer[i] - min_value) >> 1);
		
		//printf("\n buffer[%lu] = %lu, index = %lu  s(%lu) = %lu\n", i, buffer[i], index,preimage[i],buffer[i] );

		//So my theory right now is that each file is an array with blocks max_solutions large
		//So all m's such that s(n) = m are held in an array index 

	
		for (j = 0; j < max_solutions; j++){
			if (file[index + j] == 0){
				file[index + j] = preimage[i];
				break;
			}
		}

		if (j == max_solutions){
			char err[240];
			sprintf(err, "Not enough space to write all the preimages of %lu. Terminating.\n", buffer[i]);
			perror(err);
			exit(1);
		}
	}
	omp_unset_lock(&u_lock);
	munmap(file, size);
	close(fd);
}

ulong divSum(ulong num) { 
    // Final result of summation of divisors 
    ulong result = 0; 
  
    // find all divisors which divides 'num' 
    for (ulong i=2; i<=sqrt(num); i++) { 

        // if 'i' is divisor of 'num' 
        if (num%i==0) 
        { 
            // if both divisors are same then add 
            // it only once else add both 
            if (i==(num/i)) 
                result += i; 
            else
                result += (i + num/i); 
        } 
    } 
  
    // Add 1 to the result as 1 is also a divisor 
    return (result + 1); 
} 


//this function takes a buffer containing images and a preimage buffer such that s(preimage[i]) = buffer[i]
//and counts the number of preimages for each number
void writeToCharac(unsigned char * characFunc, ulong * buffer, ulong * preimage, size_t bufferSize){
	for(int i = 0; i < bufferSize; i++){
	
		characFunc[buffer[i]/2] =characFunc[buffer[i]/2] +1;
	
	}

}

void runCharac(unsigned char * characFunc){
	unsigned int statArray[max_solutions];
	for(int i = 0; i < max_solutions; i++){
		statArray[i] = 0;
	}

	for(int i = 0; i < max_bound/2; i++){
		statArray[characFunc[i]] +=1;
		//printf("%d has %d preimages \n", i*2, characFunc[i]);
	}

	for(int i = 0; i < max_solutions; i++){
		printf("Number of %d parent numbers: %d\n", i, statArray[i]);
	}
}


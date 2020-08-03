#include <fcntl.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <sys/mman.h>
#include <math.h>
#include <sys/time.h>

#include <string.h>
#include "../config/config.h"

int * getPreImageArray();
unsigned long * openFile();
int divSum();
const char dataFolderPath[] = "../data";

int main(int argc, char * argv[]){
    
    //timing
	struct timeval begin, end;
	double exec_time;
	gettimeofday(&begin, NULL);

    //read from the config file found in the data folder and assign variables
    
    int args[10];
    readConfig(6, args);
    unsigned long limit, num_files, max_preimages, bound_0, imagesPerFile, preImageSpaces;

    if(args != NULL){
        num_files = args[3]; //number of files data is split between
        max_preimages = args[4]; //the maximuim number of preimages that can be stored
        bound_0 = args[5]; //this is the total solutions that have been computed 
        imagesPerFile = bound_0 / num_files;
        preImageSpaces = imagesPerFile * max_preimages;
    }else exit(EXIT_FAILURE);
    
    //accumlates stats on number of preimages
    int accArr[max_preimages];
    for(int i = 0; i < max_preimages; i++) accArr[i] = 0;
        
    int image = 0;//The number which to loop for preimages

    //loop though all files
    for(int i = 0; i < num_files; i++){
        unsigned long * file = openFile(i, preImageSpaces);

        //loop through image "slots" in file ensuring that we dont exceed the total bound 
        for(int j = 0; j < imagesPerFile && image < bound_0; j++){
            image = (j + i * imagesPerFile)*2 +2;//This conversion is magic to me
            int position = j * max_preimages;
            
            int preImages[max_preimages];
            getPreImageArray(preImages, file, max_preimages, position);

            int index = 0;
           // if(preImages[index] == 0)printf("%d\n", image);
            //printf("The preimages of %d: [ ", image);
            while(preImages[index] != 0 && index < max_preimages){
              // printf("%d ", preImages[index]);
              // if(divSum(preImages[index]) == 1) printf("prime");
                index++;
            }
            
          // printf("]\n");
            
            accArr[index]++;
        }
    }

    for(int i = 0; i < max_preimages; i++){
            printf("Ratio %d parent numbers: %f count: %d\n", i, (float)accArr[i]/(bound_0), accArr[i]);
    }

    gettimeofday(&end, NULL);
	exec_time = (end.tv_sec*1e6 + end.tv_usec) - (begin.tv_sec*1e6 + begin.tv_usec);
	printf("\nComputed in %.3f sec.\n\n", exec_time / 1000000.0);
}

//This function returns an integer array of preimages of some image
//int * preImages: pointer to store the preimages in (size = max_preimages)
//unsigned * long file: memory mapped file from which to read
//int max_sol: maximium number of preimages that can be written to one images "slot"
//int position: the 0th index of the images "slot" for preimages
int * getPreImageArray(int * preImages, unsigned long * file, int max_sol, int position){
    int countPreimages = 0;
    
    for (int w = 0; w < max_sol; w++) preImages[w] = 0;

    for (int q = 0; q < max_sol; q++){

        if(file[position+q] != 0){
            preImages[countPreimages] = file[position+q];
            countPreimages++;
        } 
    }
}

//This function returns a memory map of the file at index in dataFolderPath
//int index: index of file to open
//int preImageSpaces: total number of preimages that could be possibly written to the file  
unsigned long * openFile(int index, int preImageSpaces){

    //construct file name and open file
    char filename[80];
    sprintf(filename,"%s/u%d", dataFolderPath, index);
   
    int fd = open(filename, O_RDWR);

    if (fd == -1){
        char err[240];
        sprintf(err, "unable to open filefor writing.\n");
        perror(err);
        exit(1);
    }

    unsigned long * file = (unsigned long *) mmap(0, preImageSpaces * sizeof(unsigned long), PROT_READ, MAP_SHARED, fd, 0);
    close(fd);
    return file;
}

int divSum(int num) { 

    if(num == 1) return 0;
    // Final result of summation of divisors 
    int result = 1; 
  
    // find all divisors which divides 'num' 
    for (int i=2; i<=sqrt(num); i++) { 

        // if 'i' is divisor of 'num' 
        if (num%i==0) { 
            // if both divisors are same then add 
            // it only once else add both 
            if (i==(num/i))  result += i; 
            else result += (i + num/i); 
        } 
    } 
    return (result); 
} 


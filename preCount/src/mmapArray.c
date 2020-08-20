#include "mmapArray.h"


#ifndef MAP_HUGETLB
#define MAP_HUGETLB 0x40
#endif



//Creates and returns a pointer to a memory mapped char array
//This function is only recommended for creating a new mmap'ed array
//Be sure to call closeByteArray() once done with the file
unsigned char * createByteArray(char * fileName, size_t arraySize){

    int fd, result;
    unsigned char *map; 

    size_t fileSize = arraySize ;

   

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

    result = write(fd, "", 1);
    if (result != 1) {
	    close(fd);
	    perror("Error writing last byte of the file");
	    exit(EXIT_FAILURE);
    }
    
    //We need MAP_SHARED to push changes to the actual file on disk
    map = mmap(NULL, fileSize , PROT_READ | PROT_WRITE, MAP_SHARED , fd, 0);
    madvise(map, fileSize, MADV_SEQUENTIAL);
    close(fd);

    if (map == MAP_FAILED) {
	    perror("Error mmapping the file");
	    exit(EXIT_FAILURE);
    }

    return map;   
}

//memory maps an existing file containing a byte array and returns the pointer
//This function will not create a new file
//Be sure to call closeByteArray() once done with the file
unsigned char * openByteArray(char * fileName, size_t arraySize){

    size_t fileSize = arraySize * sizeof(unsigned char);
    
    int fd;
    unsigned char *map;  /* mmapped array of int's */

    fd = open(fileName, O_RDWR);
    if (fd == -1) {
	    perror("Error opening file %d for writing");
        printf("%s", fileName);
	    exit(EXIT_FAILURE);
    }

    map = mmap(NULL, fileSize ,  PROT_READ | PROT_WRITE , MAP_SHARED, fd, 0);
    close(fd);

    if (map == MAP_FAILED) {
	    perror("Error mmapping the file");
	    exit(EXIT_FAILURE);
    }

    return map;
}

void closeByteArray(unsigned char *map, size_t arraySize){
    if (munmap(map, arraySize * sizeof(unsigned char)) == -1) {
	    perror("Error un-mmapping the file");
        exit(EXIT_FAILURE);
    }
    
}


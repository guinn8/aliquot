#include <stdio.h>
#include <stdlib.h>
#include <sys/mman.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <dirent.h>

#ifndef MMAP_ARRAY_H_
#define MMAP_ARRAY_H_

unsigned char * openByteArray(char * fileName, size_t arraySize);
unsigned char * createByteArray(char * fileName, size_t arraySize);
void closeByteArray(unsigned char *map, size_t arraySize);


#endif
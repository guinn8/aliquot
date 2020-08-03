py3Simple was authored by Anton Mosunov and modified by Gavin Guinn. This program uses the Pomerance-Yang algorithm (*read the notes*) to calculate all m such that s(m) = n for all even n upto
a specifed bound. These solutions are saved in a number of binary file in the 'preimageTabulator/data' folder. 

    Compliling and Running py3Simple:
        -Dependancies: https://github.com/hyperrealm/libconfig
        -py3Simple uses *nix system calls to delete files, it has only be tested on unbuntu.
        -The make file should be suffient for complilation. In py3Simple just call 'make' without additional arguments.
        -Running py3Simple requires a number of command line arguments
        - ./py3Simple [blocksize] [files_start] [files_end] [sigma_files] [max_solutions] [max_bound]

    Parameters:
        1: [files_start] - This argument is left over from a previous iteration of this program that read s(n) data from files instead of computing it on the fly. 
                           it does have an impact but I am not sure what exactly it does
        2: [files_end] - Same deal as file_start
        3: [sigma_files] - total number of  output files to produce;
        4: [max_solutions] - max number of m such that s(m) = n is recorded
        5: [max_bound] - check for all preimages below this bound 

    Example:
        ./py3Simple 0 8 10 32 1000000

        These arguments produce 10 files named respectivly u0, u1,..., u9. In each of the files you will find 1000000/10 = 100000 sets of m such that s(m)=n. Each of those sets will be 32 
        unsigned longs wide. 
        The impact of the [files_start] and [files_end] arguments are unclear but it appears that [files_end] should be >= 8.

    File Format:
        The output files are arranged as such:

        Each preimage is recorded as a unsigned long int (32bit). Each number n that we are computing the preimages of is 'allocated' a block of unsigned longs [max_soultions] in size.
        The program fills each block with preimages as they are discovered unused space is simply left as zeros. So if a number has no preimages you will find a block of unsigned longs [max_solutions]
        long containin just zeros. 

        I believe that the conversion from 'file index' to the number of preimages can be done as follows:

        image = (j + i * imagesPerFile)*2 +2;

        Where j is the number of images we are into the file. In other words j is the number of unsigned longs * [max_soultions] blocks that occured sofar in the file. (zero indexed)
        Where i is the suffix of the file being processed for instance for u1, i=1
        Where imagesPerFile = imagesPerFile = max_bound / num_files;

        An implementation of this can be found in pydata.

notes: 
    - [max_bound] / [sigma_files] has to be divisible by 4.

    - In Anton's commenting it was noted that:
    #######################################################################################
    The output is a collection of files with even values of
    s(n) < 2*[blocksize]*[sigma_files],
    where n is not a square of an odd composite number.
    In order to compute even s(n) when n is a square odd composite number,
    please use part 2 of the program.
    #######################################################################################
    It appears that this program does not implement the entire P-Y algorithm, that could account for the
    off statistics

    - For all primes p it is recorded that s(p) = p+1 THIS IS INCORRECT, However it is the case that
      s(p^2) = p+1. This problem is left in because it would require 64bit ints to repersent the square of each
      of these primes.   
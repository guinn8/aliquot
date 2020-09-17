This program counts "k-parent" aliquot numbers using the Pomerance-Yang algorithm

Denpendancies:
FLINT (actually optional undefine it in main.c)
MPFR (FLINT dependancy also optional)
OpenMP (Should be included with most systems)

To use:

1. choose desired defines in main.c
    - TIMING: ifdef the program will output information about about preimages as they are written. Also an estimation of remaining preimages
    - FLINT: ifdef the program will use FLINT to run the sum-of-proper-divisors function. FLINT appears to speed up the computation in most cases
    - ASSERT: ifdef enables a bunch of assertions checking that the program is behaving in a sane manner. Slows down runs considerably, leave undefined for production runs

2. ./make

3. ./main.exe [max_bound] OR ./main.exe [max_bound] [buffer_size]
    -max_bound: count all aliquot orphans upto this bound
        -must be even
        -must be less than 281474976700000 to avoid int wrap possibility
    -buffer_size: this controls how big of a buffer each thread gets before writing to central characterstic function. It does have an impact on runtime but I dont know how to optimize. If not used the program uses a default vale

4. find the result of the computation in counts.csv 



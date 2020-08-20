# aliquot
Numerical investigation of Aliquot seqs

This repo contains 3 programs

1. preCount implements the Pomerance-Yang algorithm to calculate the number of preimages under the sum-of-divisors function for every value upto some bound. The program does this entirely in system memory to avoid time expensive disk I/O. The program output the counts of k-parent numbers into counts.csv.

2. preEnnum implements the simple ennumeration of s(n) for all n's such that n's image under s(.) lands below the bound

3. preEst implements a generaliztion of Pollack/Pomerance conj. density for non-aliquots.

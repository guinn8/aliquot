# aliquot
Numerical investigation of Aliquot seqs

This repo contains 3 programs

1. preCount implements the Pomerance-Yang algorithm to calculate the number of preimages under the sum-of-divisors function for every value upto some bound. The program does this entirely in system memory to avoid time expensive disk I/O. The program output the counts of k-parent numbers into counts.csv.

2. preEnnum implements the simple ennumeration of s(n) for all n's such that n's image under s(.) lands below the bound

3. preEst implements a generaliztion of Pollack/Pomerance conj. density for non-aliquots.

See aliquot_parents.pdf for a (unfinished) explaination of the math for the estimator and the generalization.

See the following papers for information about aliquot orphans:

K. Chum, R. K. Guy, Jr. M. J. Jacobson, and A. S. Mosunov,Numerical and statisticalanalysis of aliquot sequences, Experimental Mathematics0(2018), no. 0, 1–12.

Paul Pollack and Carl Pomerance,Some problems of erdős on the sum-of-divisorsfunction, Transactions of the American Mathematical Society, Series B3(2016), 1–26.

CARL POMERANCE and HEE-SUNG YANG,Variant of a theorem of erdŐs onthe sum-of-proper-divisors function, Mathematics of Computation83(2014), no. 288,1903–1913.

Anton's Untouchable Counts
Bound,                  Untouchable Counts (even)
100000000000,           16988116408

200000000000,           34059307042

300000000000,           51156680232

400000000000,           68270208721

500000000000,           85395279510

600000000000,           102529360014 

700000000000,           119670797250

800000000000,           136818383893

900000000000,           153971157175

1000000000000,          171128671373

1099511627776,          188206399402

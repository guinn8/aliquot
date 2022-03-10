Aliquot Sequences: Computational Number-Theory
==========

```none
  ___ ______ __         __     _____        __ _       _ _        ___  
 |__ \____  / /    _____\ \   |_   _|      / _(_)     (_) |      |__ \ 
    ) |  / / /_   |______\ \    | |  _ __ | |_ _ _ __  _| |_ _   _  ) |
   / /  / / '_ \   ______ > >   | | | '_ \|  _| | '_ \| | __| | | |/ / 
  / /_ / /| (_) | |______/ /   _| |_| | | | | | | | | | | |_| |_| |_|  
 |____/_/  \___/        /_/   |_____|_| |_|_| |_|_| |_|_|\__|\__, (_)  
                                                              __/ |    
                                                             |___/  
```

Documentation
-------------

Documentation can be found [here](https://guinn8.github.io/aliquot/html/index.html) or `docs/index.html`

Repo can be found [here](https://github.com/guinn8/aliquot)

View data collected using the Pomerance-Yang algorithm [here](https://docs.google.com/spreadsheets/d/1VtGl4Ibcjozck1wItciN3UREl9lPy5-rlCBIPKwWrgc/edit?usp=sharing)

Introduction to Aliquot Sequences
------------

This repo contains my project in computational number-theory, studying [aliquot sequences](https://en.wikipedia.org/wiki/Aliquot_sequence) and properties of the [sum-of-proper-divisors](https://en.wikipedia.org/wiki/Divisor_function) function. 

Aliquot sequences are iterations of the sum-of-proper-divisors function, for example:

```none
s(12) = 1 + 2 + 3 + 4 + 6
s(16) = 1 + 2 + 4 + 8
s(15) = 1 + 3 + 5
s(9) = 1 + 3
s(4) = 1 + 2
s(3) = 1
```

But not all sequences are so timid, [for example](http://factordb.com/sequences.php?se=1&aq=276&action=last2):

```none
s(276) -> s(396) -> ...over 2000 iterations... -> s(646 × 10^212)
```

This sequence certainly seems to be tending towards infinity, in 1888 *Catalan-Dickson* conjectured:

> All aliquot sequences are bounded and thus must terminate or enter a cycle.

In 1965 Guy-Selfridge countered this conjecture:

> The Catalan-Dickson conjecture is false, perhaps most aliquot sequences are unbounded.

The following problem statement from the former [Dr. Richard Guy](https://en.wikipedia.org/wiki/Richard_K._Guy) motivates this work:

> Think of a number!! Say 36%, which is nice and divisible. It appears that about 36% of the even numbers are ”orphans”.
>
> Divide by 1. For about 36% of the (even) values of n there is just one positive integer m such that s(m) = n. These values of n have just one ”parent”.
>
> Divide by 2. About 18% of the even values of n have exactly two parents.
>
> Divide by 3. About 6% of the even values of n have three parents.
>
> Divide by 4. About 1.5% of the even values of n have just 4 parents.
>
> This suggests that 1/(e ·p!) of the even numbers have p parents. Experiments suggests that these values are a bit large for small values of p and a bit small for larger values of p.
>
> Can anything be proved?

**Definition:** n is a k-parent aliquot number if *s^{-1}(n)* has k solutions.

My Project
----------

[This presentation](https://github.com/guinn8/aliquot/blob/master/pdf/kparent_density_technical_presentation.pdf) details an extension of the [Pollack and Pomerance] which uses heuristics to estimate the density of k-parent numbers; [huer_model.c](https://guinn8.github.io/aliquot/html/huer__model_8c.html) computes this model.

| **Delta_k** | **k = 0** | **k = 1** | **k = 2** | **k = 3** | **k = 4** |
|-------------|-----------|-----------|-----------|-----------|-----------|
| y = 10^3    | 0.155     | 0.164     | 0.100     | 0.046     | 0.018     |
| y = 10^4    | 0.161     | 0.165     | 0.098     | 0.044     | 0.017     |
| y = 10^5    | 0.165     | 0.166     | 0.097     | 0.044     | 0.016     |
| y = 10^6    | 0.167     | 0.166     | 0.096     | 0.043     | 0.016     |
| y = 10^7    | 0.169     | 0.167     | 0.096     | 0.042     | 0.016     |

| 1/2(p! * e) | p = 0     | p = 1     | p = 2     | p = 3     | p = 4     |
|-------------|-----------|-----------|-----------|-----------|-----------|
| -           | 0.184     | 0.184     | 0.092     | 0.031     | 0.008      |

Where Delta_k is my model's prediction of k-parent numbers and 1/2(p! * e) is Dr. Guy's prediction.

[This document](https://github.com/guinn8/aliquot/blob/master/pdf/kparent_aliquot_interm_report.pdf) contains an run-down of effort to computationally verify this model. [pomyang_kparent](https://guinn8.github.io/aliquot/html/pomyang__kparent_8c.html) contains the implementation.

| **Bound** | **0-parent** | **1-parent** | **2-parent** | **3-parent** | **4-parent** | **5-parent** | **6-parent** | **7-parent** | **8-parent** |
|---|---|---|---|---|---|---|---|---|---|
| 10000 | 0.121 | 0.187 | 0.121 | 0.048 | 0.015 | 0.005 | 0.002 | 0.002 | 0.000 |
| 1000000 | 0.150 | 0.178 | 0.103 | 0.042 | 0.016 | 0.006 | 0.003 | 0.001 | 0.001 |
| 10000000 | 0.157 | 0.175 | 0.099 | 0.042 | 0.016 | 0.006 | 0.003 | 0.001 | 0.001 |
| 100000000 | 0.162 | 0.173 | 0.097 | 0.041 | 0.016 | 0.006 | 0.003 | 0.001 | 0.001 |
| 1000000000 | 0.166 | 0.171 | 0.096 | 0.040 | 0.015 | 0.006 | 0.002 | 0.001 | 0.001 |
| 10000000000 | 0.168 | 0.171 | 0.095 | 0.040 | 0.015 | 0.006 | 0.002 | 0.001 | 0.001 |
| 100000000000 | 0.170 | 0.170 | 0.094 | 0.040 | 0.015 | 0.006 | 0.002 | 0.001 | 0.001 |
| 1000000000000 | 0.171 | 0.170 | 0.094 | 0.040 | 0.015 | 0.006 | 0.002 | 0.001 | 0.001 |

Which contains the computed densities of even k-parent number less than bound.

Citations
-----------

- [Chum et al.] Chum, K., Guy, R. K., Jacobson, J. M. J., and Mosunov, A. S. (2018).
      Numerical and statistical analysis of aliquot sequences. Experimental Mathematics, 29(4):414–425.

Todo
-----

- Take additive_sieve off of doxygen exclude list.
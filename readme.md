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

Introduction
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

[This presentation](https://github.com/guinn8/aliquot/blob/master/pdf/kparent_density_technical_presentation.pdf) details an extension of the [Pollack and Pomerance] which uses heuristics to estimate the density of non-aliquots. [huer_model.c](https://guinn8.github.io/aliquot/html/huer__model_8c.html) computes this model.



[huer_model.c](https://guinn8.github.io/aliquot/html/huer__model_8c.html) computes 

Documentation can be found [here](https://guinn8.github.io/aliquot/html/index.html) or `docs/index.html`

This repo contains several project, most related to sieve the sum-of-proper divisors function. More details can be found in detailed src documentation. PDF's are my writings on the subject of k-parent aliquot numbers.

pomyang_kparent counts the occurrence of k-parent numbers, by far the most developed.
additive_sieve is a sum-of-proper divisors sieve of my own development.
huer_model computes a model for density of kparent numbers.
plot_aliquot_families is a experiment in plotting aliquot families.

[Accessible Presentation](https://raw.githubusercontent.com/guinn8/aliquot/master/pdf/kparent_density_basic_presentation.pdf)

Citations
-----------
 -  [Chum et al.] Chum, K., Guy, R. K., Jacobson, J. M. J., and Mosunov, A. S. (2018).
      Numerical and statistical analysis of aliquot sequences. Experimental Mathematics, 29(4):414–425.

Todo
-----
 - Take additive_sieve off of doxygen exclude list.
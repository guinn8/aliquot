import matplotlib.pyplot as plt
from sympy.ntheory import factorint
import numpy as np

def s(n):
    """ @brief Compute the sum-of-proper-divisors of n """
    prime_factorization = factorint(n)
    total_product = 1
    for p in prime_factorization:
        prime_exp_sum = 0
        for e in range(0, prime_factorization[p]+1): # need inclusive upper lim
            prime_exp_sum += p ** e
        total_product *= prime_exp_sum
    return total_product - n
  
def print_aliquot_seq(n, lim):
    """ @brief Prints aliquot seq starting at n, halts if s(n) > lim"""
    s_n = [n]
    while n > 0 and n < lim:
        n = s(n)
        if n in s_n:
            print("Repeats with", n)
            break
        s_n.append(n)
    return s_n

def ennum_sn(max):
    """ @brief Computes https://oeis.org/A001065 upto and including max"""
    sn_ennum = np.zeros(max)
    for i in range(1, max + 1):
        sn_ennum[i-1] = s(i)
    return sn_ennum

def get_preimages(sn_ennum, n):
    """ @brief Searches for preimages of n in an enumerated range of s(n)
        @note Pomerance-Yang algorithm could vastly speed this up"""
    pre_images = list(np.where(sn_ennum == n)[0])
    pre_images = [x+1 for x in pre_images] # gotta add 1 to offset indexes to numbers
    return pre_images

max_depth = -50
def recurse_preimages(ennum, image, preimages, ind):
    """ @brief Recursively explore preimage chains, aliquot families
        @param ennum: Enumerated range of sn
        @param image: Current image being connected to it's preimages
        @param preimages: List of preimages of image
        @param ind: Recursive depth """
    for j in preimages:
        plt.plot([ind-1, ind], [j, image])
        new_pre = get_preimages(ennum, j)
        if len(new_pre) and ind > max_depth:
            print("preimages {0} = {1}".format(j, new_pre))
            recurse_preimages(ennum, j, new_pre, ind -1)

def plot_aliquot_family(start, sn_ennum):
    """ @brief Recursively plot explore preimage chains or aliquot families"""
    seq = print_aliquot_seq(start, len(sn_ennum))
    print(seq)

    for i in range(1, len(seq)):
        if seq[i] == 1: # every prime is a preimage of 1 
            break

        pre = get_preimages(sn_ennum, seq[i])
        if pre: #non-aliquot check?
            if i in seq:
                pre.remove(seq[i - 1])

        plt.plot([i-1, i], [seq[i-1], seq[i]])
        if seq[i] > 1:
            recurse_preimages(sn_ennum, seq[i], pre, i)
    plt.show()
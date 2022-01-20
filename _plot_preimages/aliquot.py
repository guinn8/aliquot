import matplotlib.pyplot as plt
from math import sqrt
from sympy.ntheory import factorint

from numpy import asarray
from numpy import savetxt
import math



def s(n):
    prime_factorization = factorint(n)
    total_product = 1
    for p in prime_factorization:
        prime_exp_sum = 0
        for e in range(0, prime_factorization[p]+1): # need inclusive upper lim
            prime_exp_sum += p ** e
        total_product *= prime_exp_sum
    return total_product - n
  
def printAliquot(n, lim):
    s_n = [n]
    while n > 0 and n < lim:
        n = s(n)
        if n in s_n:
            print("Repeats with", n)
            break
        s_n.append(n)
    return s_n

def ennum_sn(max):
    max += 1
    s_n = [1 for i in range(max)]
    s_n[0] = 0
    s_n[1] = 0
    
    for i in range(2, math.floor(max / 2)):
        j = 2 * i
        while j < max:
            s_n[j] += i
            j += i
        
        file_size = 10**9
        if 2*i % file_size == 0 or 2*i == math.floor((max-1) / 2):
            print(str(2*i) + " 2i: "+ str(math.floor((max-1) / 2)))
            lower_slice = 2*i - file_size + 1
            file = asarray(s_n[lower_slice:2*i])
            savetxt("data/" + str(lower_slice) + '-' + str(2*i) +'.csv', file, delimiter=',', fmt='%u')
            print("data/" + str(lower_slice) + '-' + str(2*i) +'.csv')
    return s_n

# def print_abud_ratio(max):
#     for i in range(1, max):s
#         s_n = getSum(i)
#         ratio = s_n/i
#         f = '{0:.3g}'.format(ratio)
#         print(f)
#         print("num: " + str(i) + " s_n: " + str(s_n) + " ratio: " + str(ratio))





def get_preimages(sn_ennum, n):
    pre_images = []
    for j in range(len(sn_ennum)):
        if n == sn_ennum[j]:
            # if j != n:
            # if j % 2 == 0:
            pre_images.append(j)
            # plt.plot([ind-1, ind], [j, n])
            # plt.plot([i-1, i], [j, seq[i]], 'ro')
            # print(str(ind)+ ": " + str(n) + " index: "+ str(j))
    return pre_images

def recurse_pre(ennum, image, pre, ind):
    for j in pre:
        plt.plot([ind-1, ind], [j, image])
        new_pre = get_preimages(ennum, j)
        # print(j)
        # print(new_pre)
        # print()
        if len(new_pre) and ind > -50:
            recurse_pre(ennum, j, new_pre, ind -1)


# max = 10**6
# seq = printAliquot(276, max)
# print(seq)
ennum = ennum_sn(10**6)
# last_in_seq = None
# for i in range(1, len(seq)):
#     pre = get_preimages(ennum, seq[i])
#     if pre: #non-aliquot check?
#         if i in seq:
#             pre.remove(seq[i - 1])

    
#     print(seq[i])
#     # print(pre)
#     # print()

#     plt.plot([i-1, i], [seq[i-1], seq[i]])
#     if seq[i] > 1:
#         recurse_pre(ennum, seq[i], pre, i)




    
# plt.show()

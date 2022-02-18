#include "../inc/properSumDiv.h"

const int initWheel[3] = {2, 3, 5};
const int inc[8] = {4, 2, 4, 2, 4, 6, 2, 6};

uint64_t sumdiv_sigma(uint64_t n) {
    if (n == 0) {
        return 0;
    }

    // uint64_t savedN = n;
    uint64_t curr_sum = 1;
    uint64_t curr_term = 1;
    uint64_t res = 1;

    for (int i = 0; i < 3; i++) {
        if (n % initWheel[i] == 0) {
            do {
                n = n / initWheel[i];
                curr_term *= initWheel[i];
                curr_sum += curr_term;
            } while (n % initWheel[i] == 0);

            res *= curr_sum;
            curr_sum = curr_term = 1;
        }
    }

    uint64_t k = 7;
    int i = 0;
    while (k <= n) {
        if (n % k == 0) {
            do {
                n = n / k;
                curr_term *= k;
                curr_sum += curr_term;
            } while (n % k == 0);
            res *= curr_sum;
            curr_sum = curr_term = 1;
        } else {
            k += inc[i];
            if (i < 7)
                i++;
            else
                i = 0;
        }
    }
    res *= curr_sum;
    return res;
}

uint64_t sumdiv_s(uint64_t n) {
    return sumdiv_sigma(n) - n;
}

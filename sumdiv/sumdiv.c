#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "../sumdiv/sumdiv.h"
#include "../sumdiv/factor.h"

uint64_t sumdiv(uint64_t n);

// int main(int argc, char const *argv[]) {
//     if (argc != 3) {
//         assert(0);
//     }
//     uint64_t n = atol(argv[1]);
//     uint64_t sn = sumdiv(n);
//     printf("s(%ld) = %ld\n", n, sn);
//     return 0;
// }

// https://www2.math.upenn.edu/~deturck/m170/wk3/lecture/sumdiv.html
uint64_t sumdiv(uint64_t n) {
    mpz_t integ;
    mpz_init_set_ui(integ, n);
    struct factors fs = {0};
    factor(integ, &fs);
    uint64_t sn = 1;
    for (int i = 0; i < fs.nfactors; i++) {
        uint64_t exp = 1;
        uint64_t row_total = 0;
        for (size_t j = 0; j <= fs.e[i]; j++) {
            // printf("%ld,\t", exp);
            row_total += exp;
            exp *= mpz_get_ui(fs.p[i]);
        }
        sn *= row_total;
        // printf("row sum = %ld\n", row_total);
    }
    factor_clear(&fs);
    mpz_clear(integ);
    return sn;
}

uint64_t s(uint64_t n){
    return sumdiv(n) - n;
}
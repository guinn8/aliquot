#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>


#include "../pomyang_kparent/inc/sumdiv_util.h"

int main(int argc, char const *argv[])
{
    if (argc != 2) {
        assert(0);
    }
    uint64_t n = atol(argv[1]);
    uint64_t sn = sumdiv_sigma(n);
    printf("s(%ld) = %ld\n", n, sn);

    return 0;
}

/**
 * @file pomyang_kparent.c
 * @author Gavin Guinn (gavinguinn1@gmail.com)
 * @date 2021-2-17
 * 
 * @brief Counts the even numbers with k-preimages under the sum-of-proper-divisors function
 * @copyright Public Domain (Please credit me; if you find this code useful I would love to hear about your work!)
 * 
 * NOTE
 * -----
 * **The counts of non-aliquots reported by this program will be exactly one less than the published figures!**
 * This algorithm enumerates all *even* k-parent numbers, odd k-parent numbers exist but they are not considered.
 * The user is expected to adjust the count to consider *5* which is the only odd aliquot (see sect. 3, [Chum et al.]).
 * It may seem odd to only count even k-parent numbers while odd k-parents also exist,
 * this program must be seen in the context of using heuristics to approach the Guy-Selfridge conjecture,
 * see [Guy and Selfridge] and [Chum et al.] for details.
 * 
 * ALGORITHMS
 * ----------
 * This is a rewrite of Anton Mosunov's implementation of the tabulation of untouchable numbers, section 3 of [Chum et al.].
 * As described in the paper this program is capable of counting of aliquot untouchables and more generally counting numbers
 * with k preimages under the sum-of-proper-divisors (`sumdiv` / `s(n)`) function, these numbers are referred to as k-parent. The algorithm employed in
 * [Chum et al.] was originally developed in a paper of [Pomerance and Yang] where the author's describe a family of algorithms for enumerating
 * preimages, the algorithm for s(n) will be referred to as the Pomerance-Yang (`pomyang`) algorithm. The `pomyang` algorithm takes
 * the sum-of-divisors (sigma) for all odd numbers in the range (1, `bound`) as input, to produce these values a technique of [Moews and Moews]
 * for sieving ranges of sigma is used to efficiently generate these values. 
 * 
 * TECHNIQUE
 * ---------
 * The techniques used in this implementation broadly differ from those described in [Chum et al.], notably
 * no disc operations are used. The sum-of-divisors (sigma) for all odd numbers in the range (1, `bound`) are sieved
 * on-the-fly rather than ^pre-computed and stored. A buffer is required to hold the counts of preimages for all even numbers less than `bound`,
 * this can require huge amounts of memory (using 8 bit counter takes `~500gb` at `bound` of 2^40). When enumerating non-aliquots only 1-bit per
 * number is required, however `sizeof(bool) == sizeof(uint8_t)` as a byte is the smallest addressable unit of memory. A PackedArray data structure
 * is used to create a buffer that can store 1, 2, or 4 bits of information in exactly that much memory, vasty decreasing memory requirements.
 * Threads must share this buffer when running the otherwise embarrassingly parallel pomyang algorithm, this is accomplished by allocating
 * an array of locks to protect portions of the resource. Locking the entire buffer whenever a thread wants to record a image bottlenecks
 * the whole process. 
 *
 *  ^   [Chum et al.] indicates only sieving odd sigma upto `bound` / 2 but odd sigma upto `bound` is certainly a required input for the
 *      Pomerance-Yang algorithm
 * 
 * PERFORMANCE
 * ----------
 * The parameter num_locks should be set to the highest value possible given memory constraits to optimize for speed. 
 * The parameter seg_len has a far more ambiguous effect on performance, if set either to low or high it can seriously impact performance
 * see moews_moews_sieve.c for a more detailed analysis. 
 *  
 * CITATIONS
 * -----------
 *  -  [Chum et al.] Chum, K., Guy, R. K., Jacobson, J. M. J., and Mosunov, A. S. (2018).
 *      Numerical and statistical analysis of aliquot sequences. Experimental Mathematics, 29(4):414–425.
 * 
 *  -  [Pomerance and Yang] Pomerance, C. and Yang, H.-S. (2014).
 *      Variant of a theorem of Erdos on the sum-of-proper-divisors function. Mathematics of Computation, 83(288):1903–1913.
 * 
 *  -  [Moews and Moews] Moews, D. and Moews, P. C. (1991). A search for aliquot cycles below 10 10.
 *      Mathematics of Computation, 57(196):849.
 * 
 *  -  [Guy and Selfridge] Guy, R. K. and Selfridge, J. L. (1975). What drives an aliquot sequence?
 *      Mathematics of Computation, 29(129):101–101.
 *
 * NOTATION
 * ----------
 *  -  sumdiv / s(n) == sum-of-proper-divisors
 * 
 *  -  sigma == sum-of-divisors
 * 
 *  -  kparent == number with k preimages under s(n) 
 * 
 * TODO: Finish collecting data to 2^40
 * TODO: Revisit heuristic model computation code
 * TODO: Move all relevant .tex files into this repo and find a vsCode tex workflow
 * 
 */

#include "../inc/pomyang_kparent.h"

#include <assert.h>
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <alloca.h>
#include <time.h>

#include "../inc/moewsmoews_sieve.h"
#include "../inc/math_macros.h"

static inline void _record_image(uint64_t x, PackedArray *f);
static uint64_t *_tabulate_kparent(const PackedArray *f);
static void _print_config(const pomyang_config *cfg, double odd_comp_bound);

static bool quiet = false;

/* See header for documentation */
PackedArray *pomyang_algorithm(const pomyang_config *cfg) {
    assert(EVEN(cfg->bound));
    assert(cfg->bound > 0);
    assert(EVEN(cfg->seg_len));
    assert(cfg->seg_len > 0);
    assert(cfg->seg_len <= cfg->bound);
    assert(DIVIDES(cfg->seg_len, cfg->bound));
    assert(DIVIDES(cfg->preimage_count_bits, 32));
    const double odd_comp_bound = pow(cfg->bound, TWO_THIRDS);
    quiet = cfg->quiet;
    omp_set_num_threads(cfg->num_threads);
    _print_config(cfg, odd_comp_bound);

    sieve_config_t *sieve_cfg = moews_init_sieve(cfg->bound, cfg->seg_len);
    PackedArray *f = PackedArray_create(cfg->preimage_count_bits, cfg->bound / 2, cfg->num_locks);
#pragma omp parallel shared(f)
    {
        sieve_worker_t *worker = moews_init_worker(sieve_cfg);

#pragma omp for schedule(static, 1)
        for (size_t seg_start = 0; seg_start < cfg->bound; seg_start += cfg->seg_len) {
            moews_sieve_odd_standard(worker, seg_start);
            for (uint64_t m = seg_start + 1; m < seg_start + cfg->seg_len; m += 2) {
                uint64_t sigma_m = moews_lookup_sigma_m(worker, m);

                if (EVEN(sigma_m)) {
                    uint64_t t = (3 * sigma_m) - (2 * m);
                    while (t <= cfg->bound) {
                        _record_image(t, f);
                        t = (2 * t) + sigma_m;
                    }
                }

                if (moews_check_prime(worker, m)) {
                    _record_image(m + 1, f);
                }
            }

            if (seg_start < odd_comp_bound) {
                moews_sieve_odd_squared(worker, seg_start);
                uint64_t seg_bound = MIN((seg_start + cfg->seg_len), odd_comp_bound);
                for (uint64_t m = seg_start + 1; m < seg_bound; m += 2) {
                    uint64_t sumdiv_m_sq = moews_lookup_sigma_m_squared(worker, m) - SQUARE(m);

                    if (!moews_check_prime(worker, m) && (sumdiv_m_sq > 0) && (sumdiv_m_sq <= cfg->bound)) {
                        _record_image(sumdiv_m_sq, f);
                    }
                }
            }
        }

        destroy_worker(worker);
        // LOG(quiet, "thread %d completed in %.2fs\n", omp_get_thread_num(), omp_get_wtime() - thread_start);
    }
    printf("\n");
    moews_destroy_sieve(sieve_cfg);
    return f;
}

/* See header for documentation */
uint64_t *pomyang_count_kparent(const pomyang_config *cfg) {
    PackedArray *f = pomyang_algorithm(cfg);
    uint64_t *count = _tabulate_kparent(f);
    free(f);
    return count;
}

/* See header for documentation */
void print_to_file(pomyang_config *cfg, const char *filename, uint64_t *count, float runtime) {
    const size_t max_line_len = 10000;
    char *header_line = alloca(max_line_len * sizeof(char));
    memset(header_line, 0, max_line_len * sizeof(char));
    if (0 != access(filename, F_OK)) {
        snprintf(header_line, max_line_len, "timestamp, bound, segment_length, num_locks, timing, num_threads");
        for (size_t i = 0; i <= UINT8_MAX; i++) {
            snprintf(header_line + strlen(header_line), max_line_len - strlen(header_line), ", %ld", i);
        }
    }

    FILE *fp = fopen(filename, "a");
    fprintf(fp, "%s\n", header_line);
    fprintf(fp, "%ld, ", time(NULL));
    fprintf(fp, "%ld, ", cfg->bound);
    fprintf(fp, "%ld, ", cfg->seg_len);
    fprintf(fp, "%ld, ", cfg->num_locks);
    fprintf(fp, "%.2f, ", runtime);
    fprintf(fp, "%ld", cfg->num_threads);

    for (size_t i = 0; i <= UINT8_MAX; i++) {
        fprintf(fp, ", %ld", count[i]);
    }
    fclose(fp);
}

/**
 * @brief increments preimage counter at appropiate index in buffer 
 *
 * @param sumdiv_m result of s(m) = n to be recorded (ie. n has one more preimage)
 * @param f PackedArray holding preimage counters
 */
static inline void _record_image(uint64_t sumdiv_m, PackedArray *f) {
    DEBUG_ASSERT(assert(sumdiv_m > 0));
    DEBUG_ASSERT(assert(EVEN(sumdiv_m)));
    size_t offset = (sumdiv_m / 2) - 1;

    PackedArray_lock_offset(f, offset);
    uint8_t count = PackedArray_get(f, offset);
    if (count + 1 < (1 << f->bitsPerItem)) {  // 2^f->bitsPerItem, prevents overflow
        PackedArray_set(f, offset, count + 1);
    }
    PackedArray_unlock_offset(f, offset);
}

/**
 * @brief counts the occurrences of numbers with k-preimages (ie. there are n numbers with 0 preimages) 
 *
 * @param f PackedArray holding finalized preimage counts
 * @return uint64_t* buffer of length (UINT8_MAX + 1) with occurrence counts, must be free'd by caller 
 */
uint64_t *_tabulate_kparent(const PackedArray *f) {
    printf("\nTABULATION\n");
    double time_tabulate = omp_get_wtime();
    uint64_t *count = calloc(sizeof(uint64_t), UINT8_MAX + 1);  // we can have num_preimages == UINT8_MAX
    for (size_t i = 0; i < f->count; i++) {
        const uint8_t num_preimages = PackedArray_get(f, i);
        count[num_preimages]++;
    }

    printf("Time: %0.2fs\n", omp_get_wtime() - time_tabulate);
    printf("Odd k-parent count under %ld:\n", f->count * 2);
    for (size_t i = 0; i < 8; i++) {
        printf("%ld: %ld\n", i, count[i]);
    }

    return count;
}

/**
 * @brief prints config information
 *
 * @param cfg config struct (see definition)
 * @param odd_comp_bound computed bound relevant to pomyang algorithm
 */
void _print_config(const pomyang_config *cfg, double odd_comp_bound) {
    printf("\nPOMYANG CONFIG\n");
    printf("-> Using %ld bits per number, count between 0-%d preimages\n", cfg->preimage_count_bits, (1 << cfg->preimage_count_bits) - 1);
    printf("-> Bound = %ld\n", cfg->bound);
    printf("-> Segment Length = %ld\n", cfg->seg_len);
    printf("-> Number of locks = %ld\n", cfg->num_locks);
    printf("-> Number of segments = %ld\n", cfg->bound / cfg->seg_len);
    printf("-> Max number of threads = %d\n", omp_get_max_threads());
    printf("-> Bound^(2/3) = %.2f\n", odd_comp_bound);
    printf("\n");
    if (true == cfg->est_heap) {
        exit(EXIT_SUCCESS);
    }
}

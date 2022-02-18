/**
 * @file pom_yang_cli.c
 * @author Gavin Guinn (gavinguinn1@gmail.com)
 * @brief command line interface to interact with the pom_yang algorithm
 * @date 2022-02-02
 *
 * @copyright Public Domain (Please credit me; if you find this code useful I would love to hear about your work!)
 *
 */

#include <assert.h>
#include <getopt.h>
#include <locale.h>
#include <stdbool.h>
#include <stdio.h>
#include <time.h>

#include "../inc/pomyang_kparent.h"
#include "../inc/moewsmoews_sieve.h"

#define OUTPUT_FILE "dat/counts.csv"
#define BYTES_TO_GB 0.000000001

static struct option long_options[] = {
    {"bound", required_argument, NULL, 'b'},
    {"seg_len", required_argument, NULL, 's'},
    {"num_locks", required_argument, NULL, 'n'},
    {"est_heap", no_argument, NULL, 'e'},
    {"num_threads", required_argument, NULL, 't'},
    {"preimage_count_bits", required_argument, NULL, 'p'},
    {"quiet", required_argument, NULL, 'q'},
    {0, 0, 0, 0},
};

static void usage(void);
static void get_args(pomyang_config *cfg, int argc, char **argv);

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wformat"  // printf ' flag
int main(int argc, char **argv) {
    pomyang_config cfg = {0};
    get_args(&cfg, argc, argv);

    size_t total = 0;
    total += PackedArray_estimate_heap(cfg.preimage_count_bits, cfg.bound / 2, cfg.num_locks);
    total += moews_estimate_heap_usage(cfg.bound, cfg.seg_len, cfg.num_threads);
    printf("This configuration will use a minimum of: \n");
    setlocale(LC_NUMERIC, "");
    printf("\t%'ld B\n", total);
    printf("\t%.2f GB\n", total * BYTES_TO_GB);

    clock_t start = clock();
    uint64_t *count = pomyang_count_kparent(&cfg);
    print_to_file(&cfg, OUTPUT_FILE, count, (clock() - start) / CLOCKS_PER_SEC);
    free(count);
    exit(EXIT_SUCCESS);
}
#pragma GCC diagnostic pop

static void usage(void) {
    printf("\nFast and Memory effiecent implementation of the Pomerance-Yang algorithm enumerating preimages under s()\n");
    printf("Usage: [--bound=][--seg_len=][--num_locks=][(OPTIONAL)--est_heap][(OPTIONAL)--num_threads=][(OPTIONAL)--preimage_count_bits]\n");
    printf("Assertions protect this program from invalid input, if it blows an assert you need to modify your input.\n\n");
    printf("%s\n\n", bound_desc);
    printf("%s\n\n", preimage_count_bits_desc);
    printf("%s\n\n", seg_len_desc);
    printf("%s\n\n", num_locks_desc);
    printf("%s\n\n", est_heap_desc);
    printf("%s\n\n", num_threads_desc);
    printf("%s\n\n", preimage_count_bits_desc);
    printf("\n");
}

static void get_args(pomyang_config *cfg, int argc, char **argv) {
    cfg->preimage_count_bits = 8;
    cfg->est_heap = false;
    cfg->quiet = false;
    cfg->num_threads = 1;

    int opt_code, option_index;
    while (-1 != (opt_code = getopt_long(argc, argv, "", long_options, &option_index))) {
        switch (opt_code) {
            case 'b':
                cfg->bound = strtol(optarg, NULL, 10);
                break;
            case 's':
                cfg->seg_len = strtol(optarg, NULL, 10);
                break;
            case 'n':
                cfg->num_locks = strtol(optarg, NULL, 10);
                break;
            case 'e':
                cfg->est_heap = true;
                break;
            case 't':
                cfg->num_threads = strtol(optarg, NULL, 10);
                break;
            case 'p':
                cfg->preimage_count_bits = strtol(optarg, NULL, 10);
                break;
            case 'q':
                cfg->quiet = true;
                break;
            default:
                usage();
                assert(0);
        }
    }

    if (0 == cfg->bound ||
        0 == cfg->seg_len ||
        0 == cfg->num_locks) {  // required args
        usage();
        assert(0);
    }
}

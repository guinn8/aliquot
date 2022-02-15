/**
 * @file pom_yang_cli.c
 * @author Gavin Guinn (gavinguinn1@gmail.com)
 * @brief 
 * @date 2022-02-02
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#include <getopt.h>
#include <stdbool.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <assert.h>
#include <time.h>

#include "../inc/pom_yang.h"

#define OUTPUT_FILE "counts.csv"

static struct option long_options[] = {
    {"bound", required_argument, NULL, 'b'},
    {"seg_len", required_argument, NULL, 's'},
    {"num_locks", required_argument, NULL, 'n'},
    {"just_config", no_argument, NULL, 'j'},
    {"num_threads", required_argument, NULL, 't'},
    {"preimage_count_bits", required_argument, NULL, 'p'},
    {"quiet", required_argument, NULL, 'q'},
    {0, 0, 0, 0},
};

static void print_to_file(PomYang_config *cfg, uint64_t *count, float runtime);
static void usage(void);
static void get_args(PomYang_config *cfg, int argc, char **argv);

int main(int argc, char **argv) {
    PomYang_config cfg = {0};
    get_args(&cfg, argc, argv);

    clock_t start = clock();
    uint64_t *count = count_kparent_aliquot(&cfg);
    print_to_file(&cfg, count, clock() - start);
    exit(EXIT_SUCCESS);
}

static void usage(void) {
    printf("\nFast and Memory effiecent implementation of the Pomerance-Yang algorithm enumerating preimages under s()\n");
    printf("Usage: [--bound=][--seg_len=][--num_locks=][(OPTIONAL)--just_config][(OPTIONAL)--num_threads=][(OPTIONAL)--preimage_count_bits]\n\n");
}

static void print_to_file(PomYang_config *cfg, uint64_t *count, float runtime) {
    const size_t max_line_len = 10000;
    char *header_line = calloc(max_line_len, sizeof(char));
    if (0 != access(OUTPUT_FILE, F_OK)) {
        snprintf(header_line, max_line_len, "timestamp, bound, segment_length, num_locks, timing, num_threads");
        for (size_t i = 0; i <= UINT8_MAX; i++) {
            snprintf(header_line + strlen(header_line), max_line_len - strlen(header_line), ", %ld", i);
        }
    }

    FILE *fp = fopen(OUTPUT_FILE, "a");
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

static void get_args(PomYang_config *cfg, int argc, char **argv) {
    cfg->preimage_count_bits = 8;
    cfg->just_config = false;
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
            case 'j':
                cfg->just_config = true;
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

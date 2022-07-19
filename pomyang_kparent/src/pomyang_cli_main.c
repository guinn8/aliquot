/**
 * @file pomyang_cli_main.c
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
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <time.h>

#include "../inc/math_macros.h"
#include "../inc/moewsmoews_sieve.h"
#include "../inc/pomyang_kparent.h"

/** @brief CLI runs output file.*/
#define OUTPUT_FILE "dat/counts.csv"

static struct option long_options[] = {
    {"bound", required_argument, NULL, 'b'},
    {"seg_len", required_argument, NULL, 's'},
    {"num_locks", required_argument, NULL, 'n'},
    {"est_heap", no_argument, NULL, 'e'},
    {"num_threads", required_argument, NULL, 't'},
    {"preimage_count_bits", required_argument, NULL, 'p'},
    {"quiet", required_argument, NULL, 'q'},
    {"to_file", required_argument, NULL, 'f'},
    {0, 0, 0, 0},
};

static void usage(void);
static void get_args(pomyang_config_t *cfg, int argc, char **argv);
static void write_pomyang(PackedArray *arr, pomyang_config_t *cfg);

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wformat"  // printf ' flag

/** @brief Command line interface to the Pomerance-Yang algorithm. */
int main(int argc, char **argv)
{
    pomyang_config_t cfg = {0};
    get_args(&cfg, argc, argv);

    size_t total = 0;
    printf("\nESTIMATED MEMORY USAGE\n");
    total += PackedArray_estimate_heap(cfg.preimage_count_bits, cfg.bound / 2, cfg.num_locks);
    total += moews_estimate_heap_usage(cfg.bound, cfg.seg_len, cfg.num_threads);
    printf("This configuration will use a minimum of: \n");
    setlocale(LC_NUMERIC, "");
    printf("\t%'ld B\n", total);
    printf("\t%.2f GB\n", total * BYTES_TO_GB);

    clock_t start_time = clock();
    if (NULL == cfg.filename) {  // tabulate k-parent stats at specific bound
        uint64_t *count = pomyang_count_kparent(&cfg);
        print_to_file(&cfg, OUTPUT_FILE, count, (clock() - start_time) / CLOCKS_PER_SEC);
        free(count);
    } else {  // store k-parent count array to disk
        PackedArray *arr = pomyang_algorithm(&cfg);
        write_pomyang(arr, &cfg);
        free(arr);
    }
    exit(EXIT_SUCCESS);
}
#pragma GCC diagnostic pop

/**
 * @brief   Writes the PackedArray of even preimage counts to a set of files on disc
 *          The written files are named `start_end` and contain the counts of pre-images
 *          for even numbers between the (inclusive) bounds encoded as a uint8_t array.
 * 
 * @param arr   Contains preimage counts from a run of the Pomerance-Yang algorithm
 * @param cfg   Contains a filename which will contain the zipped data files and a 
 *              seg_len which determines the amount of results to be written in a single file    
 */
static void write_pomyang(PackedArray *arr, pomyang_config_t *cfg)
{
    // ! delete any directory with the name supplied
    char rmdir[128];
    snprintf(rmdir, sizeof(rmdir), "rm -rf %s", cfg->filename);
    assert(-1 != system(rmdir));
    mkdir(cfg->filename, 0777);

    // open a configuration file that will contain the name's of all data files in directory
    char conf_file[128];
    snprintf(conf_file, sizeof(conf_file), "%s/conf.txt", cfg->filename);
    FILE *conf = fopen(conf_file, "w");
    assert(conf);

    size_t num_segs = cfg->bound / cfg->seg_len;
    for (size_t i = 0; i < num_segs; i++) {
        size_t start = (i * cfg->seg_len) + 2;
        size_t end = (i + 1) * cfg->seg_len;

        // open a new data file within directory and write the down-casted array to disc
        char data_file_path[128];
        snprintf(data_file_path, sizeof(data_file_path), "%s/%ld_%ld", cfg->filename, start, end);
        fprintf(conf, "%s\n", data_file_path);

        // unpack a piece of the PackedArray and downcast to uint8_t
        uint32_t *array = calloc(cfg->seg_len / 2, sizeof(uint32_t));
        PackedArray_unpack(arr, F_OFFSET(start), array, cfg->seg_len / 2);
        uint8_t *array_downcast = calloc(cfg->seg_len / 2, sizeof(uint8_t));
        for (size_t j = 0; j < cfg->seg_len / 2; j++) {
            array_downcast[j] = (uint8_t)array[j];
        }
        free(array);

        // write the downcasted array to disc
        FILE *data_file = fopen(data_file_path, "wb");
        assert(data_file);
        assert(cfg->seg_len / 2 == fwrite(array_downcast, sizeof(uint8_t), cfg->seg_len / 2, data_file));
        free(array_downcast);
        fclose(data_file);

        // zip up files
        char zip_file_cmd[512];
        snprintf(zip_file_cmd, sizeof(zip_file_cmd), "gzip %s", data_file_path);
        assert(-1 != system(zip_file_cmd));
    }
    fclose(conf);
}

/** @brief How to use this interface.*/
static void usage(void)
{
    printf("\nINVALID INPUT!\n");
    printf("\nFast and Memory effiecent implementation of the Pomerance-Yang algorithm enumerating preimages under s()\n");

    // TODO(Gavin): autogenerate usage string based on long_options
    printf("Usage: [--bound=][--seg_len=][--num_locks=][(OPTIONAL)--est_heap][(OPTIONAL)--num_threads=][(OPTIONAL)--preimage_count_bits][(OPTIONAL)--to_file=\"filename\"]\n");

    printf("Assertions protect this program from invalid input, if it blows an assert you need to modify your input.\n\n");
    printf("See https://guinn8.github.io/aliquot/html/structpomyang__config.html for documentation.\n\n");
    printf("EXAMPLE:\nmake cli && ./bin/cli --bound=$((10**9)) --seg_len=$((10**6)) --num_locks=$((10**7)) --num_threads=12  --preimage_count_bits=8\n\n");
}

static void get_args(pomyang_config_t *cfg, int argc, char **argv)
{
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
        case 'f':
            strncpy(cfg->filename, optarg, sizeof(cfg->filename) - 1);
            break;
        default:
            usage();
            exit(0);
        }
    }

    if (0 == cfg->bound ||
        0 == cfg->seg_len ||
        0 == cfg->num_locks) {  // required args
        usage();
        exit(0);
    }
}

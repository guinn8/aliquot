/**
 * @file pomyang_readfile.c
 * @author Gavin Guinn (gavinguinn1@gmail.com)
 * @brief Reads k-parent count array file into memory
 * @date 2022-07-10
 *
 * @copyright Copyright (c) 2022
 *
 */

#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>

int main(int argc, char **argv)
{
    if (argc != 2) {
        printf("Usage: ./pomyang_readfile <filename>\n");
        exit(EXIT_FAILURE);
    }

    // open configuration file
    char *filename = argv[1];
    char conf_file[128];
    snprintf(conf_file, sizeof(conf_file), "%s/conf.txt", filename);
    FILE *conf = fopen(conf_file, "rb");
    assert(conf != NULL);

    uint64_t count[UINT8_MAX + 1] = {0};  // accumulates count of n-parent numbers
    char datafile_name[128];  // change fscanf parameters if this length changes
    while (EOF != fscanf(conf, " %128[^\n]", datafile_name)) {
        // unzip and open file
        char zip_file_cmd[512];
        snprintf(zip_file_cmd, sizeof(zip_file_cmd), "gzip -d %s", datafile_name);
        assert(-1 != system(zip_file_cmd));
        FILE *datafile = fopen(datafile_name, "rb");
        assert(datafile);

        // parse out start and end of data file
        char *tok;
        char datafile_name_parse[512];
        strncpy(datafile_name_parse, datafile_name, sizeof(datafile_name_parse));
        assert(NULL != (tok = strtok(datafile_name_parse, "/")));
        assert(NULL != (tok = strtok(NULL, "_")));
        size_t start = strtoul(tok, NULL, 10);
        assert(NULL != (tok = strtok(NULL, "_")));
        size_t end = strtoul(tok, NULL, 10);

        // read the file back into memory and re-zip file
        size_t array_count = (end - start + 2) / 2;  // count of even numbers between inclusive bounds
        uint8_t *array = calloc(array_count, sizeof(uint8_t));
        assert(array_count == fread(array, sizeof(uint8_t), array_count, datafile));
        fclose(datafile);
        snprintf(zip_file_cmd, sizeof(zip_file_cmd), "gzip %s", datafile_name);
        assert(-1 != system(zip_file_cmd));

        // accumulates count of n-parent numbers
        for (size_t i = 0; i < array_count; i++) {
            // printf("%ld: %d\n", i, array[i]);
            count[array[i]]++;
        }

        printf("\nOdd k-parent count under %ld:\n", end);
        for (size_t i = 0; i < 8; i++) {
            printf("%ld: %ld\n", i, count[i]);
        }
        free(array);
    }

    exit(EXIT_SUCCESS);
}

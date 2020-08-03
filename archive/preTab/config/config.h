#include <libconfig.h>
#include <stdlib.h>
#ifndef CONFIG_H
#define CONFIG_H
void readConfig(int count, int * args);
int writeConfig(const long unsigned args[]);
#endif /* CONFIG_H*/
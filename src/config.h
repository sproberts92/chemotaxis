#ifndef CONFIG_H
#define CONFIG_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "yajl/yajl_tree.h"

typedef struct{
	int dim;
	int arr_dim;
	int n_neigh;
	int iter;
	int age;
	int branch_window;
	double noise;

	int slow;
	int write_frames;
} config;

int check_parse_error(yajl_val node, char *errbuf, int v);
char *read_file(long int *length, const char *file_name);
int read_config(config *cf);
int get_json_int(yajl_val *node, const char *name);
double get_json_double(yajl_val *node, const char *name);

#endif

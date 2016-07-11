#include "config.h"

char *read_file(long int *length, const char *file_name)
{
	char *buffer = NULL;

	FILE *fp = fopen(file_name, "rb");
	if(!fp)
	{
		fprintf(stderr, "Could not open %s\n", file_name);
		*length = 0;
	}
	else
	{
		fseek(fp, 0, SEEK_END);
		long int string_size = ftell(fp);
		rewind(fp);

		buffer = malloc (sizeof(char) * (string_size + 1));
		long int read_size = fread(buffer, sizeof(char), string_size, fp);

		buffer[string_size] = '\0';

		*length = read_size;

		if(string_size != read_size)
		{
			fprintf(stderr, "Read failed.\n");
			free(buffer);
			buffer = NULL;
			*length = 0;
		}
	}

	fclose(fp);

	return(buffer);
}

int read_config(config *cf)
{
	char errbuf[1024];
	yajl_val node;

	long length;
	char *file_data = read_file(&length, "config/config.json");

	node = yajl_tree_parse((const char *) file_data, errbuf, sizeof(errbuf));

	if(check_parse_error(node, errbuf, 1)) return 0;

	cf->dim             = get_json_int(&node, "dim");
	cf->arr_dim         = cf->dim * cf->dim;

	cf->n_neigh         = get_json_int(&node, "n_neigh");
	cf->iter            = get_json_int(&node, "iter");
	cf->age             = get_json_int(&node, "age");
	cf->branch_window   = get_json_int(&node, "branch_window");
	cf->propagate_noise = get_json_double(&node, "propagate_noise");
	cf->select_noise    = get_json_double(&node, "select_noise");
	cf->total_act_target    = get_json_double(&node, "total_act_target");

	cf->slow            = get_json_int(&node, "slow");
	cf->write_frames    = get_json_int(&node, "write_frames");

	yajl_tree_free(node);

	free(file_data);

	return 1;
}

int get_json_int(yajl_val *node, const char *name)
{
	const char * path[] = {name, (const char *) 0 };
	yajl_val v = yajl_tree_get(*node, path, yajl_t_number);

	if(v)
		return (int)YAJL_GET_INTEGER(v);
	else
	{
		fprintf(stderr, "Can't find node %s\n", name);
		return -1;
	}
}

double get_json_double(yajl_val *node, const char *name)
{
	const char * path[] = {name, (const char *) 0 };
	yajl_val v = yajl_tree_get(*node, path, yajl_t_number);

	if(v)
		return YAJL_GET_DOUBLE(v);
	else
	{
		fprintf(stderr, "Can't find node %s\n", name);
		return -1.0f;
	}
}

int check_parse_error(yajl_val node, char *errbuf, int v)
{
	/* 0 for success, 1 for failure.
	 * makes a call like
	 *     if(check_parse_error(...)){ handle error }
	 * more natural.
	 * v argument is verbose, should error be printed to stderr */

	if (node == NULL)
	{
		if (v)
		{
			fprintf(stderr, "Parse error: ");

			if (strlen(errbuf))
				fprintf(stderr, " %s", errbuf);
			else
				fprintf(stderr, "unknown error.");

			fprintf(stderr, "\n");
		}

		return 1;
	}
	else
		return 0;
}

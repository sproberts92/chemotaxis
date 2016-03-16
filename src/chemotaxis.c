#pragma warning( disable : 4710 4711 4820 4996 )

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "mt19937ar.c"

typedef struct{
	const int dim;
	const int arr_dim;
	const int n_neigh;
	const int iter;
	const int age;
} config;

int ind(int x, int y, int d);
void get_neighbours(unsigned int *nb, int x, int y, int d, int td);
void write_array(FILE *stream, unsigned int *arr, unsigned int dim);
int modulo(int i, int n);
void propagate_1(unsigned int *lattice, unsigned int *lattice_t, unsigned int *lattice_th, unsigned int iter, config *cf);
void propagate_2(unsigned int *lattice, unsigned int *lattice_t, config *cf);
double getRandNum(void);
void initSeed(void);
void write_pcount(FILE *stream, unsigned int *arr, unsigned int dim);

int main(void)
{
	initSeed();

	// config cf = {15, 225, 6, 5000, 20};
	config cf = {50, 2500, 6, 5000, 50};

	if (cf.arr_dim != cf.dim * cf.dim)
	{
		fprintf(stderr, "Incorrect dimensions.\n");
		return EXIT_FAILURE;
	}

	unsigned int *lattice    = calloc(cf.arr_dim, sizeof(unsigned int));
	unsigned int *lattice_t  = calloc(cf.arr_dim, sizeof(unsigned int));
	unsigned int *lattice_th = calloc(cf.arr_dim, sizeof(unsigned int));

	lattice[1275] = 1;

	for (int iter = 2; iter < cf.iter + 2; ++iter)
	{
		// write_array(stdout, lattice_th, cf.dim);

		// for (int i = 0; i < 100000; ++i)
		// {
		// 	double *d = malloc(1000 * sizeof(double));
		// 	free(d);
		// }

		propagate_1(lattice, lattice_t, lattice_th, iter, &cf);
		propagate_2(lattice, lattice_t, &cf);

		char f_name[64];
		sprintf(f_name, "output/pcount_%d.dat", iter-2);
		FILE *fp = fopen(f_name, "w");

		write_pcount(fp, lattice_th, cf.dim);
		fclose(fp);
	}

	free(lattice);
	free(lattice_t);

	return EXIT_SUCCESS;
}

void propagate_1(unsigned int *lattice, unsigned int *lattice_t, unsigned int *lattice_th, unsigned int iter, config *cf)
{
	memset(lattice_t, 0, cf->arr_dim * sizeof(unsigned int));
 	unsigned int ct = 0;
	for (int i = 0; i < cf->dim; ++i)
	{
		for (int j = 0; j < cf->dim; ++j)
		{
			unsigned int index = modulo(ind(i,j,cf->dim), cf->arr_dim);
			unsigned int *elem = &lattice[index];

			if (*elem >= 1)
			{
				unsigned int *neighbours = malloc(cf->n_neigh * sizeof(unsigned int));
				get_neighbours(neighbours, i, j, cf->dim, cf->arr_dim);

				unsigned int highest = 0;
				unsigned int highest_index = 0;
				for (int n = 0; n < cf->n_neigh; ++n)
				{
					if (lattice_th[neighbours[n]] > highest && (int)lattice_th[neighbours[n]] < (int)iter - (int)cf->age)
					{
						highest = lattice_th[neighbours[n]];
						highest_index = n;
					}
				}

				if (highest > 0)
					lattice_t[neighbours[highest_index]] = 1;
				else
				{
					unsigned int rn;

					for (int ii = 0; ii < 36; ++ii)
					{
						rn = (unsigned int)(getRandNum() * 6.0f);
						if ((int)lattice_th[neighbours[rn]] < (int)iter - (int)cf->age)
							break;
					}

					lattice_t[neighbours[rn]] = 1;
				}

				lattice_th[index] = iter;

				free(neighbours);
				ct++;
			}
		}
	}
	if(ct == 0)
	{
		printf("Signal died.\n");
		exit(EXIT_FAILURE);
	}
}

void propagate_2(unsigned int *lattice, unsigned int *lattice_t, config *cf)
{
	for (int i = 0; i < cf->arr_dim; ++i)
		lattice[i] = lattice_t[i];
}

void write_array(FILE *stream, unsigned int *arr, unsigned int dim)
{
	system("cls");
	for (unsigned int i = 0; i < dim; ++i)
	{
		for (unsigned int j = 0; j < 2*((i+1) % 2); ++j)
			printf(" ");
		for (unsigned int j = 0; j < dim; ++j)
		{
			int index = i + dim * ((j + i/2) % dim);
			if (arr[index] == 0)
				fprintf(stream, "    ");
			else
				fprintf(stream, "%3d ", arr[index]);
		}

		fprintf(stream, "\n\n");
	}

	fprintf(stream, "\n");
}

void write_pcount(FILE *stream, unsigned int *arr, unsigned int dim)
{
	for (unsigned int i = 0; i < dim; ++i)
	{
		for (unsigned int j = 0; j < dim; ++j)
			fprintf(stream, "%d ", arr[i * dim + j]);
		fprintf(stream, "\n");
	}
}

int ind(int x, int y, int d)
{
	return (x * d) + y;
}

void get_neighbours(unsigned int *nb, int x, int y, int d, int td)
{
	nb[0] = modulo(ind(x + 1, y,     d), td);
	nb[1] = modulo(ind(x + 1, y + 1, d), td);
	nb[2] = modulo(ind(x,     y + 1, d), td);
	nb[3] = modulo(ind(x - 1, y,     d), td);
	nb[4] = modulo(ind(x - 1, y - 1, d), td);
	nb[5] = modulo(ind(x,     y - 1, d), td);
}

int modulo(int i, int n)
{
	if(i > -1) return i % n;
	else return(i % n + n);     // Incorrect implementation but faster, works in this case as input values are constrained such that the second % will never be needed.
//	else return(i % n + n) % n; // Proper implementation of mod function
}

void initSeed(void)
{
	init_genrand((unsigned long)time(NULL));
}

double getRandNum(void)
{
	return genrand();
}

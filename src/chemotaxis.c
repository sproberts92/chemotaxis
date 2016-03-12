#pragma warning( disable : 4710 4711 4820 4996 )

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "mt19937ar.c"

typedef struct{
	const int dim;
	const int arr_dim;
	const int n_neigh;
} config;


int ind(int x, int y, int d);
void get_neighbours(unsigned int *nb, int x, int y, int d, int td);
void write_array(FILE *stream, unsigned int *arr, unsigned int dim);
int modulo(int i, int n);
void collide(unsigned int *lattice, unsigned int *lattice_t, double *lattice_th, config *cf);
void propagate(unsigned int *lattice, unsigned int *lattice_t, config *cf);
unsigned int rotate(unsigned int word, unsigned int length, unsigned int n);
int popcountQbits(unsigned int x, int start, int end);
double getRandNum(void);
void initSeed(void);
void write_pcount(FILE *stream, unsigned int *arr, unsigned int dim);

int main(void)
{
	initSeed();

	config cf = {15, 225, 6};

	unsigned int *lattice   = calloc(cf.arr_dim, sizeof(unsigned int));
	unsigned int *lattice_t = calloc(cf.arr_dim, sizeof(unsigned int));
	double *lattice_th = calloc(cf.arr_dim, sizeof(double));

	for (int i = 0; i < cf.arr_dim; ++i)
		lattice_th[i] = 1.0f;

	lattice[42] = 63;

	// write_array(stdout, lattice, cf.dim);

	for (int iter = 0; iter < 500; ++iter)
	{
		collide(lattice, lattice_t, lattice_th, &cf);
		propagate(lattice, lattice_t, &cf);

		char f_name[64];
		sprintf(f_name, "output/pcount_%d.dat", iter);
		FILE *fp = fopen(f_name, "w");
		
		// write_array(stdout, lattice, cf.dim);
		write_pcount(fp, lattice, cf.dim);
		fclose(fp);
	}

	free(lattice);
	free(lattice_t);

	return EXIT_SUCCESS;
}

void collide(unsigned int *lattice, unsigned int *lattice_t, double *lattice_th, config *cf)
{
	memset(lattice_t, 0, cf->arr_dim * sizeof(unsigned int));

	for (int i = 0; i < cf->dim; ++i)
	{
		for (int j = 0; j < cf->dim; ++j)
		{
			unsigned int *neighbours = malloc(cf->n_neigh * sizeof(unsigned int));
			get_neighbours(neighbours, i, j, cf->dim, cf->arr_dim);

			for (int n = 0; n < cf->n_neigh; ++n)
			{
				unsigned int elem = lattice[neighbours[n]];
				unsigned int n_in = ((n + cf->n_neigh/2)%cf->n_neigh);
				if((elem >> n_in) & 1)
				{
					elem &= ~(1 << n_in);
					lattice_t[modulo(ind(i,j,cf->dim), cf->arr_dim)] |= 1 << n;
				}
			}

			unsigned int index = modulo(ind(i,j,cf->dim), cf->arr_dim);
			unsigned int elem = lattice_t[index];
			unsigned int pcount = popcountQbits(elem, 0, cf->n_neigh);

			if(pcount > 4)
			{
				int rn = (int)(getRandNum() * 6.0f);
				lattice_t[index] = 1 << rn;
			}
			else if (pcount < 1 && lattice_th[index] > 0.5 && getRandNum() > 0.99f)
			{
				for (int m = 0; m < 2; ++m)
				{
					int rn = (int)(getRandNum() * 6.0f);
					lattice_t[index] |= 1 << rn;
				}
				lattice_th[index] = 0.0;
			}
			else{
				lattice_th[index] += 0.2;
				lattice_t[index] =0;
			}

			free(neighbours);
		}
	}
}

void propagate(unsigned int *lattice, unsigned int *lattice_t, config *cf)
{
	for (int i = 0; i < cf->arr_dim; ++i)
		lattice[i] = rotate(lattice_t[i], cf->n_neigh, cf->n_neigh/2);

	for (int i = 0; i < 3000000; ++i)
	{
		double *d = malloc(1000 * sizeof(double));
		free(d);
	}
}

unsigned int rotate(unsigned int word, unsigned int length, unsigned int n)
{
	for(unsigned int i = 0; i < n; i++)
	{
		unsigned int bit = (word >> (length - 1)) & 1;
		word &= ~(1 << (length - 1));
		word <<= 1;
		if(bit) word |= 1;
		else word &= ~1;
	}
	return word;
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
			fprintf(stream, "%d ", popcountQbits(arr[i * dim + j], 0, 6));
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

int popcountQbits(unsigned int x, int start, int end)
{
	int count = 0;

	for(int i = start; i < end+1; i++)
		count += (x >> i) & 1;

	return count;
}

void initSeed(void)
{
	init_genrand((unsigned long)time(NULL));
}

double getRandNum(void)
{
	return genrand();
}

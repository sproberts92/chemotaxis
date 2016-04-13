#pragma warning( disable : 4710 4711 4820 4996 )

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "tinymt64.h"

// #define SLOW

tinymt64_t tinymt_gen;

typedef struct{
	const int dim;
	const int arr_dim;
	const int n_neigh;
	const int iter;
	const int age;
	const int target_activity;
	const double noise;
} config;

int ind(int x, int y, int d);
void get_neighbours(unsigned int *nb, int x, int y, int d, int td);
void write_array(FILE *stream, unsigned int *arr, unsigned int dim);
int modulo(int i, int n);
void propagate_1(unsigned int *lattice, unsigned int *lattice_t, unsigned int *lattice_th, unsigned int iter, config *cf);
void propagate_2(unsigned int *lattice, unsigned int *lattice_t, config *cf);
double getRandNum(void);
void initSeed(void);
void write_last_visited(FILE *stream, unsigned int *arr, unsigned int dim);
void wait_for_ms(clock_t wait_time);
int choose_site(unsigned int *neighbours, unsigned int *lattice_t, unsigned int *lattice_th, unsigned int iter, config *cf);

int main(void)
{
	initSeed();

	config cf = {250, 250*250, 6, 5000, 50, 50, 0.4f};

	if (cf.arr_dim != cf.dim * cf.dim)
	{
		fprintf(stderr, "Incorrect dimensions.\n");
		return EXIT_FAILURE;
	}

	/* lattice    - primary lattice
	 * lattice_t  - temporary propagation lattice
	 * lattice_th - stores last visited time */
	unsigned int *lattice    = calloc(cf.arr_dim, sizeof(unsigned int));
	unsigned int *lattice_t  = calloc(cf.arr_dim, sizeof(unsigned int));
	unsigned int *lattice_th = calloc(cf.arr_dim, sizeof(unsigned int));

	/* Put a signal in the centre of the lattice */
	lattice[cf.arr_dim/2 - cf.dim/2] = 1;

	for (int iter = 2; iter < cf.iter + 2; ++iter)
	{
		printf("\r%d", iter);
		// write_array(stdout, lattice_th, cf.dim);

		#ifdef SLOW
			wait_for_ms(100);
		#endif

		propagate_1(lattice, lattice_t, lattice_th, iter, &cf);
		propagate_2(lattice, lattice_t, &cf);

		/* Write out the last visited time to unique files, one per
		 * time step. Careful! Can produce a lot of output if left
		 * running for too long! */
		char f_name[64];
		sprintf(f_name, "output/pcount_%d.dat", iter-2);
		FILE *fp = fopen(f_name, "w");

		write_last_visited(fp, lattice_th, cf.dim);
		fclose(fp);
	}

	free(lattice);
	free(lattice_t);

	return EXIT_SUCCESS;
}

void propagate_1(unsigned int *lattice, unsigned int *lattice_t,
	unsigned int *lattice_th, unsigned int iter, config *cf)
{
	/* Signal is propagated into a tempoarary lattice, lattice_t. If this
	 * is not done then there can be confusion as the lattice is updated
	 * sequentially. Each signal should be propagated once per time step
	 * but without this method one signal can easily be updated many times
	 * per time step if it is propagated into a cell that has not yet been
	 * updated itself. */

	/* CLear the temporary lattice_t */
	memset(lattice_t, 0, cf->arr_dim * sizeof(unsigned int));

 	/* Used as a check to see if the signal is still alive, mainly
 	 * useful in debugging when introducing new propagation logic */
 	unsigned int ct = 0;

	for (int i = 0; i < cf->dim; ++i)
		for (int j = 0; j < cf->dim; ++j)
		{
			/* Shortcuts */
			unsigned int index = modulo(ind(i,j,cf->dim), cf->arr_dim);
			unsigned int *elem = &lattice[index];

			/* Focus on those lattice elements where there is a signal.
			 * Can probably get a big speed increase by keeping track of where the
			 * signal is and only touching that element. Oh well this is easier for now */
			if (*elem >= 1)
			{
				/* Build a vector containing indices of the surrounding elements */
				unsigned int *neighbours = malloc(cf->n_neigh * sizeof(unsigned int));
				get_neighbours(neighbours, i, j, cf->dim, cf->arr_dim);

				/* Search the surrounding elements to see which was visited the most recently,
				 * but not too recently (more than cf.age time steps ago) */
				unsigned int highest = 0;
				unsigned int highest_index = 0;
				for (int n = 0; n < cf->n_neigh; ++n)
					if (lattice_th[neighbours[n]] > highest && (int)lattice_th[neighbours[n]] < (int)iter - (int)cf->age)
					{
						highest = lattice_th[neighbours[n]];
						highest_index = n;
					}

				/* If we have been to one of the neighbours before then go there (with
				 * a bit of noise), otherwise if none of the neighbours have been visited
				 * randomly select one */
				if (highest > 0)// && getRandNum() > cf->noise)
				{
					lattice_t[neighbours[highest_index]] = 1;
					ct++;
				}
				else
					ct += choose_site(neighbours, lattice_t, lattice_th, iter, cf);

				/* Update the time last visited in the lattice_th array */
				lattice_th[index] = iter;

				free(neighbours);
			}
		}

	if(ct < cf->target_activity)
	{
		unsigned int *signals = calloc(ct, sizeof(unsigned int));

		unsigned int ap = 0;
		for (int i = 0; i < cf->arr_dim; ++i)
			if(lattice_t[i])
				signals[ap++] = i;

		unsigned int *neighbours = malloc(cf->n_neigh * sizeof(unsigned int));
		do
		{
			unsigned int k = signals[(int)(getRandNum() * ap - 1)];
			get_neighbours(neighbours, k/cf->dim, k%cf->dim, cf->dim, cf->arr_dim);

			ct += choose_site(neighbours, lattice_t, lattice_th, iter, cf);
		}
		while(getRandNum() < 0.85);

		free(neighbours);
		free(signals);
	}

	/* Check that the signal is still alive (explained further above) */
	if(ct == 0)
	{
		printf("Signal died.\n");
		exit(EXIT_FAILURE);
	}
}

int choose_site(unsigned int *neighbours, unsigned int *lattice_t, unsigned int *lattice_th, unsigned int iter, config *cf)
{
	unsigned int rn;

	/* Choose a random neighbour that hasn't been visited too recently
	 * i.e., at least cf.age time steps ago. */
	for (int ii = 0; ii < 36; ++ii)
	{
		rn = (unsigned int)(getRandNum() * (float)cf->n_neigh); /* Six neighbours */

		/* Make sure it hasn't been visited too recently */
		if ((int)lattice_th[neighbours[rn]] < abs((int)iter - (int)cf->age))
			break;
	}

	lattice_t[neighbours[rn]] = 1;

	return 1;
}

void propagate_2(unsigned int *lattice, unsigned int *lattice_t, config *cf)
{
	/* Move the signal form the temporary lattice back to the primary lattice */
	for (int i = 0; i < cf->arr_dim; ++i)
		if(getRandNum() > cf->noise)
			lattice[i] = lattice_t[i];
		else
			lattice[i] = 0;
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

void write_last_visited(FILE *stream, unsigned int *arr, unsigned int dim)
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
	/* Simple wrapper for calculating the linear index of a 2D
	 * coordinate. */
	return (x * d) + y;
}

void get_neighbours(unsigned int *nb, int x, int y, int d, int td)
{
	/* Determine the indices of the elements surrounding the current
	 * element. I would like to write this more generally, something
	 * like x + i%1 but this is a quick and easy solution for now
	 * providing the lattice coordination doesn't change. */
	nb[0] = modulo(ind(x + 1, y,     d), td);
	nb[1] = modulo(ind(x + 1, y + 1, d), td);
	nb[2] = modulo(ind(x,     y + 1, d), td);
	nb[3] = modulo(ind(x - 1, y,     d), td);
	nb[4] = modulo(ind(x - 1, y - 1, d), td);
	nb[5] = modulo(ind(x,     y - 1, d), td);
}

int modulo(int i, int n)
{
	/* The `%' operator does not behave as a proper modulo operator
	 * for negative values so a real one is implemented here.*/
	if(i > -1) return i % n;

	/* Incorrect implementation but faster, works in this case as
	 * input values are constrained such that the second % will
	 * never be needed. */
	else return(i % n + n);

	/* Proper implementation of mod function */
//	else return(i % n + n) % n;
}

void initSeed(void)
{
    tinymt64_init(&tinymt_gen, (unsigned long)time(NULL));
}

double getRandNum(void)
{
	return tinymt64_generate_double(&tinymt_gen);
}

void wait_for_ms(clock_t wait_time)
{
	clock_t ret_time = clock() + wait_time;
	while (clock() < ret_time);
}
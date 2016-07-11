#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <time.h>

#include "tinymt64.h"

#include "config.h"

/* Don't really want to make this global but its very messy otherwise. 
 * I think in mt19937ar you had no choice anyway. */
tinymt64_t tinymt_gen;

/* AoS will hopefully give better cache performance than SoA */
typedef struct {
	int v; /* Value */
	int t; /* Temp */
	int l; /* Last */
	int c; /* Count total visits */
} latt_site;

int choose_site(unsigned int *neighbours, latt_site *lattice_t, unsigned int iter, config *cf);
int propagate_1(latt_site *lattice, unsigned int *histogram, unsigned int *branches, unsigned int iter, config *cf);
void propagate_2(latt_site *lattice, config *cf);
void get_neighbours(unsigned int *nb, int x, int y, int d, int td);

int is_in_arr(unsigned int n, unsigned int *arr, int l);
int modulo(int i, int n);
int ind(int x, int y, int d);

void initSeed(void);
double getRandNum(void);

void write_array(FILE *stream, latt_site *lattice, int offset, unsigned int dim);
void write_last_visited(FILE *stream, latt_site *lattice, int offset, unsigned int dim);

void *malloc_s(size_t size);
void *calloc_s(size_t num, size_t size);

void wait_for_ms(clock_t wait_time);

int main(void)
{
	config cf;
	read_config(&cf);

	initSeed();

	clock_t start = clock();

	latt_site *lattice = calloc(cf.arr_dim, sizeof(latt_site));

	/* Distribution time since last visit */
	unsigned int *histogram  = calloc_s(cf.age * 2, sizeof(unsigned int));

	/* Points at which branching is allowed */
	unsigned int *branches   = malloc_s(0.1f * cf.arr_dim * sizeof(unsigned int));

	/* Randomly choose branching points */
	for (int i = 0; i < 0.1f * cf.arr_dim; ++i)
		branches[i] =  getRandNum() * cf.arr_dim;

	/* Put a signal in the centre of the lattice */
	lattice[cf.arr_dim/2 - cf.dim/2].v = 1;

	printf("  0%%");

	for (int iter = 2; iter < cf.iter + 2; ++iter)
	{
		if(iter % (cf.iter / 100) == 0) printf("\r%3d%%", 100 * iter / cf.iter);

		if(cf.slow > 0)
			wait_for_ms(cf.slow);

		if(!propagate_1(lattice, histogram, branches, iter, &cf))
			/* Signal has died */
			break;

		propagate_2(lattice, &cf);

		if(cf.write_frames)
		{
			/* Write out the last visited time to unique files, one per
			 * time step. Careful! Can produce a lot of output if left
			 * running for too long! */
			char f_name[64];
			sprintf(f_name, "output/pcount_%d.dat", iter-2);
			FILE *fp = fopen(f_name, "w");

			write_last_visited(fp, lattice, offsetof(latt_site, l), cf.dim);
			fclose(fp);
		}
	}

	FILE *fp = fopen("output/total_visits.dat", "w");
	write_last_visited(fp, lattice, offsetof(latt_site, c), cf.dim);
	fclose(fp);

	printf("\n\nSimulation run time: %.2lfs\n", (clock() - start)/1000.0f);

	fp = fopen("output/visited_hist.dat", "w");
	for (int i = 0; i < cf.age*2; ++i)
		fprintf(fp, "%d\n", histogram[i]);
	fclose(fp);

	free(lattice);

	return EXIT_SUCCESS;
}

int propagate_1(latt_site *lattice, unsigned int *histogram, unsigned int *branches, unsigned int iter, config *cf)
{
	/* Signal is propagated into a tempoarary lattice, lattice_t. If this
	 * is not done then there can be confusion as the lattice is updated
	 * sequentially. Each signal should be propagated once per time step
	 * but without this method one signal can easily be updated many times
	 * per time step if it is propagated into a cell that has not yet been
	 * updated itself. */

	/* CLear the temporary lattice.t */
	for (int i = 0; i < cf->arr_dim; ++i)
		lattice[i].t = 0;

 	/* Used as a check to see if the signal is still alive, mainly
 	 * useful in debugging when introducing new propagation logic */
 	unsigned int ct = 0;

	for (int i = 0; i < cf->dim; ++i)
		for (int j = 0; j < cf->dim; ++j)
		{
			/* Shortcuts */
			unsigned int index = modulo(ind(i,j,cf->dim), cf->arr_dim);

			/* Focus on those lattice elements where there is a signal.
			 * Can probably get a big speed increase by keeping track of where the
			 * signal is and only touching that element. Oh well this is easier for now */
			if (lattice[index].v >= 1)
			{
				/* Build a vector containing indices of the surrounding elements */
				unsigned int *neighbours = malloc_s(cf->n_neigh * sizeof(unsigned int));
				get_neighbours(neighbours, i, j, cf->dim, cf->arr_dim);

				/* Search the surrounding elements to see which was visited the most recently,
				 * but not too recently (more than cf.age time steps ago) */
				unsigned int highest = 0;
				unsigned int highest_index = 0;

				for (int n = 0; n < cf->n_neigh; ++n)
					if (lattice[neighbours[n]].l > highest && lattice[neighbours[n]].l < (int)iter - cf->age)
					{
						highest = lattice[neighbours[n]].l;
						highest_index = n;
					}

				/* If we have been to one of the neighbours before then go there (with
				 * a bit of noise), otherwise if none of the neighbours have been visited
				 * randomly select one */
				if (highest > 0 && getRandNum() > cf->select_noise)
				{
					lattice[neighbours[highest_index]].t = 1;
					ct++;
				}
				else
					ct += choose_site(neighbours, lattice, iter, cf);

				/* Store time since last visit in histogram */
				if(lattice[index].l > 0 && iter - lattice[index].l < cf->age * 2)
					histogram[iter - lattice[index].l]++;

				/* If we are at a branch site and the last visited time is in the appropriate window, then branch */
				if(is_in_arr(i * cf->dim + j, branches, 0.1f * cf->arr_dim) > 0 && abs(iter - lattice[index].l - cf->age) < cf->branch_window)
					ct += choose_site(neighbours, lattice, iter, cf);

				/* Total activity level, can be disabled in config.json */
				if (ct < cf->total_act_target)
					ct += choose_site(neighbours, lattice, iter, cf);

				/* Update the time last visited */
				lattice[index].l = iter;

				free(neighbours);
			}
		}

	/* Check that the signal is still alive (explained further above) */
	if(ct == 0)
	{
		printf("\n\nSignal died.\n");
		return 0;
	}
	else
		return 1;
}

int is_in_arr(unsigned int n, unsigned int *arr, int l)
{
	/* Check whether n is in the array arr of length l */
	for (int i = 0; i < l; ++i)
		if(arr[i] == n)
			return i;

	return 0;
}

int choose_site(unsigned int *neighbours, latt_site *lattice, unsigned int iter, config *cf)
{
	unsigned int rn;

	/* Choose a random neighbour that hasn't been visited too recently
	 * i.e., at least cf.age time steps ago. */
	for (int ii = 0; ii < 36; ++ii)
	{
		rn = (unsigned int)(getRandNum() * (float)cf->n_neigh); /* Six neighbours */

		/* Make sure it hasn't been visited too recently */
		if (lattice[neighbours[rn]].l < abs((int)iter - cf->age))
			break;
	}

	lattice[neighbours[rn]].t = 1;

	return 1;
}

void propagate_2(latt_site *lattice, config *cf)
{
	/* Move the signal form the temporary lattice back to the primary lattice */
	for (int i = 0; i < cf->arr_dim; ++i)
		if(getRandNum() > cf->propagate_noise)
		{
			lattice[i].v = lattice[i].t;
			if(lattice[i].v)
				lattice[i].c++;
		}
		else
			lattice[i].v = 0;
}

void write_array(FILE *stream, latt_site *lattice, int offset, unsigned int dim)
{
	/* As the lattice is hexagonal a bit of extra work goes into making
	 * the cells line up as they should. The entire array is skewed to 
	 * give hexagonal alignment and wrapped to bring it back to a square
	 * on screen. */
	 
	system("cls");
	for (unsigned int i = 0; i < dim; ++i)
	{
		for (unsigned int j = 0; j < 2*((i+1) % 2); ++j)
			printf(" ");
		for (unsigned int j = 0; j < dim; ++j)
		{
			int index = i + dim * ((j + i/2) % dim);
			if (lattice[index].l == 0)
				fprintf(stream, "    ");
			else
			{
				char* latt_elem = (char*)&lattice[index];
				fprintf(stream, "%3d ", *(int *)(latt_elem + offset));
			}
		}

		fprintf(stream, "\n\n");
	}

	fprintf(stream, "\n");
}

void write_last_visited(FILE *stream, latt_site *lattice, int offset, unsigned int dim)
{
	for (unsigned int i = 0; i < dim; ++i)
	{
		for (unsigned int j = 0; j < dim; ++j)
		{
			char* latt_elem = (char*)&lattice[i * dim + j];
			fprintf(stream, "%d ", *(int *)(latt_elem + offset));
		}
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

void *malloc_s(size_t size)
{
	void *block = NULL;
	block = malloc(size);

	if(block == NULL)
	{
		fprintf(stderr, "Failed to allocate memory.\n");
		exit(EXIT_FAILURE);
	}

	return block;
}

void *calloc_s(size_t num, size_t size)
{
	void *block = NULL;
	block = calloc(num, size);

	if(block == NULL)
	{
		fprintf(stderr, "Failed to allocate memory.\n");
		exit(EXIT_FAILURE);
	}

	return block;
}

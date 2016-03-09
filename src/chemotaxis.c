#include <stdio.h>
#include <stdlib.h>



int ind(int x, int y, int z, int xd, int yd, int zd);

void get_neighbours(int *nb, int x, int y, int z, int d);

int main(void)
{
	/* Must be even */
	const int n_neigh = 6;
	if (n_neigh % 2 != 0)
	{
		printf("n_neigh must be even. Aborting.\n");
		return EXIT_FAILURE;
	}
	
	const int dim = 20;
	const int arr_dim = pow(dim, n_neigh/2);

	int lattice;

	for (int i = 0; i < dim; ++i)
	{
		for (int j = 0; j < dim; ++j)
		{
			for (int k = 0; k < dim; ++k)
			{
				int neighbours[6];
				get_neighbours(neighbours, i, j, k, dim);
				free neighbours;
			}
		}
	}

	return EXIT_SUCCESS;
}

int ind(int x, int y, int z, int d)
{
	return (x * d * d) + (y * d) + z;
}

void get_neighbours(int *nb, int x, int y, int z, int d)
{
	nb[0] = ind(x + 1, y, z, d);
	nb[1] = ind(x, y + 1, z, d);
	nb[2] = ind(x, y, z + 1, d);
	nb[3] = ind(x - 1, y, z, d);
	nb[4] = ind(x, y - 1, z, d);
	nb[5] = ind(x, y, z - 1, d);
}

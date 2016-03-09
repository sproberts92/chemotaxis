#include <iostream>
#include <array>

int ind(int x, int y, int z, int xd, int yd, int zd);

void get_neighbours(int *nb, int x, int y, int z, int d);

int main(void)
{
	const int dim = 20;
	const int arr_dim = dim * dim * dim;

	int lattice;


	for (int i = 0; i < dim; ++i)
	{
		for (int j = 0; j < dim; ++j)
		{
			for (int k = 0; k < dim; ++k)
			{
				int neighbours[6];
				get_neighbours(i, j, k, dim);
				free neighbours;
			}
		}
	}

	return 0;
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

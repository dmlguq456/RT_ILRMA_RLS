#include "header.h"



double *CreateVectorD(int size)
{
	double *v;
	int *i;

	v = (double *)malloc((size + 1) * sizeof(double));
	i = (int *)v;
	*i = size;
	return v;
}

int VectorSizeD(double *v)
{
	int *i;

	i = (int *)v;
	return *i;
}

void ZeroVectorD(double *v)
{
	int i, n;

	n = VectorSizeD(v);
	for (i = 1; i <= n; i++) v[i] = 0.0;
}

double **Create2dVectorD(int row, int col)
{
	double **v;
	int *i, j;

	v = (double **)malloc((row + 1) * sizeof(double *));
	i = (int *)v; *i = row;

	for (j = 1; j <= row; j++)
	{
		v[j] = (double *)malloc((col + 1) * sizeof(double));
		i = (int *)v[j]; *i = col;
	}

	return v;
}

void Vector2dSizeD(double **v, int *row, int *col)
{
	int *i;

	i = (int *)v;
	*row = *i;

	i = (int *)v[1];
	*col = *i;

}
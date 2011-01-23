/* Author: Heng Li; License: MIT/X11 */

/*
   A good SSE/SSE2/SSE3 tutorial:

http://developer.apple.com/hardwaredrivers/ve/sse.html

This program demonstrates single-precision doubleing point operations
using SSE instructions.
 */

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>

double **mm_init(int n)
{
	double **m;
	int i;
	m = (double**)malloc(n * sizeof(void*));
	for (i = 0; i < n; ++i)
		m[i] = calloc(n, sizeof(double));
	return m;
}
void mm_destroy(int n, double **m)
{
	int i;
	for (i = 0; i < n; ++i) free(m[i]);
	free(m);
}
double **mm_gen(int n)
{
	double **m;
	int i, j;
	m = mm_init(n);
	for (i = 0; i < n; ++i)
		for (j = 0; j < n; ++j)
			m[i][j] = 2 * drand48() - 1.0;
	return m;
}
// better cache performance by transposing the second matrix
double **mm_mul2(int n, double *const *a, double *const *b)
{
	int i, j, k;
	double **m, **c;
	m = mm_init(n); c = mm_init(n);
	for (i = 0; i < n; ++i) // transpose
		for (j = 0; j < n; ++j)
			c[i][j] = b[j][i];
	for (i = 0; i < n; ++i) {
		double *p = a[i], *q = m[i];
		for (j = 0; j < n; ++j) {
			double t = 0.0, *r = c[j];
			for (k = 0; k < n; ++k)
				t += p[k] * r[k];
			q[j] = t;
		}
	}
	mm_destroy(n, c);
	return m;
}
int main(int argc, char *argv[])
{
	int n = 100;
	double **a, **b, **m;
	clock_t t;
	if (argc > 1) n = atoi(argv[1]);
	n = (n + 3) / 4 * 4; // for simplicity, n can be divided by 4
	srand48(11);
	a = mm_gen(n); b = mm_gen(n);

	t = clock();
	m = mm_mul2(n, a, b);
	fprintf(stderr, "cache:  %lf sec; M[%d][%d]=%f\n", (double)(clock() - t) / CLOCKS_PER_SEC, n/2, n/2, m[n/2][n/2]);

	mm_destroy(n, a); mm_destroy(n, b); mm_destroy(n, m);
	return 0;
}


/* Author: Heng Li; License: MIT/X11 */

/*
   A good SSE/SSE2/SSE3 tutorial:

http://developer.apple.com/hardwaredrivers/ve/sse.html

This program demonstrates single-precision floating point operations
using SSE instructions.
 */

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <emmintrin.h>

float **mm_init(int n)
{
	float **m;
	int i;
	m = (float**)malloc(n * sizeof(void*));
	for (i = 0; i < n; ++i)
		m[i] = calloc(n, sizeof(float));
	return m;
}
void mm_destroy(int n, float **m)
{
	int i;
	for (i = 0; i < n; ++i) free(m[i]);
	free(m);
}
float **mm_gen(int n)
{
	float **m;
	int i, j;
	m = mm_init(n);
	for (i = 0; i < n; ++i)
		for (j = 0; j < n; ++j)
			m[i][j] = 2 * drand48() - 1.0;
	return m;
}
// the most naive version
float **mm_mul0(int n, float *const *a, float *const *b)
{
	int i, j, k;
	float **m;
	m = mm_init(n); // all entries are zero
	for (i = 0; i < n; ++i)
		for (j = 0; j < n; ++j)
			for (k = 0; k < n; ++k)
				m[i][j] += a[i][k] * b[k][j];
	return m;
}
// explicitly instruct the compiler where unnecessary memory acces can be saved.
float **mm_mul1(int n, float *const *a, float *const *b)
{
	int i, j, k;
	float **m;
	m = mm_init(n);
	for (i = 0; i < n; ++i) {
		float *p = a[i], *q = m[i];
		for (j = 0; j < n; ++j) {
			float t = 0.0;
			for (k = 0; k < n; ++k)
				t += p[k] * b[k][j];
			q[j] = t;
		}
	}
	return m;
}
// better cache performance by transposing the second matrix
float **mm_mul2(int n, float *const *a, float *const *b)
{
	int i, j, k;
	float **m, **c;
	m = mm_init(n); c = mm_init(n);
	for (i = 0; i < n; ++i) // transpose
		for (j = 0; j < n; ++j)
			c[i][j] = b[j][i];
	for (i = 0; i < n; ++i) {
		float *p = a[i], *q = m[i];
		for (j = 0; j < n; ++j) {
			float t = 0.0, *r = c[j];
			for (k = 0; k < n; ++k)
				t += p[k] * r[k];
			q[j] = t;
		}
	}
	mm_destroy(n, c);
	return m;
}
// explicit SSE optimization for the inner loop
float **mm_mul3(int n, float *const *a, float *const *b)
{
	int i, j, k;
	float **m, **c, x[4];
	m = mm_init(n); c = mm_init(n);
	for (i = 0; i < n; ++i) // transpose
		for (j = 0; j < n; ++j)
			c[i][j] = b[j][i];
	for (i = 0; i < n; ++i) {
		float *p = a[i], *q = m[i];
		for (j = 0; j < n; ++j) {
			__m128 t = _mm_setzero_ps();
			float *r = c[j];
			for (k = 0; k < n; k += 4) // four operations in one CPU cycle
				t = _mm_add_ps(t, _mm_mul_ps(_mm_load_ps(p+k), _mm_load_ps(r+k)));
			_mm_store_ps(x, t);
			q[j] = x[0] + x[1] + x[2] + x[3];
		}
	}
	mm_destroy(n, c);
	return m;
}

int main(int argc, char *argv[])
{
	int n = 100;
	float **a, **b, **m;
	clock_t t;
	if (argc > 1) n = atoi(argv[1]);
	n = (n + 3) / 4 * 4; // for simplicity, n can be divided by 4
	srand48(11);
	a = mm_gen(n); b = mm_gen(n);

	t = clock();
	m = mm_mul0(n, a, b);
	fprintf(stderr, "naive:  %lf sec; M[%d][%d]=%f\n", (double)(clock() - t) / CLOCKS_PER_SEC, n/2, n/2, m[n/2][n/2]);

	t = clock();
	m = mm_mul1(n, a, b);
	fprintf(stderr, "better: %lf sec; M[%d][%d]=%f\n", (double)(clock() - t) / CLOCKS_PER_SEC, n/2, n/2, m[n/2][n/2]);

	t = clock();
	m = mm_mul2(n, a, b);
	fprintf(stderr, "cache:  %lf sec; M[%d][%d]=%f\n", (double)(clock() - t) / CLOCKS_PER_SEC, n/2, n/2, m[n/2][n/2]);

	t = clock();
	m = mm_mul3(n, a, b);
	fprintf(stderr, "SSE:    %lf sec; M[%d][%d]=%f\n", (double)(clock() - t) / CLOCKS_PER_SEC, n/2, n/2, m[n/2][n/2]);

	mm_destroy(n, a); mm_destroy(n, b); mm_destroy(n, m);
	return 0;
}


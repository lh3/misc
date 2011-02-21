#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>

double n_unequal_jk(double theta0, int g, const int *m, const double *theta, double *var)
{
	double *h, est;
	int j;
	int64_t n = 0;
	h = calloc(g, sizeof(double));
	for (j = 0; j < g; ++j) n += m[j];
	for (j = 0; j < g; ++j) h[j] = (double)n / m[j];
	for (j = 0, est = 0.0; j < g; ++j) est += (1. - (double)m[j] / n) * theta[j];
	est = g * theta0 - est;
	*var = 0;
	for (j = 0; j < g; ++j) {
		double x = (h[j] * theta0 - (h[j] - 1) * theta[j] - est);
		*var += 1. / (h[j] - 1) * x * x;
	}
	*var /= g;
	free(h);
	return est;
}

double n_unequal_jk2(int g, int *x, int *y, double *e, double *var)
{
	int64_t xsum, ysum;
	int j;
	double theta0, *theta;
	xsum = ysum = 0;
	for (j = 0; j < g; ++j) xsum += x[j], ysum += y[j];
	theta0 = (double)ysum / xsum;
	theta = calloc(g, sizeof(double));
	for (j = 0; j < g; ++j)
		theta[j] = (double)(ysum - y[j]) / (xsum - x[j]);
	*e = n_unequal_jk(theta0, g, x, theta, var);
	free(theta);
	return theta0;
}

#ifdef UNEQUAL_JK2_MAIN
#include <math.h>
#include <stdint.h>
#include <time.h>
#include <unistd.h>
#include <string.h>

int main(int argc, char *argv[])
{
	FILE *fp;
	int j, g, x, y, *xa, *ya, max, b, B;
	double theta0, e, var, bmean, bvar;

	B = 1000;
	srand48(time(0) ^ (long)getpid());

	g = max = 0;
	xa = ya = 0;
	if (argc == 1) {
		fprintf(stderr, "Usage: unequal-jk <2col.txt>\n");
		return 1;
	}
	fp = strcmp(argv[1], "-") == 0? fdopen(fileno(stdin), "r") : fopen(argv[1], "r");
	while (!feof(fp) && fscanf(fp, "%d%d", &x, &y)) {
		if (x == 0) continue;
		if (g == max) {
			max = max? max<<1 : 16;
			xa = realloc(xa, sizeof(int) * max);
			ya = realloc(ya, sizeof(int) * max);
		}
		xa[g] = x; ya[g++] = y;
	}
	theta0 = n_unequal_jk2(g, xa, ya, &e, &var);
	for (b = 0, bmean = bvar = 0.0; b < B; ++b) {
		double btheta;
		int64_t bxsum = 0, bysum = 0;
		for (j = 0; j < g; ++j) {
			int bj = (int)(g * drand48());
			bxsum += xa[bj]; bysum += ya[bj];
		}
		btheta = (double)bysum / bxsum;
		bmean += btheta; bvar += btheta * btheta;
	}
	bmean /= B;
	bvar = bvar / B - bmean * bmean;
	fclose(fp); free(xa); free(ya);
	printf("%lf\t%lf\t%lf\t%lf\t%lf\n", theta0, e, sqrt(var), bmean, sqrt(bvar));
	return 0;
}
#endif

#ifdef UNEQUAL_JK3_MAIN
#include <math.h>
#include <stdint.h>
#include <string.h>
#include "ksort.h"

typedef struct {
	int x, y;
	double z;
} data_t;

#define __lt(a, b) ((a).z < (b).z)
KSORT_INIT(d, data_t, __lt)

int main(int argc, char *argv[])
{
	FILE *fp;
	int k, j, g, x, y, max, n_bins = 10, *xa, *ya, rest, shift;
	double theta0, e0, var0, z;
	data_t *a;

	g = max = 0; a = 0;
	if (argc == 1) {
		fprintf(stderr, "Usage: unequal-jk3 <3col.txt> [n_bins]\n");
		return 1;
	}
	if (argc >= 3) n_bins = atoi(argv[2]);
	fp = strcmp(argv[1], "-") == 0? fdopen(fileno(stdin), "r") : fopen(argv[1], "r");
	while (!feof(fp) && fscanf(fp, "%lf%d%d", &z, &x, &y)) {
		if (x == 0) continue;
		if (g == max) {
			max = max? max<<1 : 16;
			a = realloc(a, sizeof(data_t) * max);
		}
		a[g].z = z; a[g].x = x; a[g++].y = y;
	}
	xa = calloc(g, sizeof(int));
	ya = calloc(g, sizeof(int));
	for (j = 0; j < g; ++j) xa[j] = a[j].x, ya[j] = a[j].y;
	theta0 = n_unequal_jk2(g, xa, ya, &e0, &var0);
	ks_introsort(d, g, a);
	rest = g; shift = 0;
	for (k = 0; k < n_bins; ++k) {
		int g2 = (int)((double)rest / (n_bins - k) + .499);
		double zmean = 0., theta, e, var;
		for (j = 0; j < g2; ++j) {
			xa[j] = a[j+shift].x, ya[j] = a[j+shift].y;
			zmean += a[j+shift].z;
		}
		shift += g2; rest -= g2;
		zmean /= g2;
		theta = n_unequal_jk2(g2, xa, ya, &e, &var);
		printf("%lf\t%lf\t%lf\t%lf\n", zmean, theta / theta0, sqrt(var) / theta0, (double)(k+1) / n_bins);
	}
	fclose(fp); free(xa); free(ya); free(a);
	fprintf(stderr, "%lf\t%lf\t%lf\n", theta0, e0, sqrt(var0));
	return 0;
}
#endif

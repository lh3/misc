#include <zlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <assert.h>
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

static int g_len = 35;
static double g_ratio = 0.5;

static inline unsigned get_cnt(int c)
{
	int c1, c2;
	c1 = (c - 63) >> 3; c2 = (c - 63) & 7;
	c1 = c1? 1<<(c1-1) : 0; c2 = c2? 1<<(c2-1) : 0;
	return c1<<16 | c2;
}

static inline int is_good(int c1, int c2)
{
	return (c1 == 1 && c2 == 0);
}

int main(int argc, char *argv[])
{
	int l, i, c;
	long long cnt[5], tot;
	kseq_t *seq;
	gzFile fp;
	cnt[0] = cnt[1] = cnt[2] = cnt[3] = cnt[4] = 0;
	while ((c = getopt(argc, argv, "l:r:")) >= 0) {
		switch (c) {
		case 'l': g_len = atoi(optarg); break;
		case 'r': g_ratio = atof(optarg); break;
		}
	}
	if (argc == optind) {
		fprintf(stderr, "Usage: gen_mask [-l %d] [-r %.2lf] <in.rawMask.fa>\n", g_len, g_ratio);
		return 1;
	}
	fp = gzopen(argv[optind], "r");
	seq = kseq_init(fp);
	while ((l = kseq_read(seq)) >= 0) {
		int n_good = 0, n_all = 0, n_mid = 0;
		printf(">%s %d %.3lf", seq->name.s, g_len, g_ratio);
		for (i = 0; i < l + g_len - 1; ++i) {
			int c1, c2;
			unsigned x = i < l? get_cnt(seq->seq.s[i]) : 0;
			c1 = x>>16; c2 = x&0xffff;
			if (c1 == 1) ++cnt[4];
			if (c1) {
				++n_all;
				if (is_good(c1, c2)) ++n_good;
				if (c1 == 1) ++n_mid;
			}
			x = i >= g_len? get_cnt(seq->seq.s[i - g_len]) : 0;
			c1 = x>>16; c2 = x&0xffff;
			if (c1) {
				--n_all;
				if (is_good(c1, c2)) --n_good;
				if (c1 == 1) --n_mid;
			}
			assert(n_all <= g_len && n_good <= n_all);
			if (i % 60 == 0) putchar('\n');
			x = n_all == 0? 0 : (double)n_good/n_all >= g_ratio? 3 : (double)n_mid/n_all >= g_ratio? 2 : 1;
			putchar(x + '0');
			cnt[x]++;
		}
		putchar('\n');
	}
	tot = cnt[1] + cnt[2] + cnt[3];
	fprintf(stderr, "%lld, %lld, %lld, %lld, %lld\n", cnt[0], cnt[1], cnt[2], cnt[3], cnt[4]);
	fprintf(stderr, "%lf, %lf, %lf\n", (double)cnt[3] / tot, (double)(cnt[2] + cnt[3]) / tot, (double)cnt[4] / tot);
	kseq_destroy(seq);
	gzclose(fp);
	return 0;
}

#include <zlib.h>
#include <stdio.h>
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

static int g_len = 35;

int main(int argc, char *argv[])
{
	int l, i, j;
	kseq_t *seq;
	gzFile fp;
	if (argc == 1) {
		fprintf(stderr, "splitfa <in.fa> [len=%d]\n", g_len);
		return 1;
	}
	if (argc >= 3) g_len = atoi(argv[2]);
	fp = gzopen(argv[1], "r");
	seq = kseq_init(fp);
	while ((l = kseq_read(seq)) >= 0) {
		if (l < g_len) continue;
		for (i = 0; i <= l - g_len; ++i) {
			printf(">%s_%d\n", seq->name.s, i+1);
			for (j = 0; j < g_len; ++j)
				putchar(seq->seq.s[i + j]);
			putchar('\n');
		}
	}
	kseq_destroy(seq);
	gzclose(fp);
	return 0;
}

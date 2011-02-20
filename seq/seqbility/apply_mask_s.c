#include <zlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <assert.h>
#include <stdint.h>
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

typedef struct {
	int ori_len;
	uint32_t *mask;
} mask32_t;

#include "khash.h"
KHASH_MAP_INIT_STR(s, mask32_t)

static khash_t(s) *load_mask(const char *fn)
{
	kseq_t *seq;
	gzFile fp;
	khash_t(s) *h;
	h = kh_init(s);
	fp = gzopen(fn, "r");
	seq = kseq_init(fp);
	while (kseq_read(seq) >= 0) {
		khint_t k;
		int ret, i;
		mask32_t *p;
		k = kh_put(s, h, strdup(seq->name.s), &ret);
		assert(ret); // duplicated name
		p = &kh_val(h, k);
		p->ori_len = seq->seq.l;
		p->mask = (uint32_t*)calloc((seq->seq.l+31)/32, 4);
		for (i = 0; i < seq->seq.l; ++i)
			if (seq->seq.s[i] == '3')
				p->mask[i/32] |= 1u<<i%32;
	}
	kseq_destroy(seq);
	gzclose(fp);
	return h;
}

int main(int argc, char *argv[])
{
	gzFile fp;
	kseq_t *seq;
	khash_t(s) *hash;
	int c, min_q = 10;

	while ((c = getopt(argc, argv, "q:")) >= 0) {
		switch (c) {
		case 'q': min_q = atoi(optarg); break;
		}
	}
	if (argc <= optind + 1) {
		fprintf(stderr, "Usage: apply_mask_s [-q %d] <in.mask.fa> <in.fa>\n", min_q);
		return 1;
	}

	hash = load_mask(argv[optind]);
	fp = strcmp(argv[optind+1], "-")? gzopen(argv[optind+1], "r") : gzdopen(fileno(stdin), "r");
	seq = kseq_init(fp);
	while (kseq_read(seq) >= 0) {
		khint_t k;
		mask32_t *p;
		int i;
		k = kh_get(s, hash, seq->name.s);
		p = &kh_val(hash, k);
		if (k == kh_end(hash)) {
			fprintf(stderr, "[apply_mask_s] skip sequence %s\n", seq->name.s);
			continue;
		}
		for (i = 0; i < seq->seq.l; ++i) {
			if (isupper(seq->seq.s[i])) {
				int do_lower = 0;
				if (seq->seq.s[i] == 'N') do_lower = 1;
				else if (seq->qual.l == seq->seq.l && seq->qual.s[i] - 33 < min_q) do_lower = 1;
				else if (i >= p->ori_len || !(p->mask[i/32]&1u<<i%32)) do_lower = 1;
				if (do_lower) seq->seq.s[i] = tolower(seq->seq.s[i]);
			}
		}
		printf(">%s", seq->name.s);
		for (i = 0; i < seq->seq.l; ++i) {
			if (i%60 == 0) putchar('\n');
			putchar(seq->seq.s[i]);
		}
		putchar('\n');
	}
	kseq_destroy(seq);
	gzclose(fp);
	// hash table is not freed...
	return 0;
}

#include <zlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
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
	kstream_t *ks;
	khash_t(s) *hash;
	mask32_t *q = 0;
	kstring_t *str;
	int i, dret, c, last = 0;

	while ((c = getopt(argc, argv, "")) >= 0) {
	}
	if (argc <= optind + 1) {
		fprintf(stderr, "Usage: uniq-dist <in.mask.fa> <in-sorted.list>\n");
		return 1;
	}

	str = (kstring_t*)calloc(1, sizeof(kstring_t));
	fprintf(stderr, "[uniq-dist] loading mask...\n");
	hash = load_mask(argv[optind]);
	fp = gzopen(argv[optind+1], "r");
	ks = ks_init(fp);
	fprintf(stderr, "[uniq-dist] calculating unique distance...\n");
	while (ks_getuntil(ks, 0, str, &dret) >= 0) {
		khint_t k;
		mask32_t *p;
		int pos;
		k = kh_get(s, hash, str->s);
		p = (k != kh_end(hash))? &kh_val(hash, k) : 0;
		ks_getuntil(ks, 0, str, &dret);
		pos = atoi(str->s) - 1;
		if (p && pos >= 0 && pos < p->ori_len) {
			if (p != q) q = p; // change of reference
			else {
				if (last >= pos) {
					fprintf(stderr, "[uniq-dist] out of order: %s:%d <= %d\n", kh_key(hash, k), pos+1, last+1);
				} else {
					for (i = last, c = 0; i < pos; ++i)
						if (p->mask[i/32] & 1u<<i%32) ++c;
					if (last > 0) printf("%s\t%d\t%d\t%d\n", kh_key(hash, k), last, pos, c);
				}
			}
			last = pos;
		}
		if (dret != '\n')
			while ((c = ks_getc(ks)) != -1 && c != '\n');
	}
	ks_destroy(ks);
	gzclose(fp);
	free(str->s); free(str);
	// hash table is not freed...
	return 0;
}

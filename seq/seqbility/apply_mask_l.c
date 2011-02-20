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
	kstring_t *str;
	int dret, c;

	while ((c = getopt(argc, argv, "")) >= 0) {
	}
	if (argc <= optind + 1) {
		fprintf(stderr, "Usage: apply_mask_l <in.mask.fa> <in.list>\n");
		return 1;
	}

	str = (kstring_t*)calloc(1, sizeof(kstring_t));
	fprintf(stderr, "[apply_mask_l] loading mask...\n");
	hash = load_mask(argv[optind]);
	fp = gzopen(argv[optind+1], "r");
	ks = ks_init(fp);
	fprintf(stderr, "[apply_mask_l] filtering list...\n");
	while (ks_getuntil(ks, 0, str, &dret) >= 0) {
		khint_t k;
		mask32_t *p;
		int pos, do_print = 0;
		k = kh_get(s, hash, str->s);
		p = (k != kh_end(hash))? &kh_val(hash, k) : 0;
		ks_getuntil(ks, 0, str, &dret);
		pos = atoi(str->s) - 1;
		if (p && pos < p->ori_len && (p->mask[pos/32]&1u<<pos%32)) do_print = 1;
		if (do_print) printf("%s\t%d", kh_key(hash, k), pos + 1);
		if (dret != '\n') {
			if (do_print) putchar('\t');
			while ((c = ks_getc(ks)) != -1 && c != '\n')
				if (do_print) putchar(c);
		}
		if (do_print) putchar('\n');
	}
	ks_destroy(ks);
	gzclose(fp);
	free(str->s); free(str);
	// hash table is not freed...
	return 0;
}

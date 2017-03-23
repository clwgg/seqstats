#include <zlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

// Heng Li klib
#include "klib/kseq.h"
#include "klib/ksort.h"

// initialize kseq
KSEQ_INIT(gzFile, gzread)

// initialize sorting and define types and sorting scheme
#define pair_lt(a, b) ((a) < (b))
KSORT_INIT(pair, long, pair_lt)
KSORT_INIT_GENERIC(long)


int main(int argc, char *argv[])
{
  gzFile fp;
  kseq_t *seq;
  int l;
  int n = 0;
  double full = 0;
  if (argc == 1) {
    fprintf(stderr, "Usage: %s <in.fasta|in.fastq>\n", argv[0]);
    return 1;
  }
  if (strcmp( argv[1], "-") == 0) {
    argv[1] = "/dev/stdin";
  }

  fp = gzopen(argv[1], "r");

  if (!fp) {
    fprintf(stderr, "Can't open input file.\n\n");
    fprintf(stderr, "Usage: %s <in.fasta|in.fastq>\n", argv[0]);
    return 1;
  }

  seq = kseq_init(fp);

  int s = 10000;

  long *a;
  a = malloc(s * sizeof *a);
  if (!a) goto memerr;

  while ((l = kseq_read(seq)) >= 0) {

    if (l == 0) continue;

    if (n == s) {
      s = s * 2;
      a = realloc(a, s * sizeof *a);
      if (!a) goto memerr;
    }

    a[n] = l;

    full = full + l;
    n++;
  }
  if (n == 0) goto seqerr;

  ks_mergesort(pair, n, a, 0);

  double adding = 0;
  int j = -1;

  while (adding < full/2) {

    j++;
    adding = adding + a[j];
    
  }

  printf("Total n:\t%i\n", n);
  printf("Total seq:\t%.0f bp\n", full);
  printf("Avg. seq:\t%.2f bp\n", full/n);
  printf("Median seq:\t%.2f bp\n", (a[(int)ceil((n-1)/2.0)]+a[(int)floor((n-1)/2.0)])/2.0);
  printf("N 50:\t\t%ld bp\n", a[j]);
  printf("Min seq:\t%ld bp\n", a[0]);
  printf("Max seq:\t%ld bp\n", a[n-1]);

  kseq_destroy(seq);
  free(a);
  gzclose(fp);
  return 0;

  memerr:
    printf("\n\tCould not allocate sufficient continuous memory.\n\n");
    kseq_destroy(seq);
    gzclose(fp);
    return 1;

  seqerr:
    printf("\n\tSequence file is empty.\n\n");
    kseq_destroy(seq);
    gzclose(fp);
    return 1;

}


#include <zlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <stdint.h>

// Heng Li klib
#include "klib/kseq.h"
#include "klib/ksort.h"
#include "klib/khash.h"
#include "klib/kvec.h"

// initialize kseq
KSEQ_INIT(gzFile, gzread)

// initialize sorting and define types and sorting scheme
KSORT_INIT_GENERIC(long)

// map: sequence length -> count
KHASH_MAP_INIT_INT64(len_count, uint64_t)


int main(int argc, char *argv[])
{
  gzFile fp;
  kseq_t *seq;
  int l;
  uint64_t n = 0;
  double full = 0;
  if (argc == 1) {
    printf("Usage: %s <in.fasta|in.fastq>\n", argv[0]);
    return 1;
  }
  else if (argc > 2) {
    printf("\n\tPlease use\n\n\tcat in1.fa/q[.gz] in2.fa/q[.gz] | %s -", argv[0]);
    printf("\n\n\t          to specify more than one input.\n\n");
    return 1;
  }

  if (strcmp( argv[1], "-") == 0) {
    argv[1] = "/dev/stdin";
  }

  fp = gzopen(argv[1], "r");

  if (!fp) {
    printf("Can't open input file.\n\n");
    printf("Usage: %s <in.fasta|in.fastq>\n", argv[0]);
    return 1;
  }

  seq = kseq_init(fp);

  khash_t(len_count) *counts = kh_init(len_count);
  if (counts == NULL) {
    printf("\n\tCould not allocate memory.\n\n");
    kseq_destroy(seq);
    gzclose(fp);
    return 1;
  }

  long min_len = LONG_MAX;
  long max_len = LONG_MIN;

  while ((l = kseq_read(seq)) >= 0) {

    if (l == 0) continue;

    if ((long)l < min_len) min_len = (long)l;
    if ((long)l > max_len) max_len = (long)l;

    // increment count for length l
    khint_t k;
    int ret;
    k = kh_put(len_count, counts, (int64_t)l, &ret);
    if (ret == 0) {
      kh_value(counts, k) += 1;
    } else {
      kh_value(counts, k) = 1;
    }

    full = full + l;
    n++;
  }
  if (n == 0) goto seqerr;

  // collect unique lengths and sort
  kvec_t(long) uniq; kv_init(uniq);
  for (khint_t i_k = kh_begin(counts); i_k != kh_end(counts); ++i_k) {
    if (!kh_exist(counts, i_k)) continue;
    kv_push(long, uniq, (long)kh_key(counts, i_k));
  }
  if (kv_size(uniq) == 0) goto seqerr;
  ks_introsort_long(kv_size(uniq), uniq.a);

  double adding = 0;
  long n50_val = 0;

  // compute N50: traverse from largest length downwards
  for (long idx = (long)kv_size(uniq) - 1; idx >= 0; --idx) {
    long len = uniq.a[idx];
    khint_t k = kh_get(len_count, counts, (int64_t)len);
    if (k == kh_end(counts)) continue; // shouldn't happen
    uint64_t c = kh_value(counts, k);
    adding += (double)len * (double)c;
    if (adding >= full / 2.0) { n50_val = len; break; }
  }

  // compute median from counts
  uint64_t mid1 = (n - 1) / 2;  // 0-based
  uint64_t mid2 = n / 2;        // 0-based
  uint64_t cum = 0;
  long med1_val = 0, med2_val = 0;
  for (size_t i = 0; i < kv_size(uniq); ++i) {
    long len = uniq.a[i];
    khint_t k = kh_get(len_count, counts, (int64_t)len);
    if (k == kh_end(counts)) continue;
    uint64_t c = kh_value(counts, k);
    uint64_t prev = cum;
    cum += c;
    if (prev <= mid1 && mid1 < cum) med1_val = len;
    if (prev <= mid2 && mid2 < cum) { med2_val = len; break; }
  }

  printf("Total n:\t%llu\n", (unsigned long long)n);
  printf("Total seq:\t%.0f bp\n", full);
  printf("Avg. seq:\t%.2f bp\n", full/n);
  printf("Median seq:\t%.2f bp\n", (med1_val + med2_val) / 2.0);
  printf("N 50:\t\t%ld bp\n", n50_val);
  printf("Min seq:\t%ld bp\n", min_len);
  printf("Max seq:\t%ld bp\n", max_len);

  kseq_destroy(seq);
  kv_destroy(uniq);
  kh_destroy(len_count, counts);
  gzclose(fp);
  return 0;

  seqerr:
    printf("\n\tSequence file is empty.\n\n");
    kseq_destroy(seq);
    if (counts) kh_destroy(len_count, counts);
    gzclose(fp);
    return 1;

}
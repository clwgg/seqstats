seqstats: seqstats.c klib/kseq.h klib/ksort.h
	gcc $< -Wall -O3 -lz -o $@

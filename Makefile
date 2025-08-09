CC = gcc
CFLAGS ?= -O3 -Wall
LDFLAGS ?= -lz

SRC := seqstats.c
DEPS := klib/kseq.h klib/ksort.h klib/khash.h klib/kvec.h

seqstats: $(SRC) $(DEPS)
	$(CC) $(CFLAGS) $(SRC) $(LDFLAGS) -o $@

.PHONY: test
test: seqstats
	./tests/run_tests.sh

.PHONY: clean
clean:
	rm -f seqstats

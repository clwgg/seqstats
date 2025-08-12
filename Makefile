CC = gcc
CFLAGS ?= -O3 -Wall
LDFLAGS ?= -lz

SRC := seqstats.c
DEPS := klib/kseq.h klib/ksort.h klib/khash.h klib/kvec.h klib/ketopt.h

seqstats: $(SRC) $(DEPS)
	$(CC) $(CFLAGS) $(SRC) $(LDFLAGS) -o $@

.PHONY: debug
debug: $(SRC) $(DEPS)
	$(CC) -g -O0 -Wall $(SRC) $(LDFLAGS) -o seqstats

.PHONY: test
test: seqstats
	./tests/run_tests.sh

.PHONY: clean
clean:
	rm -f seqstats

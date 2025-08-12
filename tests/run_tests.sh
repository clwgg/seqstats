#!/usr/bin/env bash
set -euo pipefail

ROOT=$(cd "$(dirname "$0")/.." && pwd)
BIN="$ROOT/seqstats"
DATA="$ROOT/tests/data"
EXP="$ROOT/tests/expected"

if [ ! -x "$BIN" ]; then
  echo "Building seqstats..." >&2
  make -C "$ROOT"
fi

fail=0

run_case() {
  local input="$1" expected="$2"
  local out
  out=$("$BIN" "$input") || return 1
  # Normalize both outputs: strip trailing spaces and trailing blank lines
  if diff -u \
    <(printf "%s\n" "$out" | sed -e 's/[[:space:]]\+$//' -e :a -e '/^$/{$d;N;ba' -e '}') \
    <(sed -e 's/[[:space:]]\+$//' -e :a -e '/^$/{$d;N;ba' -e '}' "$expected"); then
    echo "PASS $input"
  else
    echo "FAIL $input" >&2
    fail=1
  fi
}

run_error_case() {
  local input="$1" expected_error="$2"
  local out err
  if out=$("$BIN" "$input" 2>&1); then
    echo "FAIL $input (expected error but got success)" >&2
    fail=1
    return
  fi
  if echo "$out" | grep -q "$expected_error"; then
    echo "PASS $input (error case)"
  else
    echo "FAIL $input (unexpected error output)" >&2
    echo "Expected: $expected_error" >&2
    echo "Got: $out" >&2
    fail=1
  fi
}

# Basic functionality tests
run_case "$DATA/test1.fa" "$EXP/test1.out"
run_case "$DATA/test2.fa" "$EXP/test2.out"
run_case "$DATA/test3.fq" "$EXP/test3.out"

# Edge case tests
run_case "$DATA/single.fa" "$EXP/single.out"
run_case "$DATA/identical.fa" "$EXP/identical.out"
run_case "$DATA/mixed.fa" "$EXP/mixed.out"
run_case "$DATA/multiline.fa" "$EXP/multiline.out"

# Error case tests
run_error_case "$DATA/empty.fa" "Sequence file is empty"
run_error_case "$DATA/zero_length.fa" "Sequence file is empty"

# Malformed input tests
run_case "$DATA/malformed2.fa" "$EXP/malformed2.out"
run_error_case "$DATA/malformed3.fq" "Incomplete FASTQ quartet"
run_error_case "$DATA/malformed4.fq" "FASTQ quality string length"
run_case "$DATA/malformed5.fa" "$EXP/malformed5.out"

# gzipped case generated on the fly
tmp_gz="$DATA/tmp_gz.fa.gz"
printf ">x\nACGT\n>y\nA\n" | gzip -c > "$tmp_gz"
exp_tmp=$(mktemp)
cat > "$exp_tmp" <<EOF
Total n:	2
Total seq:	5 bp
Avg. seq:	2.50 bp
Median seq:	2.50 bp
N 50:		4 bp
Min seq:	1 bp
Max seq:	4 bp
EOF
run_case "$tmp_gz" "$exp_tmp"
rm -f "$tmp_gz" "$exp_tmp"

if [ "$fail" -eq 0 ]; then
  echo "All tests passed."
else
  echo "Some tests failed." >&2
  exit 1
fi


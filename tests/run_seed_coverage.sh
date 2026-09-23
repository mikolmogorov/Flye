#!/usr/bin/env bash
set -euo pipefail
repo=$(cd "$(dirname "$0")/.." && pwd)
test_build=$(mktemp -d "${TMPDIR:-/tmp}/flye-seed-coverage.XXXXXXXX")
trap 'rm -f "$test_build/seed_coverage"; rmdir "$test_build"' EXIT
"${CXX:-g++}" -O0 -std=c++11 -pthread \
    -I"$repo/src" -I"$repo/lib/libcuckoo" \
    -I"$repo/lib/interval_tree" -I"$repo/lib/minimap2" \
    "$repo/tests/seed_coverage.cpp" "$repo/src/assemble/chimera.o" \
    "$repo"/src/sequence/*.o \
    -L"$repo/lib/minimap2" -lminimap2 -lz -o "$test_build/seed_coverage"
"$test_build/seed_coverage" "$repo/flye/config/bin_cfg/asm_raw_reads.cfg" partial-first
"$test_build/seed_coverage" "$repo/flye/config/bin_cfg/asm_raw_reads.cfg" full-first

#!/usr/bin/env bash
set -euo pipefail
repo=$(cd "$(dirname "$0")/.." && pwd)
test_build=$(mktemp -d "${TMPDIR:-/tmp}/flye-tip-pruning.XXXXXXXX")
trap 'rm -f "$test_build/tip_pruning"; rmdir "$test_build"' EXIT

# Build the production objects first (make THREADS=1); do not link main.o.
"${CXX:-g++}" -O0 -std=c++11 -pthread \
    -I"$repo/src" -I"$repo/lib/libcuckoo" \
    -I"$repo/lib/interval_tree" -I"$repo/lib/minimap2" \
    "$repo/tests/tip_pruning.cpp" \
    "$repo"/src/repeat_graph/*.o "$repo"/src/sequence/*.o \
    -L"$repo/lib/minimap2" -lminimap2 -lz -o "$test_build/tip_pruning"
"$test_build/tip_pruning" "$repo/flye/config/bin_cfg/asm_raw_reads.cfg"

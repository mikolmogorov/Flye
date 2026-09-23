#!/usr/bin/env bash
set -euo pipefail
repo=$(cd "$(dirname "$0")/.." && pwd)
test_build=$(mktemp -d "${TMPDIR:-/tmp}/flye-contained.XXXXXXXX")
trap 'rm -f "$test_build/contained"; rmdir "$test_build"' EXIT
"${CXX:-g++}" -O0 -std=c++11 -pthread \
    -I"$repo/src" -I"$repo/lib/libcuckoo" \
    -I"$repo/lib/interval_tree" -I"$repo/lib/minimap2" \
    "$repo/tests/contained_disjointigs.cpp" "$repo"/src/assemble/*.o \
    "$repo"/src/sequence/*.o \
    -L"$repo/lib/minimap2" -lminimap2 -lz -o "$test_build/contained"
"$test_build/contained" "$repo/flye/config/bin_cfg/asm_raw_reads.cfg"

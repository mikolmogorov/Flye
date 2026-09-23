#!/usr/bin/env bash
set -euo pipefail
repo=$(cd "$(dirname "$0")/.." && pwd)
test_build=$(mktemp -d "${TMPDIR:-/tmp}/flye-chimera-strands.XXXXXXXX")
trap 'rm -f "$test_build/chimera_strands"; rmdir "$test_build"' EXIT
"${CXX:-g++}" -O0 -std=c++11 -pthread \
    -I"$repo/src" -I"$repo/lib/libcuckoo" \
    -I"$repo/lib/interval_tree" -I"$repo/lib/minimap2" \
    "$repo/tests/chimera_strands.cpp" "$repo/src/assemble/chimera.o" \
    "$repo"/src/sequence/*.o \
    -L"$repo/lib/minimap2" -lminimap2 -lz -o "$test_build/chimera_strands"
"$test_build/chimera_strands" "$repo/flye/config/bin_cfg/asm_raw_reads.cfg"

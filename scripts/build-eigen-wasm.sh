#!/bin/sh
set -eu

ROOT="$(CDPATH= cd -- "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

EM_LLVM_ROOT="${EM_LLVM_ROOT:-/opt/local/libexec/llvm-22/bin}"
OUTPUT="src/eigen_solver.generated.mjs"

mkdir -p src

EM_LLVM_ROOT="$EM_LLVM_ROOT" emcc wasm/eigen_solver.cpp \
  -O3 \
  -std=c++17 \
  -I/opt/local/include/eigen3 \
  -sMODULARIZE=1 \
  -sEXPORT_ES6=1 \
  -sSINGLE_FILE=1 \
  -sWASM_ASYNC_COMPILATION=0 \
  -sALLOW_MEMORY_GROWTH=1 \
  -sENVIRONMENT=web,node \
  -sEXPORTED_FUNCTIONS=_solve_hermitian_eigen,_malloc,_free \
  -sEXPORTED_RUNTIME_METHODS=cwrap,getValue,setValue \
  -o "$OUTPUT"

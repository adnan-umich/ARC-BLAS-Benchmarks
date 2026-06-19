#!/usr/bin/env bash

set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

die() {
    printf 'error: %s\n' "$*" >&2
    exit 1
}

log() {
    printf '\n==> %s\n' "$*"
}

require_cmd() {
    command -v "$1" >/dev/null 2>&1 || die "$2"
}

require_loaded_paths() {
    local name="$1"
    local include_dir="$2"
    local lib_dir="$3"

    [[ -d "$include_dir" ]] || die "$name include directory not found: $include_dir"
    [[ -d "$lib_dir" ]] || die "$name library directory not found: $lib_dir"
}

require_hpro_source_tree() {
    local hpro_dir="$ROOT_DIR/hlr/hpro"
    local missing=()
    local required=(
        "configure"
        "SConstruct.in"
        "src/SConscript"
        "include/hpro/config.h.in"
        "bin/hpro-config.in"
        "include/hpro/base/types.hh"
        "include/hpro/base/config.hh"
        "include/hpro/blas/Algebra.hh"
        "include/hpro/matrix/TMatrix.hh"
    )

    for path in "${required[@]}"; do
        [[ -e "$hpro_dir/$path" ]] || missing+=("$path")
    done

    if ((${#missing[@]} > 0)); then
        printf 'error: bundled HPro/HLIBcore source tree is incomplete.\n' >&2
        printf 'Missing files under %s:\n' "$hpro_dir" >&2
        printf '  %s\n' "${missing[@]}" >&2
        printf '\nRestore/populate the full HPro/HLIBcore source tree before running this script.\n' >&2
        exit 1
    fi
}

prompt_prefix() {
    local prefix="${BUILD_PREFIX:-}"

    while [[ -z "$prefix" ]]; do
        read -e -r -p "Install prefix for BLAS-Benchmarks and HLR: " prefix
        if [[ "$prefix" == *[$'\001'-$'\037'$'\177']* ]]; then
            printf 'error: install prefix contains control characters; please enter a normal path.\n' >&2
            prefix=""
        fi
    done

    prefix="${prefix/#\~/$HOME}"
    mkdir -p "$prefix"
    cd "$prefix"
    pwd -P
}

find_boost_root() {
    if [[ -n "${BOOST_ROOT:-}" ]]; then
        printf '%s\n' "$BOOST_ROOT"
    elif [[ -n "${BOOST_DIR:-}" ]]; then
        printf '%s\n' "$BOOST_DIR"
    else
        die "Boost is not loaded. Set BOOST_ROOT or BOOST_DIR before running this script."
    fi
}

copy_hlr_outputs() {
    local prefix="$1"
    local bindir="$prefix/hlr/bin"
    local libdir="$prefix/hlr/lib"
    local includedir="$prefix/hlr/include"
    local sharedir="$prefix/hlr/share"

    mkdir -p "$bindir" "$libdir" "$includedir" "$sharedir"

    [[ -f "$ROOT_DIR/hlr/libhlr.a" ]] || die "HLR library was not produced at hlr/libhlr.a"
    cp "$ROOT_DIR/hlr/libhlr.a" "$libdir/"
    cp -R "$ROOT_DIR/hlr/include/hlr" "$includedir/"
    cp "$ROOT_DIR/hlr/"*.conf "$sharedir/" 2>/dev/null || true
    cp -R "$ROOT_DIR/hlr/scripts" "$sharedir/"

    while IFS= read -r exe; do
        cp "$exe" "$bindir/"
    done < <(
        find "$ROOT_DIR/hlr/programs" -type f -perm -111 \
            ! -name '*.cc' ! -name '*.hh' ! -name '*.o' \
            -print
    )

    if ! find "$bindir" -type f -perm -111 | grep -q .; then
        die "No HLR benchmark executables were found to install."
    fi
}

log "Checking required loaded software"
require_cmd gcc "GCC is not available. Load a GCC compiler module first."
require_cmd g++ "G++ is not available. Load a GCC compiler module first."
require_cmd gfortran "gfortran is not available. Load a GCC compiler module that includes Fortran support."
require_cmd pkg-config "pkg-config is required to locate FlexiBLAS."
require_cmd python "Python is not available. Load a Python module first."
require_cmd scons "SCons is not available. Load a Python module that includes the scons command."

require_hpro_source_tree

python - <<'PY' || die "Python cannot import SCons. Load a Python module that includes scons, or install scons into it."
import SCons
PY

pkg-config --exists flexiblas || die "FlexiBLAS is not available through pkg-config. Load FlexiBLAS and ensure PKG_CONFIG_PATH is set."

BOOST_ROOT="$(find_boost_root)"

require_loaded_paths "Boost" "$BOOST_ROOT/include" "$BOOST_ROOT/lib"

BUILD_PREFIX="$(prompt_prefix)"
VENV_DIR="$BUILD_PREFIX/venv"
JOBS="${JOBS:-$(getconf _NPROCESSORS_ONLN 2>/dev/null || printf '1')}"

if [[ -d "$VENV_DIR" ]]; then
    [[ -f "$VENV_DIR/bin/activate" ]] || die "Existing venv is missing bin/activate: $VENV_DIR"
    log "Reusing existing BLAS-Benchmarks Python environment at $VENV_DIR"
else
    log "Creating BLAS-Benchmarks Python environment at $VENV_DIR"
    python -m venv "$VENV_DIR"
fi
# shellcheck disable=SC1091
source "$VENV_DIR/bin/activate"
python -m pip install --upgrade pip
python -m pip install "$ROOT_DIR"

log "Configuring and building bundled HPro/HLIBcore"
cd "$ROOT_DIR/hlr/hpro"
./configure \
    --prefix="$ROOT_DIR/hlr/hpro" \
    --cc=gcc \
    --cxx=g++ \
    --fc=gfortran \
    --with-boost="$BOOST_ROOT" \
    --without-mkl \
    --without-zlib \
    --without-metis \
    --without-scotch \
    --without-mongoose \
    --without-gsl \
    --without-netcdf \
    --without-cgal \
    --without-amdlibm \
    --without-acml \
    --without-cuda \
    --without-hip
scons -j "$JOBS"

log "Building HLR benchmarks with benchmark.sh"
cd "$ROOT_DIR/hlr"
./scripts/benchmark.sh

log "Installing HLR artifacts into $BUILD_PREFIX/hlr"
copy_hlr_outputs "$BUILD_PREFIX"

cat <<EOF

Build complete.

BLAS-Benchmarks Python environment:
  $VENV_DIR

HLR install tree:
  $BUILD_PREFIX/hlr

To use the Python package:
  source "$VENV_DIR/bin/activate"

HLR benchmark executables are in:
  $BUILD_PREFIX/hlr/bin
EOF

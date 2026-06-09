#!/bin/bash
#
# CPU benchmarks using various HLR programs (mm, lu, etc.)
#

######################################################################
##
## common settings
##
######################################################################

# set of FlexiBLAS backends to benchmark (e.g. "OpenBLAS MKL BLIS")
FLEXIBLAS_BACKENDS="default"

# set of frameworks to use
FRAMEWORKS="seq tbb"

# set of programs
PROGRAMS="approx-mm approx-lu uniform-mm uniform-lu"

# set of applications
APPLICATIONS="laplace materncov"

# if 1, just build benchmark programs
JUST_BUILD=0

# if 1, only simulate benchmarking
DRY_RUN=0

# scons command for HLR
SCONS='scons -Q -D -j 8'

# directories
TOP_DIR=$( pwd )
CONTRIB_DIR=$TOP_DIR/contrib
TBB_DIR=/opt/local/tbb/2021.8

# determine CPU platform
ARCH=$( gcc -march=native -Q --help=target | grep march | head -n 1 | sed "s/ *-march=[ \t]*//" )

######################################################################
##
## parse command line
##
######################################################################

POSITIONAL_ARGS=()

while [[ $# -gt 0 ]]; do
    case $1 in
        --apps)
            APPLICATIONS="$2"
            shift # past argument
            shift # past value
            ;;
        --blas)
            FLEXIBLAS_BACKENDS="$2"
            shift # past argument
            shift # past value
            ;;
        --build)
            JUST_BUILD=1
            shift # past argument
            ;;
        --dry-run)
            DRY_RUN=1
            shift # past argument
            ;;
        --fworks)
            FRAMEWORKS="$2"
            shift # past argument
            shift # past value
            ;;
        --progs)
            PROGRAMS="$2"
            shift # past argument
            shift # past value
            ;;
        -h|--help)
            echo "usage: benchmark [options]"
            echo "with"
            echo "    --apps \"...\"   : set of applications (laplace,materncov)"
            echo "    --blas \"...\"   : set of FlexiBLAS backends to benchmark (e.g. \"OpenBLAS MKL BLIS\")"
            echo "    --build        : only build benchmark programs"
            echo "    --dry-run      : do not actually do anything"
            echo "    --fworks \"...\" : set of frameworks to use (seq,tbb)"
            echo "    --progs \"...\"  : set of programs to use (approx-mm,approx-lu,uniform-mm,uniform-lu)"
            echo "    -h/--help      : print this usage info"
            echo
            echo "example:"
            echo "    benchmark --blas \"mkl blis\" --frameworks seq --programs \"approx-lu uniform-lu\""
            exit 0
            shift # past argument
            ;;
        -*|--*)
            echo "Unknown option $1"
            exit 1
            ;;
        *)
            POSITIONAL_ARGS+=("$1") # save positional arg
            shift # past argument
            ;;
    esac
done

set -- "${POSITIONAL_ARGS[@]}" # restore positional parameters

#
# just make sure, the contrib directory exists
#
if [ $DRY_RUN == 0 ]; then
    mkdir -p $CONTRIB_DIR

    if [ "$?" != "0" ]; then
        echo "failed creating directory $CONTRIB_DIR"
        exit 1
    fi
fi

######################################################################
##
## install HPROcore alternatives
##
######################################################################

HPRO_DIR=/clusterhome/rok/programming/hpro/core

## TODO

######################################################################
##
## detect FlexiBLAS
##
######################################################################

FLEXIBLAS_LFLAGS=""
FLEXIBLAS_CFLAGS=""

if [[ $DRY_RUN == 0 ]]; then
    if pkg-config --exists flexiblas 2>/dev/null; then
        FLEXIBLAS_LFLAGS=$( pkg-config --libs flexiblas )
        FLEXIBLAS_CFLAGS=$( pkg-config --cflags flexiblas )
        echo
        echo "found FlexiBLAS via pkg-config: $FLEXIBLAS_LFLAGS"
        echo
    else
        for dir in /usr /usr/local /opt/flexiblas; do
            for libdir in lib lib64 lib/x86_64-linux-gnu; do
                if [ -f $dir/$libdir/libflexiblas.so ]; then
                    FLEXIBLAS_LFLAGS="-L$dir/$libdir -lflexiblas -llapack"
                    FLEXIBLAS_CFLAGS="-I$dir/include"
                    break 2
                fi
            done
        done

        if [ "$FLEXIBLAS_LFLAGS" == "" ]; then
            echo "FlexiBLAS not found; set PKG_CONFIG_PATH or install FlexiBLAS"
            exit 1
        else
            echo
            echo "found FlexiBLAS in $dir/$libdir"
            echo
        fi
    fi
fi

######################################################################
##
## generate configurations
##
######################################################################

if [[ $DRY_RUN == 0 ]]; then
    cat > laplace.conf <<EOF
[app]   
appl = laplaceslp
adm = std
cluster = h

[arith]
ntile = 64
eps = 1e-6
EOF

    cat > materncov.conf <<EOF
[app]
appl = materncov
adm = weak
cluster = h

[arith]
ntile = 64
eps = 1e-6
EOF

fi

######################################################################
##
## build programs
##
######################################################################

#
# update HLR
#

if [[ $DRY_RUN == 0 ]]; then

    echo
    echo "updating HLR"
    echo
    
    git pull

    if [ "$?" != "0" ]; then
        echo "error while updating HLR"
        exit 1
    fi
fi

#
# define common settings without building anything
#
SCONS="$SCONS programs=approx-mm,approx-lu,uniform-mm,uniform-lu frameworks=seq,tbb buildtype=release hpro=$HPRO_DIR tbb=$TBB_DIR"

build () {
    NAME=$1
    LFLAGS=$2

    $SCONS lapack=user lapackflags="$LFLAGS"

    if [ "$?" != "0" ]; then
        echo "error while building HLR with $LFLAGS"
        exit 1
    fi

    for fwork in seq tbb; do
        for prog in approx/approx-mm-${fwork} approx/approx-lu-${fwork} uniform/uniform-mm-${fwork} uniform/uniform-lu-${fwork} ; do

            # wait till it becomes available (distributed FS issues)
            while ! [ -x programs/$prog ]; do sleep 1 ; done
            
            bprog=`basename programs/$prog`
            mv programs/$prog ${bprog}-${NAME}
        done
    done
}

if [[ $DRY_RUN == 0 ]]; then

    echo
    echo "building programs against FlexiBLAS"
    echo

    build "flexiblas" "$FLEXIBLAS_LFLAGS"

fi

# stop if only building required
[[ $JUST_BUILD == 1 ]] && ( exit 0 )

######################################################################
##
## benchmarks
##
######################################################################

# ARITH=std
# APPROX=svd

# optional numactl call
NUMACTL=""

run_benchmark () {
    prog=$1
    app=$2
    args="$3"

    logbase=$app--$( basename $prog )

    for backend in $FLEXIBLAS_BACKENDS ; do
        if [[ $DRY_RUN == 0 ]]; then
            echo "    FLEXIBLAS=$backend $prog-flexiblas $args"
            FLEXIBLAS=$backend $NUMACTL ${prog}-flexiblas $args > ${logbase}--${backend}.log
        else
            echo "    FLEXIBLAS=$backend $NUMACTL ${prog}-flexiblas $args > ${logbase}--${backend}.log"
        fi
    done
}

########################################
#
# sequential
#
########################################

if [[ $FRAMEWORKS =~ seq ]]; then
    echo
    echo "running sequential benchmarks"
    echo

    COMMON_ARGS="-e 1e-4 --nbench 10 --tbench 1e10 --arith std --approx svd"

    if [[ $APPLICATIONS =~ laplace ]]; then
        echo "  Laplace"
        ARGS="--config laplace.conf --grid sphere-6 $COMMON_ARGS"
        [[ $PROGRAMS =~ approx-mm  ]] && ( run_benchmark ./approx-mm-seq laplace "$ARGS" )
        [[ $PROGRAMS =~ approx-lu  ]] && ( run_benchmark ./approx-lu-seq laplace "$ARGS" )
        [[ $PROGRAMS =~ uniform-mm ]] && ( run_benchmark ./uniform-mm-seq laplace "$ARGS" )
        [[ $PROGRAMS =~ uniform-lu ]] && ( run_benchmark ./uniform-lu-seq laplace "$ARGS" )
    fi

    if [[ $APPLICATIONS =~ materncov ]]; then
        echo "  MaternCov"
        ARGS="--config materncov.conf --grid randcube-32768 $COMMON_ARGS"
        [[ $PROGRAMS =~ approx-mm  ]] && ( run_benchmark ./approx-mm-seq materncov "$ARGS" )
        [[ $PROGRAMS =~ approx-lu  ]] && ( run_benchmark ./approx-lu-seq materncov "$ARGS" )
        [[ $PROGRAMS =~ uniform-mm ]] && ( run_benchmark ./uniform-mm-seq materncov "$ARGS" )
        [[ $PROGRAMS =~ uniform-lu ]] && ( run_benchmark ./uniform-lu-seq materncov "$ARGS" )
    fi
fi

########################################
#
# parallel
#
########################################

if [[ $FRAMEWORKS =~ tbb ]]; then
    NUMACTL="numactl -i all"

    echo
    echo "running parallel benchmarks"
    echo

    COMMON_ARGS="-e 1e-6 --nbench 10 --tbench 1e10 --approx svd"

    if [[ $APPLICATIONS =~ laplace ]]; then
        echo "  Laplace"
        ARGS="--config laplace.conf --grid sphere-8 $COMMON_ARGS"
        [[ $PROGRAMS =~ approx-mm  ]] && ( run_benchmark ./approx-mm-tbb laplace "$ARGS --arith std" )
        [[ $PROGRAMS =~ approx-lu  ]] && ( run_benchmark ./approx-lu-tbb laplace "$ARGS --arith dagstd" )
        [[ $PROGRAMS =~ uniform-mm ]] && ( run_benchmark ./uniform-mm-tbb laplace "$ARGS --arith std" )
    fi

    if [[ $APPLICATIONS =~ materncov ]]; then
        echo "  MaternCov"
        ARGS="--config materncov.conf --grid randcube-131072 $COMMON_ARGS"
        [[ $PROGRAMS =~ approx-mm  ]] && ( run_benchmark ./approx-mm-tbb materncov "$ARGS --arith std" )
        [[ $PROGRAMS =~ approx-lu  ]] && ( run_benchmark ./approx-lu-tbb materncov "$ARGS --arith dagstd" )
        [[ $PROGRAMS =~ uniform-mm ]] && ( run_benchmark ./uniform-mm-tbb materncov "$ARGS --arith std" )
    fi
fi

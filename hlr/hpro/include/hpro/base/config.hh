#ifndef __HPRO_BASE_CONFIG_HH
#define __HPRO_BASE_CONFIG_HH
//
// Project     : HLIBpro
// File        : config.hh
// Description : global variables
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include <string>

#include "hpro/base/types.hh"
#include "hpro/cluster/types.hh"
#include "hpro/blas/types.hh"
#include "hpro/algebra/types.hh"

///////////////////////////////////////////////////////////////////
//
// configuration management of HLIBpro
//
///////////////////////////////////////////////////////////////////

namespace Hpro
{

namespace CFG
{

////////////////////////////////////////////////
//
// initialisation: set up default values
//
////////////////////////////////////////////////

void init ();

////////////////////////////////////////////////
//
// finalisation: reset values
//
////////////////////////////////////////////////

void done ();

////////////////////////////////////////////////
//
// threads
//
////////////////////////////////////////////////

// return number of threads used by HLIBpro
uint  nthreads ();

// set (default) number of threads to be used by HLIBpro
void  set_nthreads ( const uint  n );
    
////////////////////////////////////////////////
//
// misc.
//
////////////////////////////////////////////////

// version information
uint         major_version ();
uint         minor_version ();
uint         patch_level   ();
std::string  version       ();

// verbosity level
extern uint  verbosity;

// define verbosity level
void set_verbosity ( const uint level );

// return true if given verbosity level is reached
inline bool verbose ( const uint level ) { return ( level <= verbosity ); }

// print list of all parameters with values
void print_parameters ();

////////////////////////////////////////////////
//
// parameters for BLAS interface
//
////////////////////////////////////////////////

namespace BLAS
{

// check for and remove zero rows/columns in input matrices of SVD (default: false)
extern bool                  check_zeroes;

// check for INF/NAN values in input data (default: false)
extern bool                  check_inf_nan;

// upper matrix size limit for using *gesvd (*gesdd for larger matrices)
extern idx_t                 gesvd_limit;

// use gesvj instead of gesvd/gesdd
extern bool                  use_gesvj;

// use double precision SVD for single precision types
extern bool                  use_double_prec;

// low-rank approximation method to use (default: SVD)
extern Hpro::approx_t        approx_method;

// low-rank truncation method to use (default: SVD)
extern Hpro::approx_t        trunc_method;

// number of power iteration steps in randomized SVD
extern uint                  power_steps;

// number of samples per step in adaptive randomized SVD
extern uint                  sample_size;

// size of oversampling in randomized SVD
extern uint                  oversampling;

}// namespace BLAS

////////////////////////////////////////////////
//
// parameters for clustering
//
////////////////////////////////////////////////

namespace Cluster
{

// default n_min (default: 60)
extern uint                  nmin;

// default split mode during geometrical clusterung (default: adaptive_split_axis)
extern split_axis_mode_t     split_mode;

// sort sub clusters according to size, e.g. larger clusters first (default: false)
extern bool                  sort_wrt_size;

// synchronise depth of interface clusters with domain clusters in ND case (default: true)
extern bool                  sync_interface_depth;

// adjust bounding box during clustering to set of indices and not as defined by parent
// partitioning (default: true)
extern bool                  adjust_bvol;

// during block cluster tree construction: permit clusters of different level or not
extern cluster_level_mode_t  cluster_level_mode;

// if true, build SCCs before graph partitioning (default: on)
extern bool                  build_scc;

// if true, METIS uses random seed for RNG (default: false)
extern bool                  METIS_random;

}// namespace Cluster

////////////////////////////////////////////////
//
// parameters for matrix building
//
////////////////////////////////////////////////

namespace Build
{

// default flags for recompression and coarsening
extern bool            recompress;
extern bool            coarsen;

// if true, low-rank matrices with large rank (≥ min(nrows,ncols)/2)
// will be converted to dense matrices
extern bool            to_dense;

// switching point from low-rank to dense: rank >= min(rows,cols) * ratio
extern double          to_dense_ratio;

// of true, low-rank matrices will always be converted to dense
extern bool            pure_dense;

// if true, TSparseMatrix is used during matrix construction
extern bool            use_sparsemat;

// if true, TZeroMatrix is used for domain-domain couplings
extern bool            use_zeromat;

// if true, ghost matrices are used during construction
extern bool            use_ghostmat;

// indicate checking of return values of callback functions
extern bool            check_cb_ret;

// symmetrise matrices after building
extern bool            symmetrise;

// maximal rank ratio (lowrank/fullrank) before stopping iteration (default: 0.25)
extern double          aca_max_ratio;

}// namespace Build

////////////////////////////////////////////////
//
// parameters for arithmetics
//
////////////////////////////////////////////////

namespace Arith
{

// recompression of low-rank matrices after non-optimal approximation
extern bool            recompress;

// default absolute error bound for low-rank SVD
extern double          abs_eps;
    
// default flags for coarsening
extern bool            coarsen;

// if true, low-rank matrices with large rank (≥ min(nrows,ncols)/2)
// will be converted to dense matrices
extern bool            to_dense;

// switching point from low-rank to dense: rank >= min(rows,cols) * ratio
extern double          to_dense_ratio;

// default evaluation type of triangular/diagonal matrices
extern eval_type_t     eval_type;

// default storage type for diagonal block during factorization
extern storage_type_t  storage_type;

// maximal size for sequential mode
extern size_t          max_seq_size;

// maximal size for sequential mode for vector operations
extern size_t          max_seq_size_vec;

// use DAG based functions instead of recursion (default: true)
extern bool            use_dag;

// version of DAG system to use (1: old, 2: new; default: 1)
extern uint            dag_version;

// if true, unnecessary nodes are removed from DAGs
extern bool            dag_optimise;

// default maximal matrix size for additional checks
// during arithmetics
extern uint            max_check_size;

// try to fix singular blocks
extern bool            fix_singular;

// try to fix blocks with bad condition
extern bool            fix_bad_cond;

// use accuracy based pseudo inversion instead of real inversion
extern bool            pseudo_inversion;

// use accumulator based arithmetic (default: off)
extern bool            use_accu;

// use lazy evaluation in arithmetic (default: eager)
extern bool            lazy_eval;

// sort updates based on norm before summation (only for lazy evaluation)
extern bool            sort_updates;

// use sum approximation for summing up direct/parent updates in lazy mode
extern bool            sum_approx;

// approximation method to use for sums
extern Hpro::approx_t  sum_apx_type;

// split leaf matrices during updates as long as destination is blocked
extern bool            split_update;

// use accumulator also for dense blocks (default: false)
extern bool            dense_accu;

// symmetrise matrices after factorisation
extern bool            symmetrise;

// enable/disable truncation of low-rank matrices if one summand has zero rank
extern bool            zero_sum_trunc;

// algorithm for triangular vector solves (0: auto, 1: rec, 2: global, 3: dag)
extern uint            vector_solve_method;

}// namespace Arith

////////////////////////////////////////////////
//
// GPU parameters
//
////////////////////////////////////////////////

namespace GPU
{

// lower boundary of matrix size for GPU based QR factorization
extern size_t          fp32_qr_size;
extern size_t          fp64_qr_size;

// lower boundary of matrix size for GPU based SVD factorization
extern size_t          fp32_svd_size;
extern size_t          fp64_svd_size;

}// namespace GPU

////////////////////////////////////////////////
//
// solver parameters
//
////////////////////////////////////////////////

namespace Solver
{

// maximal number of iterations
extern uint            max_iter;

// relative residual reduction, e.g., stop if |r_n| / |r_0| < ε
extern double          rel_res_red;

// absolute residual reduction, e.g., stop if |r_n| < ε
extern double          abs_res_red;

// relative residual growth (divergence), e.g., stop if |r_n| / |r_0| > ε
extern double          rel_res_growth;

// default restart for GMRES iteration
extern uint            gmres_restart;

// initialise start value before iteration (default: true)
extern bool            init_start_value;

// compute exact residual during iteration (default: false)
extern bool            use_exact_residual;

}// namespace Solver

////////////////////////////////////////////////
//
// machine parameters
//
////////////////////////////////////////////////

namespace Mach
{

// indicate availability of SSE2/SSE3/AVX instruction set
bool    has_sse2    ();
bool    has_sse3    ();
bool    has_avx     ();
bool    has_avx2    ();
bool    has_mic     ();
bool    has_avx512f ();
bool    has_vsx     ();
bool    has_neon    ();

// for vector sizes when using vector instructions: defines preferred multiple
// size_t  vec_size_mult ();

// yields size with correct padding w.r.t. SIMD operations
template < typename value_t >
size_t  simd_padded_size ( const size_t  n );

}// namespace Mach

////////////////////////////////////////////////
//
// machine parameters
//
////////////////////////////////////////////////

namespace BEM
{

// default quadrature order
extern uint  quad_order;

// use distance adaptive quadrature order
extern bool  adaptive_quad_order;

// use special vector functions if available (default: on)
extern bool  use_simd;

// use special SSE3 vector functions if available (default: on)
extern bool  use_simd_sse3;

// use special AVX vector functions if available (default: on)
extern bool  use_simd_avx;

// use special AVX2 vector functions if available (default: on)
extern bool  use_simd_avx2;

// use special MIC vector functions if available (default: on)
extern bool  use_simd_mic;

// use special AVX512F vector functions if available (default: on)
extern bool  use_simd_avx512f;

// use special VSX vector functions if available (default: on)
extern bool  use_simd_vsx;

// use special NEON vector functions if available (default: on)
extern bool  use_simd_neon;

}// BEM

////////////////////////////////////////////////
//
// I/O parameters
//
////////////////////////////////////////////////

namespace IO
{

// use Matlab syntax for printing complex numbers, vectors and
// matrices to stdio (default: off)
extern bool            use_matlab_syntax;

// use color in terminal output
extern term_color_t    color_mode;

// use unicode in terminal output
extern term_charset_t  charset_mode;

// permute H-matrix before saving as dense matrix
extern bool            permute_save;

}// namespace IO

}// namespace CFG

//
// for compatibility
//
inline bool verbose ( const uint lvl ) { return CFG::verbose( lvl ); }

}// namespace

#endif // __HPRO_BASE_CONFIG_HH

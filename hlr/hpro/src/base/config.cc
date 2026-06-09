//
// Project     : HLIBpro
// File        : config.cc
// Description : global variables
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include <hpro/config.h>

#include <iostream>

#include <boost/filesystem.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <boost/lexical_cast.hpp>

#include "hpro/base/config.hh"
#include "hpro/base/traits.hh"
#include "hpro/base/error.hh"
#include "hpro/base/System.hh"

namespace Hpro
{

namespace
{

///////////////////////////////////////////////////////////////////
//
// local functions
//
///////////////////////////////////////////////////////////////////

//
// run CPUID call
//
void
cpuid ( int     lvl,
        uint &  eax,
        uint &  ebx,
        uint &  ecx,
        uint &  edx )
{
    #if ( HPRO_HAS_CPUID == 1 ) &&                                   \
        (( defined(__GNUC__) || defined(__INTEL_COMPILER) ) &&  \
         ( defined(LINUX)    || defined(DARWIN) ))

    int  cnt = 0;
    
    __asm__ ( "cpuid"
              : "=a" (eax),
                "=b" (ebx),
                "=c" (ecx),
                "=d" (edx)
              : "a"  (lvl),
                "c"  (cnt) );
    
    // asm volatile ( "movl %%ebx, %%edi;"  // 32bit PIC: don't clobber ebx
    //                "cpuid;"
    //                "movl %%ebx, %%esi;"
    //                "movl %%edi, %%ebx;"
    //                :"+a" (eax), "=S" (ebx), "=c" (ecx), "=d" (edx)
    //                : :"edi" );

    #else

    eax = 0;
    ebx = 0;
    ecx = 0;
    edx = 0;
    
    #endif
}

//
// inquire environment variable and return true if 
// available or <false> otherwise
//
bool
get_environ ( const std::string &  envname,
              std::string &        content )
{
    char  * val = nullptr;
    
#if defined(WINDOWS) || defined(_WIN32) || defined(_WIN64)

    size_t  val_len  = 0;

    if ( ! _dupenv_s( & val, & val_len, envname.c_str() ) && ( val != nullptr ))
    {
        content = val;
        free( val );
    }// if
#else
    if (( val = ::getenv( envname.c_str() ) ) != nullptr )
        content = val;
#endif 
    else
    {
        content = "";
        return false;
    }// else

    return true;
}

//
// read configuration file
//
void
read_config ()
{
    std::string  cfg_filename = ".hlibpro.conf";
    std::string  buffer;

    if ( get_environ( "HPRO_CONFIG", buffer ) )
        cfg_filename = buffer;

    if ( ! boost::filesystem::exists( cfg_filename ) )
        return;
    
    boost::property_tree::ptree  cfg;
    
    boost::property_tree::ini_parser::read_ini( cfg_filename, cfg );

    #define EVAL_CFG( str, var )        var = cfg.get( str, var );
    #define EVAL_CF2( str, var, type )  var = decltype( var )( cfg.get( str, type( var ) ) );
    
    EVAL_CFG( "BLAS.check_zeroes",    CFG::BLAS::check_zeroes );
    EVAL_CFG( "BLAS.check_inf_nan",   CFG::BLAS::check_inf_nan );
    EVAL_CFG( "BLAS.gesvd_limit",     CFG::BLAS::gesvd_limit );
    EVAL_CFG( "BLAS.use_gesvj",       CFG::BLAS::use_gesvj );
    EVAL_CFG( "BLAS.use_double_prec", CFG::BLAS::use_double_prec );
    EVAL_CF2( "BLAS.approx_method",   CFG::BLAS::approx_method, int );
    EVAL_CF2( "BLAS.trunc_method",    CFG::BLAS::trunc_method, int );
    EVAL_CFG( "BLAS.power_steps",     CFG::BLAS::power_steps );
    EVAL_CFG( "BLAS.sample_size",     CFG::BLAS::sample_size );
    EVAL_CFG( "BLAS.oversampling",    CFG::BLAS::oversampling );

    EVAL_CFG( "Cluster.nmin",                 CFG::Cluster::nmin );
    EVAL_CFG( "Cluster.sort_wrt_size",        CFG::Cluster::sort_wrt_size );
    EVAL_CFG( "Cluster.sync_interface_depth", CFG::Cluster::sync_interface_depth );
    EVAL_CFG( "Cluster.adjust_bvol",          CFG::Cluster::adjust_bvol );
    EVAL_CF2( "Cluster.same_cluster_level",   CFG::Cluster::cluster_level_mode, bool ) 
    EVAL_CFG( "Cluster.build_scc",            CFG::Cluster::build_scc );
    EVAL_CFG( "Cluster.METIS_random",         CFG::Cluster::METIS_random );

    EVAL_CFG( "Build.recompress",        CFG::Build::recompress );
    EVAL_CFG( "Build.coarsen_build",     CFG::Build::coarsen );
    EVAL_CFG( "Build.to_dense_build",    CFG::Build::to_dense );
    EVAL_CFG( "Build.to_dense_ratio",    CFG::Build::to_dense_ratio );
    EVAL_CFG( "Build.pure_dense",        CFG::Build::pure_dense );
    EVAL_CFG( "Build.use_sparsemat",     CFG::Build::use_sparsemat );
    EVAL_CFG( "Build.use_zeromat",       CFG::Build::use_zeromat );
    EVAL_CFG( "Build.use_ghostmat",      CFG::Build::use_ghostmat );
    EVAL_CFG( "Build.check_cb_ret",      CFG::Build::check_cb_ret );
    EVAL_CFG( "Build.symmetrise",        CFG::Build::symmetrise );
    EVAL_CFG( "Build.aca_max_ratio",     CFG::Build::aca_max_ratio );

    EVAL_CFG( "Arith.recompress",        CFG::Arith::recompress );
    EVAL_CFG( "Arith.abs_eps",           CFG::Arith::abs_eps );
    EVAL_CFG( "Arith.coarsen",           CFG::Arith::coarsen );
    EVAL_CFG( "Arith.to_dense",          CFG::Arith::to_dense );
    EVAL_CFG( "Arith.to_dense_ratio",    CFG::Arith::to_dense_ratio );
    EVAL_CFG( "Arith.max_seq_size",      CFG::Arith::max_seq_size );
    EVAL_CFG( "Arith.max_seq_size_vec",  CFG::Arith::max_seq_size_vec );
    EVAL_CFG( "Arith.use_dag",           CFG::Arith::use_dag );
    EVAL_CFG( "Arith.dag_version",       CFG::Arith::dag_version );
    EVAL_CFG( "Arith.dag_optimise",      CFG::Arith::dag_optimise );
    EVAL_CFG( "Arith.max_check_size",    CFG::Arith::max_check_size );
    EVAL_CFG( "Arith.pseudo_inversion",  CFG::Arith::pseudo_inversion );
    EVAL_CFG( "Arith.use_accu",          CFG::Arith::use_accu );
    EVAL_CFG( "Arith.lazy_eval",         CFG::Arith::lazy_eval );
    EVAL_CFG( "Arith.sum_approx",        CFG::Arith::sum_approx );
    EVAL_CF2( "Arith.sum_apx_type",      CFG::Arith::sum_apx_type, int );
    EVAL_CFG( "Arith.split_update",      CFG::Arith::split_update );
    EVAL_CFG( "Arith.sort_updates",      CFG::Arith::sort_updates );
    EVAL_CFG( "Arith.dense_accu",        CFG::Arith::dense_accu );
    EVAL_CFG( "Arith.symmetrise",        CFG::Arith::symmetrise );
    EVAL_CFG( "Arith.zero_sum_trunc",    CFG::Arith::zero_sum_trunc );
    EVAL_CFG( "Arith.vector_solve_method", CFG::Arith::vector_solve_method );

    EVAL_CFG( "GPU.fp32_qr_size",          CFG::GPU::fp32_qr_size );
    EVAL_CFG( "GPU.fp64_qr_size",          CFG::GPU::fp64_qr_size );
    EVAL_CFG( "GPU.fp32_svd_size",         CFG::GPU::fp32_svd_size );
    EVAL_CFG( "GPU.fp64_svd_size",         CFG::GPU::fp64_svd_size );

    EVAL_CFG( "Solver.max_iter",           CFG::Solver::max_iter );
    EVAL_CFG( "Solver.rel_res_red",        CFG::Solver::rel_res_red );
    EVAL_CFG( "Solver.abs_res_red",        CFG::Solver::abs_res_red );
    EVAL_CFG( "Solver.rel_res_growth",     CFG::Solver::rel_res_growth );
    EVAL_CFG( "Solver.gmres_restart",      CFG::Solver::gmres_restart );
    EVAL_CFG( "Solver.init_start_value",   CFG::Solver::init_start_value );
    EVAL_CFG( "Solver.use_exact_residual", CFG::Solver::use_exact_residual );

    EVAL_CFG( "BEM.quad_order",          CFG::BEM::quad_order );
    EVAL_CFG( "BEM.adaptive_quad_order", CFG::BEM::adaptive_quad_order );
    EVAL_CFG( "BEM.use_simd",            CFG::BEM::use_simd );
    EVAL_CFG( "BEM.use_simd_sse3",       CFG::BEM::use_simd_sse3 );
    EVAL_CFG( "BEM.use_simd_avx",        CFG::BEM::use_simd_avx );
    EVAL_CFG( "BEM.use_simd_avx2",       CFG::BEM::use_simd_avx2 );
    EVAL_CFG( "BEM.use_simd_mic",        CFG::BEM::use_simd_mic );
    EVAL_CFG( "BEM.use_simd_avx512f",    CFG::BEM::use_simd_avx512f );
    EVAL_CFG( "BEM.use_simd_vsx",        CFG::BEM::use_simd_vsx );
    EVAL_CFG( "BEM.use_simd_neon",       CFG::BEM::use_simd_neon );

    EVAL_CFG( "IO.use_matlab_syntax",  CFG::IO::use_matlab_syntax );
    EVAL_CFG( "IO.permute_save",       CFG::IO::permute_save );

    switch ( cfg.get< int >( "IO.color_mode",  CFG::IO::color_mode ) )
    {
        case 0 : CFG::IO::color_mode = no_color;   break;
        case 1 : CFG::IO::color_mode = use_color;  break;
        default:
        case 2 : CFG::IO::color_mode = auto_color; break;
    }// switch

    switch ( cfg.get< int >( "IO.charset_mode",  CFG::IO::charset_mode ) )
    {
        case 0 : CFG::IO::charset_mode = ascii_charset;   break;
        case 1 : CFG::IO::charset_mode = unicode_charset; break;
        default:
        case 2 : CFG::IO::charset_mode = auto_charset;    break;
    }// switch
    
}

///////////////////////////////////////////////////////////////////
//
// local variables
//
///////////////////////////////////////////////////////////////////

struct local_variables_t
{
    bool    CPUID_SSE2;
    bool    CPUID_SSE3;
    bool    CPUID_AVX;
    bool    CPUID_AVX2;
    bool    CPUID_AVX512F;

    bool    MIC_ARCH;

    bool    VSX;
    bool    NEON;

    // size of SIMD registers in bytes
    size_t  SIMD_SIZE;

    //
    // initialise local variables
    //
    local_variables_t ()
    {
        //
        // default values
        //
    
        CPUID_SSE2     = false;
        CPUID_SSE3     = false;
        CPUID_AVX      = false;
        CPUID_AVX2     = false;
        CPUID_AVX512F  = false;
        
        MIC_ARCH       = false;

        VSX            = false;
        NEON           = false;
        
        SIMD_SIZE      = 8;
        
        //
        // check CPU features
        //
    
        #if HPRO_HAS_CPUID == 1 && ( defined(__GNUC__) || defined(__INTEL_COMPILER) )
        {
            uint  lvl, eax, ebx, ecx, edx;

            // check cpuid lvl
            cpuid( 0, eax, ebx, ecx, edx );

            lvl = eax;
            
            if ( lvl >= 1 )
            {
                cpuid( 1, eax, ebx, ecx, edx );

                #if HPRO_HAS_SSE3 == 1
                
                if (( CPUID_SSE3 = (( ecx & ( 1 <<  0 ) ) != 0 ) ))
                    SIMD_SIZE = 128 / 8;
                
                #endif
                
                #if HPRO_HAS_AVX == 1
                
                if (( CPUID_AVX = (( ecx & ( 1 << 28 ) ) != 0 ) ))
                    SIMD_SIZE = 256 / 8;
                
                #endif
                
                #if HPRO_HAS_AVX2 == 1
                
                if (( CPUID_AVX2 = (( ecx & ( 1 << 12 ) ) != 0 ) ))
                    SIMD_SIZE = 256 / 8;
                
                #endif
            }// if

            if ( lvl >= 7 )
            {
                cpuid( 7, eax, ebx, ecx, edx );
                
                #if HPRO_HAS_AVX512F == 1
                
                if (( CPUID_AVX512F = (( ebx & ( 1 << 16 ) ) != 0 ) ))
                    SIMD_SIZE = 512 / 8;
                
                #endif
            }// if
        }
        #endif

        #if HPRO_HAS_MIC == 1
        
        MIC_ARCH  = true;
        SIMD_SIZE = 512 / 8;
        
        #endif

        #if HPRO_HAS_VSX == 1

        VSX       = true;
        SIMD_SIZE = 128 / 8;
        
        #endif

        #if HPRO_HAS_NEON == 1

        NEON      = true;
        SIMD_SIZE = 128 / 8;
        
        #endif

        if ( verbose( LOG_INFO ) )
        {
            std::cout << "architecture:" << std::endl
                      << "  SSE3   = " << int( CPUID_SSE3 ) << std::endl
                      << "  AVX    = " << int( CPUID_AVX ) << std::endl
                      << "  AVX2   = " << int( CPUID_AVX2 )  << std::endl
                      << "  AVX512 = " << int( CPUID_AVX512F ) << std::endl
                      << "  MIC    = " << int( MIC_ARCH ) << std::endl
                      << "  VSX    = " << int( VSX ) << std::endl
                      << "  NEON   = " << int( NEON ) << std::endl;
        }// if
    }
};

//
// return single instance
//
local_variables_t &
get_local_variables ()
{
    static  local_variables_t  _local_variables;

    return  _local_variables;
}

}// namespace anonymous

///////////////////////////////////////////////////////////////////
//
// global variables
//
///////////////////////////////////////////////////////////////////

namespace
{

const uint    HPRO_MAJOR_VERSION = 3;
const uint    HPRO_MINOR_VERSION = 1;
const uint    HPRO_PATCH_LEVEL   = 2;
const char *  HPRO_VERSION       = "3.1.2";

}// namespace anonymous

///////////////////////////////////////////////////////////////////
//
// configuration management of HLIBpro
//
///////////////////////////////////////////////////////////////////

namespace CFG
{

////////////////////////////////////////////////
//
// initialisation: set up default values
//
////////////////////////////////////////////////

void
init ()
{
    //
    // read config from config file
    //

    read_config();
    
    //
    // set config variables with value of environment variables
    //

    std::string  content;

    #define GET_ENV(  str )             get_environ( str, content )
    #define EVAL_ENV( str, var, type )  { if ( GET_ENV( str ) ) var = boost::lexical_cast< type >( content ); }
    #define EVAL_EN2( str, var, type )  { if ( GET_ENV( str ) ) var = decltype( var )( boost::lexical_cast< type >( content ) ); }
    
    EVAL_ENV( "HPRO_BLAS_check_zeroes",    BLAS::check_zeroes,    bool );
    EVAL_ENV( "HPRO_BLAS_check_inf_nan",   BLAS::check_inf_nan,   bool );
    EVAL_ENV( "HPRO_BLAS_gesvd_limit",     BLAS::gesvd_limit,     uint );
    EVAL_ENV( "HPRO_BLAS_use_gesvj",       BLAS::use_gesvj,       bool );
    EVAL_ENV( "HPRO_BLAS_use_double_prec", BLAS::use_double_prec, bool );
    EVAL_EN2( "HPRO_BLAS_approx_method",   BLAS::approx_method,   int );
    EVAL_EN2( "HPRO_BLAS_trunc_method",    BLAS::trunc_method,    int );
    EVAL_ENV( "HPRO_BLAS_power_steps",     BLAS::power_steps,     uint );
    EVAL_ENV( "HPRO_BLAS_sample_size",     BLAS::sample_size,     uint );
    EVAL_ENV( "HPRO_BLAS_oversampling",    BLAS::oversampling,    uint );

    EVAL_ENV( "HPRO_Cluster_nmin",                 Cluster::nmin,                 uint );
    EVAL_ENV( "HPRO_Cluster_sort_wrt_size",        Cluster::sort_wrt_size,        bool );
    EVAL_ENV( "HPRO_Cluster_sync_interface_depth", Cluster::sync_interface_depth, bool );
    EVAL_ENV( "HPRO_Cluster_adjust_bvol",          Cluster::adjust_bvol,          bool );
    EVAL_EN2( "HPRO_Cluster_same_cluster_level",   Cluster::cluster_level_mode,   bool );
    EVAL_ENV( "HPRO_Cluster_build_scc",            Cluster::build_scc,            bool );
    EVAL_ENV( "HPRO_Cluster_METIS_random",         Cluster::METIS_random,         bool );

    EVAL_ENV( "HPRO_Build_recompress",     Build::recompress,     bool );
    EVAL_ENV( "HPRO_Build_coarsen_build",  Build::coarsen,        bool );
    EVAL_ENV( "HPRO_Build_to_dense_build", Build::to_dense,       bool );
    EVAL_ENV( "HPRO_Build_to_dense_ratio", Build::to_dense_ratio, double );
    EVAL_ENV( "HPRO_Build_pure_dense",     Build::pure_dense,     bool );
    EVAL_ENV( "HPRO_Build_use_sparsemat",  Build::use_sparsemat,  bool );
    EVAL_ENV( "HPRO_Build_use_zeromat",    Build::use_zeromat,    bool );
    EVAL_ENV( "HPRO_Build_use_ghostmat",   Build::use_ghostmat,   bool );
    EVAL_ENV( "HPRO_Build_check_cb_ret",   Build::check_cb_ret,   bool );
    EVAL_ENV( "HPRO_Build_symmetrise",     Build::symmetrise,     bool );
    EVAL_ENV( "HPRO_Build_aca_max_ratio",  Build::aca_max_ratio,  double );
    
    EVAL_ENV( "HPRO_Arith_recompress",       Arith::recompress, bool );
    EVAL_ENV( "HPRO_Arith_abs_eps",          Arith::abs_eps, double );
    EVAL_ENV( "HPRO_Arith_coarsen",          Arith::coarsen, bool );
    EVAL_ENV( "HPRO_Arith_to_dense",         Arith::to_dense, bool );
    EVAL_ENV( "HPRO_Arith_to_dense_ratio",   Arith::to_dense_ratio, double );
    EVAL_ENV( "HPRO_Arith_max_seq_size",     Arith::max_seq_size, uint );
    EVAL_ENV( "HPRO_Arith_max_seq_size_vec", Arith::max_seq_size_vec, uint );
    EVAL_ENV( "HPRO_Arith_use_dag",          Arith::use_dag, bool );
    EVAL_ENV( "HPRO_Arith_dag_version",      Arith::dag_version, uint );
    EVAL_ENV( "HPRO_Arith_dag_optimise",     Arith::dag_optimise, bool );
    EVAL_ENV( "HPRO_Arith_max_check_size",   Arith::max_check_size, uint );
    EVAL_ENV( "HPRO_Arith_pseudo_inversion", Arith::pseudo_inversion, bool );
    EVAL_ENV( "HPRO_Arith_use_accu",         Arith::use_accu, bool );
    EVAL_ENV( "HPRO_Arith_lazy_eval",        Arith::lazy_eval, bool );
    EVAL_ENV( "HPRO_Arith_sum_approx",       Arith::sum_approx, bool );
    EVAL_EN2( "HPRO_Arith_sum_apx_type",     Arith::sum_apx_type, int );
    EVAL_ENV( "HPRO_Arith_split_update",     Arith::split_update, bool );
    EVAL_ENV( "HPRO_Arith_sort_updates",     Arith::sort_updates, bool );
    EVAL_ENV( "HPRO_Arith_dense_accu",       Arith::dense_accu, bool );
    EVAL_ENV( "HPRO_Arith_symmetrise",       Arith::dense_accu, bool );
    EVAL_ENV( "HPRO_Arith_zero_sum_trunc",   Arith::zero_sum_trunc, bool );
    EVAL_ENV( "HPRO_Arith_vector_solve_method", Arith::vector_solve_method, uint );

    EVAL_ENV( "HPRO_GPU_fp32_qr_size",          GPU::fp32_qr_size, size_t );
    EVAL_ENV( "HPRO_GPU_fp64_qr_size",          GPU::fp64_qr_size, size_t );
    EVAL_ENV( "HPRO_GPU_fp32_svd_size",         GPU::fp32_svd_size, size_t );
    EVAL_ENV( "HPRO_GPU_fp64_svd_size",         GPU::fp64_svd_size, size_t );

    EVAL_ENV( "HPRO_Solver_max_iter",            Solver::max_iter,           uint   );
    EVAL_ENV( "HPRO_Solver_rel_res_red",         Solver::rel_res_red,        double );
    EVAL_ENV( "HPRO_Solver_abs_res_red",         Solver::abs_res_red,        double );
    EVAL_ENV( "HPRO_Solver_rel_res_growth",      Solver::rel_res_growth,     double );
    EVAL_ENV( "HPRO_Solver_gmres_restart",       Solver::gmres_restart,      uint   );
    EVAL_ENV( "HPRO_Solver_init_start_value",    Solver::init_start_value,   bool );
    EVAL_ENV( "HPRO_Solver_use_exact_residual",  Solver::use_exact_residual, bool );

    EVAL_ENV( "HPRO_BEM_quad_order",          BEM::quad_order,          uint );
    EVAL_ENV( "HPRO_BEM_adaptive_quad_order", BEM::adaptive_quad_order, bool );
    EVAL_ENV( "HPRO_BEM_use_simd",            BEM::use_simd,            bool );
    EVAL_ENV( "HPRO_BEM_use_simd_sse3",       BEM::use_simd_sse3,       bool );
    EVAL_ENV( "HPRO_BEM_use_simd_avx",        BEM::use_simd_avx,        bool );
    EVAL_ENV( "HPRO_BEM_use_simd_avx2",       BEM::use_simd_avx2,       bool );
    EVAL_ENV( "HPRO_BEM_use_simd_mic",        BEM::use_simd_mic,        bool );
    EVAL_ENV( "HPRO_BEM_use_simd_avx512f",    BEM::use_simd_avx512f,    bool );
    EVAL_ENV( "HPRO_BEM_use_simd_vsx",        BEM::use_simd_vsx,        bool );
    EVAL_ENV( "HPRO_BEM_use_simd_neon",       BEM::use_simd_neon,       bool );

    EVAL_ENV( "HPRO_IO_use_matlab_syntax",  IO::use_matlab_syntax, bool   );
    EVAL_ENV( "HPRO_IO_permute_save",       IO::permute_save,      bool   );
    
    if ( GET_ENV( "HPRO_IO_color_mode" ) )
    {
        switch ( boost::lexical_cast< int   >( content ) )
        {
            case 0 : IO::color_mode = no_color;   break;
            case 1 : IO::color_mode = use_color;  break;
            default:
            case 2 : IO::color_mode = auto_color; break;
        }// switch
    }// if
    if ( GET_ENV( "HPRO_IO_charset_mode" ) )
    {
        switch ( boost::lexical_cast< int   >( content ) )
        {
            case 0 : IO::charset_mode = ascii_charset;   break;
            case 1 : IO::charset_mode = unicode_charset; break;
            default:
            case 2 : IO::charset_mode = auto_charset;    break;
        }// switch
    }// if

    EVAL_ENV( "HPRO_verbosity", verbosity, uint );


    // for compatibility
    EVAL_ENV( "HLIB_BLAS_check_zeroes",    BLAS::check_zeroes,    bool );
    EVAL_ENV( "HLIB_BLAS_check_inf_nan",   BLAS::check_inf_nan,   bool );
    EVAL_ENV( "HLIB_BLAS_gesvd_limit",     BLAS::gesvd_limit,     uint );
    EVAL_ENV( "HLIB_BLAS_use_gesvj",       BLAS::use_gesvj,       bool );
    EVAL_ENV( "HLIB_BLAS_use_double_prec", BLAS::use_double_prec, bool );
    EVAL_EN2( "HLIB_BLAS_approx_method",   BLAS::approx_method,   int );
    EVAL_EN2( "HLIB_BLAS_trunc_method",    BLAS::trunc_method,    int );
    EVAL_ENV( "HLIB_BLAS_power_steps",     BLAS::power_steps,     uint );
    EVAL_ENV( "HLIB_BLAS_sample_size",     BLAS::sample_size,     uint );
    EVAL_ENV( "HLIB_BLAS_oversampling",    BLAS::oversampling,    uint );

    EVAL_ENV( "HLIB_Cluster_nmin",                 Cluster::nmin,                 uint );
    EVAL_ENV( "HLIB_Cluster_sort_wrt_size",        Cluster::sort_wrt_size,        bool );
    EVAL_ENV( "HLIB_Cluster_sync_interface_depth", Cluster::sync_interface_depth, bool );
    EVAL_ENV( "HLIB_Cluster_adjust_bvol",          Cluster::adjust_bvol,          bool );
    EVAL_EN2( "HLIB_Cluster_same_cluster_level",   Cluster::cluster_level_mode,   bool );
    EVAL_ENV( "HLIB_Cluster_build_scc",            Cluster::build_scc,            bool );
    EVAL_ENV( "HLIB_Cluster_METIS_random",         Cluster::METIS_random,         bool );

    EVAL_ENV( "HLIB_Build_recompress",     Build::recompress,     bool );
    EVAL_ENV( "HLIB_Build_coarsen_build",  Build::coarsen,        bool );
    EVAL_ENV( "HLIB_Build_to_dense_build", Build::to_dense,       bool );
    EVAL_ENV( "HLIB_Build_to_dense_ratio", Build::to_dense_ratio, double );
    EVAL_ENV( "HLIB_Build_pure_dense",     Build::pure_dense,     bool );
    EVAL_ENV( "HLIB_Build_use_sparsemat",  Build::use_sparsemat,  bool );
    EVAL_ENV( "HLIB_Build_use_zeromat",    Build::use_zeromat,    bool );
    EVAL_ENV( "HLIB_Build_use_ghostmat",   Build::use_ghostmat,   bool );
    EVAL_ENV( "HLIB_Build_check_cb_ret",   Build::check_cb_ret,   bool );
    EVAL_ENV( "HLIB_Build_symmetrise",     Build::symmetrise,     bool );
    EVAL_ENV( "HLIB_Build_aca_max_ratio",  Build::aca_max_ratio,  double );
    
    EVAL_ENV( "HLIB_Arith_recompress",       Arith::recompress, bool );
    EVAL_ENV( "HLIB_Arith_abs_eps",          Arith::abs_eps, double );
    EVAL_ENV( "HLIB_Arith_coarsen",          Arith::coarsen, bool );
    EVAL_ENV( "HLIB_Arith_to_dense",         Arith::to_dense, bool );
    EVAL_ENV( "HLIB_Arith_to_dense_ratio",   Arith::to_dense_ratio, double );
    EVAL_ENV( "HLIB_Arith_max_seq_size",     Arith::max_seq_size, uint );
    EVAL_ENV( "HLIB_Arith_max_seq_size_vec", Arith::max_seq_size_vec, uint );
    EVAL_ENV( "HLIB_Arith_use_dag",          Arith::use_dag, bool );
    EVAL_ENV( "HLIB_Arith_dag_version",      Arith::dag_version, uint );
    EVAL_ENV( "HLIB_Arith_dag_optimise",     Arith::dag_optimise, bool );
    EVAL_ENV( "HLIB_Arith_max_check_size",   Arith::max_check_size, uint );
    EVAL_ENV( "HLIB_Arith_pseudo_inversion", Arith::pseudo_inversion, bool );
    EVAL_ENV( "HLIB_Arith_use_accu",         Arith::use_accu, bool );
    EVAL_ENV( "HLIB_Arith_lazy_eval",        Arith::lazy_eval, bool );
    EVAL_ENV( "HLIB_Arith_sum_approx",       Arith::sum_approx, bool );
    EVAL_EN2( "HLIB_Arith_sum_apx_type",     Arith::sum_apx_type, int );
    EVAL_ENV( "HLIB_Arith_split_update",     Arith::split_update, bool );
    EVAL_ENV( "HLIB_Arith_sort_updates",     Arith::sort_updates, bool );
    EVAL_ENV( "HLIB_Arith_dense_accu",       Arith::dense_accu, bool );
    EVAL_ENV( "HLIB_Arith_symmetrise",       Arith::dense_accu, bool );
    EVAL_ENV( "HLIB_Arith_zero_sum_trunc",   Arith::zero_sum_trunc, bool );
    EVAL_ENV( "HLIB_Arith_vector_solve_method", Arith::vector_solve_method, uint );

    EVAL_ENV( "HLIB_GPU_fp32_qr_size",          GPU::fp32_qr_size, size_t );
    EVAL_ENV( "HLIB_GPU_fp64_qr_size",          GPU::fp64_qr_size, size_t );
    EVAL_ENV( "HLIB_GPU_fp32_svd_size",         GPU::fp32_svd_size, size_t );
    EVAL_ENV( "HLIB_GPU_fp64_svd_size",         GPU::fp64_svd_size, size_t );

    EVAL_ENV( "HLIB_Solver_max_iter",            Solver::max_iter,           uint   );
    EVAL_ENV( "HLIB_Solver_rel_res_red",         Solver::rel_res_red,        double );
    EVAL_ENV( "HLIB_Solver_abs_res_red",         Solver::abs_res_red,        double );
    EVAL_ENV( "HLIB_Solver_rel_res_growth",      Solver::rel_res_growth,     double );
    EVAL_ENV( "HLIB_Solver_gmres_restart",       Solver::gmres_restart,      uint   );
    EVAL_ENV( "HLIB_Solver_init_start_value",    Solver::init_start_value,   bool );
    EVAL_ENV( "HLIB_Solver_use_exact_residual",  Solver::use_exact_residual, bool );

    EVAL_ENV( "HLIB_BEM_quad_order",          BEM::quad_order,          uint );
    EVAL_ENV( "HLIB_BEM_adaptive_quad_order", BEM::adaptive_quad_order, bool );
    EVAL_ENV( "HLIB_BEM_use_simd",            BEM::use_simd,            bool );
    EVAL_ENV( "HLIB_BEM_use_simd_sse3",       BEM::use_simd_sse3,       bool );
    EVAL_ENV( "HLIB_BEM_use_simd_avx",        BEM::use_simd_avx,        bool );
    EVAL_ENV( "HLIB_BEM_use_simd_avx2",       BEM::use_simd_avx2,       bool );
    EVAL_ENV( "HLIB_BEM_use_simd_mic",        BEM::use_simd_mic,        bool );
    EVAL_ENV( "HLIB_BEM_use_simd_avx512f",    BEM::use_simd_avx512f,    bool );
    EVAL_ENV( "HLIB_BEM_use_simd_vsx",        BEM::use_simd_vsx,        bool );
    EVAL_ENV( "HLIB_BEM_use_simd_neon",       BEM::use_simd_neon,       bool );

    EVAL_ENV( "HLIB_IO_use_matlab_syntax",  IO::use_matlab_syntax, bool   );
    EVAL_ENV( "HLIB_IO_permute_save",       IO::permute_save,      bool   );
    
    if ( GET_ENV( "HLIB_IO_color_mode" ) )
    {
        switch ( boost::lexical_cast< int   >( content ) )
        {
            case 0 : IO::color_mode = no_color;   break;
            case 1 : IO::color_mode = use_color;  break;
            default:
            case 2 : IO::color_mode = auto_color; break;
        }// switch
    }// if
    if ( GET_ENV( "HLIB_IO_charset_mode" ) )
    {
        switch ( boost::lexical_cast< int   >( content ) )
        {
            case 0 : IO::charset_mode = ascii_charset;   break;
            case 1 : IO::charset_mode = unicode_charset; break;
            default:
            case 2 : IO::charset_mode = auto_charset;    break;
        }// switch
    }// if

    EVAL_ENV( "HLIB_verbosity", verbosity, uint );
}

////////////////////////////////////////////////
//
// finalisation: reset values
//
////////////////////////////////////////////////

void
done ()
{
}

////////////////////////////////////////////////
//
// misc.
//
////////////////////////////////////////////////

// version information
uint         major_version () { return HPRO_MAJOR_VERSION; }
uint         minor_version () { return HPRO_MINOR_VERSION; }
std::string  version       () { return HPRO_VERSION; }

// verbosity level
uint  verbosity = 0;

// define verbosity level
void set_verbosity ( const uint level )
{
    verbosity = level;
}

//
// print list of all parameters with values
//
void
print_parameters ()
{
    std::cout << "HPRO_BLAS_check_zeroes             = " << BLAS::check_zeroes << std::endl
              << "HPRO_BLAS_check_inf_nan            = " << BLAS::check_inf_nan << std::endl
              << "HPRO_BLAS_gesvd_limit              = " << BLAS::gesvd_limit << std::endl
              << "HPRO_BLAS_use_gesvj                = " << BLAS::use_gesvj << std::endl
              << "HPRO_BLAS_use_double_prec          = " << BLAS::use_double_prec << std::endl
              << "HPRO_BLAS_approx_method            = " << BLAS::approx_method << std::endl
              << "HPRO_BLAS_trunc_method             = " << BLAS::trunc_method << std::endl
              << "HPRO_BLAS_power_steps              = " << BLAS::power_steps << std::endl
              << "HPRO_BLAS_sample_size              = " << BLAS::sample_size << std::endl
              << "HPRO_BLAS_oversampling             = " << BLAS::oversampling << std::endl
              << std::endl;

    std::cout << "HPRO_Cluster_nmin                  = " << Cluster::nmin << std::endl
              << "HPRO_Cluster_sort_wrt_size         = " << Cluster::sort_wrt_size << std::endl
              << "HPRO_Cluster_sync_interface_depth  = " << Cluster::sync_interface_depth << std::endl
              << "HPRO_Cluster_adjust_bvol           = " << Cluster::adjust_bvol << std::endl
              << "HPRO_Cluster_same_cluster_level    = " << (Cluster::cluster_level_mode != 0 ? true : false) << std::endl
              << "HPRO_Cluster_build_scc             = " << Cluster::build_scc << std::endl
              << "HPRO_Cluster_METIS_random          = " << Cluster::METIS_random << std::endl
              << std::endl;

    std::cout << "HPRO_Build_recompress              = " << Build::recompress << std::endl
              << "HPRO_Build_coarsen                 = " << Build::coarsen << std::endl
              << "HPRO_Build_to_dense                = " << Build::to_dense << std::endl
              << "HPRO_Build_to_dense_ratio          = " << Build::to_dense_ratio << std::endl
              << "HPRO_Build_pure_dense              = " << Build::pure_dense << std::endl
              << "HPRO_Build_use_sparsemat           = " << Build::use_sparsemat << std::endl
              << "HPRO_Build_use_zeromat             = " << Build::use_zeromat << std::endl
              << "HPRO_Build_use_ghostmat            = " << Build::use_ghostmat << std::endl
              << "HPRO_Build_check_cb_ret            = " << Build::check_cb_ret << std::endl
              << "HPRO_Build_symmetrise              = " << Build::symmetrise << std::endl
              << "HPRO_Build_aca_max_ratio           = " << Build::aca_max_ratio << std::endl
              << std::endl;
    
    std::cout << "HPRO_Arith_recompress              = " << Arith::recompress << std::endl
              << "HPRO_Arith_abs_eps                 = " << Arith::abs_eps << std::endl
              << "HPRO_Arith_coarsen                 = " << Arith::coarsen << std::endl
              << "HPRO_Arith_to_dense                = " << Arith::to_dense << std::endl
              << "HPRO_Arith_to_dense_ratio          = " << Arith::to_dense_ratio << std::endl
              << "HPRO_Arith_max_seq_size            = " << Arith::max_seq_size << std::endl
              << "HPRO_Arith_max_seq_size_vec        = " << Arith::max_seq_size_vec << std::endl
              << "HPRO_Arith_use_dag                 = " << Arith::use_dag << std::endl
              << "HPRO_Arith_dag_version             = " << Arith::dag_version << std::endl
              << "HPRO_Arith_dag_optimise            = " << Arith::dag_optimise << std::endl
              << "HPRO_Arith_max_check_size          = " << Arith::max_check_size << std::endl
              << "HPRO_Arith_pseudo_inversion        = " << Arith::pseudo_inversion << std::endl
              << "HPRO_Arith_use_accu                = " << Arith::use_accu << std::endl
              << "HPRO_Arith_lazy_eval               = " << Arith::lazy_eval << std::endl
              << "HPRO_Arith_sum_approx              = " << Arith::sum_approx << std::endl
              << "HPRO_Arith_sum_apx_type            = " << Arith::sum_apx_type << std::endl
              << "HPRO_Arith_split_update            = " << Arith::split_update << std::endl
              << "HPRO_Arith_sort_updates            = " << Arith::sort_updates << std::endl
              << "HPRO_Arith_dense_accu              = " << Arith::dense_accu << std::endl
              << "HPRO_Arith_symmetrise              = " << Arith::symmetrise << std::endl
              << "HPRO_Arith_zero_sum_trunc          = " << Arith::zero_sum_trunc << std::endl
              << "HPRO_Arith_vector_solve_method     = " << Arith::vector_solve_method << std::endl
              << std::endl;

    std::cout << "HPRO_GPU_fp32_qr_size              = " << GPU::fp32_qr_size << std::endl
              << "HPRO_GPU_fp64_qr_size              = " << GPU::fp64_qr_size << std::endl
              << "HPRO_GPU_fp32_svd_size             = " << GPU::fp32_svd_size << std::endl
              << "HPRO_GPU_fp64_svd_size             = " << GPU::fp64_svd_size << std::endl
              << std::endl;
        
    std::cout << "HPRO_Solver_max_iter               = " << Solver::max_iter << std::endl
              << "HPRO_Solver_rel_res_red            = " << Solver::rel_res_red << std::endl
              << "HPRO_Solver_abs_res_red            = " << Solver::abs_res_red << std::endl
              << "HPRO_Solver_rel_res_growth         = " << Solver::rel_res_growth << std::endl
              << "HPRO_Solver_gmres_restart          = " << Solver::gmres_restart << std::endl
              << "HPRO_Solver_init_start_value       = " << Solver::init_start_value << std::endl
              << "HPRO_Solver_use_exact_residual     = " << Solver::use_exact_residual << std::endl
              << std::endl;

    std::cout << "HPRO_BEM_quad_order                = " << BEM::quad_order << std::endl
              << "HPRO_BEM_adaptive_quad_order       = " << BEM::adaptive_quad_order << std::endl
              << "HPRO_BEM_use_simd                  = " << BEM::use_simd << std::endl
              << "HPRO_BEM_use_simd_sse3             = " << BEM::use_simd_sse3 << std::endl
              << "HPRO_BEM_use_simd_avx              = " << BEM::use_simd_avx << std::endl
              << "HPRO_BEM_use_simd_avx2             = " << BEM::use_simd_avx2 << std::endl
              << "HPRO_BEM_use_simd_mic              = " << BEM::use_simd_mic << std::endl
              << "HPRO_BEM_use_simd_avx512f          = " << BEM::use_simd_avx512f << std::endl
              << "HPRO_BEM_use_simd_vsx              = " << BEM::use_simd_vsx << std::endl
              << "HPRO_BEM_use_simd_neon             = " << BEM::use_simd_neon << std::endl
              << std::endl;

    std::cout << "HPRO_IO_use_matlab_syntax          = " << IO::use_matlab_syntax << std::endl
              << "HPRO_IO_color_mode                 = " << IO::color_mode << std::endl
              << "HPRO_IO_charset_mode               = " << IO::charset_mode << std::endl
              << "HPRO_IO_permute_save               = " << IO::permute_save << std::endl
              << std::endl;
}

////////////////////////////////////////////////
//
// parameters for BLAS interface
//
////////////////////////////////////////////////

namespace BLAS
{

// check for and remove zero rows/columns in input matrices of SVD
bool   check_zeroes    = false;

// check for INF/NAN values in input data
bool   check_inf_nan   = false;

// upper matrix size limit for using *gesvd (*gesdd for larger matrices)
idx_t  gesvd_limit     = 1000;

// use gesvj instead of gesvd/gesdd
bool   use_gesvj       = false;

// use double precision SVD for single precision types
bool   use_double_prec = false;

// low-rank approximation method
Hpro::approx_t   approx_method = Hpro::use_svd;

// low-rank truncation method
Hpro::approx_t   trunc_method  = Hpro::use_svd;

// number of power iteration steps in randomized SVD
uint   power_steps     = 0;

// number of samples per step in adaptive randomized SVD
uint   sample_size     = 8;

// size of oversampling in randomized SVD
uint   oversampling    = 0;

}// namespace BLAS

////////////////////////////////////////////////
//
// parameters for clustering
//
////////////////////////////////////////////////

namespace Cluster
{

// default n_min
uint                  nmin                 = 60;

// default split mode during geometrical clusterung (default: adaptive_split_axis)
split_axis_mode_t     split_mode           = adaptive_split_axis;

// sort sub clusters according to size, e.g. larger clusters first
bool                  sort_wrt_size        = false;

// synchronise depth of interface clusters with domain clusters in ND case (default: true)
bool                  sync_interface_depth = true;

// adjust bounding volume during clustering to set of indices and not as defined by parent
// partitioning (default: true)
bool                  adjust_bvol          = true;

// during block cluster tree construction: permit clusters of different level or not
cluster_level_mode_t  cluster_level_mode   = cluster_level_same;

// if true, build SCCs before graph partitioning (default: on)
bool                  build_scc            = true;

// if true, METIS uses random seed for RNG (default: false)
bool                  METIS_random         = false;

}// namespace

////////////////////////////////////////////////
//
// parameters for matrix building
//
////////////////////////////////////////////////
    
namespace Build
{

// default flags for recompression and coarsening
bool            recompress          = true;
bool            coarsen             = false;

// if true, low-rank matrices with large rank (≥ min(nrows,ncols)/2)
// will be converted to dense matrices
bool            to_dense            = true;

// switching point from low-rank to dense: rank >= min(rows,cols) * ratio
double          to_dense_ratio      = 0.5;

// of true, low-rank matrices will always be converted to dense
bool            pure_dense          = false;

// if true, TSparseMatrix is used during matrix construction
bool            use_sparsemat       = false;

// if true, TZeroMatrix is used for domain-domain couplings
bool            use_zeromat         = true;

// if true, ghost matrices are used during construction
bool            use_ghostmat        = false;

// indicate checking of return values of callback functions
bool            check_cb_ret        = true;

// symmetrise matrices after building
bool            symmetrise          = true;

// maximal rank ratio (lowrank/fullrank) before stopping iteration (default: 0.25)
double          aca_max_ratio       = 0.25;

}// namespace Build

////////////////////////////////////////////////
//
// parameters for arithmetics
//
////////////////////////////////////////////////
    
namespace Arith
{

// recompression of low-rank matrices after non-optimal approximation
bool            recompress          = true;

// default absolute error bound for low-rank SVD
double          abs_eps             = 0; // Math::square( Limits::epsilon<real>() );
    
// default flags for coarsening
bool            coarsen             = false;

// if true, low-rank matrices with large rank (≥ min(nrows,ncols)/2)
// will be converted to dense matrices
bool            to_dense            = true;

// switching point from low-rank to dense: rank >= min(rows,cols) * ratio
double          to_dense_ratio      = 0.5;

// use blocked LU factorisation
eval_type_t     eval_type           = block_wise;

// storage type for diagonal block during factorisation
storage_type_t  storage_type        = store_inverse;

// maximal size for sequential mode
size_t          max_seq_size        = 100;

// maximal size for sequential mode for vector operations
size_t          max_seq_size_vec    = 250;

// use DAG based methods instead of recursion
bool            use_dag             = true;

// version of DAG system to use (1: old, 2: new; default: 1)
uint            dag_version         = 2;

// remove unnecessary nodes from DAG
bool            dag_optimise        = false;

// default maximal matrix size for additional checks
// during arithmetics
uint            max_check_size      = 0;

// try to fix singular blocks
bool            fix_singular        = false;

// try to fix blocks with bad condition
bool            fix_bad_cond        = false;

// use accuracy based pseudo inversion instead of real inversion
bool            pseudo_inversion    = false;

// use accumulator based arithmetic
bool            use_accu            = false;

// use lazy version of factorisation
bool            lazy_eval           = false;

// use sum approximation for summing up direct/parent updates in lazy mode
bool            sum_approx          = false;

// approximation method to use for sums
Hpro::approx_t  sum_apx_type        = use_aca;

// split leaf matrices during updates as long as destination is blocked
bool            split_update        = false;

// sort updates based on norm before summation (only for lazy evaluation)
bool            sort_updates        = false;

// use accumulator based arithmetic also for dense matrices
bool            dense_accu          = false;

// symmetrise matrices after factorisation
bool            symmetrise          = true;

// enable/disable truncation of low-rank matrices if one summand has zero rank
bool            zero_sum_trunc      = true;

// algorithm for triangular vector solves (0: auto, 1: rec, 2: global, 3: dag)
uint            vector_solve_method = 0;

}// namespace Arith

////////////////////////////////////////////////
//
// GPU parameters
//
////////////////////////////////////////////////

namespace GPU
{

// lower boundary of matrix size for GPU based QR factorization
size_t          fp32_qr_size  = 2048;
size_t          fp64_qr_size  = 2048;

// lower boundary of matrix size for GPU based SVD factorization
size_t          fp32_svd_size = 512;
size_t          fp64_svd_size = 512;

}// namespace GPU

////////////////////////////////////////////////
//
// solver parameters
//
////////////////////////////////////////////////

namespace Solver
{

// maximal number of iterations
uint    max_iter       = 100;

// relative residual reduction, e.g., stop if |r_n| / |r_0| < ε
double  rel_res_red    = 1e-8;

// absolute residual reduction, e.g., stop if |r_n| < ε
double  abs_res_red    = 1e-14;

// relative residual growth (divergence), e.g., stop if |r_n| / |r_0| > ε
double  rel_res_growth = 1e6;

// default restart for GMRES iteration
uint    gmres_restart  = 20;

// initialise start value before iteration (default: true)
bool    init_start_value = true;

// compute exact residual during iteration (default: false)
bool    use_exact_residual = false;

}// namespace Solver

////////////////////////////////////////////////
//
// machine parameters
//
////////////////////////////////////////////////
    
namespace Mach
{

bool
has_sse2 ()
{
    return get_local_variables().CPUID_SSE2;
}

bool
has_sse3 ()
{
    return get_local_variables().CPUID_SSE3;
}

bool
has_avx ()
{
    return get_local_variables().CPUID_AVX;
}

bool
has_avx2 ()
{
    return get_local_variables().CPUID_AVX2;
}

bool
has_mic ()
{
    return get_local_variables().MIC_ARCH;
}

bool
has_avx512f ()
{
    return get_local_variables().CPUID_AVX512F;
}

bool
has_vsx ()
{
    return get_local_variables().VSX;
}

bool
has_neon ()
{
    return get_local_variables().NEON;
}

// yields size with correct padding w.r.t. SIMD operations
template < typename value_t >
size_t
simd_padded_size ( const size_t  n )
{
    using  real_t = real_type_t< value_t >;
    
    const size_t  mult = get_local_variables().SIMD_SIZE / sizeof(real_t);
    
    return ( n + ( mult - n % mult ) );
}

template size_t simd_padded_size< float  > ( const size_t );
template size_t simd_padded_size< double > ( const size_t );
template size_t simd_padded_size< std::complex< float >  > ( const size_t );
template size_t simd_padded_size< std::complex< double > > ( const size_t );

}// namespace

////////////////////////////////////////////////
//
// parameters for BEM
//
////////////////////////////////////////////////
    
namespace BEM
{

// default quadrature order
uint  quad_order          = 4;

// use distance adaptive quadrature order
bool  adaptive_quad_order = true;

// use special vector functions for SSE2/AVX/MIC instruction sets
bool  use_simd         = true;
bool  use_simd_sse3    = true;
bool  use_simd_avx     = true;
bool  use_simd_avx2    = true;
bool  use_simd_mic     = true;
bool  use_simd_avx512f = true;
bool  use_simd_vsx     = true;
bool  use_simd_neon    = true;

}// namespace BEM

////////////////////////////////////////////////
//
// I/O parameters
//
////////////////////////////////////////////////

namespace IO
{

// use Matlab syntax for stdio
bool            use_matlab_syntax = false;

// use color in terminal output
term_color_t    color_mode        = auto_color;

// use unicode in terminal output
term_charset_t  charset_mode      = auto_charset;

// permute H-matrix before saving as dense matrix
bool            permute_save      = true;

}// namespace IO

}// namespace CFG

}// namespace

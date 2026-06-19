#ifndef __HPRO_TPROGRESSBAR_HH
#define __HPRO_TPROGRESSBAR_HH
//
// Project     : HLIBpro
// File        : TProgressBar.hh
// Description : class for showing the progress in various formats
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include <ostream>

#include "hpro/base/types.hh"
#include "hpro/parallel/TMutex.hh"

namespace Hpro
{

//!
//! \class  TProgressBar
//! \brief  base class defining interface
//!
class TProgressBar
{
protected:
    // minimal, maximal and current value of progress
    double  _min, _max, _curr;

    // for thread-saftety
    TMutex  _mutex;
    
    // is progress bar initialised
    bool    _is_init;

    // signal cancelation request
    bool    _do_cancel;
    
public:
    //
    // constructor and destructor
    //
    
    TProgressBar ();

    TProgressBar ( const double   amin,
                   const double   amax,
                   const double   acurr );

    virtual ~TProgressBar ();

    //
    // get status
    //

    //! return true if progress meter is initialised
    virtual bool is_initialised () const { return _is_init; }

    
    //! signal cancelation request
    virtual void cancel    ()       { _do_cancel = true; }

    //! request cancelation
    virtual bool do_cancel () const { return _do_cancel; }

    
    //! return progress interval (minimum value)
    virtual double min () const { return _min; }

    //! return progress interval (maximum value)
    virtual double max () const { return _max; }

    //! return current value in progress interval
    virtual double val () const { return _curr; }

    
    //! return percentage of progress in interval [0,1]
    virtual double percentage () const { return ( _curr - _min ) / ( _max - _min ); }
    
    //
    // change status
    //

    //! initialise status of progress bar, e.g. [\a min, \a max]
    //! with current value \a curr
    virtual void init ( const double  min,
                        const double  max,
                        const double  curr );

    //! reset status, e.g. set new values without intialisation
    virtual void reset ( const double  min,
                         const double  max,
                         const double  curr );

    //! advance progress by \a f
    virtual void advance ( const double f  = 1.0 );

    //! finish progress bar
    virtual void finish ();

protected:
    //
    // react upon change in status
    //
    virtual void update () = 0;

    //
    // return mutex for locking
    //
    TMutex &  mutex () { return _mutex; }
};

//!
//! \class  TConsoleProgressBar
//! \brief  class for a progress bar printing progress on standard console (or via IO streams)
//!   - the format string follows printf style with
//!     - %b : print horizontal percentage bar with optional length, e.g. %30b for 30 characters
//!     - %B : print vertical percentage bar
//!     - %c : print percentage circle
//!     - %i : print progress indicator
//!     - %p : print progress percentage
//!     - %e : print time to finish (ETA)
//!     - %t : print elapsed time
//!     - %m : print current memory usage (if supported)
//!   - if format == "", a default definition is used
//!
class TConsoleProgressBar : public TProgressBar
{
private:
    // io stream to use
    std::ostream *        _out;
    
    // format of the output
    std::string           _format;
    
    // length of text-line
    uint                  _text_len;

    // last value in output
    double                _old_val;

    // progress indicator status
    int                   _indicator_status;

    // character set to be used
    const term_charset_t  _charset;
    
    // character set to be used
    const term_color_t    _color_mode;
    
public:
    //
    // constructor and destructor
    //

    //! construct progress bar using std::cout output stream
    //! and standard output format
    TConsoleProgressBar ( const term_charset_t  acharset = CFG::IO::charset_mode,
                          const term_color_t    acolor   = CFG::IO::color_mode );

    //! construct progress bar using std::cout output stream 
    //! and output format as defined by \a format
    TConsoleProgressBar ( const std::string &   format,
                          const term_charset_t  acharset = CFG::IO::charset_mode,
                          const term_color_t    acolor   = CFG::IO::color_mode );

    //! construct progress bar with output stream \a out_stream
    //! and standard output format
    TConsoleProgressBar ( std::ostream &        out_stream,
                          const term_charset_t  acharset = CFG::IO::charset_mode,
                          const term_color_t    acolor   = CFG::IO::color_mode );

    //! construct progress bar with output stream \a out_stream
    //! and output format as defined by \a format
    TConsoleProgressBar ( std::ostream &        out_stream,
                          const std::string &   format,
                          const term_charset_t  acharset = CFG::IO::charset_mode,
                          const term_color_t    acolor   = CFG::IO::color_mode );

    //! dtor
    virtual ~TConsoleProgressBar ();

    //
    // change status
    //

    // initialise status
    virtual void init ( const double  min,
                        const double  max,
                        const double  curr );

    // finish output
    virtual void finish ();

protected:
    //
    // react upon change in status
    //
    virtual void update ();

    //
    // clear current line in output
    //
    void clear_line ();
};

// alias for old name
using TStreamProgressBar = TConsoleProgressBar;

//!
//! \class  TCBProgressBar
//! \brief  class for progress bar calling user defined call back function
//!         upon change in status
//!
class TCBProgressBar : public TProgressBar
{
public:
    //!
    //! callback function type
    //! \a values is an array with values as defined by enum below
    //! \a cancel indicates cancelation request by user
    //! \a arg    user arguments to call back function
    using  func_t = void (*) ( const double * values,
                               int *          cancel,
                               void *         arg );

    //!
    //! fields in "values" array of progress bar callback function
    //!
    enum { MIN_VAL     = 0,   // minimal value of progress bar
           MAX_VAL     = 1,   // maximal value of progress bar
           CURRENT_VAL = 2    // current value of progress bar
    };
    
protected:
    // callback function pointer
    const func_t  _cb;

    // argument to call back function
    void *        _cb_arg;
    
public:
    //
    // constructor and destructor
    //

    //! construct progress bar with \a cb as call back function
    //! and \a cb_arg as optional arguments provided by each call
    //! to \a cb
    TCBProgressBar ( const func_t   cb,
                     void *         cb_arg = NULL )
            : _cb(cb), _cb_arg(cb_arg)
    {}

    virtual ~TCBProgressBar () {}

protected:
    //
    // react upon change in status: call external function
    //
    virtual void update ()
    {
        int           acancel   = 0;
        const double  values[3] = { min(), max(), val() };
        
        _cb( values, & acancel, _cb_arg );

        if ( acancel != 0 )
            cancel();
    }

    DISABLE_COPY_OP( TCBProgressBar );
};

//!
//! \class  TPartProgressBar
//! \brief  handles progress in part of a parent progress bar
//!
class TPartProgressBar : public TProgressBar
{
protected:

    //! @cond
    
    // parent progress bar
    TProgressBar *  _parent;

    // local part of the total progress in parent
    const double    _total_part;
    
public:
    //
    // constructor and destructor
    //

    TPartProgressBar ( TProgressBar  *  parent,
                       const double     tot_part );

    virtual ~TPartProgressBar () {}

    //! signal cancelation request
    virtual void cancel    ()       { _parent->cancel(); }

    //! request cancelation
    virtual bool do_cancel () const { return _parent->do_cancel(); }
    
    //! advance progress by \a f
    virtual void advance ( const double f );

protected:
    // react upon change in status
    virtual void update ();

    DISABLE_COPY_OP( TPartProgressBar );
};

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

//
// class for automatically finishing progress meters upon leaving a block
//
struct TProgressBarAutoFinish
{
    // handled progress bar
    TProgressBar *  progress;

    // ctor
    TProgressBarAutoFinish ( TProgressBar *  aprogress )
            : progress( aprogress )
    {}

    // dtor
    ~TProgressBarAutoFinish ()
    {
        if ( progress != nullptr )
            progress->finish();
    }
};

}// namespace Hpro

#endif  // __HPRO_TPROGRESSBAR_HH

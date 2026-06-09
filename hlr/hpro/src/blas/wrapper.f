*
* wrapper for complex dot product to ensure correct
* return mechanism
*

****************************************************************
*
* single precision
*
****************************************************************

      subroutine xcdotc ( n, zx, incx, zy, incy, retval )

      complex        cdotc, zx(*), zy(*), retval
      integer        n, incx, incy

      external       cdotc

      retval = cdotc( n, zx, incx, zy, incy )

      return
      end



      subroutine xcdotu ( n, zx, incx, zy, incy, retval )

      complex        cdotu, zx(*), zy(*), retval
      integer        n, incx, incy

      external       cdotu

      retval = cdotu( n, zx, incx, zy, incy )

      return
      end

****************************************************************
*
* double precision
*
****************************************************************

      subroutine xzdotc ( n, zx, incx, zy, incy, retval )

      double complex zdotc, zx(*), zy(*), retval
      integer        n, incx, incy

      external       zdotc

      retval = zdotc( n, zx, incx, zy, incy )

      return
      end



      subroutine xzdotu ( n, zx, incx, zy, incy, retval )

      double complex zdotu, zx(*), zy(*), retval
      integer        n, incx, incy

      external       zdotu

      retval = zdotu( n, zx, incx, zy, incy )

      return
      end

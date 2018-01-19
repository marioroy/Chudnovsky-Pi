
/*
 * https://sourceforge.net/p/predef/wiki/Architectures/
 */

#if defined( __amd64 ) || defined( __x86_64 ) || \
    defined( _M_IX86 ) || defined( __i386   ) || \
    defined( _M_X64  )

  #if GMP_NUMB_BITS == 32
  # include "x86/longlong.h"
  #endif

  #if GMP_NUMB_BITS == 64
  # include "x86_64/longlong.h"
  #endif

#elif defined( __arm__ ) || defined( _M_ARM )
  #include "arm/longlong.h"

#else
  Error, error, this data is for arm, x86 and x86_64

#endif


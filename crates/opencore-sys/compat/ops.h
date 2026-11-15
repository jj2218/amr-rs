/* Our own shim (not opencore code): opencore 0.1.6 calls `add` (r_fft.cpp, unbuilt upstream)
 * but declares only add_16. C99 forbids the implicit declaration, and GCC >= 14 and Clang >= 16
 * reject it, so this prototype is force-included into the NB build; compat/ops.c defines it. */
#include "typedef.h"
Word16 add(Word16 var1, Word16 var2, Flag *pOverflow);

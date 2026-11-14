/* Our own shim (not opencore code): opencore 0.1.6 ships no l_deposit_h.h, which its unbuilt
 * r_fft.cpp includes. The function itself is common/src/l_deposit_h.cpp. */
#include "typedef.h"
Word32 L_deposit_h(Word16 var1);

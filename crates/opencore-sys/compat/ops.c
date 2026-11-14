/* Our own shim (not opencore code): opencore's unbuilt vad2.cpp and r_fft.cpp call these
 * operators as functions, but opencore 0.1.6 defines them only as `static inline` (basic_op.h)
 * or as add_16. Forward to opencore's own implementations (shim.c exports the inline ones). */
#include "typedef.h"

Word16 add_16(Word16 var1, Word16 var2, Flag *pOverflow);
Word16 shim_shl(Word16 a, Word16 b, Flag *o);
Word32 shim_L_shl(Word32 a, Word16 b, Flag *o);
Word32 shim_L_shr(Word32 a, Word16 b, Flag *o);
Word16 shim_abs_s(Word16 a);

Word16 add(Word16 a, Word16 b, Flag *o) { return add_16(a, b, o); }
Word16 shl(Word16 a, Word16 b, Flag *o) { return shim_shl(a, b, o); }
Word32 L_shl(Word32 a, Word16 b, Flag *o) { return shim_L_shl(a, b, o); }
Word32 L_shr(Word32 a, Word16 b, Flag *o) { return shim_L_shr(a, b, o); }
Word16 abs_s(Word16 a) { return shim_abs_s(a); }

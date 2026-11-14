/* Our own shim (not opencore code): exports opencore's `static inline` basic operators so
 * the Rust port can be differentially tested against them. */
#include "basic_op.h"

Word16 shim_shl(Word16 a, Word16 b, Flag *o) { return shl(a, b, o); }
Word16 shim_mult(Word16 a, Word16 b, Flag *o) { return mult(a, b, o); }
Word32 shim_L_mult(Word16 a, Word16 b, Flag *o) { return L_mult(a, b, o); }
Word32 shim_L_add(Word32 a, Word32 b, Flag *o) { return L_add(a, b, o); }
Word32 shim_L_sub(Word32 a, Word32 b, Flag *o) { return L_sub(a, b, o); }
Word32 shim_L_mac(Word32 a, Word16 b, Word16 c, Flag *o) { return L_mac(a, b, c, o); }
Word32 shim_L_msu(Word32 a, Word16 b, Word16 c, Flag *o) { return L_msu(a, b, c, o); }
Word32 shim_L_shl(Word32 a, Word16 b, Flag *o) { return L_shl(a, b, o); }
Word32 shim_L_shr(Word32 a, Word16 b, Flag *o) { return L_shr(a, b, o); }
Word16 shim_abs_s(Word16 a) { return abs_s(a); }
Word32 shim_Mpy_32(Word16 a, Word16 b, Word16 c, Word16 d, Flag *o) { return Mpy_32(a, b, c, d, o); }
Word32 shim_Mpy_32_16(Word16 a, Word16 b, Word16 c, Flag *o) { return Mpy_32_16(a, b, c, o); }
Word32 shim_Mac_32(Word32 x, Word16 a, Word16 b, Word16 c, Word16 d, Flag *o) { return Mac_32(x, a, b, c, d, o); }
Word32 shim_Mac_32_16(Word32 x, Word16 a, Word16 b, Word16 c, Flag *o) { return Mac_32_16(x, a, b, c, o); }

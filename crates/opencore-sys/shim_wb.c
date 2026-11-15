/* Our own shim (not opencore code): opencore's AMR-WB frame loop (`amrwb/wrapper.cpp`,
 * `D_IF_decode`) on the serial parameter array instead of RFC 4867 octets, so the 3GPP TS
 * 26.174 vectors (one word per bit) feed it directly: NO_DATA/lost frames keep the last mode,
 * homing per TS 26.173, the 2 LSBs cleared (14-bit output). */
#include <stdlib.h>
#include <string.h>
#include "pvamrwbdecoder_api.h"
#include "pvamrwbdecoder.h"
#include "pvamrwbdecoder_cnst.h"
#include "dtx.h"

typedef struct {
    void *st;
    unsigned char *pt_st;
    int16 *scratch;
    int16 mode_old;
    int16 reset_flag;
    int16 reset_flag_old;
} OcWb;

void *oc_wb_init(void) {
    OcWb *s = (OcWb *) calloc(1, sizeof(OcWb));
    s->reset_flag_old = 1;
    s->pt_st = (unsigned char *) malloc(pvDecoder_AmrWbMemRequirements());
    pvDecoder_AmrWb_Init(&s->st, s->pt_st, &s->scratch);
    return s;
}

void oc_wb_exit(void *p) {
    OcWb *s = (OcWb *) p;
    free(s->pt_st);
    free(s);
}

/* `frame_type`: RX_* (dtx.h); `mode`: 0..8, ignored for NO_DATA/lost; `prms`: the frame's
 * bits as BIT_0/BIT_1 words (SID frames: 35); `out`: 320 samples. */
void oc_wb_decode(void *p, int16 mode, int16 frame_type, int16 *prms, int16 *out) {
    OcWb *s = (OcWb *) p;
    int16 i, frame_length;
    if (frame_type == RX_NO_DATA || frame_type == RX_SPEECH_LOST) {
        mode = s->mode_old;
        s->reset_flag = 0;
    } else {
        s->mode_old = mode;
        if (s->reset_flag_old == 1)
            s->reset_flag = pvDecoder_AmrWb_homing_frame_test_first(prms, mode);
    }
    if (s->reset_flag != 0 && s->reset_flag_old != 0) {
        for (i = 0; i < AMR_WB_PCM_FRAME; i++)
            out[i] = EHF_MASK;
    } else {
        pvDecoder_AmrWb(mode, prms, out, &frame_length, s->st, frame_type, s->scratch);
    }
    for (i = 0; i < AMR_WB_PCM_FRAME; i++)
        out[i] &= 0xfffC;
    if (s->reset_flag_old == 0)
        s->reset_flag = pvDecoder_AmrWb_homing_frame_test(prms, mode);
    if (s->reset_flag != 0)
        pvDecoder_AmrWb_Reset(s->st, 1);
    s->reset_flag_old = s->reset_flag;
}

/* Our own shim (not vo-amrwbenc code): vo-amrwbenc's encoder API in its default output format
 * (`VOAMRWB_DEFAULT`): per frame 0x6b21, the TX frame type, the mode, then one word per bit
 * (±127), as the 3GPP TS 26.174 `.cod` vectors store it. The input is copied: vo clears the 2
 * LSBs of the caller's buffer in place. */
#include <stdlib.h>
#include <string.h>
#include "voAMRWB.h"
#include "cmnMemory.h"

typedef struct {
    VO_AUDIO_CODECAPI api;
    VO_HANDLE handle;
    VO_MEM_OPERATOR mem;
    VO_CODEC_INIT_USERDATA user;
} VoWbEnc;

void *vo_wbenc_init(void) {
    VoWbEnc *s = (VoWbEnc *) calloc(1, sizeof(VoWbEnc));
    int format = VOAMRWB_DEFAULT;
    voGetAMRWBEncAPI(&s->api);
    s->mem.Alloc = cmnMemAlloc;
    s->mem.Copy = cmnMemCopy;
    s->mem.Free = cmnMemFree;
    s->mem.Set = cmnMemSet;
    s->mem.Check = cmnMemCheck;
    s->user.memflag = VO_IMF_USERMEMOPERATOR;
    s->user.memData = (VO_PTR) &s->mem;
    s->api.Init(&s->handle, VO_AUDIO_CodingAMRWB, &s->user);
    s->api.SetParam(s->handle, VO_PID_AMRWB_FRAMETYPE, &format);
    return s;
}

void vo_wbenc_exit(void *p) {
    VoWbEnc *s = (VoWbEnc *) p;
    s->api.Uninit(s->handle);
    free(s);
}

/* One 320-sample frame at `mode` (0..8); `out` takes up to 3 + 477 words. Returns the words
 * written. */
int vo_wbenc_encode(void *p, int mode, int dtx, const short *speech, short *out) {
    VoWbEnc *s = (VoWbEnc *) p;
    short pcm[320];
    VO_CODECBUFFER in, o;
    VO_AUDIO_OUTPUTINFO info;
    memcpy(pcm, speech, sizeof pcm);
    s->api.SetParam(s->handle, VO_PID_AMRWB_MODE, &mode);
    s->api.SetParam(s->handle, VO_PID_AMRWB_DTX, &dtx);
    in.Buffer = (unsigned char *) pcm;
    in.Length = sizeof pcm;
    o.Buffer = (unsigned char *) out;
    s->api.SetInputData(s->handle, &in);
    s->api.GetOutputData(s->handle, &o, &info);
    return (int) (o.Length / 2);
}

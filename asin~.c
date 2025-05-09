/*
asin~ - Signal rate arcsine external for Pure Data
2025, Ben Wesch
Functionality:
- Takes one signal inlet (input should be in range -1..1)
- Outputs the arcsine of the input as a signal
- Optional 'turn' argument to output in range 0..1 instead of -π/2..π/2
- Supports multichannel input/output
*/

#include "m_pd.h"
#include <math.h>

static t_class *asin_tilde_class;

typedef struct _asin_tilde {
    t_object x_obj;
    t_sample f;
    float scale;     // 1/2π for turns, 1.0 for radians, 180/π for degrees
} t_asin_tilde;

static t_int *asin_tilde_perform(t_int *w) {
    t_asin_tilde *x = (t_asin_tilde *)(w[1]);
    t_sample *in = (t_sample *)(w[2]);
    t_sample *out = (t_sample *)(w[3]);
    int n = (int)(w[4]);
    int nchans = (int)(w[5]);
    
    while (n--) {
        for (int ch = 0; ch < nchans; ch++) {
            float angle = asinf(in[ch]);
            out[ch] = angle * x->scale;
        }
        in += nchans;
        out += nchans;
    }
    
    return (w+6);
}

static void asin_tilde_dsp(t_asin_tilde *x, t_signal **sp) {
    signal_setmultiout(&sp[1], sp[0]->s_nchans);
    dsp_add(asin_tilde_perform, 5, x,
            sp[0]->s_vec, sp[1]->s_vec, sp[0]->s_n, sp[0]->s_nchans);
}

static void *asin_tilde_new(t_symbol *s, int argc, t_atom *argv) {
    (void)s;
    
    t_asin_tilde *x = (t_asin_tilde *)pd_new(asin_tilde_class);
    outlet_new(&x->x_obj, gensym("signal"));
    
    x->scale = 0.5f / M_PI;  // Default to turns (0..1)
    
    if (argc > 0 && argv[0].a_type == A_SYMBOL) {
        t_symbol *mode = atom_getsymbol(&argv[0]);
        if (mode == gensym("rad")) {
            x->scale = 1.0f;            // Output in radians
        } else if (mode == gensym("deg")) {
            x->scale = 180.0f / M_PI;   // Output in degrees
        }
    }
    
    return (void *)x;
}

void asin_tilde_setup(void) {
    asin_tilde_class = class_new(gensym("asin~"),
        (t_newmethod)asin_tilde_new,
        0, sizeof(t_asin_tilde),
        CLASS_MULTICHANNEL,
        A_GIMME, 0);
    
    class_addmethod(asin_tilde_class, (t_method)asin_tilde_dsp, gensym("dsp"), A_CANT, 0);
    CLASS_MAINSIGNALIN(asin_tilde_class, t_asin_tilde, f);
} 
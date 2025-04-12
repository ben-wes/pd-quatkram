/*
tan~ - Signal rate tangent external for Pure Data
2025, Ben Wesch
Functionality:
- Takes one signal inlet
- Outputs the tangent of the input as a signal
- Default output in turns (0..1)
- Optional 'rad' argument for input in radians
- Optional 'deg' argument for input in degrees
- Supports multichannel input/output
*/

#include "m_pd.h"
#include <math.h>

static t_class *tan_tilde_class;

typedef struct _tan_tilde {
    t_object x_obj;
    t_sample f;
    float scale;     // 2π for turns, 1.0 for radians, π/180 for degrees
} t_tan_tilde;

static t_int *tan_tilde_perform(t_int *w) {
    t_tan_tilde *x = (t_tan_tilde *)(w[1]);
    t_sample *in = (t_sample *)(w[2]);
    t_sample *out = (t_sample *)(w[3]);
    int n = (int)(w[4]);
    int nchans = (int)(w[5]);
    
    while (n--) {
        for (int ch = 0; ch < nchans; ch++) {
            float angle = in[ch] * x->scale;
            out[ch] = tanf(angle);
        }
        in += nchans;
        out += nchans;
    }
    
    return (w+6);
}

static void tan_tilde_dsp(t_tan_tilde *x, t_signal **sp) {
    signal_setmultiout(&sp[1], sp[0]->s_nchans);
    dsp_add(tan_tilde_perform, 5, x,
            sp[0]->s_vec, sp[1]->s_vec, sp[0]->s_n, sp[0]->s_nchans);
}

static void *tan_tilde_new(t_symbol *s, int argc, t_atom *argv) {
    (void)s;
    
    t_tan_tilde *x = (t_tan_tilde *)pd_new(tan_tilde_class);
    outlet_new(&x->x_obj, gensym("signal"));
    
    x->scale = 2.0f * M_PI;  // Default to turns (0..1)
    
    if (argc > 0 && argv[0].a_type == A_SYMBOL) {
        t_symbol *mode = atom_getsymbol(&argv[0]);
        if (mode == gensym("rad")) {
            x->scale = 1.0f;            // Input in radians
        } else if (mode == gensym("deg")) {
            x->scale = M_PI / 180.0f;   // Input in degrees
        }
    }
    
    return (void *)x;
}

void tan_tilde_setup(void) {
    tan_tilde_class = class_new(gensym("tan~"),
        (t_newmethod)tan_tilde_new,
        0, sizeof(t_tan_tilde),
        CLASS_MULTICHANNEL,
        A_GIMME, 0);
    
    class_addmethod(tan_tilde_class, (t_method)tan_tilde_dsp, gensym("dsp"), A_CANT, 0);
    CLASS_MAINSIGNALIN(tan_tilde_class, t_tan_tilde, f);
} 
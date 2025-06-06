/*
atan2~ - Signal rate atan2 external for Pure Data
2024, Assistant
Functionality:
- Takes two signal inlets for y and x
- Outputs the atan2 of these values as a signal
- Optional 'turn' argument to output in range 0..1 instead of 0..2π
- Supports multichannel input/output
- Left inlet (y) determines channel count, right inlet (x) is wrapped if needed
Usage:
1. Create object: [atan2~] or [atan2~ turn]
2. Connect signal inputs to the left (y) and right (x) inlets
3. Use the outlet for the atan2 result
Note: This code was developed with assistance from an AI language model.
*/

#include "m_pd.h"
#include <math.h>

static t_class *atan2_tilde_class;

typedef struct _atan2_tilde {
    t_object x_obj;
    t_sample f_dummy;
    t_inlet *x_in2;          // Inlet for the x value
    float scale;             // 1/2π for turns, 1.0 for radians, 180/π for degrees
} t_atan2_tilde;

static t_int *atan2_tilde_perform(t_int *w) {
    t_atan2_tilde *x = (t_atan2_tilde *)(w[1]);
    t_sample *y_in = (t_sample *)(w[2]);  // y input
    t_sample *x_in = (t_sample *)(w[3]);  // x input
    t_sample *out = (t_sample *)(w[4]);
    int n = (int)(w[5]);
    int nchans = (int)(w[6]);
    int x_nchans = (int)(w[7]);
    
    while (n--) {
        for (int ch = 0; ch < nchans; ch++) {
            // Wrap the x input channel if needed
            int x_ch = ch % x_nchans;
            float angle = atan2f(y_in[ch], x_in[x_ch]);
            if (angle < 0) angle += 2 * M_PI;
            out[ch] = angle * x->scale;
        }
        y_in += nchans;
        x_in += x_nchans;
        out += nchans;
    }
    
    return (w+8);
}

static void atan2_tilde_dsp(t_atan2_tilde *x, t_signal **sp) {
    int y_nchans = sp[0]->s_nchans;
    int x_nchans = sp[1]->s_nchans;
    
    // Left inlet (y) determines output channel count
    signal_setmultiout(&sp[2], y_nchans);
    
    dsp_add(atan2_tilde_perform, 7, x,
            sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec, 
            sp[0]->s_n, y_nchans, x_nchans);
}

static void *atan2_tilde_new(t_symbol *s, int argc, t_atom *argv) {
    (void)s;
    
    t_atan2_tilde *x = (t_atan2_tilde *)pd_new(atan2_tilde_class);
    x->x_in2 = inlet_new(&x->x_obj, &x->x_obj.ob_pd, &s_signal, &s_signal);
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

void atan2_tilde_setup(void) {
    atan2_tilde_class = class_new(gensym("atan2~"),
        (t_newmethod)atan2_tilde_new,
        0, 
        sizeof(t_atan2_tilde),
        CLASS_MULTICHANNEL,
        A_GIMME, 0);
    
    class_addmethod(atan2_tilde_class, (t_method)atan2_tilde_dsp, gensym("dsp"), A_CANT, 0);
    CLASS_MAINSIGNALIN(atan2_tilde_class, t_atan2_tilde, f_dummy);
}
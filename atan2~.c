/*
atan2~ - Signal rate atan2 external for Pure Data

2024, Assistant

Functionality:
- Takes two signal inlets for y and x
- Outputs the atan2 of these values as a signal

Usage:
1. Connect signal inputs to the left (y) and right (x) inlets
2. Use the outlet for the atan2 result

Note: This code was developed with assistance from an AI language model.
*/

#include "m_pd.h"
#include <math.h>

static t_class *atan2_tilde_class;

typedef struct _atan2_tilde {
    t_object x_obj;
    t_sample f_dummy;
    t_inlet *x_in2;    // Inlet for the x value
    t_outlet *x_out;
} t_atan2_tilde;

static t_int *atan2_tilde_perform(t_int *w) {
    t_atan2_tilde *x = (t_atan2_tilde *)(w[1]);
    t_sample *in1 = (t_sample *)(w[2]);  // y input
    t_sample *in2 = (t_sample *)(w[3]);  // x input
    t_sample *out = (t_sample *)(w[4]);
    int n = (int)(w[5]);
    
    while (n--) {
        *out++ = atan2f(*in1++, *in2++);
    }
    
    return (w+6);
}

static void atan2_tilde_dsp(t_atan2_tilde *x, t_signal **sp) {
    dsp_add(atan2_tilde_perform, 5, x,
            sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec, sp[0]->s_n);
}

static void *atan2_tilde_new(void) {
    t_atan2_tilde *x = (t_atan2_tilde *)pd_new(atan2_tilde_class);
    x->x_in2 = inlet_new(&x->x_obj, &x->x_obj.ob_pd, &s_signal, &s_signal);
    x->x_out = outlet_new(&x->x_obj, &s_signal);
    return (void *)x;
}

void atan2_tilde_setup(void) {
    atan2_tilde_class = class_new(gensym("atan2~"),
        (t_newmethod)atan2_tilde_new,
        0, sizeof(t_atan2_tilde),
        CLASS_DEFAULT,
        0);
    
    class_addmethod(atan2_tilde_class, (t_method)atan2_tilde_dsp, gensym("dsp"), A_CANT, 0);
    CLASS_MAINSIGNALIN(atan2_tilde_class, t_atan2_tilde, f_dummy);
}
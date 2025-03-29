/*
atan~ - Signal rate atan external for Pure Data
2024, Assistant
Functionality:
- Takes one signal inlet
- Outputs the arctangent of the input as a signal
- Default output in turns (0..1)
- Optional 'rad' argument for radians (-π/2..π/2)
- Optional 'deg' argument for degrees (-90..90)
Usage:
1. Create object: [atan~] for turns, [atan~ rad] for radians, or [atan~ deg] for degrees
2. Connect signal input
3. Use the outlet for the atan result
*/

#include "m_pd.h"
#include <math.h>

static t_class *atan_tilde_class;

typedef struct _atan_tilde {
    t_object x_obj;
    t_sample f;
    t_outlet *x_out;
    float scale;     // 1/(2π) for turns, 1.0 for radians, 180/π for degrees
} t_atan_tilde;

static t_int *atan_tilde_perform(t_int *w) {
    t_atan_tilde *x = (t_atan_tilde *)(w[1]);
    t_sample *in = (t_sample *)(w[2]);
    t_sample *out = (t_sample *)(w[3]);
    int n = (int)(w[4]);
    
    while (n--) {
        *out++ = atanf(*in++) * x->scale;
    }
    
    return (w+5);
}

static void atan_tilde_dsp(t_atan_tilde *x, t_signal **sp) {
    dsp_add(atan_tilde_perform, 4, x,
            sp[0]->s_vec, sp[1]->s_vec, sp[0]->s_n);
}

static void *atan_tilde_new(t_symbol *s, int argc, t_atom *argv) {
    (void)s;
    
    t_atan_tilde *x = (t_atan_tilde *)pd_new(atan_tilde_class);
    x->x_out = outlet_new(&x->x_obj, &s_signal);
    
    x->scale = 1.0f / (2.0f * M_PI);  // Default to turns (0..1)
    
    // Check for mode argument
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

void atan_tilde_setup(void) {
    atan_tilde_class = class_new(gensym("atan~"),
        (t_newmethod)atan_tilde_new,
        0, sizeof(t_atan_tilde),
        CLASS_DEFAULT,
        A_GIMME, 0);
    
    class_addmethod(atan_tilde_class, (t_method)atan_tilde_dsp, gensym("dsp"), A_CANT, 0);
    CLASS_MAINSIGNALIN(atan_tilde_class, t_atan_tilde, f);
} 
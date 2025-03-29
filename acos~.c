/*
acos~ - Signal rate arccosine external for Pure Data
2025, Ben Wesch
Functionality:
- Takes one signal inlet (input should be in range -1..1)
- Outputs the arccosine of the input as a signal
- Optional 'turn' argument to output in range 0..1 instead of 0..π
*/

#include "m_pd.h"
#include <math.h>

static t_class *acos_tilde_class;

typedef struct _acos_tilde {
    t_object x_obj;
    t_sample f;
    t_outlet *x_out;
    float scale;     // 1/2π for turns, 1.0 for radians, 180/π for degrees
} t_acos_tilde;

static t_int *acos_tilde_perform(t_int *w) {
    t_acos_tilde *x = (t_acos_tilde *)(w[1]);
    t_sample *in = (t_sample *)(w[2]);
    t_sample *out = (t_sample *)(w[3]);
    int n = (int)(w[4]);
    
    while (n--) {
        float angle = acosf(*in++);
        *out++ = angle * x->scale;
    }
    
    return (w+5);
}

static void acos_tilde_dsp(t_acos_tilde *x, t_signal **sp) {
    dsp_add(acos_tilde_perform, 4, x,
            sp[0]->s_vec, sp[1]->s_vec, sp[0]->s_n);
}

static void *acos_tilde_new(t_symbol *s, int argc, t_atom *argv) {
    (void)s;
    
    t_acos_tilde *x = (t_acos_tilde *)pd_new(acos_tilde_class);
    x->x_out = outlet_new(&x->x_obj, &s_signal);
    
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

void acos_tilde_setup(void) {
    acos_tilde_class = class_new(gensym("acos~"),
        (t_newmethod)acos_tilde_new,
        0, sizeof(t_acos_tilde),
        CLASS_DEFAULT,
        A_GIMME, 0);
    
    class_addmethod(acos_tilde_class, (t_method)acos_tilde_dsp, gensym("dsp"), A_CANT, 0);
    CLASS_MAINSIGNALIN(acos_tilde_class, t_acos_tilde, f);
} 
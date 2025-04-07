/*
tanh~ - Signal rate hyperbolic tangent external for Pure Data
2025, Ben Wesch
Functionality:
- Takes one signal inlet
- Outputs the hyperbolic tangent of the input as a signal
- Output is always in range -1..1
*/

#include "m_pd.h"
#include <math.h>

static t_class *tanh_tilde_class;

typedef struct _tanh_tilde {
    t_object x_obj;
    t_float x_f;
} t_tanh_tilde;

static void *tanh_tilde_new(void)
{
    t_tanh_tilde *x = (t_tanh_tilde *)pd_new(tanh_tilde_class);
    outlet_new(&x->x_obj, gensym("signal"));
    return (x);
}

static t_int *tanh_tilde_perform(t_int *w)
{
    t_sample *in = (t_sample *)(w[2]);
    t_sample *out = (t_sample *)(w[3]);
    int n = (int)(w[4]);
    
    while (n--) {
        *out++ = tanhf(*in++);
    }
    
    return (w+5);
}

static void tanh_tilde_dsp(t_tanh_tilde *x, t_signal **sp)
{
    dsp_add(tanh_tilde_perform, 4, x,
            sp[0]->s_vec, sp[1]->s_vec, sp[0]->s_n);
}

void tanh_tilde_setup(void)
{
    tanh_tilde_class = class_new(gensym("tanh~"), 
        (t_newmethod)tanh_tilde_new, 0,
        sizeof(t_tanh_tilde), CLASS_DEFAULT, 0);
    class_addmethod(tanh_tilde_class, (t_method)tanh_tilde_dsp,
        gensym("dsp"), A_CANT, 0);
    CLASS_MAINSIGNALIN(tanh_tilde_class, t_tanh_tilde, x_f);
} 
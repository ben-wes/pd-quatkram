/*
faccbounce~ - Float accumulator and bounce external for Pure Data

2024, Ben Wesch

Functionality:
- Accumulates incoming float values on a sample-by-sample basis
- Bounces the accumulated value between -1 and 1
- Maintains an internal float state that represents the accumulated value

Usage:
1. Send a signal to the inlet (float values to accumulate)
2. Receive the resulting signal from the outlet (accumulated and bounced value)
3. Use 'set <value>' message to reset the internal state to a specific value

Note: This code was developed with assistance from the Anthropic Claude AI language model.
*/

#include "m_pd.h"
#include <math.h>

static t_class *faccbounce_tilde_class;

typedef struct _faccbounce_tilde {
    t_object  x_obj;
    t_sample f_dummy;  // dummy float for signal inlet
    t_float f_accum;   // Accumulated float value
    t_outlet *x_out;   // Signal outlet
} t_faccbounce_tilde;

static float bounce(float value) {
    value = fmodf(value + 1.0f, 4.0f) - 1.0f;
    return (value > 1.0f) ? (2.0f - value) : 
           (value < -1.0f) ? (-2.0f - value) : value;
}

t_int *faccbounce_tilde_perform(t_int *w)
{
    t_faccbounce_tilde *x = (t_faccbounce_tilde *)(w[1]);
    t_sample *in = (t_sample *)(w[2]);
    t_sample *out = (t_sample *)(w[3]);
    int n = (int)(w[4]);
    t_float accum = x->f_accum;

    while (n--) {
        // Accumulate the input
        accum += *in++;
        // Bounce between -1 and 1
        accum = bounce(accum);
        // Output the bounced result
        *out++ = accum;
    }

    // Store the final accumulated state for the next DSP cycle
    x->f_accum = accum;
    return (w+5);
}

void faccbounce_tilde_dsp(t_faccbounce_tilde *x, t_signal **sp)
{
    dsp_add(faccbounce_tilde_perform, 4, x, sp[0]->s_vec, sp[1]->s_vec, sp[0]->s_n);
}

void faccbounce_tilde_set(t_faccbounce_tilde *x, t_floatarg f)
{
    x->f_accum = bounce(f);
}

void *faccbounce_tilde_new(void)
{
    t_faccbounce_tilde *x = (t_faccbounce_tilde *)pd_new(faccbounce_tilde_class);
    // Initialize accumulated value to 0
    x->f_accum = 0;
    x->x_out = outlet_new(&x->x_obj, &s_signal);
    return (void *)x;
}

void faccbounce_tilde_setup(void)
{
    faccbounce_tilde_class = class_new(gensym("faccbounce~"),
        (t_newmethod)faccbounce_tilde_new,
        0, sizeof(t_faccbounce_tilde),
        CLASS_DEFAULT,
        0);
    class_addmethod(faccbounce_tilde_class, (t_method)faccbounce_tilde_dsp, gensym("dsp"), A_CANT, 0);
    class_addmethod(faccbounce_tilde_class, (t_method)faccbounce_tilde_set, gensym("set"), A_FLOAT, 0);
    CLASS_MAINSIGNALIN(faccbounce_tilde_class, t_faccbounce_tilde, f_dummy);
}
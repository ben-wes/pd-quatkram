/*
nchans~ - Multichannel counter external for Pure Data
2024, Ben Wesch
Functionality:
- Outputs the number of input channels on DSP activation
- Outputs the number of input channels when receiving a bang
Usage:
1. Create object with [nchans~]
2. Connect signal inputs
3. Receive channel count from outlet on DSP activation or bang
Note: This code was developed with assistance from the Anthropic Claude AI language model.
*/

#include "m_pd.h"

static t_class *nchans_tilde_class;

typedef struct _nchans_tilde {
    t_object x_obj;
    t_sample f_dummy;     // dummy float for signal inlet
    t_outlet *x_outlet;   // outlet for channel count
    int n_channels;       // stored channel count
} t_nchans_tilde;

static t_int *nchans_tilde_perform(t_int *w)
{
    // Minimal perform routine - we don't process audio
    return (w + 2);
}

static void nchans_tilde_dsp(t_nchans_tilde *x, t_signal **sp)
{
    // Store and output the number of channels
    x->n_channels = sp[0]->s_nchans;
    outlet_float(x->x_outlet, (t_float)x->n_channels);
    
    // Add perform routine
    dsp_add(nchans_tilde_perform, 1, x);
}

static void nchans_tilde_bang(t_nchans_tilde *x)
{
    // Output current channel count on bang
    outlet_float(x->x_outlet, (t_float)x->n_channels);
}

static void *nchans_tilde_new(void)
{
    t_nchans_tilde *x = (t_nchans_tilde *)pd_new(nchans_tilde_class);
    
    x->n_channels = 0;
    x->x_outlet = outlet_new(&x->x_obj, &s_float);
    
    return (void *)x;
}

void nchans_tilde_setup(void)
{
    nchans_tilde_class = class_new(gensym("nchans~"),
        (t_newmethod)nchans_tilde_new,
        0,  // no special free needed
        sizeof(t_nchans_tilde),
        CLASS_DEFAULT | CLASS_MULTICHANNEL,
        0); // no creation arguments
    
    class_addmethod(nchans_tilde_class,
        (t_method)nchans_tilde_dsp,
        gensym("dsp"),
        A_CANT,
        0);
        
    class_addmethod(nchans_tilde_class,
        (t_method)nchans_tilde_bang,
        gensym("bang"),
        0);
        
    CLASS_MAINSIGNALIN(nchans_tilde_class, t_nchans_tilde, f_dummy);
}
/*
zc~ - Zero crossing detector external for Pure Data
2024, Ben Wesch
Functionality:
- Detects zero crossings in incoming signals on a sample-by-sample basis
- Outputs impulses for both upward and downward zero crossings
- Supports multiple channels
Usage:
1. Send signals to the inlet
2. Receive impulses from two outlets:
   - Left outlet: upward zero crossings
   - Right outlet: downward zero crossings
Note: This code was developed with assistance from the Anthropic Claude AI language model.
*/
#include "m_pd.h"
#include <math.h>

static t_class *zc_tilde_class;

typedef struct _zc_tilde {
    t_object  x_obj;
    t_sample f_dummy;  // dummy float for signal inlet
    t_sample **in;     // Array of pointers to input channels
    t_sample **out_up;    // Array of pointers to output channels (upward)
    t_sample **out_down;  // Array of pointers to output channels (downward)
    t_float *prev_sample; // Array to store previous sample for each channel
    int n_channels;    // Number of channels
    t_outlet *x_out_down; // Outlet for downward zero crossings
} t_zc_tilde;

t_int *zc_tilde_perform(t_int *w)
{
    t_zc_tilde *x = (t_zc_tilde *)(w[1]);
    int n = (int)(w[2]);
    int i, j;
    
    for (i = 0; i < x->n_channels; i++) {
        t_sample *in = x->in[i];
        t_sample *out_up = x->out_up[i];
        t_sample *out_down = x->out_down[i];
        t_float prev = x->prev_sample[i];
        
        for (j = 0; j < n; j++) {
            t_float curr = in[j];
            if (prev <= 0 && curr > 0) {
                out_up[j] = 1.0f;
            } else {
                out_up[j] = 0.0f;
            }
            if (prev >= 0 && curr < 0) {
                out_down[j] = 1.0f;
            } else {
                out_down[j] = 0.0f;
            }
            prev = curr;
        }
        
        x->prev_sample[i] = prev;
    }
    
    return (w + 3);
}

void zc_tilde_dsp(t_zc_tilde *x, t_signal **sp)
{
    int i;
    int n_channels = sp[0]->s_nchans;
    int vec_size = sp[0]->s_n;

    // Reallocate memory if number of channels has changed
    if (n_channels != x->n_channels) {
        x->in = (t_sample **)resizebytes(x->in, x->n_channels * sizeof(t_sample *), n_channels * sizeof(t_sample *));
        x->out_up = (t_sample **)resizebytes(x->out_up, x->n_channels * sizeof(t_sample *), n_channels * sizeof(t_sample *));
        x->out_down = (t_sample **)resizebytes(x->out_down, x->n_channels * sizeof(t_sample *), n_channels * sizeof(t_sample *));
        x->prev_sample = (t_float *)resizebytes(x->prev_sample, x->n_channels * sizeof(t_float), n_channels * sizeof(t_float));

        // Initialize new channels
        for (i = x->n_channels; i < n_channels; i++) {
            x->prev_sample[i] = 0;
        }

        x->n_channels = n_channels;
    }

    // Assign signal vectors
    for (i = 0; i < n_channels; i++) {
        x->in[i] = sp[0]->s_vec + vec_size * i;
    }

    // Set up output channels
    signal_setmultiout(&sp[1], n_channels);
    signal_setmultiout(&sp[2], n_channels);
    for (i = 0; i < n_channels; i++) {
        x->out_up[i] = sp[1]->s_vec + sp[1]->s_n * i;
        x->out_down[i] = sp[2]->s_vec + sp[2]->s_n * i;
    }

    dsp_add(zc_tilde_perform, 2, x, vec_size);
}

void *zc_tilde_new(void)
{
    t_zc_tilde *x = (t_zc_tilde *)pd_new(zc_tilde_class);
    
    x->n_channels = 0;
    x->in = NULL;
    x->out_up = NULL;
    x->out_down = NULL;
    x->prev_sample = NULL;
    
    outlet_new(&x->x_obj, &s_signal); // Outlet for upward zero crossings
    x->x_out_down = outlet_new(&x->x_obj, &s_signal); // Outlet for downward zero crossings
    return (void *)x;
}

void zc_tilde_free(t_zc_tilde *x)
{
    if (x->in) freebytes(x->in, x->n_channels * sizeof(t_sample *));
    if (x->out_up) freebytes(x->out_up, x->n_channels * sizeof(t_sample *));
    if (x->out_down) freebytes(x->out_down, x->n_channels * sizeof(t_sample *));
    if (x->prev_sample) freebytes(x->prev_sample, x->n_channels * sizeof(t_float));
}

void zc_tilde_setup(void)
{
    zc_tilde_class = class_new(gensym("zc~"),
        (t_newmethod)zc_tilde_new,
        (t_method)zc_tilde_free,
        sizeof(t_zc_tilde),
        CLASS_DEFAULT | CLASS_MULTICHANNEL,
        0);
    
    class_addmethod(zc_tilde_class, (t_method)zc_tilde_dsp, gensym("dsp"), A_CANT, 0);
    CLASS_MAINSIGNALIN(zc_tilde_class, t_zc_tilde, f_dummy);
}
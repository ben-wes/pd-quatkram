#include "m_pd.h"
#include <string.h>

static t_class *mc_route_tilde_class;

typedef struct _mc_route_tilde {
    t_object  x_obj;
    t_sample f;         // Dummy float for main signal inlet
    int n_in;           // Number of input channels
    int n_out;          // Number of output channels (same as n_in)
    int n_route;        // Number of routing channels
    t_sample *input_buffer; // Pre-allocated buffer for input samples
    int buffer_size;    // Size of the input buffer
} t_mc_route_tilde;

static t_int *mc_route_tilde_perform(t_int *w) {
    t_mc_route_tilde *x = (t_mc_route_tilde *)(w[1]);
    int n = (int)(w[2]);
    t_sample *in = (t_sample *)(w[3]);
    t_sample *route = (t_sample *)(w[4]);
    t_sample *out = (t_sample *)(w[5]);
    
    // Copy input samples to our pre-allocated buffer
    memcpy(x->input_buffer, in, x->n_in * n * sizeof(t_sample));
    
    // Perform routing
    for (int i = 0; i < x->n_out; i++) {
        int route_index = i % x->n_route;
        int route_value = (int)(route[route_index*n]) - 1;
        
        if (route_value >= 0 && route_value < x->n_in) {
            memcpy(out + i*n, x->input_buffer + route_value*n, n * sizeof(t_sample));
        } else {
            memset(out + i*n, 0, n * sizeof(t_sample));
        }
    }

    return (w + 6);
}

static void mc_route_tilde_dsp(t_mc_route_tilde *x, t_signal **sp) {
    x->n_in = sp[0]->s_nchans;
    x->n_route = sp[1]->s_nchans;
    x->n_out = x->n_in;

    // Reallocate input buffer if necessary
    int new_buffer_size = x->n_in * sp[0]->s_n;
    if (new_buffer_size > x->buffer_size) {
        x->input_buffer = (t_sample *)resizebytes(x->input_buffer, 
                                                  x->buffer_size * sizeof(t_sample),
                                                  new_buffer_size * sizeof(t_sample));
        x->buffer_size = new_buffer_size;
    }

    // Set up output channels
    signal_setmultiout(&sp[2], x->n_out);

    dsp_add(mc_route_tilde_perform, 5, x, sp[0]->s_n, sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec);
}

static void mc_route_tilde_free(t_mc_route_tilde *x) {
    if (x->input_buffer) freebytes(x->input_buffer, x->buffer_size * sizeof(t_sample));
}

static void *mc_route_tilde_new(void) {
    t_mc_route_tilde *x = (t_mc_route_tilde *)pd_new(mc_route_tilde_class);
    x->n_in = x->n_out = x->n_route = 1;
    x->buffer_size = 64;  // Initial buffer size
    x->input_buffer = (t_sample *)getbytes(x->buffer_size * sizeof(t_sample));
    inlet_new(&x->x_obj, &x->x_obj.ob_pd, &s_signal, &s_signal);
    outlet_new(&x->x_obj, &s_signal);
    return (void *)x;
}

void mc_route_tilde_setup(void) {
    mc_route_tilde_class = class_new(gensym("mc_route~"),
        (t_newmethod)mc_route_tilde_new,
        (t_method)mc_route_tilde_free,
        sizeof(t_mc_route_tilde),
        CLASS_MULTICHANNEL,
        0);
    class_addmethod(mc_route_tilde_class, (t_method)mc_route_tilde_dsp, gensym("dsp"), A_CANT, 0);
    CLASS_MAINSIGNALIN(mc_route_tilde_class, t_mc_route_tilde, f);
}
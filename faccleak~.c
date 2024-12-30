/*
faccleak~ - Leaky integrator external for Pure Data
2024
Functionality:
- Implements a leaky integrator (exponential moving average)
- Signal rate processing for both input and leak rate
- Configurable initial leak rate
- Supports multichannel processing
Usage:
1. Creation args: [leak_rate] (default: 0)
2. Signal inlets: Input signal, leak rate (multichannel)
3. Outlet: Integrated signal (multichannel)
Messages:
- reset [float]: Reset all channels' integrator to specified value
*/

#include "m_pd.h"
#include <math.h>

static t_class *faccleak_tilde_class;

typedef struct _faccleak_tilde {
    t_object x_obj;
    t_sample f_dummy;      // Dummy float for main signal inlet
    t_float *state;        // Array of states for each channel
    t_sample **in;         // Array of pointers to input channels
    t_sample **leak;       // Array of pointers to leak rate channels
    t_sample **out;        // Array of pointers to output channels
    int n_channels;        // Number of channels
    int leak_nchans;       // Actual number of leak channels
    t_float init_leak;     // Initial leak rate value
} t_faccleak_tilde;

// Perform routine - processes audio at signal rate
static t_int *faccleak_tilde_perform(t_int *w)
{
    t_faccleak_tilde *x = (t_faccleak_tilde *)(w[1]);
    int n = (int)(w[2]);
    int i, j;
    int leak_nchans = x->leak_nchans;  // Number of actual leak rate channels
    
    for (i = 0; i < x->n_channels; i++) {
        t_sample *in = x->in[i];
        t_sample *leak = x->leak[i % leak_nchans];  // Wrap leak channel index
        t_sample *out = x->out[i];
        t_float state = x->state[i];
        
        // Process each sample
        for (j = 0; j < n; j++) {
            t_float current_leak = leak[j];
            // Constrain leak rate
            if (current_leak < 0.0f) current_leak = 0.0f;
            if (current_leak > 1.0f) current_leak = 1.0f;
            
            // Update state: new_state = (1-leak)*old_state + input
            state = (1.0f - current_leak) * state + in[j];
            out[j] = state;
        }
        
        x->state[i] = state;  // Store state for next block
    }
    
    return (w + 3);
}

static void faccleak_tilde_dsp(t_faccleak_tilde *x, t_signal **sp)
{
    int i;
    int n_channels = sp[0]->s_nchans;
    int leak_nchans = sp[1]->s_nchans;  // Get number of leak channels
    int vec_size = sp[0]->s_n;
    
    // Reallocate memory if number of channels has changed
    if (n_channels != x->n_channels) {
        x->in = (t_sample **)resizebytes(x->in, 
            x->n_channels * sizeof(t_sample *), 
            n_channels * sizeof(t_sample *));
        x->leak = (t_sample **)resizebytes(x->leak, 
            x->leak_nchans * sizeof(t_sample *), 
            leak_nchans * sizeof(t_sample *));
        x->out = (t_sample **)resizebytes(x->out, 
            x->n_channels * sizeof(t_sample *), 
            n_channels * sizeof(t_sample *));
        x->state = (t_float *)resizebytes(x->state, 
            x->n_channels * sizeof(t_float), 
            n_channels * sizeof(t_float));
        
        // Initialize new channels
        for (i = x->n_channels; i < n_channels; i++) {
            x->state[i] = 0.0f;
        }
        
        x->n_channels = n_channels;
    }
    
    // Update leak channels count
    x->leak_nchans = leak_nchans;
    
    // Assign signal vectors
    for (i = 0; i < n_channels; i++) {
        x->in[i] = sp[0]->s_vec + vec_size * i;
    }
    for (i = 0; i < leak_nchans; i++) {
        x->leak[i] = sp[1]->s_vec + vec_size * i;
    }
    
    // If leak signal is unconnected, use the init value
    if (!sp[1]->s_vec) {
        // Create a constant signal with init_leak value
        for (i = 0; i < n_channels; i++) {
            x->leak[i] = (t_sample *)getbytes(vec_size * sizeof(t_sample));
            if (x->leak[i]) {
                for (int j = 0; j < vec_size; j++) {
                    x->leak[i][j] = x->init_leak;
                }
            }
        }
    }
    
    // Set up output channels
    signal_setmultiout(&sp[2], n_channels);
    for (i = 0; i < n_channels; i++) {
        x->out[i] = sp[2]->s_vec + sp[2]->s_n * i;
    }
    
    dsp_add(faccleak_tilde_perform, 2, x, vec_size);
}

// Message handler for 'reset' message
static void faccleak_tilde_reset(t_faccleak_tilde *x, t_floatarg f)
{
    for (int i = 0; i < x->n_channels; i++) {
        x->state[i] = f;
    }
}

static void *faccleak_tilde_new(t_floatarg leak_rate)
{
    t_faccleak_tilde *x = (t_faccleak_tilde *)pd_new(faccleak_tilde_class);
    
    // Store initial leak rate, with bounds checking
    x->init_leak = leak_rate;
    if (x->init_leak < 0.0f) x->init_leak = 0.0f;
    if (x->init_leak > 1.0f) x->init_leak = 1.0f;
    
    // Initialize
    x->n_channels = 0;
    x->leak_nchans = 1;
    
    x->in = (t_sample **)getbytes(sizeof(t_sample *));
    x->leak = (t_sample **)getbytes(sizeof(t_sample *));
    x->out = (t_sample **)getbytes(sizeof(t_sample *));
    x->state = (t_float *)getbytes(sizeof(t_float));
    
    // Default state for the initial allocation is 0
    if (x->state) x->state[0] = 0.0f;
    
    // Create signal inlet for leak rate with initial value
    signalinlet_new(&x->x_obj, leak_rate);
    
    // Create signal outlet
    outlet_new(&x->x_obj, &s_signal);
    
    return (void *)x;
}

static void faccleak_tilde_free(t_faccleak_tilde *x)
{
    if (x->in) freebytes(x->in, x->n_channels * sizeof(t_sample *));
    if (x->leak) freebytes(x->leak, x->n_channels * sizeof(t_sample *));
    if (x->out) freebytes(x->out, x->n_channels * sizeof(t_sample *));
    if (x->state) freebytes(x->state, x->n_channels * sizeof(t_float));
}

void faccleak_tilde_setup(void)
{
    faccleak_tilde_class = class_new(gensym("faccleak~"),
        (t_newmethod)faccleak_tilde_new,
        (t_method)faccleak_tilde_free,
        sizeof(t_faccleak_tilde),
        CLASS_MULTICHANNEL,
        A_DEFFLOAT, 0);
    
    class_addmethod(faccleak_tilde_class,
        (t_method)faccleak_tilde_dsp, gensym("dsp"), A_CANT, 0);
    class_addmethod(faccleak_tilde_class,
        (t_method)faccleak_tilde_reset, gensym("reset"), A_FLOAT, 0);
    
    CLASS_MAINSIGNALIN(faccleak_tilde_class, t_faccleak_tilde, f_dummy);
}

#include "m_pd.h"
#include <string.h>

static t_class *sampdel_tilde_class;

typedef struct _sampdel_tilde {
    t_object  x_obj;
    t_sample f;             // Dummy float for signal inlet
    t_sample **delay_buffers; // Array of delay buffers, one per channel
    int *write_indices;     // Array of write positions
    int buffer_size;        // Size of each delay buffer
    int n_in;              // Number of input channels
    int n_delays;          // Number of delay control channels
} t_sampdel_tilde;

static t_int *sampdel_tilde_perform(t_int *w) {
    t_sampdel_tilde *x = (t_sampdel_tilde *)(w[1]);
    int n = (int)(w[2]);
    t_sample *in = (t_sample *)(w[3]);
    t_sample *delays = (t_sample *)(w[4]);
    t_sample *out = (t_sample *)(w[5]);
    
    // Determine number of output channels (max of input and delay channels)
    int n_out = x->n_in > x->n_delays ? x->n_in : x->n_delays;
    
    // Process each output channel
    for (int chan = 0; chan < n_out; chan++) {
        // Get corresponding input and delay channels (with wrapping)
        int in_chan = chan % x->n_in;
        int delay_chan = chan % x->n_delays;
        
        // Process each sample in the channel
        for (int i = 0; i < n; i++) {
            // Get delay time from control signal (with bounds checking)
            int delay_samples = (int)delays[delay_chan * n + i];
            if (delay_samples < 0) delay_samples = 0;
            if (delay_samples >= x->buffer_size) delay_samples = x->buffer_size - 1;
            
            // Store input in buffer (wrapped input channel)
            x->delay_buffers[chan][x->write_indices[chan]] = in[in_chan * n + i];
            
            // Calculate read position
            int read_index = x->write_indices[chan] - delay_samples;
            if (read_index < 0) read_index += x->buffer_size;
            
            // Output delayed sample
            out[chan * n + i] = x->delay_buffers[chan][read_index];
            
            // Update write position
            x->write_indices[chan] = (x->write_indices[chan] + 1) % x->buffer_size;
        }
    }
    
    return (w + 6);
}

static void sampdel_tilde_dsp(t_sampdel_tilde *x, t_signal **sp) {
    // Get channel counts from signals
    x->n_in = sp[0]->s_nchans;
    x->n_delays = sp[1]->s_nchans;
    
    // Set up output signal with max number of channels
    int n_out = x->n_in > x->n_delays ? x->n_in : x->n_delays;
    signal_setmultiout(&sp[2], n_out);
    
    // Reallocate buffers for max number of channels
    x->delay_buffers = (t_sample **)resizebytes(x->delay_buffers,
                                               0,
                                               n_out * sizeof(t_sample *));
    x->write_indices = (int *)resizebytes(x->write_indices,
                                         0,
                                         n_out * sizeof(int));
    
    // Initialize new channels if necessary
    for (int i = 0; i < n_out; i++) {
        x->delay_buffers[i] = (t_sample *)getbytes(x->buffer_size * sizeof(t_sample));
        memset(x->delay_buffers[i], 0, x->buffer_size * sizeof(t_sample));
        x->write_indices[i] = 0;
    }
    
    dsp_add(sampdel_tilde_perform, 5,
            x,
            sp[0]->s_n,      // Block size
            sp[0]->s_vec,    // Input signals
            sp[1]->s_vec,    // Delay control signals
            sp[2]->s_vec);   // Output signals
}

static void sampdel_tilde_free(t_sampdel_tilde *x) {
    // Free all delay buffers
    if (x->delay_buffers) {
        int n_channels = x->n_in > x->n_delays ? x->n_in : x->n_delays;
        for (int i = 0; i < n_channels; i++) {
            if (x->delay_buffers[i]) freebytes(x->delay_buffers[i], x->buffer_size * sizeof(t_sample));
        }
        freebytes(x->delay_buffers, n_channels * sizeof(t_sample *));
    }
    
    // Free write indices array
    if (x->write_indices) {
        int n_channels = x->n_in > x->n_delays ? x->n_in : x->n_delays;
        freebytes(x->write_indices, n_channels * sizeof(int));
    }
}

static void *sampdel_tilde_new(t_floatarg max_samples) {
    t_sampdel_tilde *x = (t_sampdel_tilde *)pd_new(sampdel_tilde_class);
    
    // Set buffer size based on max samples (default 44100 if not specified)
    x->buffer_size = max_samples > 0 ? (int)max_samples : 44100;
    
    // Initialize with single channel
    x->n_in = 1;
    x->n_delays = 1;
    
    // Allocate initial buffers
    x->delay_buffers = (t_sample **)getbytes(sizeof(t_sample *));
    x->delay_buffers[0] = (t_sample *)getbytes(x->buffer_size * sizeof(t_sample));
    memset(x->delay_buffers[0], 0, x->buffer_size * sizeof(t_sample));
    
    // Initialize write indices
    x->write_indices = (int *)getbytes(sizeof(int));
    x->write_indices[0] = 0;
    
    // Create signal inlets and outlet
    inlet_new(&x->x_obj, &x->x_obj.ob_pd, &s_signal, &s_signal);
    outlet_new(&x->x_obj, &s_signal);
    
    return (void *)x;
}

void sampdel_tilde_setup(void) {
    sampdel_tilde_class = class_new(gensym("sampdel~"),
        (t_newmethod)sampdel_tilde_new,
        (t_method)sampdel_tilde_free,
        sizeof(t_sampdel_tilde),
        CLASS_MULTICHANNEL, // Enable multichannel support
        A_DEFFLOAT, // Optional argument for maximum samples
        0);
    
    class_addmethod(sampdel_tilde_class,
        (t_method)sampdel_tilde_dsp,
        gensym("dsp"),
        A_CANT,
        0);
    
    CLASS_MAINSIGNALIN(sampdel_tilde_class, t_sampdel_tilde, f);
}

#include "m_pd.h"

static t_class *zcflip_tilde_class;

typedef struct _zcflip_tilde {
    t_object x_obj;
    t_sample f;
    
    // Signal buffers
    t_sample *buffer1;
    t_sample *buffer2;
    int buffer_size;
    
    // Buffer indices
    int write_pos1;
    int write_pos2;
    
    // State tracking
    int current_in;          // which input is active (1 or 2)
    int switch_pending;      // flag for pending switch
    int silent_zc_found;     // flag for found zero crossing in silent signal
    int delay_samples;       // current delay in samples
    int stored_delay;        // stored delay for silent ZC detection
    
    t_outlet *out1;
    t_outlet *out2;
    t_outlet *active_chan;
} t_zcflip_tilde;

static t_int *zcflip_tilde_perform(t_int *w) {
    t_zcflip_tilde *x = (t_zcflip_tilde *)(w[1]);
    t_sample *in1 = (t_sample *)(w[2]);
    t_sample *in2 = (t_sample *)(w[3]);
    t_sample *trig = (t_sample *)(w[4]);
    t_sample *out1 = (t_sample *)(w[5]);
    t_sample *out2 = (t_sample *)(w[6]);
    t_sample *active = (t_sample *)(w[7]);
    int n = (int)(w[8]);
    
    // Process each sample
    for (int i = 0; i < n; i++) {
        // First check for trigger
        if (trig[i] > 0 && !x->switch_pending) {
            x->switch_pending = 1;
            x->silent_zc_found = 0;
            x->stored_delay = 0;
        }
        
        // Write new samples to buffers
        x->buffer1[x->write_pos1] = in1[i];
        x->buffer2[x->write_pos2] = in2[i];
        
        // Monitor zero crossings in the silent signal
        int silent_buffer = (x->current_in == 1) ? 2 : 1;
        t_sample silent_prev = (silent_buffer == 1) ? 
            x->buffer1[(x->write_pos1 - 1 + x->buffer_size) % x->buffer_size] :
            x->buffer2[(x->write_pos2 - 1 + x->buffer_size) % x->buffer_size];
        t_sample silent_curr = (silent_buffer == 1) ? in1[i] : in2[i];
        
        // Track upward zero crossing in silent signal
        if (!x->silent_zc_found && silent_prev <= 0 && silent_curr > 0) {
            x->silent_zc_found = 1;
            x->stored_delay = -1;  // Start counting from -1 to compensate
        }
        
        // If we found a silent zc, increment its delay count
        if (x->silent_zc_found && x->stored_delay < x->buffer_size) {
            x->stored_delay++;
        }
        
        // Check for switching condition
        if (x->switch_pending && x->silent_zc_found) {
            // Check for zero crossing in delayed active signal
            int active_read_pos = (x->current_in == 1) ? 
                x->write_pos1 : x->write_pos2;
            active_read_pos -= x->delay_samples;  // Apply current delay
            if (active_read_pos < 0) active_read_pos += x->buffer_size;
            
            // Get delayed samples for zero crossing check
            t_sample active_prev = (x->current_in == 1) ? 
                x->buffer1[(active_read_pos - 1 + x->buffer_size) % x->buffer_size] :
                x->buffer2[(active_read_pos - 1 + x->buffer_size) % x->buffer_size];
            t_sample active_curr = (x->current_in == 1) ? 
                x->buffer1[active_read_pos] :
                x->buffer2[active_read_pos];
            
            // Switch at active signal's upward zero crossing
            if (active_prev <= 0 && active_curr > 0) {
                x->current_in = (x->current_in == 1) ? 2 : 1;
                x->switch_pending = 0;
                x->delay_samples = x->stored_delay;
            }
        }
        
        // Calculate read positions based on current delay
        int read_pos1 = x->write_pos1 - x->delay_samples;
        int read_pos2 = x->write_pos2 - x->delay_samples;
        if (read_pos1 < 0) read_pos1 += x->buffer_size;
        if (read_pos2 < 0) read_pos2 += x->buffer_size;
        
        // Output based on current active input
        if (x->current_in == 1) {
            out1[i] = x->buffer1[read_pos1];
            out2[i] = 0;
            active[i] = 0;  // Output 0 for channel 1
        } else {
            out1[i] = 0;
            out2[i] = x->buffer2[read_pos2];
            active[i] = 1;  // Output 1 for channel 2
        }
        
        // Update write positions
        x->write_pos1 = (x->write_pos1 + 1) % x->buffer_size;
        x->write_pos2 = (x->write_pos2 + 1) % x->buffer_size;
    }
    
    return (w+9);
}

// Constructor
static void *zcflip_tilde_new(void) {
    t_zcflip_tilde *x = (t_zcflip_tilde *)pd_new(zcflip_tilde_class);
    
    x->buffer_size = 44100;  // 1 second at 44.1kHz
    x->buffer1 = (t_sample *)getbytes(x->buffer_size * sizeof(t_sample));
    x->buffer2 = (t_sample *)getbytes(x->buffer_size * sizeof(t_sample));
    
    x->write_pos1 = 0;
    x->write_pos2 = 0;
    x->current_in = 1;
    x->switch_pending = 0;
    x->silent_zc_found = 0;
    x->delay_samples = 0;
    x->stored_delay = 0;
    
    // Clear buffers
    for (int i = 0; i < x->buffer_size; i++) {
        x->buffer1[i] = 0;
        x->buffer2[i] = 0;
    }
    
    inlet_new(&x->x_obj, &x->x_obj.ob_pd, &s_signal, &s_signal);  // in2
    inlet_new(&x->x_obj, &x->x_obj.ob_pd, &s_signal, &s_signal);  // trigger
    x->out1 = outlet_new(&x->x_obj, &s_signal);
    x->out2 = outlet_new(&x->x_obj, &s_signal);
    x->active_chan = outlet_new(&x->x_obj, &s_signal);
    
    return (x);
}

static void zcflip_tilde_free(t_zcflip_tilde *x) {
    if (x->buffer1) freebytes(x->buffer1, x->buffer_size * sizeof(t_sample));
    if (x->buffer2) freebytes(x->buffer2, x->buffer_size * sizeof(t_sample));
    outlet_free(x->out1);
    outlet_free(x->out2);
    outlet_free(x->active_chan);
}

static void zcflip_tilde_dsp(t_zcflip_tilde *x, t_signal **sp) {
    dsp_add(zcflip_tilde_perform, 8,
            x, sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec,
            sp[3]->s_vec, sp[4]->s_vec, sp[5]->s_vec, sp[0]->s_n);
}

void zcflip_tilde_setup(void) {
    zcflip_tilde_class = class_new(gensym("zcflip~"),
        (t_newmethod)zcflip_tilde_new,
        (t_method)zcflip_tilde_free,
        sizeof(t_zcflip_tilde),
        CLASS_DEFAULT,
        0);
    
    class_addmethod(zcflip_tilde_class,
        (t_method)zcflip_tilde_dsp,
        gensym("dsp"), A_CANT, 0);
    CLASS_MAINSIGNALIN(zcflip_tilde_class, t_zcflip_tilde, f);
}
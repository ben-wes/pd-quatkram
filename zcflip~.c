#include "m_pd.h"
#include <math.h>
#include <string.h>

static t_class *zcflip_tilde_class;

typedef struct _zcflip_tilde {
    t_object x_obj;
    t_sample f;
    
    t_sample *buffers[2];
    int write_pos[2];
    int buffer_size;
    
    int stored_delays[2];
    
    int active_chan;         // which input is active (0 or 1)
    int switch_pending;      // flag for pending switch
    int delay_samples;       // current delay in samples
    
    int has_delay_outs;      // flag for -d parameter
    int reset_silent_delay;  // flag for -r parameter
    float max_step;          // maximum allowed step size at zero-crossing (-s parameter)
    
    int cooldown_samples;    // number of samples to wait after switching (-c parameter)
    int cooldown_counter;    // counts down after a switch
    
    t_outlet *outs[2];
    t_outlet *active_chan_out;
    t_outlet *delay_outs[2];
} t_zcflip_tilde;

static t_int *zcflip_tilde_perform(t_int *w) {
    t_zcflip_tilde *x = (t_zcflip_tilde *)(w[1]);
    t_sample *inputs[2] = {(t_sample *)(w[2]), (t_sample *)(w[3])};
    t_sample *trig = (t_sample *)(w[4]);
    t_sample *outputs[2] = {(t_sample *)(w[5]), (t_sample *)(w[6])};
    t_sample *active = (t_sample *)(w[7]);
    t_sample *delays[2] = {
        x->has_delay_outs ? (t_sample *)(w[8]) : NULL,
        x->has_delay_outs ? (t_sample *)(w[9]) : NULL
    };
    int n = (int)(w[x->has_delay_outs ? 10 : 8]);
    
    for (int i = 0; i < n; i++) {
        // Store samples in circular buffers
        for (int chan = 0; chan < 2; chan++) {
            x->buffers[chan][x->write_pos[chan]] = inputs[chan][i];
            
            // Track zero-crossings
            int prev_pos = (x->write_pos[chan] - 1 + x->buffer_size) % x->buffer_size;
            t_sample prev = x->buffers[chan][prev_pos];
            
            // Check for zero-crossings and reset delays
            if (x->cooldown_counter == 0) {
                t_sample step = fabs(inputs[chan][i] - prev);
                if (prev <= 0 && inputs[chan][i] >= 0 && 
                    (x->max_step <= 0 || step <= x->max_step)) {
                    x->stored_delays[chan] = 0;
                }
            }
            
            // Increment and check delays
            if (x->stored_delays[chan] >= 0) {
                x->stored_delays[chan]++;
                if (x->stored_delays[chan] >= x->buffer_size)
                    x->stored_delays[chan] = -1;
            }
        }
        
        // Handle switching
        if (trig[i] > 0 && !x->switch_pending) x->switch_pending = 1;
        
        // Check for zero-crossing in active channel to complete switch
        if (x->switch_pending) {
            // Check if the channel we would switch to has a valid delay
            int new_delay = x->stored_delays[!x->active_chan] - 1;
            if (new_delay >= 0)
            {
                int active_read_pos = x->write_pos[x->active_chan] - x->delay_samples;
                if (active_read_pos < 0) active_read_pos += x->buffer_size;
                
                int prev_pos = (active_read_pos - 1 + x->buffer_size) % x->buffer_size;
                t_sample active_prev = x->buffers[x->active_chan][prev_pos];
                t_sample active_curr = x->buffers[x->active_chan][active_read_pos];
                
                // Calculate step size for the active channel
                t_sample active_step = fabs(active_curr - active_prev);
                
                // Only switch on zero-crossings with step below the maximum (if max_step > 0)
                if (active_prev <= 0 && active_curr >= 0 && 
                    (x->max_step <= 0 || active_step <= x->max_step))
                {
                    // Switch active channel
                    x->active_chan = !x->active_chan;
                    x->delay_samples = new_delay;
                    x->switch_pending = 0;
                    
                    // Start cooldown
                    x->cooldown_counter = x->cooldown_samples;
                    
                    // If reset_silent_delay is enabled, reset the delay counter for silent channel
                    if (x->reset_silent_delay)
                        x->stored_delays[!x->active_chan] = -1;
                }
            }
        }
        
        // Output
        for (int chan = 0; chan < 2; chan++) {
            int read_pos = (x->write_pos[chan] - x->delay_samples + x->buffer_size) % x->buffer_size;
            outputs[chan][i] = (chan == x->active_chan) ? x->buffers[chan][read_pos] : 0;
            if (x->has_delay_outs) {
                delays[chan][i] = (chan == x->active_chan) ? x->delay_samples : 0;
            }
            x->write_pos[chan] = (x->write_pos[chan] + 1) % x->buffer_size;
        }
        
        active[i] = x->active_chan;
    }
    
    if (x->cooldown_counter > 0) x->cooldown_counter--;
    
    return (w + (x->has_delay_outs ? 11 : 9));
}

static void zcflip_tilde_dsp(t_zcflip_tilde *x, t_signal **sp) {
    if (x->has_delay_outs) {
        dsp_add(zcflip_tilde_perform, 10,
                x, sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec,
                sp[3]->s_vec, sp[4]->s_vec, sp[5]->s_vec,
                sp[6]->s_vec, sp[7]->s_vec, sp[0]->s_n);
    } else {
        dsp_add(zcflip_tilde_perform, 8,
                x, sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec,
                sp[3]->s_vec, sp[4]->s_vec, sp[5]->s_vec, sp[0]->s_n);
    }
}

static void zcflip_tilde_maxstep(t_zcflip_tilde *x, t_floatarg f) {
    x->max_step = f;
}

static void zcflip_tilde_cooldown(t_zcflip_tilde *x, t_floatarg f) {
    // Convert ms to samples, minimum 2
    x->cooldown_samples = (int)(sys_getsr() * f * 0.001f);
    if (x->cooldown_samples < 2) x->cooldown_samples = 2;
}

static void *zcflip_tilde_new(t_symbol *s, int argc, t_atom *argv) {
    (void)s;
    t_zcflip_tilde *x = (t_zcflip_tilde *)pd_new(zcflip_tilde_class);
    
    // Parse arguments
    float buffer_ms = 1000.0f;
    x->has_delay_outs = 0;
    x->reset_silent_delay = 0;
    x->max_step = 0;  // Default: no step limiting
    
    x->cooldown_samples = 2;  // default minimum cooldown
    x->cooldown_counter = 0;
    
    while (argc > 0) {
        if (argv->a_type == A_SYMBOL) {
            t_symbol *arg = atom_getsymbol(argv);
            if (arg == gensym("-d")) {
                x->has_delay_outs = 1;
            }
            else if (arg == gensym("-r")) {
                x->reset_silent_delay = 1;
            }
            else if (arg == gensym("-s") && argc > 1 && argv[1].a_type == A_FLOAT) {
                zcflip_tilde_maxstep(x, atom_getfloat(argv + 1));
                argc--; argv++;
            }
            if (arg == gensym("-c") && argc > 1 && argv[1].a_type == A_FLOAT) {
                zcflip_tilde_cooldown(x, atom_getfloat(argv + 1));
                argc--; argv++;
            }
        }
        else if (argv->a_type == A_FLOAT) {
            buffer_ms = atom_getfloat(argv);
        }
        argc--; argv++;  // Always consume one argument
    }
    
    // Convert ms to samples
    x->buffer_size = (int)(sys_getsr() * buffer_ms * 0.001f);
    
    // Buffer allocation
    for (int i = 0; i < 2; i++) {
        x->buffers[i] = (t_sample *)getbytes(x->buffer_size * sizeof(t_sample));
        memset(x->buffers[i], 0, x->buffer_size * sizeof(t_sample));
        x->write_pos[i] = 0;
        x->stored_delays[i] = -1;
        
        // Create outlets
        x->outs[i] = outlet_new(&x->x_obj, &s_signal);
        if (x->has_delay_outs) {
            x->delay_outs[i] = outlet_new(&x->x_obj, &s_signal);
        }
    }
    
    x->active_chan = 0;
    x->switch_pending = 0;
    x->delay_samples = 0;
    
    inlet_new(&x->x_obj, &x->x_obj.ob_pd, &s_signal, &s_signal);
    inlet_new(&x->x_obj, &x->x_obj.ob_pd, &s_signal, &s_signal);
    x->active_chan_out = outlet_new(&x->x_obj, &s_signal);
    
    return (x);
}

static void zcflip_tilde_free(t_zcflip_tilde *x) {
    for (int i = 0; i < 2; i++) {
        if (x->buffers[i]) freebytes(x->buffers[i], x->buffer_size * sizeof(t_sample));
        outlet_free(x->outs[i]);
        if (x->has_delay_outs) {
            outlet_free(x->delay_outs[i]);
        }
    }
    outlet_free(x->active_chan_out);
}

void zcflip_tilde_setup(void) {
    zcflip_tilde_class = class_new(gensym("zcflip~"),
        (t_newmethod)zcflip_tilde_new,
        (t_method)zcflip_tilde_free,
        sizeof(t_zcflip_tilde),
        CLASS_DEFAULT,
        A_GIMME, 0);

    class_addmethod(zcflip_tilde_class,
        (t_method)zcflip_tilde_dsp,
        gensym("dsp"), A_CANT, 0);

    class_addmethod(zcflip_tilde_class,
        (t_method)zcflip_tilde_maxstep,
        gensym("maxstep"), A_FLOAT, 0);

    class_addmethod(zcflip_tilde_class,
        (t_method)zcflip_tilde_cooldown,
        gensym("cooldown"), A_FLOAT, 0);
        
    CLASS_MAINSIGNALIN(zcflip_tilde_class, t_zcflip_tilde, f);
}
#include "m_pd.h"

static t_class *zcflip_tilde_class;

typedef struct _zcflip_tilde {
    t_object x_obj;
    t_sample f;
    
    t_sample *buffer0;
    t_sample *buffer1; 
    int buffer_size;
    
    int write_pos0;
    int write_pos1;
    
    // State tracking
    int active_chan;         // which input is active (0 or 1)
    int switch_pending;      // flag for pending switch
    int delay_samples;       // current delay in samples
    int stored_delay0;       // stored delay for channel 0
    int stored_delay1;       // stored delay for channel 1
    
    // Flag for delay outlets
    int has_delay_outs;
    
    t_outlet *out0;
    t_outlet *out1;
    t_outlet *active_chan_out;
    t_outlet *delay_out0;
    t_outlet *delay_out1;
} t_zcflip_tilde;

static t_int *zcflip_tilde_perform(t_int *w) {
    t_zcflip_tilde *x = (t_zcflip_tilde *)(w[1]);
    t_sample *in0 = (t_sample *)(w[2]);
    t_sample *in1 = (t_sample *)(w[3]);
    t_sample *trig = (t_sample *)(w[4]);
    t_sample *out0 = (t_sample *)(w[5]);
    t_sample *out1 = (t_sample *)(w[6]);
    t_sample *active = (t_sample *)(w[7]);
    t_sample *delay0 = x->has_delay_outs ? (t_sample *)(w[8]) : NULL;
    t_sample *delay1 = x->has_delay_outs ? (t_sample *)(w[9]) : NULL;
    int n = (int)(w[x->has_delay_outs ? 10 : 8]);
    
    for (int i = 0; i < n; i++) {
        // Store current samples
        x->buffer0[x->write_pos0] = in0[i];
        x->buffer1[x->write_pos1] = in1[i];
        
        // Track zero-crossings for both channels
        t_sample prev0 = x->buffer0[(x->write_pos0 - 1 + x->buffer_size) % x->buffer_size];
        t_sample prev1 = x->buffer1[(x->write_pos1 - 1 + x->buffer_size) % x->buffer_size];
        
        // Check for zero-crossings and reset respective delays
        if (prev0 <= 0 && in0[i] > 0) x->stored_delay0 = 0;
        if (prev1 <= 0 && in1[i] > 0) x->stored_delay1 = 0;
        
        // Increment delays for both channels
        if (x->stored_delay0 >= 0) x->stored_delay0++;
        if (x->stored_delay1 >= 0) x->stored_delay1++;
        
        // Handle switching
        if (trig[i] > 0 && !x->switch_pending) x->switch_pending = 1;
        
        // Check for zero-crossing in active channel to complete switch
        if (x->switch_pending) {
            int active_read_pos = (x->active_chan == 0) ? 
                x->write_pos0 : x->write_pos1;
            active_read_pos -= x->delay_samples;
            if (active_read_pos < 0) active_read_pos += x->buffer_size;
            
            t_sample active_prev = (x->active_chan == 0) ? 
                x->buffer0[(active_read_pos - 1 + x->buffer_size) % x->buffer_size] :
                x->buffer1[(active_read_pos - 1 + x->buffer_size) % x->buffer_size];
            t_sample active_curr = (x->active_chan == 0) ? 
                x->buffer0[active_read_pos] :
                x->buffer1[active_read_pos];
            
            if (active_prev <= 0 && active_curr > 0) {
                x->active_chan = !x->active_chan;
                // Use the stored delay of the channel we're switching to
                x->delay_samples = (x->active_chan ? x->stored_delay1 : x->stored_delay0) - 1;
                x->switch_pending = 0;
            }
        }
        
        // Output using current delay
        int read_pos0 = x->write_pos0 - x->delay_samples;
        int read_pos1 = x->write_pos1 - x->delay_samples;
        if (read_pos0 < 0) read_pos0 += x->buffer_size;
        if (read_pos1 < 0) read_pos1 += x->buffer_size;
        
        active[i] = x->active_chan;
        out0[i] = x->active_chan ? 0 : x->buffer0[read_pos0];
        out1[i] = x->active_chan ? x->buffer1[read_pos1] : 0;
        
        if (x->has_delay_outs) {
            delay0[i] = x->active_chan ? 0 : x->delay_samples;
            delay1[i] = x->active_chan ? x->delay_samples : 0;
        }
        
        x->write_pos0 = (x->write_pos0 + 1) % x->buffer_size;
        x->write_pos1 = (x->write_pos1 + 1) % x->buffer_size;
    }
    
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

static void *zcflip_tilde_new(t_symbol *s, int argc, t_atom *argv) {
    t_zcflip_tilde *x = (t_zcflip_tilde *)pd_new(zcflip_tilde_class);
    
    // Parse arguments
    float buffer_ms = 1000.0f;
    x->has_delay_outs = 0;
    
    while (argc--) {
        if (argv->a_type == A_SYMBOL && 
            gensym("-d") == atom_getsymbol(argv)) {
            x->has_delay_outs = 1;
        }
        else if (argv->a_type == A_FLOAT) {
            buffer_ms = atom_getfloat(argv);
        }
        argv++;
    }
    
    // Convert ms to samples
    x->buffer_size = (int)(sys_getsr() * buffer_ms * 0.001f);
    
    x->buffer0 = (t_sample *)getbytes(x->buffer_size * sizeof(t_sample));
    x->buffer1 = (t_sample *)getbytes(x->buffer_size * sizeof(t_sample));
    
    x->write_pos0 = 0;
    x->write_pos1 = 0;
    x->active_chan = 0;
    x->switch_pending = 0;
    x->delay_samples = 0;
    x->stored_delay0 = -1;
    x->stored_delay1 = -1;
    
    for (int i = 0; i < x->buffer_size; i++) {
        x->buffer0[i] = 0;
        x->buffer1[i] = 0;
    }
    
    inlet_new(&x->x_obj, &x->x_obj.ob_pd, &s_signal, &s_signal);
    inlet_new(&x->x_obj, &x->x_obj.ob_pd, &s_signal, &s_signal);
    x->out0 = outlet_new(&x->x_obj, &s_signal);
    x->out1 = outlet_new(&x->x_obj, &s_signal);
    x->active_chan_out = outlet_new(&x->x_obj, &s_signal);
    
    // Create delay outlets if flag is set
    if (x->has_delay_outs) {
        x->delay_out0 = outlet_new(&x->x_obj, &s_signal);
        x->delay_out1 = outlet_new(&x->x_obj, &s_signal);
    }
    
    return (x);
}

static void zcflip_tilde_free(t_zcflip_tilde *x) {
    if (x->buffer0) freebytes(x->buffer0, x->buffer_size * sizeof(t_sample));
    if (x->buffer1) freebytes(x->buffer1, x->buffer_size * sizeof(t_sample));
    outlet_free(x->out0);
    outlet_free(x->out1);
    outlet_free(x->active_chan_out);
    if (x->has_delay_outs) {
        outlet_free(x->delay_out0);
        outlet_free(x->delay_out1);
    }
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
    CLASS_MAINSIGNALIN(zcflip_tilde_class, t_zcflip_tilde, f);
}
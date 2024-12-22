#include "m_pd.h"
#include <math.h>
#include <stdlib.h>

static t_class *tabsmear_tilde_class;

typedef struct _tabsmear_tilde {
    t_object x_obj;
    t_symbol *x_arrayname;
    t_float x_f;
    t_inlet *x_indexlet;
    t_inlet *x_trigglet;
    
    int is_writing;
    int loop_enabled;
    int mix_mode;
    
    // For averaging at current position
    double current_sum;       // Sum of values at current position
    int current_count;        // Number of values at current position
    int current_pos;          // Current integer position
    
    // For interpolation between positions
    double last_written_val;  // Value written to last position
    int last_written_pos;     // Last position we wrote to
} t_tabsmear_tilde;

static void write_between_positions(t_word *buf, int bufsize, t_tabsmear_tilde *x,
    double pos, double val)
{
    // Handle wrapping for looping
    if (x->loop_enabled) {
        while (pos >= bufsize) pos -= bufsize;
        while (pos < 0) pos += bufsize;
    } else {
        if (pos >= bufsize || pos < 0) return;
        pos = fmin(fmax(pos, 0), bufsize-1);
    }

    int write_pos = (int)floor(pos);
    
    if (!x->is_writing) {
        // First write - just start accumulating
        x->current_pos = write_pos;
        x->current_sum = val;
        x->current_count = 1;
        x->is_writing = 1;
        x->last_written_pos = -1;
        x->last_written_val = 0;
        return;
    }
    
    if (write_pos != x->current_pos) {
        // Write position changed - write accumulated value
        double avg_val = x->current_sum / x->current_count;
        
        // If mix mode, add to existing value
        if (x->mix_mode) {
            buf[x->current_pos].w_float += avg_val;
        } else {
            buf[x->current_pos].w_float = avg_val;
        }
        
        // Determine direction
        int direction = (write_pos > x->current_pos) ? 1 : -1;
        
        // If we skipped positions, interpolate
        if (abs(write_pos - x->current_pos) > 1) {
            double val1 = avg_val;
            double val2 = val;
            
            // Handle wrap-around for looping
            int start = x->current_pos;
            int end = write_pos;
            
            if (x->loop_enabled && direction < 0 && (start - end) > bufsize/2) {
                // Going backward through the loop point
                // First interpolate from current to end of buffer
                for (int i = start - 1; i >= 0; i--) {
                    double alpha = (double)(start - i) / (double)(start + (bufsize - end));
                    double interp_val = val1 * (1.0 - alpha) + val2 * alpha;
                    if (x->mix_mode) {
                        buf[i].w_float += interp_val;
                    } else {
                        buf[i].w_float = interp_val;
                    }
                }
                // Then from start of buffer to target
                for (int i = bufsize - 1; i > end; i--) {
                    double alpha = (double)(start + (bufsize - i)) / (double)(start + (bufsize - end));
                    double interp_val = val1 * (1.0 - alpha) + val2 * alpha;
                    if (x->mix_mode) {
                        buf[i].w_float += interp_val;
                    } else {
                        buf[i].w_float = interp_val;
                    }
                }
            } else {
                // Normal interpolation (including forward wrap-around)
                int step = direction;
                for (int i = start + step; i != end; i += step) {
                    // Handle array bounds for looping
                    int idx = i;
                    if (x->loop_enabled) {
                        while (idx >= bufsize) idx -= bufsize;
                        while (idx < 0) idx += bufsize;
                    }
                    
                    double alpha = direction > 0 ?
                        (double)(i - start) / (end - start) :
                        (double)(start - i) / (start - end);
                    
                    double interp_val = val1 * (1.0 - alpha) + val2 * alpha;
                    if (x->mix_mode) {
                        buf[idx].w_float += interp_val;
                    } else {
                        buf[idx].w_float = interp_val;
                    }
                }
            }
        }
        
        // Remember what we wrote
        x->last_written_pos = x->current_pos;
        x->last_written_val = avg_val;
        
        // Start accumulating at new position
        x->current_pos = write_pos;
        x->current_sum = val;
        x->current_count = 1;
    } else {
        // Still at same position - accumulate
        x->current_sum += val;
        x->current_count++;
    }
}

static void *tabsmear_tilde_new(t_symbol *s, int argc, t_atom *argv)
{
    (void)s; // Silence unused parameter warning
    t_tabsmear_tilde *x = (t_tabsmear_tilde *)pd_new(tabsmear_tilde_class);
    x->x_arrayname = (argc && argv->a_type == A_SYMBOL) ? 
        argv->a_w.w_symbol : gensym("array1");
    
    x->x_indexlet = inlet_new(&x->x_obj, &x->x_obj.ob_pd, &s_signal, &s_signal);
    x->x_trigglet = inlet_new(&x->x_obj, &x->x_obj.ob_pd, &s_signal, &s_signal);
    x->x_f = 0;
    
    x->is_writing = 0;
    x->loop_enabled = 1;
    x->mix_mode = 0;
    
    x->current_sum = 0;
    x->current_count = 0;
    x->current_pos = 0;
    x->last_written_pos = -1;
    x->last_written_val = 0;
    
    return (x);
}

static t_int *tabsmear_tilde_perform(t_int *w)
{
    t_tabsmear_tilde *x = (t_tabsmear_tilde *)(w[1]);
    t_sample *in = (t_sample *)(w[2]);
    t_sample *index = (t_sample *)(w[3]);
    t_sample *trigger = (t_sample *)(w[4]);
    int n = (int)(w[5]);
    t_garray *array;
    t_word *buf;
    int arraysize;

    if (!(array = (t_garray *)pd_findbyclass(x->x_arrayname, garray_class))) {
        pd_error(x, "%s: no such array", x->x_arrayname->s_name);
        return (w+6);
    }
    if (!garray_getfloatwords(array, &arraysize, &buf)) {
        pd_error(x, "%s: bad template", x->x_arrayname->s_name);
        return (w+6);
    }

    while (n--) {
        float idx = *index++;
        float trig = *trigger++;
        float curr_sample = *in++;
        
        // Handle trigger
        if (trig <= 0 || idx < 0) {
            if (x->is_writing) {
                // Write final accumulated value
                if (x->current_count > 0) {
                    double avg_val = x->current_sum / x->current_count;
                    if (x->mix_mode) {
                        buf[x->current_pos].w_float += avg_val;
                    } else {
                        buf[x->current_pos].w_float = avg_val;
                    }
                }
                x->is_writing = 0;
                x->current_count = 0;
                x->current_sum = 0;
            }
            continue;
        }
        
        write_between_positions(buf, arraysize, x, idx, curr_sample);
    }
    
    garray_redraw(array);
    return (w+6);
}

static void tabsmear_tilde_set(t_tabsmear_tilde *x, t_symbol *s)
{
    x->x_arrayname = s;
}

static void tabsmear_tilde_loop(t_tabsmear_tilde *x, t_floatarg f)
{
    x->loop_enabled = (f != 0);
}

static void tabsmear_tilde_mix(t_tabsmear_tilde *x, t_floatarg f)
{
    x->mix_mode = (f != 0);
}

static void tabsmear_tilde_dsp(t_tabsmear_tilde *x, t_signal **sp)
{
    dsp_add(tabsmear_tilde_perform, 5, x,
        sp[0]->s_vec,          // audio input
        sp[1]->s_vec,          // index input
        sp[2]->s_vec,          // trigger input
        (t_int)sp[0]->s_length);
}

static void tabsmear_tilde_free(t_tabsmear_tilde *x)
{
    inlet_free(x->x_indexlet);
    inlet_free(x->x_trigglet);
}

void tabsmear_tilde_setup(void)
{
    tabsmear_tilde_class = class_new(gensym("tabsmear~"),
        (t_newmethod)tabsmear_tilde_new,
        (t_method)tabsmear_tilde_free,
        sizeof(t_tabsmear_tilde),
        CLASS_DEFAULT,
        A_GIMME, 0);
        
    class_addmethod(tabsmear_tilde_class, (t_method)tabsmear_tilde_dsp,
        gensym("dsp"), A_CANT, 0);
    class_addmethod(tabsmear_tilde_class, (t_method)tabsmear_tilde_set,
        gensym("set"), A_SYMBOL, 0);
    class_addmethod(tabsmear_tilde_class, (t_method)tabsmear_tilde_loop,
        gensym("loop"), A_FLOAT, 0);
    class_addmethod(tabsmear_tilde_class, (t_method)tabsmear_tilde_mix,
        gensym("mix"), A_FLOAT, 0);
    CLASS_MAINSIGNALIN(tabsmear_tilde_class, t_tabsmear_tilde, x_f);
}

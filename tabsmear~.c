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
    int add_mode;
    
    // For averaging at current position
    double current_sum;       // Sum of values at current position
    int current_count;        // Number of values at current position
    int current_pos;          // Current integer position
    
    // For interpolation between positions
    double last_written_val;  // Value written to last position
    int last_written_pos;     // Last position we wrote to
} t_tabsmear_tilde;

static inline double wrap_position(double pos, int bufsize, int loop_enabled) {
    if (!loop_enabled) {
        return fmin(fmax(pos, 0), bufsize-1);
    }
    double wrapped = fmod(pos, bufsize);
    if (wrapped < 0) wrapped += bufsize;
    return wrapped;
}

static inline void write_value(t_word *buf, int pos, double val, int add_mode, float trig) {
    if (add_mode) {
        buf[pos].w_float += (trig * val);
    } else {
        buf[pos].w_float = (1.0f - trig) * buf[pos].w_float + trig * val;
    }
}

static void write_between_positions(t_tabsmear_tilde *x, t_word *buf, int bufsize,
    double pos, double val, float trig)
{
    pos = wrap_position(pos, bufsize, x->loop_enabled);
    if (!x->loop_enabled && (pos >= bufsize || pos < 0)) return;

    int write_pos = (int)floor(pos);
    
    if (!x->is_writing) {
        x->current_pos = write_pos;
        x->current_sum = val;
        x->current_count = 1;
        x->is_writing = 1;
        x->last_written_pos = -1;
        x->last_written_val = 0;
        return;
    }
    
    if (write_pos != x->current_pos) {
        // Write accumulated average
        double avg_val = x->current_sum / x->current_count;
        write_value(buf, x->current_pos, avg_val, x->add_mode, trig);
        
        // Check if we probably wrapped around
        int direction;
        if (abs(write_pos - x->current_pos) > bufsize/2) {
            // We wrapped - don't interpolate
            direction = 0;
        } else {
            direction = (write_pos > x->current_pos) ? 1 : -1;
        }
        
        // Only interpolate if we didn't wrap
        if (direction != 0 && abs(write_pos - x->current_pos) > 1) {
            int start = x->current_pos;
            int end = write_pos;
            double val1 = avg_val;
            double val2 = val;
            
            // Handle interpolation between positions
            for (int i = start + direction; i != end; i += direction) {
                int idx = (int)wrap_position(i, bufsize, x->loop_enabled);
                double alpha = direction > 0 ?
                    (double)(i - start) / (end - start) :
                    (double)(start - i) / (start - end);
                
                // Interpolate both value and trigger
                float interp_trig = trig * ((direction > 0) ? alpha : (1.0 - alpha));
                write_value(buf, idx, 
                          val1 * (1.0 - alpha) + val2 * alpha,
                          x->add_mode, interp_trig);
            }
        }
        
        // Update state
        x->last_written_pos = x->current_pos;
        x->last_written_val = avg_val;
        x->current_pos = write_pos;
        x->current_sum = val;
        x->current_count = 1;
    } else {
        // Accumulate at current position
        x->current_sum += val;
        x->current_count++;
    }
}

static void *tabsmear_tilde_new(t_symbol *s, int argc, t_atom *argv)
{
    (void)s;
    t_tabsmear_tilde *x = (t_tabsmear_tilde *)pd_new(tabsmear_tilde_class);
    x->x_arrayname = (argc && argv->a_type == A_SYMBOL) ? 
        argv->a_w.w_symbol : gensym("array1");
    
    x->x_indexlet = inlet_new(&x->x_obj, &x->x_obj.ob_pd, &s_signal, &s_signal);
    x->x_trigglet = inlet_new(&x->x_obj, &x->x_obj.ob_pd, &s_signal, &s_signal);
    x->x_f = 0;
    
    // Initialize state
    x->is_writing = 0;
    x->loop_enabled = 1;
    x->add_mode = 0;
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
        
        if (trig > 0 && idx >= 0) {
            write_between_positions(x, buf, arraysize, idx, curr_sample, trig);
        } else if (x->is_writing && x->current_count) {
            // Writing stops - save final value and redraw
            double avg_val = x->current_sum / x->current_count;
            write_value(buf, x->current_pos, avg_val, x->add_mode, trig);
            x->is_writing = 0;
            x->current_count = 0;
            x->current_sum = 0;
            garray_redraw(array);
        }
    }
    
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

static void tabsmear_tilde_add(t_tabsmear_tilde *x, t_floatarg f)
{
    x->add_mode = (f != 0);
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
    class_addmethod(tabsmear_tilde_class, (t_method)tabsmear_tilde_add,
        gensym("add"), A_FLOAT, 0);
    CLASS_MAINSIGNALIN(tabsmear_tilde_class, t_tabsmear_tilde, x_f);
}

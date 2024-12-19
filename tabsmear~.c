#include "m_pd.h"
#include <math.h>

#define RATE_EPSILON 0.01

static t_class *tabsmear_tilde_class;

typedef struct _tabsmear_tilde {
    t_object x_obj;
    t_symbol *x_arrayname;
    t_float x_f;
    t_inlet *x_indexlet;
    t_inlet *x_trigglet;
    
    // History buffer for interpolation
    #define HISTORY_SIZE 4
    double hist_pos[HISTORY_SIZE];   // Position history
    double hist_val[HISTORY_SIZE];   // Value history
    int hist_count;                  // Number of points in history
    int hist_write;                  // Write position in circular history buffer
    
    long sample_count;
    int is_writing;
    int loop_enabled;
    int mix_mode;        // 0 = replace, 1 = add/accumulate
    
    // For accumulating slow writes
    double accum_val;      // Value we're accumulating
    int accum_pos;        // Position we're currently accumulating for
    int is_accumulating;  // Whether we're in the middle of accumulating
    double last_frac_pos; // Last fractional position we wrote from
} t_tabsmear_tilde;

static void write_from_history(t_word *buf, int bufsize, t_tabsmear_tilde *x,
    double start_pos, double end_pos)
{
    // Handle wrapping for looping
    if (x->loop_enabled) {
        while (start_pos >= bufsize) start_pos -= bufsize;
        while (start_pos < 0) start_pos += bufsize;
        while (end_pos >= bufsize) end_pos -= bufsize;
        while (end_pos < 0) end_pos += bufsize;
    } else {
        if (start_pos >= bufsize || end_pos < 0) return;
        start_pos = fmin(fmax(start_pos, 0), bufsize-1);
        end_pos = fmin(fmax(end_pos, 0), bufsize-1);
    }

    if (x->mix_mode) {
        // First determine if we have a direction
        if (!x->is_accumulating) {
            // Start tracking
            x->last_frac_pos = start_pos;
            x->is_accumulating = 1;
            x->accum_pos = (int)floor(start_pos);
        } else {
            // We have a previous position, check which integers we've crossed
            double prev_pos = x->last_frac_pos;
            double curr_pos = end_pos;
            int direction = (curr_pos > prev_pos) ? 1 : -1;
            
            // Find integer positions we've crossed
            int prev_int = (int)floor(prev_pos);
            int curr_int = (int)floor(curr_pos);
            
            // If we've crossed integer boundaries
            if (prev_int != curr_int) {
                // Write to crossed positions using interpolation
                if (direction > 0) {
                    for (int i = prev_int + 1; i <= curr_int; i++) {
                        // Only write when we've definitively passed this position
                        if (curr_pos > (double)(i + 0.5)) {
                            // Do interpolation through this position using history
                            double new_val = 0;
                            int count = 0;
                            
                            // Look through history for values that crossed this position
                            for (int h = 0; h < x->hist_count - 1; h++) {
                                int idx1 = (x->hist_write - h - 1 + HISTORY_SIZE) % HISTORY_SIZE;
                                int idx2 = (x->hist_write - h - 2 + HISTORY_SIZE) % HISTORY_SIZE;
                                double pos1 = x->hist_pos[idx1];
                                double pos2 = x->hist_pos[idx2];
                                
                                // If these history points crossed our target position
                                if ((pos1 <= i && pos2 > i) || (pos1 > i && pos2 <= i)) {
                                    double alpha = ((double)i - pos1) / (pos2 - pos1);
                                    new_val += x->hist_val[idx1] * (1.0 - alpha) + 
                                             x->hist_val[idx2] * alpha;
                                    count++;
                                }
                            }
                            
                            if (count > 0) {
                                buf[i].w_float += new_val / count;
                            }
                        }
                    }
                } else {
                    // Similar code for backward direction
                    for (int i = prev_int - 1; i >= curr_int; i--) {
                        if (curr_pos < (double)(i - 0.5)) {
                            // Same interpolation code as above
                        }
                    }
                }
            }
            
            x->last_frac_pos = curr_pos;
        }
    } else {
        // Non-mix mode - simple interpolation between points
        int i_start = (int)floor(start_pos);
        int i_end = (int)ceil(end_pos);
        
        for (int i = i_start; i <= i_end; i++) {
            if (i < 0 || i >= bufsize) continue;
            
            double pos = (double)i;
            int before = -1, after = -1;
            
            for (int h = 0; h < x->hist_count; h++) {
                int idx = (x->hist_write - h - 1 + HISTORY_SIZE) % HISTORY_SIZE;
                double hist_pos = x->hist_pos[idx];
                
                if (hist_pos <= pos && (before == -1 || hist_pos > x->hist_pos[before])) {
                    before = idx;
                }
                if (hist_pos > pos && (after == -1 || hist_pos < x->hist_pos[after])) {
                    after = idx;
                }
            }
            
            if (before >= 0 && after >= 0) {
                double pos_before = x->hist_pos[before];
                double pos_after = x->hist_pos[after];
                double val_before = x->hist_val[before];
                double val_after = x->hist_val[after];
                
                double alpha = (pos - pos_before) / (pos_after - pos_before);
                double new_val = val_before + alpha * (val_after - val_before);
                
                buf[i].w_float = new_val;
            }
        }
    }
}

static void add_to_history(t_tabsmear_tilde *x, double pos, double val) 
{
    x->hist_pos[x->hist_write] = pos;
    x->hist_val[x->hist_write] = val;
    x->hist_write = (x->hist_write + 1) % HISTORY_SIZE;
    if (x->hist_count < HISTORY_SIZE) x->hist_count++;
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
    x->sample_count = 0;
    x->is_writing = 0;
    x->loop_enabled = 1;
    x->hist_count = 0;
    x->hist_write = 0;
    x->mix_mode = 0;     // Start in replace mode
    x->accum_val = 0;
    x->accum_pos = 0;
    x->is_accumulating = 0;
    x->last_frac_pos = 0;
    
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
                x->is_writing = 0;
                x->sample_count = 0;
                x->hist_count = 0;  // Clear history when stopping
            }
            continue;
        }
        
        // Handle looping
        if (x->loop_enabled) {
            while (idx >= arraysize) idx -= arraysize;
            while (idx < 0) idx += arraysize;
        } else if (idx >= arraysize) {
            continue;
        }
        
        // If we just started writing, initialize state
        if (!x->is_writing) {
            x->is_writing = 1;
            x->hist_count = 0;
            x->hist_write = 0;
        }
        
        // Add current point to history
        add_to_history(x, idx, curr_sample);
        
        // Once we have enough history, start interpolating
        if (x->hist_count >= 2) {
            // Get previous position from one slot back in history
            int prev_idx = (x->hist_write - 2 + HISTORY_SIZE) % HISTORY_SIZE;
            double prev_pos = x->hist_pos[prev_idx];
            write_from_history(buf, arraysize, x, prev_pos, idx);
        }
        
        x->sample_count++;
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

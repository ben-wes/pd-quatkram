#include "m_pd.h"
#include <math.h>

static t_class *tabloop_tilde_class;

typedef struct _tabloop_tilde {
    t_object x_obj;
    t_symbol *x_arrayname;
    t_float x_f;
    t_inlet *x_startlet;
    t_inlet *x_endlet;
    t_word *x_vec;
    int x_npoints;
    float x_start;
    float x_end;
    int x_absolute_mode;
    
    // New fields for improved interpolation
    float x_last_pos;     // Track last position for rate calculation
    float x_last_output;  // For simple lowpass when needed
} t_tabloop_tilde;

static inline float wrap_position(float pos, float start, float end, int npoints, int absolute_mode)
{
    if (absolute_mode) {
        // Improved absolute mode wrapping
        float array_pos = fmod(pos, npoints);
        if (array_pos < 0) array_pos += npoints;
        
        if (start == end) return start;
        
        if (start <= end) {
            if (array_pos >= start && array_pos <= end) {
                return array_pos;
            }
            float range_size = end - start + 1;
            return start + fmod(array_pos - start, range_size);
        } else {
            if (array_pos >= start || array_pos <= end) {
                return array_pos;
            }
            float range_size = npoints - start + end + 1;
            float offset = array_pos > end ? array_pos - end - 1 : array_pos + (npoints - end - 1);
            return start + fmod(offset, range_size);
        }
    } else {
        float loop_length = (end >= start) ? (end - start + 1) : (npoints - start + end + 1);
        if (loop_length <= 0) loop_length = npoints;
        
        float cycle_pos = fmod(pos, loop_length);
        if (cycle_pos < 0) cycle_pos += loop_length;
        
        float result = start + cycle_pos;
        if (result >= npoints) result -= npoints;
        
        return result;
    }
}

static inline float interpolate_catmull_rom(t_word *buf, float pos, int npoints)
{
    int idx1 = (int)floor(pos);
    int idx0 = idx1 - 1;
    int idx2 = idx1 + 1;
    int idx3 = idx1 + 2;
    
    // Improved index wrapping
    idx0 = (idx0 + npoints) % npoints;
    idx1 = (idx1 + npoints) % npoints;
    idx2 = idx2 % npoints;
    idx3 = idx3 % npoints;
    
    float a = buf[idx0].w_float;
    float b = buf[idx1].w_float;
    float c = buf[idx2].w_float;
    float d = buf[idx3].w_float;
    
    float frac = pos - floor(pos);
    float tension = 0.5f;
    
    // Catmull-Rom interpolation with tension parameter
    float a1 = tension * (c - a);
    float a2 = tension * (d - b);
    float t2 = frac * frac;
    float t3 = t2 * frac;
    
    return b + 0.5f * (
        (a1 + a2) * frac + 
        (2.0f * a - 5.0f * b + 4.0f * c - d) * t2 + 
        (3.0f * b - 3.0f * c + d - a) * t3
    );
}

// Simple one-pole lowpass filter for anti-aliasing
static inline float lowpass_filter(float input, float last, float coeff)
{
    return input * coeff + last * (1.0f - coeff);
}

static t_int *tabloop_tilde_perform(t_int *w)
{
    t_tabloop_tilde *x = (t_tabloop_tilde *)(w[1]);
    t_sample *pos_in = (t_sample *)(w[2]);
    t_sample *start_in = (t_sample *)(w[3]);
    t_sample *end_in = (t_sample *)(w[4]);
    t_sample *out = (t_sample *)(w[5]);
    int n = (int)(w[6]);
    
    t_word *buf = x->x_vec;
    int npoints = x->x_npoints;
    float last_output = x->x_last_output;
    float last_pos = x->x_last_pos;
    
    if (!buf || npoints < 4) goto zero;
    
    while (n--) {
        float start = *start_in++;
        float end = *end_in++;
        float pos = *pos_in++;
        
        // Bound start
        start = fmax(0, fmin(start, npoints - 1));
        // Handle end position, allowing -1 for full array
        end = (end < 0) ? npoints - 1 : fmin(end, npoints - 1);
        
        // Calculate playback rate
        float rate = fabsf(pos - last_pos);
        last_pos = pos;
        
        // Get wrapped position
        float wrapped_pos = wrap_position(pos, start, end, npoints, x->x_absolute_mode);
        
        // Get interpolated value
        float output = interpolate_catmull_rom(buf, wrapped_pos, npoints);
        
        // Apply lowpass filtering for high playback rates
        if (rate > 1.0f) {
            // Adjust filter coefficient based on playback rate
            float filter_coeff = 1.0f / rate;
            output = lowpass_filter(output, last_output, filter_coeff);
        }
        
        *out++ = output;
        last_output = output;
    }
    
    x->x_last_pos = last_pos;
    x->x_last_output = last_output;
    return (w+7);
    
zero:
    while (n--) *out++ = 0;
    x->x_last_pos = 0;
    x->x_last_output = 0;
    return (w+7);
}

static void tabloop_tilde_set(t_tabloop_tilde *x, t_symbol *s)
{
    t_garray *a;
    
    x->x_arrayname = s;
    if (!(a = (t_garray *)pd_findbyclass(x->x_arrayname, garray_class))) {
        if (*s->s_name) pd_error(x, "tabloop~: %s: no such array", x->x_arrayname->s_name);
        x->x_vec = 0;
    }
    else if (!garray_getfloatwords(a, &x->x_npoints, &x->x_vec)) {
        pd_error(x, "%s: bad template for tabloop~", x->x_arrayname->s_name);
        x->x_vec = 0;
    }
    else {
        // Set default end if not set
        if (x->x_end < 0) x->x_end = x->x_npoints - 1;
        // Ensure start is valid
        if (x->x_start < 0) x->x_start = 0;
        if (x->x_start >= x->x_npoints) x->x_start = 0;
        // Ensure end is valid
        if (x->x_end >= x->x_npoints) x->x_end = x->x_npoints - 1;
        
        garray_usedindsp(a);
    }
}

static void tabloop_tilde_dsp(t_tabloop_tilde *x, t_signal **sp)
{
    tabloop_tilde_set(x, x->x_arrayname);
    
    dsp_add(tabloop_tilde_perform, 6, x,
        sp[0]->s_vec,      // read position
        sp[1]->s_vec,      // start
        sp[2]->s_vec,      // end
        sp[3]->s_vec,      // output
        sp[0]->s_n);       // block size
}

static void *tabloop_tilde_new(t_symbol *s, int argc, t_atom *argv)
{
    (void)s;
    t_tabloop_tilde *x = (t_tabloop_tilde *)pd_new(tabloop_tilde_class);
    
    // Check for -a flag first
    x->x_absolute_mode = 0;
    if (argc > 0 && argv[0].a_type == A_SYMBOL && 
        atom_getsymbol(argv) == gensym("-a")) {
        x->x_absolute_mode = 1;
        argc--; argv++; // Skip flag for remaining arg parsing
    }
     
    // First argument should be array name
    x->x_arrayname = (argc > 0 && argv[0].a_type == A_SYMBOL) ? 
        argv[0].a_w.w_symbol : gensym("array1");
    
    // Next arguments are start and end positions
    float start_pos = (argc > 1) ? atom_getfloatarg(1, argc, argv) : 0;
    float end_pos = (argc > 2) ? atom_getfloatarg(2, argc, argv) : -1;
    
    // Main signal inlet (read position) handled by CLASS_MAINSIGNALIN
    
    // Create signal inlets with default values
    x->x_startlet = signalinlet_new(&x->x_obj, start_pos);
    x->x_endlet = signalinlet_new(&x->x_obj, end_pos);
    
    outlet_new(&x->x_obj, &s_signal);
    
    x->x_f = 0;
    x->x_start = start_pos;
    x->x_end = end_pos;
    
    // Initialize array data
    x->x_vec = 0;
    x->x_npoints = 0;
    
    // Do initial array lookup
    tabloop_tilde_set(x, x->x_arrayname);
    
    return (x);
}

static void tabloop_tilde_abs(t_tabloop_tilde *x, t_floatarg f)
{
    x->x_absolute_mode = (f != 0);
}

static void tabloop_tilde_free(t_tabloop_tilde *x)
{
    inlet_free(x->x_startlet);
    inlet_free(x->x_endlet);
}

void tabloop_tilde_setup(void)
{
    tabloop_tilde_class = class_new(gensym("tabloop~"),
        (t_newmethod)tabloop_tilde_new,
        (t_method)tabloop_tilde_free,
        sizeof(t_tabloop_tilde),
        CLASS_DEFAULT,
        A_GIMME, 0);
    
    class_addmethod(tabloop_tilde_class, (t_method)tabloop_tilde_dsp,
        gensym("dsp"), A_CANT, 0);
    class_addmethod(tabloop_tilde_class, (t_method)tabloop_tilde_set,
        gensym("set"), A_SYMBOL, 0);
    class_addmethod(tabloop_tilde_class, (t_method)tabloop_tilde_abs,
        gensym("abs"), A_FLOAT, 0);
    CLASS_MAINSIGNALIN(tabloop_tilde_class, t_tabloop_tilde, x_f);
}

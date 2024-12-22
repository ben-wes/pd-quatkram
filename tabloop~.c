// interpolation copied from nusmuk-audio's tabread4c~

#include "m_pd.h"
#include <math.h>

static t_class *tabloop_tilde_class;

typedef struct _tabloop_tilde {
    t_object x_obj;
    t_symbol *x_arrayname;
    t_float x_f;  // dummy float for main signal inlet
    
    // Additional signal inlets
    t_inlet *x_startlet;
    t_inlet *x_endlet;
    
    // Array info
    t_word *x_vec;
    int x_npoints;
    
    // Loop parameters
    float x_start;
    float x_end;

    int x_absolute_mode;    // 0 = relative to start, 1 = absolute array position
} t_tabloop_tilde;

static inline float wrap_position(float pos, float start, float end, int npoints, int absolute_mode)
{
    if (absolute_mode) {
        // First wrap to array boundaries
        float array_pos = fmod(pos, npoints);
        if (array_pos < 0) array_pos += npoints;
        
        // Special case: start == end
        if (start == end) return start;
        
        // Check if position is in active range
        if (start <= end) {
            if (array_pos >= start && array_pos <= end) {
                return array_pos;
            }
            // Not in range - offset to start
            return start + fmod(array_pos - start, end - start + 1);
        } else {
            if (array_pos >= start || array_pos <= end) {
                return array_pos;
            }
            // Not in range - wrap around end of array
            float range_size = npoints - start + end + 1;
            float offset = array_pos > end ? array_pos - end : array_pos + (npoints - end);
            return start + fmod(offset, range_size);
        }
    } else {
        // Relative mode unchanged
        float loop_length = (end >= start) ? (end - start + 1) : (npoints - start + end + 1);
        if (loop_length <= 0) loop_length = npoints;
        
        float cycle = floor(pos / loop_length);
        float cycle_pos = pos - (cycle * loop_length);
        float result = start + cycle_pos;
        if (result >= npoints) result -= npoints;
        
        return result;
    }
}

static inline float interpolate_hermite(t_word *buf, float pos, int npoints)
{
    // Get integer indices for the 4 points
    int idx1 = (int)floor(pos);
    int idx0 = idx1 - 1;
    int idx2 = idx1 + 1;
    int idx3 = idx1 + 2;
    
    // Wrap all indices
    idx0 = (idx0 + npoints) % npoints;
    idx1 = (idx1 + npoints) % npoints;
    idx2 = idx2 % npoints;
    idx3 = idx3 % npoints;
    
    // Get sample values
    float a = buf[idx0].w_float;
    float b = buf[idx1].w_float;
    float c = buf[idx2].w_float;
    float d = buf[idx3].w_float;
    
    // Get fractional part
    float frac = pos - floor(pos);
    
    // 4-point, 3rd-order Hermite
    float a1 = 0.5f * (c - a);
    float a2 = a - 2.5f * b + 2.f * c - 0.5f * d;
    float a3 = 0.5f * (d - a) + 1.5f * (b - c);
    
    return ((a3 * frac + a2) * frac + a1) * frac + b;
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
    
    if (!buf || npoints < 4) goto zero;
    
    while (n--) {
        float start = *start_in++;
        float end = *end_in++;
        float pos = *pos_in++;
        
        // Bound start and end to array size
        if (start < 0) start = 0;
        if (start >= npoints) start = npoints - 1;
        if (end < 0) end = npoints - 1;
        if (end >= npoints) end = npoints - 1;
        
        // Get wrapped position
        float wrapped_pos = wrap_position(pos, start, end, npoints, x->x_absolute_mode);

        // Get interpolated value
        *out++ = interpolate_hermite(buf, wrapped_pos, npoints);
    }
    return (w+7);
    
zero:
    while (n--) *out++ = 0;
    return (w+7);
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
    x->x_start = (argc > 1) ? atom_getfloatarg(1, argc, argv) : 0;
    x->x_end = (argc > 2) ? atom_getfloatarg(2, argc, argv) : -1;
    
    // Create additional inlets for start and end positions
    x->x_startlet = inlet_new(&x->x_obj, &x->x_obj.ob_pd, &s_signal, &s_signal);
    x->x_endlet = inlet_new(&x->x_obj, &x->x_obj.ob_pd, &s_signal, &s_signal);
    
    outlet_new(&x->x_obj, &s_signal);
    
    x->x_f = 0;
    
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

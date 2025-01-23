/*
faccbounce~ - Float accumulator and bounce external for Pure Data
2024, Ben Wesch
Functionality:
- Accumulates incoming float values on a sample-by-sample basis
- Bounces the accumulated value between user-defined lower and upper bounds (default: -1 and 1)
- Maintains an internal float state that represents the accumulated value
Usage:
1. Send a signal to the inlet (float values to accumulate)
2. Receive the resulting signal from the outlet (accumulated and bounced value)
3. Use 'set <value>' message to reset the internal state to a specific value
4. Use 'range <lower> <upper>' message to set new lower and upper bounds
Creation arguments: [lower bound] [upper bound] (optional, default: -1 1)
Note: This code was developed with assistance from the Anthropic Claude AI language model.
*/

#include <string.h>
#include "m_pd.h"
#include <math.h>

static t_class *faccbounce_tilde_class;

typedef struct _faccbounce_tilde {
    t_object  x_obj;
    t_sample f_dummy;    // dummy float for signal inlet
    t_float *f_accum;    // Pointer to array of accumulated float values
    t_float *f_dir;      // Array of directions (1 or -1) for each channel
    t_float f_lower;     // Lower bound
    t_float f_upper;     // Upper bound
    t_int f_mode;        // Mode: 0 = wrap/bounce, 1 = reflect
    t_sample **in;       // Array of pointers to input channels
    t_sample **out;      // Array of pointers to output channels
    int n_channels;      // Number of channels
} t_faccbounce_tilde;

static float bounce_wrap(float value, float lower, float upper) {
    float range = upper - lower;
    value = fmodf(value - lower, 2 * range);
    if (value < 0) value += 2 * range;
    return (value > range) ? (2 * range - value) + lower : value + lower;
}

static float bounce_reflect(float value, float lower, float upper, float *direction) {
    // If we hit a boundary, reflect the direction
    if (value > upper) {
        *direction = -1.0f;
        value = upper - (value - upper);
    }
    else if (value < lower) {
        *direction = 1.0f;
        value = lower + (lower - value);
    }
    return value;
}

t_int *faccbounce_tilde_perform(t_int *w)
{
    t_faccbounce_tilde *x = (t_faccbounce_tilde *)(w[1]);
    int n = (int)(w[2]);
    int i, j;
    
    for (i = 0; i < x->n_channels; i++) {
        t_sample *in = x->in[i];
        t_sample *out = x->out[i];
        t_float accum = x->f_accum[i];
        t_float dir = x->f_dir[i];
        t_float lower = x->f_lower;
        t_float upper = x->f_upper;
        
        for (j = 0; j < n; j++) {
            if (x->f_mode == 0) {
                accum += in[j];
                accum = bounce_wrap(accum, lower, upper);
            } else {
                accum += in[j] * dir;
                accum = bounce_reflect(accum, lower, upper, &dir);
            }
            out[j] = accum;
        }
        
        x->f_accum[i] = accum;
        x->f_dir[i] = dir;
    }
    
    return (w + 3);
}

void faccbounce_tilde_dsp(t_faccbounce_tilde *x, t_signal **sp)
{
    int i;
    int n_channels = sp[0]->s_nchans;
    int vec_size = sp[0]->s_n;

    // Reallocate memory if number of channels has changed
    if (n_channels != x->n_channels) {
        x->in = (t_sample **)resizebytes(x->in, x->n_channels * sizeof(t_sample *), n_channels * sizeof(t_sample *));
        x->out = (t_sample **)resizebytes(x->out, x->n_channels * sizeof(t_sample *), n_channels * sizeof(t_sample *));
        x->f_accum = (t_float *)resizebytes(x->f_accum, x->n_channels * sizeof(t_float), n_channels * sizeof(t_float));
        x->f_dir = (t_float *)resizebytes(x->f_dir, x->n_channels * sizeof(t_float), n_channels * sizeof(t_float));

        // Initialize new channels
        for (i = x->n_channels; i < n_channels; i++) {
            x->f_accum[i] = 0;
            x->f_dir[i] = 1.0f;
        }

        x->n_channels = n_channels;
    }

    // Assign signal vectors
    for (i = 0; i < n_channels; i++) {
        x->in[i] = sp[0]->s_vec + vec_size * i;
    }

    // Set up output channels
    signal_setmultiout(&sp[1], n_channels);
    for (i = 0; i < n_channels; i++) {
        x->out[i] = sp[1]->s_vec + sp[1]->s_n * i;
    }

    dsp_add(faccbounce_tilde_perform, 2, x, vec_size);
}

void faccbounce_tilde_mode(t_faccbounce_tilde *x, t_floatarg f)
{
    x->f_mode = (f != 0);
}

void faccbounce_tilde_reset(t_faccbounce_tilde *x, t_floatarg f)
{
    // Apply the appropriate bounce mode to ensure the reset value is within bounds
    for (int i = 0; i < x->n_channels; i++) {
        if (x->f_mode == 0) {
            x->f_accum[i] = bounce_wrap(f, x->f_lower, x->f_upper);
        } else {
            x->f_accum[i] = bounce_reflect(f, x->f_lower, x->f_upper, &x->f_dir[i]);
        }
        x->f_dir[i] = 1.0f;
    }
}

void faccbounce_tilde_range(t_faccbounce_tilde *x, t_floatarg lower, t_floatarg upper)
{
    if (lower < upper) {
        x->f_lower = lower;
        x->f_upper = upper;
        // Adjust current values to new range
        for (int i = 0; i < x->n_channels; i++) {
            if (x->f_mode == 0) {
                x->f_accum[i] = bounce_wrap(x->f_accum[i], lower, upper);
            } else {
                x->f_accum[i] = bounce_reflect(x->f_accum[i], lower, upper, &x->f_dir[i]);
            }
        }
    } else {
        pd_error(x, "faccbounce~: lower bound must be less than upper bound");
    }
}

static void *faccbounce_tilde_new(t_symbol *s, int argc, t_atom *argv)
{
    t_faccbounce_tilde *x = (t_faccbounce_tilde *)pd_new(faccbounce_tilde_class);
    x->f_lower = -1, x->f_upper = 1;  // Default values
    x->f_mode = 0;  // Default to wrap/bounce mode
    
    (void)s;  // Unused parameter
    
    // Initialize with defaults
    x->n_channels = 0;
    x->in = (t_sample **)getbytes(sizeof(t_sample *));
    x->out = (t_sample **)getbytes(sizeof(t_sample *));
    x->f_accum = NULL;
    x->f_dir = NULL;
    
    // Parse creation arguments - flags first
    while (argc && argv->a_type == A_SYMBOL) {
        if (atom_getsymbol(argv) == gensym("-r")) x->f_mode = 1;
        else pd_error(x, "faccbounce~: invalid flag %s", atom_getsymbol(argv)->s_name);
        argc--, argv++;
    }
    
    // Then parse numerical arguments (bounds)
    if (argc && argv->a_type == A_FLOAT) {
        x->f_lower = atom_getfloat(argv);
        argc--, argv++;
        
        if (argc && argv->a_type == A_FLOAT) {
            x->f_upper = atom_getfloat(argv);
        }
        if (x->f_lower >= x->f_upper) {
            pd_error(x, "faccbounce~: lower bound must be less than upper bound, using defaults");
            x->f_lower = -1, x->f_upper = 1;
        }
    }    
    
    outlet_new(&x->x_obj, &s_signal);
    return (void *)x;
}

void faccbounce_tilde_free(t_faccbounce_tilde *x)
{
    if (x->in) freebytes(x->in, x->n_channels * sizeof(t_sample *));
    if (x->out) freebytes(x->out, x->n_channels * sizeof(t_sample *));
    if (x->f_accum) freebytes(x->f_accum, x->n_channels * sizeof(t_float));
    if (x->f_dir) freebytes(x->f_dir, x->n_channels * sizeof(t_float));
}

void faccbounce_tilde_setup(void)
{
    faccbounce_tilde_class = class_new(gensym("faccbounce~"),
        (t_newmethod)faccbounce_tilde_new,
        (t_method)faccbounce_tilde_free,
        sizeof(t_faccbounce_tilde),
        CLASS_DEFAULT | CLASS_MULTICHANNEL,
        A_GIMME, 0);
    
    class_addmethod(faccbounce_tilde_class, (t_method)faccbounce_tilde_dsp, gensym("dsp"), A_CANT, 0);
    class_addmethod(faccbounce_tilde_class, (t_method)faccbounce_tilde_reset, gensym("reset"), A_FLOAT, 0);
    class_addmethod(faccbounce_tilde_class, (t_method)faccbounce_tilde_range, gensym("range"), A_FLOAT, A_FLOAT, 0);
    class_addmethod(faccbounce_tilde_class, (t_method)faccbounce_tilde_mode, gensym("mode"), A_FLOAT, 0);
    CLASS_MAINSIGNALIN(faccbounce_tilde_class, t_faccbounce_tilde, f_dummy);
}

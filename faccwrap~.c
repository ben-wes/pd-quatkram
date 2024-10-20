/*
faccwrap~ - Float accumulator and wrapper external for Pure Data
2024, Ben Wesch
Functionality:
- Accumulates incoming float values on a sample-by-sample basis
- Wraps the accumulated value between user-defined lower and upper bounds (default: -1 and 1)
- Maintains an internal float state that represents the accumulated value
Usage:
1. Send a signal to the inlet (float values to accumulate)
2. Receive the resulting signal from the outlet (accumulated and wrapped value)
3. Use 'set <value>' message to reset the internal state to a specific value
4. Use 'range <lower> <upper>' message to set new lower and upper bounds
Creation arguments: [lower bound] [upper bound] (optional, default: -1 1)
Note: This code was developed with assistance from the Anthropic Claude AI language model.
*/
#include "m_pd.h"
#include <math.h>

static t_class *faccwrap_tilde_class;

typedef struct _faccwrap_tilde {
    t_object  x_obj;
    t_sample f_dummy;  // dummy float for signal inlet
    t_float *f_accum;   // Pointer to array of accumulated float values
    t_float f_lower;   // Lower bound
    t_float f_upper;   // Upper bound
    t_sample **in;     // Array of pointers to input channels
    t_sample **out;    // Array of pointers to output channels
    int n_channels;    // Number of channels
} t_faccwrap_tilde;

static float wrap(float value, float lower, float upper) {
    float range = upper - lower;
    float offset = value - lower;
    
    offset -= range * floorf(offset / range);
    
    return lower + offset;
}

t_int *faccwrap_tilde_perform(t_int *w)
{
    t_faccwrap_tilde *x = (t_faccwrap_tilde *)(w[1]);
    int n = (int)(w[2]);
    int i, j;
    
    for (i = 0; i < x->n_channels; i++) {
        t_sample *in = x->in[i];
        t_sample *out = x->out[i];
        t_float accum = x->f_accum[i];
        t_float lower = x->f_lower;
        t_float upper = x->f_upper;
        
        for (j = 0; j < n; j++) {
            accum += in[j];
            accum = wrap(accum, lower, upper);
            out[j] = accum;
        }
        
        x->f_accum[i] = accum;
    }
    
    return (w + 3);
}

void faccwrap_tilde_dsp(t_faccwrap_tilde *x, t_signal **sp)
{
    int i;
    int n_channels = sp[0]->s_nchans;
    int vec_size = sp[0]->s_n;

    // Reallocate memory if number of channels has changed
    if (n_channels != x->n_channels) {
        x->in = (t_sample **)resizebytes(x->in, x->n_channels * sizeof(t_sample *), n_channels * sizeof(t_sample *));
        x->out = (t_sample **)resizebytes(x->out, x->n_channels * sizeof(t_sample *), n_channels * sizeof(t_sample *));
        x->f_accum = (t_float *)resizebytes(x->f_accum, x->n_channels * sizeof(t_float), n_channels * sizeof(t_float));

        // Initialize new channels
        for (i = x->n_channels; i < n_channels; i++) {
            x->f_accum[i] = 0;
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

    dsp_add(faccwrap_tilde_perform, 2, x, vec_size);
}

void faccwrap_tilde_set(t_faccwrap_tilde *x, t_floatarg f)
{
    for (int i = 0; i < x->n_channels; i++) {
        x->f_accum[i] = wrap(f, x->f_lower, x->f_upper);
    }
}

void faccwrap_tilde_reset(t_faccwrap_tilde *x)
{
    for (int i = 0; i < x->n_channels; i++) {
        x->f_accum[i] = 0;
    }
}

void faccwrap_tilde_range(t_faccwrap_tilde *x, t_floatarg lower, t_floatarg upper)
{
    if (lower < upper) {
        x->f_lower = lower;
        x->f_upper = upper;
        for (int i = 0; i < x->n_channels; i++) {
            x->f_accum[i] = wrap(x->f_accum[i], lower, upper);
        }
    } else {
        pd_error(x, "faccwrap~: lower bound must be less than upper bound");
    }
}

void *faccwrap_tilde_new(t_floatarg lower, t_floatarg upper)
{
    t_faccwrap_tilde *x = (t_faccwrap_tilde *)pd_new(faccwrap_tilde_class);
    
    x->f_lower = (lower == 0 && upper == 0) ? -1 : lower;
    x->f_upper = (lower == 0 && upper == 0) ? 1 : upper;
    
    if (x->f_lower >= x->f_upper) {
        pd_error(x, "faccwrap~: lower bound must be less than upper bound, using defaults");
        x->f_lower = -1;
        x->f_upper = 1;
    }
    
    x->n_channels = 0;
    x->in = (t_sample **)getbytes(sizeof(t_sample *));
    x->out = (t_sample **)getbytes(sizeof(t_sample *));
    x->f_accum = NULL;
    
    outlet_new(&x->x_obj, &s_signal);
    return (void *)x;
}

void faccwrap_tilde_free(t_faccwrap_tilde *x)
{
    if (x->in) freebytes(x->in, x->n_channels * sizeof(t_sample *));
    if (x->out) freebytes(x->out, x->n_channels * sizeof(t_sample *));
    if (x->f_accum) freebytes(x->f_accum, x->n_channels * sizeof(t_float));
}

void faccwrap_tilde_setup(void)
{
    faccwrap_tilde_class = class_new(gensym("faccwrap~"),
        (t_newmethod)faccwrap_tilde_new,
        (t_method)faccwrap_tilde_free,
        sizeof(t_faccwrap_tilde),
        CLASS_DEFAULT | CLASS_MULTICHANNEL,
        A_DEFFLOAT, A_DEFFLOAT, 0);
    
    class_addmethod(faccwrap_tilde_class, (t_method)faccwrap_tilde_dsp, gensym("dsp"), A_CANT, 0);
    class_addmethod(faccwrap_tilde_class, (t_method)faccwrap_tilde_set, gensym("set"), A_FLOAT, 0);
    class_addmethod(faccwrap_tilde_class, (t_method)faccwrap_tilde_reset, gensym("reset"), 0);
    class_addmethod(faccwrap_tilde_class, (t_method)faccwrap_tilde_range, gensym("range"), A_FLOAT, A_FLOAT, 0);
    CLASS_MAINSIGNALIN(faccwrap_tilde_class, t_faccwrap_tilde, f_dummy);
}

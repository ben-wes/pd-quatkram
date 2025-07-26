/*
minhold~ - Minimum hold (valley hold) external for Pure Data
2024
- Outputs min(previous, current input) * coefficient
- Coefficient in 0..1 (default 1), signal-rate or arg
- State decays to 0 if coefficient < 1
- Multichannel, supports 'reset' message
*/

#include "m_pd.h"
#include <math.h>

static t_class *minhold_tilde_class;

typedef struct _minhold_tilde {
    t_object x_obj;
    t_sample f_dummy;
    t_float *state;
    t_sample **in;
    t_sample **coef;
    t_sample **out;
    int n_channels;
    int coef_nchans;
    t_float init_coef;
} t_minhold_tilde;

static t_int *minhold_tilde_perform(t_int *w)
{
    t_minhold_tilde *x = (t_minhold_tilde *)(w[1]);
    int n = (int)(w[2]);
    int i, j;
    int coef_nchans = x->coef_nchans;
    for (i = 0; i < x->n_channels; i++) {
        t_sample *in = x->in[i];
        t_sample *coef = x->coef[i % coef_nchans];
        t_sample *out = x->out[i];
        t_float state = x->state[i];
        for (j = 0; j < n; j++) {
            t_float c = coef[j];
            state = fminf(state, in[j]) * c;
            out[j] = state;
        }
        x->state[i] = state;
    }
    return (w + 3);
}

static void minhold_tilde_dsp(t_minhold_tilde *x, t_signal **sp)
{
    int i;
    int n_channels = sp[0]->s_nchans;
    int coef_nchans = sp[1]->s_nchans;
    int vec_size = sp[0]->s_n;
    if (n_channels != x->n_channels) {
        x->in = (t_sample **)resizebytes(x->in, x->n_channels * sizeof(t_sample *), n_channels * sizeof(t_sample *));
        x->coef = (t_sample **)resizebytes(x->coef, x->coef_nchans * sizeof(t_sample *), coef_nchans * sizeof(t_sample *));
        x->out = (t_sample **)resizebytes(x->out, x->n_channels * sizeof(t_sample *), n_channels * sizeof(t_sample *));
        x->state = (t_float *)resizebytes(x->state, x->n_channels * sizeof(t_float), n_channels * sizeof(t_float));
        for (i = x->n_channels; i < n_channels; i++) x->state[i] = 0.0f;
        x->n_channels = n_channels;
    }
    x->coef_nchans = coef_nchans;
    for (i = 0; i < n_channels; i++) x->in[i] = sp[0]->s_vec + vec_size * i;
    for (i = 0; i < coef_nchans; i++) x->coef[i] = sp[1]->s_vec + vec_size * i;
    if (!sp[1]->s_vec) {
        for (i = 0; i < n_channels; i++) {
            x->coef[i] = (t_sample *)getbytes(vec_size * sizeof(t_sample));
            if (x->coef[i]) for (int j = 0; j < vec_size; j++) x->coef[i][j] = x->init_coef;
        }
    }
    signal_setmultiout(&sp[2], n_channels);
    for (i = 0; i < n_channels; i++) x->out[i] = sp[2]->s_vec + sp[2]->s_n * i;
    dsp_add(minhold_tilde_perform, 2, x, vec_size);
}

static void minhold_tilde_reset(t_minhold_tilde *x, t_floatarg f)
{
    for (int i = 0; i < x->n_channels; i++) x->state[i] = f;
}

static void *minhold_tilde_new(t_symbol *s, int argc, t_atom *argv)
{
    t_minhold_tilde *x = (t_minhold_tilde *)pd_new(minhold_tilde_class);
    t_float coef = 1;
    if (argc > 0 && argv[0].a_type == A_FLOAT)
        coef = atom_getfloat(argv);
    x->init_coef = coef;
    x->n_channels = 0;
    x->coef_nchans = 1;
    x->in = (t_sample **)getbytes(sizeof(t_sample *));
    x->coef = (t_sample **)getbytes(sizeof(t_sample *));
    x->out = (t_sample **)getbytes(sizeof(t_sample *));
    x->state = (t_float *)getbytes(sizeof(t_float));
    if (x->state) x->state[0] = 0.0f;
    signalinlet_new(&x->x_obj, coef);
    outlet_new(&x->x_obj, &s_signal);
    return (void *)x;
}

static void minhold_tilde_free(t_minhold_tilde *x)
{
    if (x->in) freebytes(x->in, x->n_channels * sizeof(t_sample *));
    if (x->coef) freebytes(x->coef, x->n_channels * sizeof(t_sample *));
    if (x->out) freebytes(x->out, x->n_channels * sizeof(t_sample *));
    if (x->state) freebytes(x->state, x->n_channels * sizeof(t_float));
}

void minhold_tilde_setup(void)
{
    minhold_tilde_class = class_new(gensym("minhold~"),
        (t_newmethod)minhold_tilde_new,
        (t_method)minhold_tilde_free,
        sizeof(t_minhold_tilde),
        CLASS_MULTICHANNEL,
        A_GIMME, 0);
    class_addmethod(minhold_tilde_class, (t_method)minhold_tilde_dsp, gensym("dsp"), A_CANT, 0);
    class_addmethod(minhold_tilde_class, (t_method)minhold_tilde_reset, gensym("reset"), A_FLOAT, 0);
    CLASS_MAINSIGNALIN(minhold_tilde_class, t_minhold_tilde, f_dummy);
} 
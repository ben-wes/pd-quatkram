#include "m_pd.h"
#include <string.h>
#include <stdlib.h>

static t_class *mc_conv_tilde_class;

typedef struct _mc_conv_tilde {
    t_object x_obj;
    t_sample f;
    t_sample **in;
    t_sample **out;
    t_float *kernel;
    int kernel_size;
    int is_cyclic;     // 1 for cyclic, 0 for non-cyclic
    int is_normalized; // 1 for normalized, 0 for non-normalized
    int n_chans;
    int vec_size;
    t_sample *in_buffer; // Preallocated input buffer
} t_mc_conv_tilde;

static void mc_conv_tilde_normalize_kernel(t_mc_conv_tilde *x)
{
    if (x->is_normalized && x->kernel) {
        float sum = 0;
        int i;
        for (i = 0; i < x->kernel_size; i++) {
            sum += x->kernel[i];
        }
        if (sum != 0) {
            for (i = 0; i < x->kernel_size; i++) {
                x->kernel[i] /= sum;
            }
        }
    }
}

static t_int *mc_conv_tilde_perform(t_int *w)
{
    t_mc_conv_tilde *x = (t_mc_conv_tilde *)(w[1]);
    int n = (int)(w[2]);
    int i, j, k;
    t_sample sum;
    int half_kernel = x->kernel_size / 2;

    // If no kernel is set, just copy input to output
    if (!x->kernel || x->kernel_size == 0) {
        for (i = 0; i < x->n_chans; i++) {
            for (j = 0; j < n; j++) {
                x->out[i][j] = x->in[i][j];
            }
        }
        return (w+3);
    }

    for (i = 0; i < n; i++) {
        // Copy input to buffer
        for (j = 0; j < x->n_chans; j++) {
            x->in_buffer[j] = x->in[j][i];
        }

        for (j = 0; j < x->n_chans; j++) {
            sum = 0;
            for (k = 0; k < x->kernel_size; k++) {
                int idx = j - half_kernel + k;
                if (x->is_cyclic) {
                    idx = (idx + x->n_chans) % x->n_chans;
                    if (idx < 0) idx += x->n_chans;
                } else {
                    if (idx < 0) idx = 0;
                    if (idx >= x->n_chans) idx = x->n_chans - 1;
                }
                sum += x->in_buffer[idx] * x->kernel[k];
            }
            x->out[j][i] = sum;
        }
    }

    return (w+3);
}

static void mc_conv_tilde_dsp(t_mc_conv_tilde *x, t_signal **sp)
{
    int i;
    int vec_size = sp[0]->s_n;
    int n_chans = sp[0]->s_nchans;

    // Free previous memory if it exists
    if (x->in) freebytes(x->in, x->n_chans * sizeof(t_sample *));
    if (x->out) freebytes(x->out, x->n_chans * sizeof(t_sample *));
    if (x->in_buffer) freebytes(x->in_buffer, x->n_chans * sizeof(t_sample));
    
    x->n_chans = n_chans;
    x->vec_size = vec_size;

    // Allocate memory for input and output channel pointers
    x->in = (t_sample **)getbytes(n_chans * sizeof(t_sample *));
    x->out = (t_sample **)getbytes(n_chans * sizeof(t_sample *));
    x->in_buffer = (t_sample *)getbytes(n_chans * sizeof(t_sample));

    if (!x->in || !x->out || !x->in_buffer) {
        pd_error(x, "mc_conv~: out of memory");
        return;
    }

    // Assign input channels
    for (i = 0; i < n_chans; i++) {
        x->in[i] = sp[0]->s_vec + vec_size * i;
    }

    // Assign output channels
    signal_setmultiout(&sp[1], n_chans);
    for (i = 0; i < n_chans; i++) {
        x->out[i] = sp[1]->s_vec + sp[1]->s_n * i;
    }

    dsp_add(mc_conv_tilde_perform, 2, x, sp[0]->s_n);
}

static void mc_conv_tilde_set_kernel(t_mc_conv_tilde *x, t_symbol *s, int argc, t_atom *argv)
{
    (void)s; // Unused parameter

    int i;

    if (argc % 2 == 0) {
        pd_error(x, "mc_conv~: kernel must have an odd number of elements");
        return;
    }

    // Free the old kernel if it exists
    if (x->kernel) {
        freebytes(x->kernel, x->kernel_size * sizeof(t_float));
    }

    // Allocate memory for the new kernel
    x->kernel_size = argc;
    x->kernel = (t_float *)getbytes(x->kernel_size * sizeof(t_float));

    if (!x->kernel) {
        pd_error(x, "mc_conv~: out of memory");
        x->kernel_size = 0;
        return;
    }

    // Set the new kernel values
    for (i = 0; i < x->kernel_size; i++) {
        x->kernel[i] = atom_getfloat(argv + i);
    }

    // Normalize if required
    mc_conv_tilde_normalize_kernel(x);
}

static void mc_conv_tilde_set_cyclic(t_mc_conv_tilde *x, t_floatarg f)
{
    x->is_cyclic = (f != 0);
}

static void mc_conv_tilde_set_normalize(t_mc_conv_tilde *x, t_floatarg f)
{
    x->is_normalized = (f != 0);
    mc_conv_tilde_normalize_kernel(x);
}

static void *mc_conv_tilde_new(t_symbol *s, int argc, t_atom *argv)
{
    (void)s;

    t_mc_conv_tilde *x = (t_mc_conv_tilde *)pd_new(mc_conv_tilde_class);

    x->kernel = NULL;
    x->kernel_size = 0;
    x->is_cyclic = 0;
    x->is_normalized = 0;
    x->in = NULL;
    x->out = NULL;
    x->n_chans = 0;
    x->vec_size = 0;
    x->in_buffer = NULL;

    // Parse creation arguments
    while (argc && argv->a_type == A_SYMBOL) {
        if (atom_getsymbol(argv) == gensym("-c")) x->is_cyclic = 1;
        else if (atom_getsymbol(argv) == gensym("-n")) x->is_normalized = 1;
        else pd_error(x, "mc_conv~: invalid argument");
        argc--, argv++;
    }
    if (argc) mc_conv_tilde_set_kernel(x, &s_, argc, argv); // read kernel from args if present

    outlet_new(&x->x_obj, &s_signal);

    return (void *)x;
}

static void mc_conv_tilde_free(t_mc_conv_tilde *x)
{
    if (x->kernel) {
        freebytes(x->kernel, x->kernel_size * sizeof(t_float));
    }
    if (x->in) {
        freebytes(x->in, x->n_chans * sizeof(t_sample *));
    }
    if (x->out) {
        freebytes(x->out, x->n_chans * sizeof(t_sample *));
    }
    if (x->in_buffer) {
        freebytes(x->in_buffer, x->n_chans * sizeof(t_sample));
    }
}

void mc_conv_tilde_setup(void)
{
    mc_conv_tilde_class = class_new(gensym("mc_conv~"),
        (t_newmethod)mc_conv_tilde_new,
        (t_method)mc_conv_tilde_free,
        sizeof(t_mc_conv_tilde),
        CLASS_MULTICHANNEL,
        A_GIMME, 0);

    class_addmethod(mc_conv_tilde_class, (t_method)mc_conv_tilde_dsp, gensym("dsp"), A_CANT, 0);
    class_addmethod(mc_conv_tilde_class, (t_method)mc_conv_tilde_set_kernel, gensym("kernel"), A_GIMME, 0);
    class_addmethod(mc_conv_tilde_class, (t_method)mc_conv_tilde_set_cyclic, gensym("cyclic"), A_FLOAT, 0);
    class_addmethod(mc_conv_tilde_class, (t_method)mc_conv_tilde_set_normalize, gensym("normalize"), A_FLOAT, 0);
    CLASS_MAINSIGNALIN(mc_conv_tilde_class, t_mc_conv_tilde, f);
}
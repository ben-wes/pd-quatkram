#include "m_pd.h"
#include <string.h>
#include <math.h>
#include <stdlib.h>

static t_class *mc_conv_tilde_class;

typedef struct _mc_conv_tilde {
    t_object x_obj;
    t_sample f;
    t_sample **in;
    t_sample **out;
    t_sample **kernel_in;
    t_float *msg_kernel;  // Message-defined kernel
    int msg_kernel_size;  // Size of message-defined kernel
    t_float *active_kernel;  // Pointer to the currently active kernel
    int kernel_size;  // Size of the active kernel
    int input_chans;
    int is_cyclic;     // 1 for cyclic, 0 for non-cyclic
    int is_normalized; // 1 for normalized, 0 for non-normalized
    int vec_size;
    t_sample *in_buffer; // Preallocated input buffer
    int use_signal_kernel; // Flag to indicate if signal kernel is being used
} t_mc_conv_tilde;

// Normalize the kernel
static void mc_conv_tilde_normalize_kernel(t_float *kernel, int size)
{
    float sum = 0;
    int i;
    for (i = 0; i < size; i++) {
        sum += fabsf(kernel[i]);
    }
    if (sum != 0) {
        for (i = 0; i < size; i++) {
            kernel[i] /= sum;
        }
    }
}

// Perform the convolution
static t_int *mc_conv_tilde_perform(t_int *w)
{
    t_mc_conv_tilde *x = (t_mc_conv_tilde *)(w[1]);
    int n = (int)(w[2]);
    int i, j, k;
    t_sample sum;
    int half_kernel = x->kernel_size / 2;
    t_float normalized_kernel[x->kernel_size];

    for (i = 0; i < n; i++) {
        // Copy input to buffer
        for (j = 0; j < x->input_chans; j++) {
            x->in_buffer[j] = x->in[j][i];
        }

        // Use signal kernel if connected
        if (x->use_signal_kernel) {
            for (k = 0; k < x->kernel_size; k++) {
                x->active_kernel[k] = x->kernel_in[k][i];
            }
        }

        // Normalize kernel if required
        t_float *current_kernel = x->active_kernel;
        if (x->is_normalized) {
            memcpy(normalized_kernel, x->active_kernel, x->kernel_size * sizeof(t_float));
            mc_conv_tilde_normalize_kernel(normalized_kernel, x->kernel_size);
            current_kernel = normalized_kernel;
        }

        // Perform convolution
        for (j = 0; j < x->input_chans; j++) {
            sum = 0;
            for (k = 0; k < x->kernel_size; k++) {
                int idx = j - half_kernel + k;
                if (x->is_cyclic) {
                    idx = (idx + x->input_chans) % x->input_chans;
                    if (idx < 0) idx += x->input_chans;
                } else {
                    if (idx < 0) idx = 0;
                    if (idx >= x->input_chans) idx = x->input_chans - 1;
                }
                sum += x->in_buffer[idx] * current_kernel[k];
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
    int input_chans = sp[0]->s_nchans;
    int kernel_chans = sp[1]->s_nchans;

    // Free previous memory
    if (x->in) freebytes(x->in, x->input_chans * sizeof(t_sample *));
    if (x->out) freebytes(x->out, x->input_chans * sizeof(t_sample *));
    if (x->in_buffer) freebytes(x->in_buffer, x->input_chans * sizeof(t_sample));
    if (x->kernel_in) freebytes(x->kernel_in, x->kernel_size * sizeof(t_sample *));
    
    x->input_chans = input_chans;
    x->vec_size = vec_size;

    // Determine if we're using a signal kernel (more than 1 channel or non-zero single channel)
    x->use_signal_kernel = (kernel_chans > 1) || (kernel_chans == 1 && sp[1]->s_vec[0] != 0);
    
    if (x->use_signal_kernel) {
        x->kernel_size = kernel_chans;
        x->active_kernel = (t_float *)getbytes(x->kernel_size * sizeof(t_float));
    } else {
        x->kernel_size = x->msg_kernel_size;
        x->active_kernel = x->msg_kernel;
    }

    // Allocate memory
    x->in = (t_sample **)getbytes(input_chans * sizeof(t_sample *));
    x->out = (t_sample **)getbytes(input_chans * sizeof(t_sample *));
    x->in_buffer = (t_sample *)getbytes(input_chans * sizeof(t_sample));
    x->kernel_in = (t_sample **)getbytes(x->kernel_size * sizeof(t_sample *));

    if (!x->in || !x->out || !x->in_buffer || !x->kernel_in || (x->use_signal_kernel && !x->active_kernel)) {
        pd_error(x, "mc_conv~: out of memory");
        return;
    }

    // Assign input channels
    for (i = 0; i < input_chans; i++) {
        x->in[i] = sp[0]->s_vec + vec_size * i;
    }

    // Assign kernel input channels
    for (i = 0; i < x->kernel_size; i++) {
        x->kernel_in[i] = (i < kernel_chans) ? sp[1]->s_vec + vec_size * i : sp[1]->s_vec;
    }

    // Assign output channels
    signal_setmultiout(&sp[2], input_chans);
    for (i = 0; i < input_chans; i++) {
        x->out[i] = sp[2]->s_vec + sp[2]->s_n * i;
    }

    dsp_add(mc_conv_tilde_perform, 2, x, sp[0]->s_n);
}

// Set the kernel via message
static void mc_conv_tilde_set_kernel(t_mc_conv_tilde *x, t_symbol *s, int argc, t_atom *argv)
{
    (void)s; // Unused parameter

    if (argc == 0) {
        pd_error(x, "mc_conv~: kernel must have at least one element");
        return;
    }

    // Ensure odd kernel size
    int new_kernel_size = argc;
    if (new_kernel_size % 2 == 0) {
        new_kernel_size++;
        post("mc_conv~: kernel size must be odd, appending 0");
    }

    // Reallocate msg_kernel if size changed
    if (new_kernel_size != x->msg_kernel_size) {
        x->msg_kernel = (t_float *)resizebytes(x->msg_kernel, x->msg_kernel_size * sizeof(t_float), new_kernel_size * sizeof(t_float));
        if (!x->msg_kernel) {
            pd_error(x, "mc_conv~: out of memory");
            x->msg_kernel_size = 0;
            return;
        }
        x->msg_kernel_size = new_kernel_size;
    }

    // Set the new kernel values
    for (int i = 0; i < argc; i++) {
        x->msg_kernel[i] = atom_getfloat(argv + i);
    }
    if (new_kernel_size > argc) {
        x->msg_kernel[argc] = 0; // Set the added element to 0 if we padded
    }

    // If not using signal kernel, update active kernel
    if (!x->use_signal_kernel) {
        x->active_kernel = x->msg_kernel;
        x->kernel_size = x->msg_kernel_size;
    }
}

static void mc_conv_tilde_set_cyclic(t_mc_conv_tilde *x, t_floatarg f)
{
    x->is_cyclic = (f != 0);
}

static void mc_conv_tilde_set_normalize(t_mc_conv_tilde *x, t_floatarg f)
{
    x->is_normalized = (f != 0);
}

static void *mc_conv_tilde_new(t_symbol *s, int argc, t_atom *argv)
{
    (void)s;

    t_mc_conv_tilde *x = (t_mc_conv_tilde *)pd_new(mc_conv_tilde_class);

    x->msg_kernel = NULL;
    x->msg_kernel_size = 0;
    x->active_kernel = NULL;
    x->kernel_size = 0;
    x->is_cyclic = 0;
    x->is_normalized = 0;
    x->in = NULL;
    x->out = NULL;
    x->kernel_in = NULL;
    x->input_chans = 0;
    x->vec_size = 0;
    x->in_buffer = NULL;
    x->use_signal_kernel = 0;

    // Parse creation arguments
    while (argc && argv->a_type == A_SYMBOL) {
        if (atom_getsymbol(argv) == gensym("-c")) x->is_cyclic = 1;
        else if (atom_getsymbol(argv) == gensym("-n")) x->is_normalized = 1;
        else pd_error(x, "mc_conv~: invalid flags");
        argc--, argv++;
    }
    if (argc) {
        mc_conv_tilde_set_kernel(x, &s_, argc, argv); // read kernel from args if present
    } else {
        // If no kernel is provided, initialize with a default kernel
        t_atom default_kernel[1];
        SETFLOAT(default_kernel, 1.0f);
        mc_conv_tilde_set_kernel(x, &s_, 1, default_kernel);
    }

    inlet_new(&x->x_obj, &x->x_obj.ob_pd, &s_signal, &s_signal); // Add second inlet for kernel
    outlet_new(&x->x_obj, &s_signal);

    return (void *)x;
}

static void mc_conv_tilde_free(t_mc_conv_tilde *x)
{
    if (x->in) freebytes(x->in, x->input_chans * sizeof(t_sample *));
    if (x->out) freebytes(x->out, x->input_chans * sizeof(t_sample *));
    if (x->in_buffer) freebytes(x->in_buffer, x->input_chans * sizeof(t_sample));
    if (x->msg_kernel) freebytes(x->msg_kernel, x->msg_kernel_size * sizeof(t_float));
    if (x->kernel_in) freebytes(x->kernel_in, x->kernel_size * sizeof(t_sample *));
    if (x->use_signal_kernel && x->active_kernel) freebytes(x->active_kernel, x->kernel_size * sizeof(t_float));
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
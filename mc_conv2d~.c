#include "m_pd.h"
#include <string.h>
#include <math.h>
#include <stdlib.h>

static t_class *mc_conv2d_tilde_class;

typedef struct _mc_conv2d_tilde {
    t_object x_obj;
    t_sample f;
    t_sample **in;
    t_sample **out;
    t_sample **kernel_in;
    t_float *msg_kernel;  // Message-defined kernel
    int msg_kernel_size;  // Size of message-defined kernel (should be square)
    t_float *active_kernel;  // Pointer to the currently active kernel
    int kernel_size;  // Size of the active kernel (should be square)
    int input_size;  // Size of the input (should be square)
    int input_dim;   // Dimension of the input (sqrt of input_size)
    int is_cyclic;     // 1 for cyclic, 0 for non-cyclic
    int is_normalized; // 1 for normalized, 0 for non-normalized
    int vec_size;
    t_sample *in_buffer; // Preallocated input buffer
    int use_signal_kernel; // Flag to indicate if signal kernel is being used
} t_mc_conv2d_tilde;

// Normalize the kernel
static void mc_conv2d_tilde_normalize_kernel(t_float *kernel, int size)
{
    float sum = 0;
    int i;
    for (i = 0; i < size * size; i++) {
        sum += fabsf(kernel[i]);
    }
    if (sum != 0) {
        for (i = 0; i < size * size; i++) {
            kernel[i] /= sum;
        }
    }
}

// Perform the 2D convolution
static t_int *mc_conv2d_tilde_perform(t_int *w)
{
    t_mc_conv2d_tilde *x = (t_mc_conv2d_tilde *)(w[1]);
    int n = (int)(w[2]);
    int i, j, k, m, row, col;
    t_sample sum;
    int half_kernel = x->kernel_size / 2;
    t_float normalized_kernel[x->kernel_size * x->kernel_size];

    for (i = 0; i < n; i++) {
        // Copy input to buffer
        for (j = 0; j < x->input_size; j++) {
            x->in_buffer[j] = x->in[j][i];
        }

        // Use signal kernel if connected
        if (x->use_signal_kernel) {
            for (k = 0; k < x->kernel_size * x->kernel_size; k++) {
                x->active_kernel[k] = x->kernel_in[k][i];
            }
        }

        // Normalize kernel if required
        t_float *current_kernel = x->active_kernel;
        if (x->is_normalized) {
            memcpy(normalized_kernel, x->active_kernel, x->kernel_size * x->kernel_size * sizeof(t_float));
            mc_conv2d_tilde_normalize_kernel(normalized_kernel, x->kernel_size);
            current_kernel = normalized_kernel;
        }

        // Perform 2D convolution
        for (j = 0; j < x->input_size; j++) {
            row = j / x->input_dim;
            col = j % x->input_dim;
            sum = 0;
            for (k = 0; k < x->kernel_size; k++) {
                for (m = 0; m < x->kernel_size; m++) {
                    int input_row, input_col;
                    
                    if (x->is_cyclic) {
                        input_row = (row - half_kernel + k + x->input_dim) % x->input_dim;
                        input_col = (col - half_kernel + m + x->input_dim) % x->input_dim;
                    } else {
                        input_row = row - half_kernel + k;
                        input_col = col - half_kernel + m;
                        // Clamp to edge for non-cyclic mode
                        if (input_row < 0) input_row = 0;
                        if (input_row >= x->input_dim) input_row = x->input_dim - 1;
                        if (input_col < 0) input_col = 0;
                        if (input_col >= x->input_dim) input_col = x->input_dim - 1;
                    }
                    
                    int idx = input_row * x->input_dim + input_col;
                    sum += x->in_buffer[idx] * current_kernel[k * x->kernel_size + m];
                }
            }
            x->out[j][i] = sum;
        }
    }

    return (w+3);
}

static void mc_conv2d_tilde_dsp(t_mc_conv2d_tilde *x, t_signal **sp)
{
    int i;
    int vec_size = sp[0]->s_n;
    int input_size = sp[0]->s_nchans;
    int kernel_size = sp[1]->s_nchans;

    // Check if input size is a perfect square
    int input_dim = (int)sqrt(input_size);
    if (input_dim * input_dim != input_size) {
        pd_error(x, "mc_conv2d~: input size must be a perfect square");
        return;
    }

    // Free previous memory
    if (x->in) freebytes(x->in, x->input_size * sizeof(t_sample *));
    if (x->out) freebytes(x->out, x->input_size * sizeof(t_sample *));
    if (x->in_buffer) freebytes(x->in_buffer, x->input_size * sizeof(t_sample));
    if (x->kernel_in) freebytes(x->kernel_in, x->kernel_size * x->kernel_size * sizeof(t_sample *));
    
    x->input_size = input_size;
    x->input_dim = input_dim;
    x->vec_size = vec_size;

    // Determine if we're using a signal kernel
    x->use_signal_kernel = (kernel_size > 1) || (kernel_size == 1 && sp[1]->s_vec[0] != 0);
    
    if (x->use_signal_kernel) {
        x->kernel_size = (int)sqrt(kernel_size);
        if (x->kernel_size * x->kernel_size != kernel_size) {
            pd_error(x, "mc_conv2d~: kernel size must be a perfect square");
            return;
        }
        x->active_kernel = (t_float *)getbytes(x->kernel_size * x->kernel_size * sizeof(t_float));
        post("mc_conv2d~: Using signal kernel mode. Kernel size: %dx%d", x->kernel_size, x->kernel_size);
    } else {
        x->kernel_size = (int)sqrt(x->msg_kernel_size);
        x->active_kernel = x->msg_kernel;
        post("mc_conv2d~: Using message kernel mode. Kernel size: %dx%d", x->kernel_size, x->kernel_size);
    }

    // Allocate memory
    x->in = (t_sample **)getbytes(input_size * sizeof(t_sample *));
    x->out = (t_sample **)getbytes(input_size * sizeof(t_sample *));
    x->in_buffer = (t_sample *)getbytes(input_size * sizeof(t_sample));
    x->kernel_in = (t_sample **)getbytes(x->kernel_size * x->kernel_size * sizeof(t_sample *));

    if (!x->in || !x->out || !x->in_buffer || !x->kernel_in || (x->use_signal_kernel && !x->active_kernel)) {
        pd_error(x, "mc_conv2d~: out of memory");
        return;
    }

    // Assign input channels
    for (i = 0; i < input_size; i++) {
        x->in[i] = sp[0]->s_vec + vec_size * i;
    }

    // Assign kernel input channels
    for (i = 0; i < x->kernel_size * x->kernel_size; i++) {
        x->kernel_in[i] = (i < kernel_size) ? sp[1]->s_vec + vec_size * i : sp[1]->s_vec;
    }

    // Assign output channels
    signal_setmultiout(&sp[2], input_size);
    for (i = 0; i < input_size; i++) {
        x->out[i] = sp[2]->s_vec + sp[2]->s_n * i;
    }

    post("mc_conv2d~: Input size: %dx%d", x->input_dim, x->input_dim);
    post("mc_conv2d~: Convolution mode: %s", x->is_cyclic ? "Cyclic" : "Non-cyclic");
    post("mc_conv2d~: Normalization: %s", x->is_normalized ? "On" : "Off");

    dsp_add(mc_conv2d_tilde_perform, 2, x, sp[0]->s_n);
}

static void mc_conv2d_tilde_set_kernel(t_mc_conv2d_tilde *x, t_symbol *s, int argc, t_atom *argv)
{
    (void)s; // Unused parameter

    if (argc == 0) {
        pd_error(x, "mc_conv2d~: kernel must have at least one element");
        return;
    }

    // Check if kernel size is a perfect square
    int new_kernel_dim = (int)sqrt(argc);
    if (new_kernel_dim * new_kernel_dim != argc) {
        pd_error(x, "mc_conv2d~: kernel size must be a perfect square");
        return;
    }

    // Reallocate msg_kernel if size changed
    if (argc != x->msg_kernel_size) {
        x->msg_kernel = (t_float *)resizebytes(x->msg_kernel, x->msg_kernel_size * sizeof(t_float), argc * sizeof(t_float));
        if (!x->msg_kernel) {
            pd_error(x, "mc_conv2d~: out of memory");
            x->msg_kernel_size = 0;
            return;
        }
        x->msg_kernel_size = argc;
    }

    // Set the new kernel values
    for (int i = 0; i < argc; i++) {
        x->msg_kernel[i] = atom_getfloat(argv + i);
    }

    // If not using signal kernel, update active kernel
    if (!x->use_signal_kernel) {
        x->active_kernel = x->msg_kernel;
        x->kernel_size = new_kernel_dim;
    }

    // Post a message about the kernel size and center handling
    if (x->kernel_size % 2 == 0) {
        post("mc_conv2d~: Even kernel size (%dx%d). Using top-left center approach.", x->kernel_size, x->kernel_size);
    } else {
        post("mc_conv2d~: Odd kernel size (%dx%d). Using standard center approach.", x->kernel_size, x->kernel_size);
    }
}

static void mc_conv2d_tilde_set_cyclic(t_mc_conv2d_tilde *x, t_floatarg f)
{
    x->is_cyclic = (f != 0);
}

static void mc_conv2d_tilde_set_normalize(t_mc_conv2d_tilde *x, t_floatarg f)
{
    x->is_normalized = (f != 0);
}

static void *mc_conv2d_tilde_new(t_symbol *s, int argc, t_atom *argv)
{
    (void)s;

    t_mc_conv2d_tilde *x = (t_mc_conv2d_tilde *)pd_new(mc_conv2d_tilde_class);

    x->msg_kernel = NULL;
    x->msg_kernel_size = 0;
    x->active_kernel = NULL;
    x->kernel_size = 0;
    x->is_cyclic = 0;
    x->is_normalized = 0;
    x->in = NULL;
    x->out = NULL;
    x->kernel_in = NULL;
    x->input_size = 0;
    x->input_dim = 0;
    x->vec_size = 0;
    x->in_buffer = NULL;
    x->use_signal_kernel = 0;

    // Parse creation arguments
    while (argc && argv->a_type == A_SYMBOL) {
        if (atom_getsymbol(argv) == gensym("-c")) x->is_cyclic = 1;
        else if (atom_getsymbol(argv) == gensym("-n")) x->is_normalized = 1;
        else pd_error(x, "mc_conv2d~: invalid flags");
        argc--, argv++;
    }
    if (argc) {
        mc_conv2d_tilde_set_kernel(x, &s_, argc, argv); // read kernel from args if present
    } else {
        // If no kernel is provided, initialize with a default kernel
        t_atom default_kernel[1];
        SETFLOAT(default_kernel, 1.0f);
        mc_conv2d_tilde_set_kernel(x, &s_, 1, default_kernel);
    }

    inlet_new(&x->x_obj, &x->x_obj.ob_pd, &s_signal, &s_signal); // Add second inlet for kernel
    outlet_new(&x->x_obj, &s_signal);

    return (void *)x;
}

static void mc_conv2d_tilde_free(t_mc_conv2d_tilde *x)
{
    if (x->in) freebytes(x->in, x->input_size * sizeof(t_sample *));
    if (x->out) freebytes(x->out, x->input_size * sizeof(t_sample *));
    if (x->in_buffer) freebytes(x->in_buffer, x->input_size * sizeof(t_sample));
    if (x->msg_kernel) freebytes(x->msg_kernel, x->msg_kernel_size * sizeof(t_float));
    if (x->kernel_in) freebytes(x->kernel_in, x->kernel_size * x->kernel_size * sizeof(t_sample *));
    if (x->use_signal_kernel && x->active_kernel) freebytes(x->active_kernel, x->kernel_size * x->kernel_size * sizeof(t_float));
}

void mc_conv2d_tilde_setup(void)
{
    mc_conv2d_tilde_class = class_new(gensym("mc_conv2d~"),
        (t_newmethod)mc_conv2d_tilde_new,
        (t_method)mc_conv2d_tilde_free,
        sizeof(t_mc_conv2d_tilde),
        CLASS_MULTICHANNEL,
        A_GIMME, 0);

    class_addmethod(mc_conv2d_tilde_class, (t_method)mc_conv2d_tilde_dsp, gensym("dsp"), A_CANT, 0);
    class_addmethod(mc_conv2d_tilde_class, (t_method)mc_conv2d_tilde_set_kernel, gensym("kernel"), A_GIMME, 0);
    class_addmethod(mc_conv2d_tilde_class, (t_method)mc_conv2d_tilde_set_cyclic, gensym("cyclic"), A_FLOAT, 0);
    class_addmethod(mc_conv2d_tilde_class, (t_method)mc_conv2d_tilde_set_normalize, gensym("normalize"), A_FLOAT, 0);
    CLASS_MAINSIGNALIN(mc_conv2d_tilde_class, t_mc_conv2d_tilde, f);
}
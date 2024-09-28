/*
qdiv~ - Quaternion division external for Pure Data

2024, Adapted from qmul~ by Ben Wesch

Functionality:
- Performs quaternion division on a sample-by-sample basis
- Two 4-channel signal inlets for input quaternions
- One 4-channel signal outlet for the resulting quaternion
- Allows initialization of the right fallback quaternion with 4 float arguments (default: 1 0 0 0)
- Uses fallback quaternion (1 0 0 0) if less than 4 channels are connected to the left inlet
- Uses custom fallback quaternion if less than 4 channels are connected to the right inlet

Usage:
1. Create the object with [qdiv~] or [qdiv~ w x y z] to set right fallback quaternion
2. Send two 4-channel signals to the inlets (w, x, y, z components of quaternions)
3. Receive the resulting 4-channel signal from the outlet

Note: This code was adapted with assistance from the Anthropic Claude AI language model.
*/

#include "m_pd.h"
#include <string.h>
#include <math.h>

static t_class *qdiv_tilde_class;

typedef struct _qdiv_tilde {
    t_object  x_obj;
    t_sample f;
    t_sample *quat_in1[4];  // Left quaternion input (4 channels)
    t_sample *quat_in2[4];  // Right quaternion input (4 channels)
    t_sample *quat_out[4];  // Output quaternion (4 channels)
    t_float fallback_quat[4];  // Right fallback quaternion (default: 1 0 0 0)
    t_sample *zero_buffer;  // Dynamic zero buffer
    int buffer_size;        // Size of the zero buffer
    int left_channel_count;  // Number of channels connected to the left inlet
    int right_channel_count; // Number of channels connected to the right inlet
} t_qdiv_tilde;

t_int *qdiv_tilde_perform(t_int *w)
{
    t_qdiv_tilde *x = (t_qdiv_tilde *)(w[1]);
    int n = (int)(w[2]);

    t_sample *in1w = x->quat_in1[0], *in1x = x->quat_in1[1], *in1y = x->quat_in1[2], *in1z = x->quat_in1[3];
    t_sample *in2w = x->quat_in2[0], *in2x = x->quat_in2[1], *in2y = x->quat_in2[2], *in2z = x->quat_in2[3];
    t_sample *ow = x->quat_out[0], *ox = x->quat_out[1], *oy = x->quat_out[2], *oz = x->quat_out[3];

    while (n--) {
        t_float q1w, q1x, q1y, q1z, q2w, q2x, q2y, q2z;

        // Use actual values or fallback to 1 0 0 0 for left input
        if (x->left_channel_count == 4) {
            q1w = *in1w++;
            q1x = *in1x++;
            q1y = *in1y++;
            q1z = *in1z++;
        } else {
            q1w = 1;
            q1x = q1y = q1z = 0;
        }

        // Use fallback quaternion if less than 4 channels are connected to the right inlet
        if (x->right_channel_count == 4) {
            q2w = *in2w++;
            q2x = *in2x++;
            q2y = *in2y++;
            q2z = *in2z++;
        } else {
            q2w = x->fallback_quat[0];
            q2x = x->fallback_quat[1];
            q2y = x->fallback_quat[2];
            q2z = x->fallback_quat[3];
        }

        // Calculate the magnitude of the second quaternion
        t_float mag = q2w*q2w + q2x*q2x + q2y*q2y + q2z*q2z;

        // Check if magnitude is close to zero to avoid division by zero
        if (mag < 1e-6) {
            *ow++ = *ox++ = *oy++ = *oz++ = 0;
        } else {
            // Calculate the conjugate of the second quaternion
            t_float conj_w = q2w / mag;
            t_float conj_x = -q2x / mag;
            t_float conj_y = -q2y / mag;
            t_float conj_z = -q2z / mag;

            // Perform quaternion division (multiplication by inverse)
            *ow++ = q1w * conj_w - q1x * conj_x - q1y * conj_y - q1z * conj_z;
            *ox++ = q1w * conj_x + q1x * conj_w + q1y * conj_z - q1z * conj_y;
            *oy++ = q1w * conj_y - q1x * conj_z + q1y * conj_w + q1z * conj_x;
            *oz++ = q1w * conj_z + q1x * conj_y - q1y * conj_x + q1z * conj_w;
        }
    }

    return (w+3);
}

void qdiv_tilde_dsp(t_qdiv_tilde *x, t_signal **sp)
{
    int vec_size = sp[0]->s_n;

    // Reallocate zero buffer if necessary
    if (x->buffer_size != vec_size) {
        if (x->zero_buffer) {
            freebytes(x->zero_buffer, x->buffer_size * sizeof(t_sample));
        }
        x->zero_buffer = (t_sample *)getbytes(vec_size * sizeof(t_sample));
        if (!x->zero_buffer) {
            pd_error(x, "qdiv~: out of memory");
            return;
        }
        x->buffer_size = vec_size;
        memset(x->zero_buffer, 0, vec_size * sizeof(t_sample));
    }

    // Assign input channels or zero buffer for left inlet
    x->left_channel_count = sp[0]->s_nchans;
    for (int i = 0; i < 4; i++) {
        x->quat_in1[i] = (i < x->left_channel_count) ? sp[0]->s_vec + vec_size * i : x->zero_buffer;
    }

    // Assign input channels or zero buffer for right inlet
    x->right_channel_count = sp[1]->s_nchans;
    for (int i = 0; i < 4; i++) {
        x->quat_in2[i] = (i < x->right_channel_count) ? sp[1]->s_vec + vec_size * i : x->zero_buffer;
    }

    // Set up output channels
    signal_setmultiout(&sp[2], 4);
    for (int i = 0; i < 4; i++) {
        x->quat_out[i] = sp[2]->s_vec + sp[2]->s_n * i;
    }

    dsp_add(qdiv_tilde_perform, 2, x, sp[0]->s_n);
}

void *qdiv_tilde_new(t_symbol *s, int argc, t_atom *argv)
{
    t_qdiv_tilde *x = (t_qdiv_tilde *)pd_new(qdiv_tilde_class);

    // Initialize fallback quaternion to identity (1 0 0 0)
    x->fallback_quat[0] = 1;
    x->fallback_quat[1] = 0;
    x->fallback_quat[2] = 0;
    x->fallback_quat[3] = 0;

    // If arguments are provided, use them to set the fallback quaternion
    if (argc > 0) {
        for (int i = 0; i < 4 && i < argc; i++) {
            x->fallback_quat[i] = atom_getfloat(argv + i);
        }
    }

    // Create inlets
    inlet_new(&x->x_obj, &x->x_obj.ob_pd, &s_signal, &s_signal);

    // Create outlet
    outlet_new(&x->x_obj, &s_signal);

    // Initialize zero buffer and channel counts
    x->zero_buffer = NULL;
    x->buffer_size = 0;
    x->left_channel_count = 0;
    x->right_channel_count = 0;

    return (void *)x;
    (void)s; // Unused argument
}

void qdiv_tilde_free(t_qdiv_tilde *x)
{
    if (x->zero_buffer) {
        freebytes(x->zero_buffer, x->buffer_size * sizeof(t_sample));
    }
}

void qdiv_tilde_setup(void)
{
    qdiv_tilde_class = class_new(gensym("qdiv~"),
        (t_newmethod)qdiv_tilde_new,
        (t_method)qdiv_tilde_free,
        sizeof(t_qdiv_tilde),
        CLASS_MULTICHANNEL,
        A_GIMME,
        0);

    class_addmethod(qdiv_tilde_class, (t_method)qdiv_tilde_dsp, gensym("dsp"), A_CANT, 0);
    CLASS_MAINSIGNALIN(qdiv_tilde_class, t_qdiv_tilde, f);
}
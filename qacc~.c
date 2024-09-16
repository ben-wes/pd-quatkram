/*
qacc~ - Quaternion accumulator external for Pure Data

2024, Ben Wesch

Functionality:
- Performs cumulative quaternion multiplication (accumulation) on a sample-by-sample basis
- Maintains an internal quaternion state that represents the accumulated transformation
- Does not normalize input quaternions, allowing for general transformations
- Normalizes the output quaternion to maintain stability

Usage:
1. Send a 4-channel signal to the inlet (w, x, y, z components of a quaternion)
2. Receive the resulting 4-channel signal from the outlet
3. Use [set 1 0 0 0( message to reset the internal state to identity (or other) quaternion

Note: This code was developed with assistance from the Anthropic Claude AI language model.
*/

#include "m_pd.h"
#include <math.h>

static t_class *qacc_tilde_class;

typedef struct _qacc_tilde {
    t_object  x_obj;
    t_sample f;
    t_sample **quat_vec;   // 4-channel quaternion input
    t_sample **out_vec;    // 4-channel quaternion output
    t_float accum_quat[4]; // Accumulated quaternion state
} t_qacc_tilde;

t_int *qacc_tilde_perform(t_int *w)
{
    t_qacc_tilde *x = (t_qacc_tilde *)(w[1]);
    int n = (int)(w[2]);

    t_sample *qw = x->quat_vec[0], *qx = x->quat_vec[1], *qy = x->quat_vec[2], *qz = x->quat_vec[3];
    t_sample *ow = x->out_vec[0], *ox = x->out_vec[1], *oy = x->out_vec[2], *oz = x->out_vec[3];

    t_float aw = x->accum_quat[0], ax = x->accum_quat[1], ay = x->accum_quat[2], az = x->accum_quat[3];

    while (n--) {
        // quaternion multiplication (accumulation step)
        t_float rw = aw * (*qw) - ax * (*qx) - ay * (*qy) - az * (*qz);
        t_float rx = aw * (*qx) + ax * (*qw) + ay * (*qz) - az * (*qy);
        t_float ry = aw * (*qy) - ax * (*qz) + ay * (*qw) + az * (*qx);
        t_float rz = aw * (*qz) + ax * (*qy) - ay * (*qx) + az * (*qw);

        // normalize result
        t_float rmag = sqrt(rw*rw + rx*rx + ry*ry + rz*rz);
        aw = rw / rmag;
        ax = rx / rmag;
        ay = ry / rmag;
        az = rz / rmag;

        // output normalized result
        *ow++ = aw;
        *ox++ = ax;
        *oy++ = ay;
        *oz++ = az;

        qw++; qx++; qy++; qz++;
    }

    // store final accumulated state for the next DSP cycle
    x->accum_quat[0] = aw;
    x->accum_quat[1] = ax;
    x->accum_quat[2] = ay;
    x->accum_quat[3] = az;

    return (w+3);
}

void qacc_tilde_dsp(t_qacc_tilde *x, t_signal **sp)
{
    int quat_channels = (int)sp[0]->s_nchans;

    // wrap quaternion channels if less than 4 given
    for (int i = 0; i < 4; i++)
        x->quat_vec[i] = sp[0]->s_vec + sp[0]->s_n * (i % quat_channels);
    
    signal_setmultiout(&sp[1], 4);
    for (int i = 0; i < 4; i++)
        x->out_vec[i] = sp[1]->s_vec + sp[1]->s_n * i;

    dsp_add(qacc_tilde_perform, 2, x, sp[0]->s_n);
}

void qacc_tilde_set(t_qacc_tilde *x, t_symbol *s, int argc, t_atom *argv)
{
    if (argc >= 4) {
        x->accum_quat[0] = atom_getfloat(argv);
        x->accum_quat[1] = atom_getfloat(argv+1);
        x->accum_quat[2] = atom_getfloat(argv+2);
        x->accum_quat[3] = atom_getfloat(argv+3);
    } else {
        pd_error(x, "qacc~: 'set' needs 4 float arguments for initial quaternion");
    }
    (void)s;
}

void *qacc_tilde_new(void)
{
    t_qacc_tilde *x = (t_qacc_tilde *)pd_new(qacc_tilde_class);

    x->quat_vec = (t_sample **)getbytes(4 * sizeof(t_sample *));
    x->out_vec = (t_sample **)getbytes(4 * sizeof(t_sample *));

    // initialize accumulated quaternion to identity
    x->accum_quat[0] = 1;
    x->accum_quat[1] = 0;
    x->accum_quat[2] = 0;
    x->accum_quat[3] = 0;

    outlet_new(&x->x_obj, &s_signal);

    return (void *)x;
}

void qacc_tilde_free(t_qacc_tilde *x)
{
    freebytes(x->quat_vec, 4 * sizeof(t_sample *));
    freebytes(x->out_vec, 4 * sizeof(t_sample *));
}

void qacc_tilde_setup(void)
{
    qacc_tilde_class = class_new(gensym("qacc~"),
        (t_newmethod)qacc_tilde_new,
        (t_method)qacc_tilde_free,
        sizeof(t_qacc_tilde),
        CLASS_MULTICHANNEL,
        0);

    class_addmethod(qacc_tilde_class, (t_method)qacc_tilde_dsp, gensym("dsp"), A_CANT, 0);
    class_addmethod(qacc_tilde_class, (t_method)qacc_tilde_set, gensym("set"), A_GIMME, 0);
    CLASS_MAINSIGNALIN(qacc_tilde_class, t_qacc_tilde, f);
}
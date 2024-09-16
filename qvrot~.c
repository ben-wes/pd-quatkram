/*
qvrot~ - Quaternion-based vector rotation external for Pure Data
2024, Ben Wesch

Functionality:
- Transforms a 3D vector using quaternion input
Usage:
1. Send a 3-channel signal to the left inlet (x, y, z components of vector to be rotated)
2. Send a 4-channel signal to the right inlet (w, x, y, z components of normalized rotation quaternion)
3. Receive the resulting 3-channel signal of the rotated vector from the outlet
*/

#include "m_pd.h"
#include <math.h>

static t_class *qvrot_tilde_class;

typedef struct _qvrot_tilde {
    t_object  x_obj;
    t_sample f;
    t_sample **vec_in;   // 3-channel vector input
    t_sample **quat_vec; // 4-channel quaternion input
    t_sample **out_vec;  // 3-channel vector output
} t_qvrot_tilde;

t_int *qvrot_tilde_perform(t_int *w)
{
    t_qvrot_tilde *x = (t_qvrot_tilde *)(w[1]);
    int n = (int)(w[2]);
    t_sample *vx = x->vec_in[0], *vy = x->vec_in[1], *vz = x->vec_in[2];
    t_sample *qw = x->quat_vec[0], *qx = x->quat_vec[1], *qy = x->quat_vec[2], *qz = x->quat_vec[3];
    t_sample *ox = x->out_vec[0], *oy = x->out_vec[1], *oz = x->out_vec[2];

    while (n--) {
        t_float rx = (*vx)*((*qw)*(*qw) + (*qx)*(*qx) - (*qy)*(*qy) - (*qz)*(*qz)) + 
                     2*(*vy)*((*qx)*(*qy) - (*qw)*(*qz)) + 
                     2*(*vz)*((*qx)*(*qz) + (*qw)*(*qy));
        t_float ry = 2*(*vx)*((*qx)*(*qy) + (*qw)*(*qz)) + 
                     (*vy)*((*qw)*(*qw) - (*qx)*(*qx) + (*qy)*(*qy) - (*qz)*(*qz)) + 
                     2*(*vz)*((*qy)*(*qz) - (*qw)*(*qx));
        t_float rz = 2*(*vx)*((*qx)*(*qz) - (*qw)*(*qy)) + 
                     2*(*vy)*((*qy)*(*qz) + (*qw)*(*qx)) + 
                     (*vz)*((*qw)*(*qw) - (*qx)*(*qx) - (*qy)*(*qy) + (*qz)*(*qz));

        *ox++ = rx;
        *oy++ = ry;
        *oz++ = rz;
        
        vx++; vy++; vz++;
        qw++; qx++; qy++; qz++;
    }
    return (w+3);
}

void qvrot_tilde_dsp(t_qvrot_tilde *x, t_signal **sp)
{
    int vec_channels = (int)sp[0]->s_nchans;
    int quat_channels = (int)sp[1]->s_nchans;

    // wrap vector channels if less than 3 given
    for (int i = 0; i < 3; i++)
        x->vec_in[i] = sp[0]->s_vec + sp[0]->s_n * (i % vec_channels);

    // wrap quaternion channels if less than 4 given
    for (int i = 0; i < 4; i++)
        x->quat_vec[i] = sp[1]->s_vec + sp[1]->s_n * (i % quat_channels);
    
    signal_setmultiout(&sp[2], 3);
    for (int i = 0; i < 3; i++)
        x->out_vec[i] = sp[2]->s_vec + sp[2]->s_n * i;

    dsp_add(qvrot_tilde_perform, 2, x, sp[0]->s_n);
}

void *qvrot_tilde_new(void)
{
    t_qvrot_tilde *x = (t_qvrot_tilde *)pd_new(qvrot_tilde_class);
    x->vec_in = (t_sample **)getbytes(3 * sizeof(t_sample *));
    x->quat_vec = (t_sample **)getbytes(4 * sizeof(t_sample *));
    x->out_vec = (t_sample **)getbytes(3 * sizeof(t_sample *));
    inlet_new(&x->x_obj, &x->x_obj.ob_pd, &s_signal, &s_signal);
    outlet_new(&x->x_obj, &s_signal);
    return (void *)x;
}

void qvrot_tilde_free(t_qvrot_tilde *x)
{
    freebytes(x->vec_in, 3 * sizeof(t_sample *));
    freebytes(x->quat_vec, 4 * sizeof(t_sample *));
    freebytes(x->out_vec, 3 * sizeof(t_sample *));
}

void qvrot_tilde_setup(void)
{
    qvrot_tilde_class = class_new(gensym("qvrot~"),
        (t_newmethod)qvrot_tilde_new,
        (t_method)qvrot_tilde_free,
        sizeof(t_qvrot_tilde),
        CLASS_MULTICHANNEL,
        0);
    class_addmethod(qvrot_tilde_class, (t_method)qvrot_tilde_dsp, gensym("dsp"), A_CANT, 0);
    CLASS_MAINSIGNALIN(qvrot_tilde_class, t_qvrot_tilde, f);
}
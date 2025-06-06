#include "m_pd.h"
#ifdef NT
#pragma warning( disable : 4244 )
#pragma warning( disable : 4305 )
#endif

/* ------------------------ qrsqrt~ ----------------------------- */



static t_class *qrsqrt_class;

typedef struct _qrsqrt
{
    t_object x_obj; 	/* obligatory header */
    t_float x_f;    	/* place to hold inlet's value if it's set by message */
} t_qrsqrt;

    /* this is the actual performance routine which acts on the samples.
    It's called with a single pointer "w" which is our location in the
    DSP call list.  We return a new "w" which will point to the next item
    after us.  Meanwhile, w[0] is just a pointer to dsp-perform itself
    (no use to us), w[1] and w[2] are the input and output vector locations,
    and w[3] is the number of points to calculate. */
static t_int *qrsqrt_perform(t_int *w)
{
    t_float *in = (t_float *)(w[1]);
    t_float *out = (t_float *)(w[2]);

    int n = (int)(w[3]);
    while (n--)
    {
        float number = *(in++);
        union {
            float    f;
            uint32_t i;
        } conv = { .f = number };
        conv.i  = 0x5f3759df - (conv.i >> 1);
        conv.f *= 1.5F - (number * 0.5F * conv.f * conv.f);
        *out++ = conv.f;
    }
    return (w+4);
}

    /* called to start DSP.  Here we call Pd back to add our perform
    routine to a linear callback list which Pd in turn calls to grind
    out the samples. */
static void qrsqrt_dsp(t_qrsqrt *x, t_signal **sp)
{
    (void)x;
    dsp_add(qrsqrt_perform, 3, sp[0]->s_vec, sp[1]->s_vec, (t_int)sp[0]->s_n);
}

static void *qrsqrt_new(void)
{
    t_qrsqrt *x = (t_qrsqrt *)pd_new(qrsqrt_class);
    outlet_new(&x->x_obj, gensym("signal"));
    x->x_f = 0;
    return (x);
}

    /* this routine, which must have exactly this name (with the "~" replaced
    by "_tilde) is called when the code is first loaded, and tells Pd how
    to build the "class". */
void qrsqrt_tilde_setup(void)
{
    qrsqrt_class = class_new(gensym("qrsqrt~"), (t_newmethod)qrsqrt_new, 0,
    	sizeof(t_qrsqrt), 0, A_DEFFLOAT, 0);
	    /* this is magic to declare that the leftmost, "main" inlet
	    takes signals; other signal inlets are done differently... */
    CLASS_MAINSIGNALIN(qrsqrt_class, t_qrsqrt, x_f);
    	/* here we tell Pd about the "dsp" method, which is called back
	when DSP is turned on. */
    class_addmethod(qrsqrt_class, (t_method)qrsqrt_dsp, gensym("dsp"), 0);
}

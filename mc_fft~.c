#include "m_pd.h"
#include <string.h>
#include <math.h>
#include <fftw3.h>

static t_class *mc_fft_tilde_class;

typedef struct _mc_fft_tilde {
    t_object  x_obj;
    t_sample  f;
    int       n_channels;
    t_sample  **in;
    t_sample  **out_real;
    t_sample  **out_imag;
    fftwf_plan plan;
    float     *fft_in;
    fftwf_complex *fft_out;
    int       vec_size;
    int       fft_size;
    t_outlet  *out_real_outlet;
    t_outlet  *out_imag_outlet;
} t_mc_fft_tilde;

static void *mc_fft_tilde_new(void)
{
    t_mc_fft_tilde *x = (t_mc_fft_tilde *)pd_new(mc_fft_tilde_class);
    
    x->n_channels = 0;
    x->in = NULL;
    x->out_real = NULL;
    x->out_imag = NULL;
    x->fft_in = NULL;
    x->fft_out = NULL;
    x->plan = NULL;
    x->fft_size = 0;
    
    x->out_real_outlet = outlet_new(&x->x_obj, &s_signal);
    x->out_imag_outlet = outlet_new(&x->x_obj, &s_signal);
    
    return (void *)x;
}

static void mc_fft_tilde_free(t_mc_fft_tilde *x)
{
    if (x->plan) fftwf_destroy_plan(x->plan);
    if (x->fft_in) fftwf_free(x->fft_in);
    if (x->fft_out) fftwf_free(x->fft_out);
    if (x->in) freebytes(x->in, x->n_channels * sizeof(t_sample *));
    if (x->out_real) freebytes(x->out_real, x->n_channels * sizeof(t_sample *));
    if (x->out_imag) freebytes(x->out_imag, x->n_channels * sizeof(t_sample *));
}

static t_int *mc_fft_tilde_perform(t_int *w)
{
    t_mc_fft_tilde *x = (t_mc_fft_tilde *)(w[1]);
    int n = (int)(w[2]);
    
    for (int i = 0; i < n; i++) {
        // Copy input channels to fft_in buffer
        for (int ch = 0; ch < x->n_channels; ch++) {
            x->fft_in[ch] = x->in[ch][i];
        }
        
        // Perform FFT
        fftwf_execute(x->plan);
        
        // Copy FFT result to output channels
        for (int ch = 0; ch < x->n_channels; ch++) {
            if (ch < x->n_channels / 2 + 1) {
                x->out_real[ch][i] = x->fft_out[ch][0];  // Real part
                x->out_imag[ch][i] = x->fft_out[ch][1];  // Imaginary part
            } else {
                // For frequencies above Nyquist, use complex conjugate
                int mirror_ch = x->n_channels - ch;
                x->out_real[ch][i] = x->fft_out[mirror_ch][0];    // Real part (same)
                x->out_imag[ch][i] = -x->fft_out[mirror_ch][1];   // Imaginary part (negated)
            }
        }
    }
    
    return (w+3);
}

static void mc_fft_tilde_dsp(t_mc_fft_tilde *x, t_signal **sp)
{
    x->vec_size = sp[0]->s_n;
    int new_n_channels = sp[0]->s_nchans;
    
    // Reallocate and reinitialize if necessary
    if (new_n_channels != x->n_channels) {
        x->n_channels = new_n_channels;
        
        if (x->in) freebytes(x->in, x->n_channels * sizeof(t_sample *));
        if (x->out_real) freebytes(x->out_real, x->n_channels * sizeof(t_sample *));
        if (x->out_imag) freebytes(x->out_imag, x->n_channels * sizeof(t_sample *));
        
        x->in = (t_sample **)getbytes(x->n_channels * sizeof(t_sample *));
        x->out_real = (t_sample **)getbytes(x->n_channels * sizeof(t_sample *));
        x->out_imag = (t_sample **)getbytes(x->n_channels * sizeof(t_sample *));
        
        if (x->fft_in) fftwf_free(x->fft_in);
        if (x->fft_out) fftwf_free(x->fft_out);
        if (x->plan) fftwf_destroy_plan(x->plan);
        
        x->fft_in = (float *)fftwf_malloc(sizeof(float) * x->n_channels);
        x->fft_out = (fftwf_complex *)fftwf_malloc(sizeof(fftwf_complex) * (x->n_channels / 2 + 1));
        x->plan = fftwf_plan_dft_r2c_1d(x->n_channels, x->fft_in, x->fft_out, FFTW_ESTIMATE);
    }
    
    // Assign input channels
    for (int i = 0; i < x->n_channels; i++) {
        x->in[i] = sp[0]->s_vec + x->vec_size * i;
    }
    
    // Set up output channels
    signal_setmultiout(&sp[1], x->n_channels);
    for (int i = 0; i < x->n_channels; i++) {
        x->out_real[i] = sp[1]->s_vec + sp[1]->s_n * i;
    }
    signal_setmultiout(&sp[2], x->n_channels);
    for (int i = 0; i < x->n_channels; i++) {
        x->out_imag[i] = sp[2]->s_vec + sp[2]->s_n * i;
    }
    
    dsp_add(mc_fft_tilde_perform, 2, x, sp[0]->s_n);
}

void mc_fft_tilde_setup(void)
{
    mc_fft_tilde_class = class_new(gensym("mc_fft~"),
        (t_newmethod)mc_fft_tilde_new,
        (t_method)mc_fft_tilde_free,
        sizeof(t_mc_fft_tilde),
        CLASS_MULTICHANNEL,
        0);
    
    class_addmethod(mc_fft_tilde_class, (t_method)mc_fft_tilde_dsp, gensym("dsp"), A_CANT, 0);
    CLASS_MAINSIGNALIN(mc_fft_tilde_class, t_mc_fft_tilde, f);
}
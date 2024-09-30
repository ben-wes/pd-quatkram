#include "m_pd.h"
#include <string.h>
#include <math.h>
#include <fftw3.h>

static t_class *mc_ifft_tilde_class;

typedef struct _mc_ifft_tilde {
    t_object  x_obj;
    t_sample  f;
    int       n_channels;
    t_sample  **in_real;
    t_sample  **in_imag;
    t_sample  **out;
    fftwf_plan plan;
    fftwf_complex *fft_in;
    float     *fft_out;
    int       vec_size;
    int       fft_size;
    t_inlet   *in_imag_inlet;
} t_mc_ifft_tilde;

static void *mc_ifft_tilde_new(void)
{
    t_mc_ifft_tilde *x = (t_mc_ifft_tilde *)pd_new(mc_ifft_tilde_class);
    
    x->n_channels = 0;
    x->in_real = NULL;
    x->in_imag = NULL;
    x->out = NULL;
    x->fft_in = NULL;
    x->fft_out = NULL;
    x->plan = NULL;
    x->fft_size = 0;
    
    x->in_imag_inlet = inlet_new(&x->x_obj, &x->x_obj.ob_pd, &s_signal, &s_signal);
    outlet_new(&x->x_obj, &s_signal);
    
    return (void *)x;
}

static void mc_ifft_tilde_free(t_mc_ifft_tilde *x)
{
    if (x->plan) fftwf_destroy_plan(x->plan);
    if (x->fft_in) fftwf_free(x->fft_in);
    if (x->fft_out) fftwf_free(x->fft_out);
    if (x->in_real) freebytes(x->in_real, x->n_channels * sizeof(t_sample *));
    if (x->in_imag) freebytes(x->in_imag, x->n_channels * sizeof(t_sample *));
    if (x->out) freebytes(x->out, x->n_channels * sizeof(t_sample *));
    inlet_free(x->in_imag_inlet);
}

static t_int *mc_ifft_tilde_perform(t_int *w)
{
    t_mc_ifft_tilde *x = (t_mc_ifft_tilde *)(w[1]);
    int n = (int)(w[2]);
    
    for (int i = 0; i < n; i++) {
        // Copy input channels to fft_in buffer
        for (int ch = 0; ch < x->n_channels; ch++) {
            if (ch < x->n_channels / 2 + 1) {
                x->fft_in[ch][0] = x->in_real[ch][i];  // Real part
                x->fft_in[ch][1] = x->in_imag[ch][i];  // Imaginary part
            } else {
                // For frequencies above Nyquist, use complex conjugate
                int mirror_ch = x->n_channels - ch;
                x->fft_in[ch][0] = x->in_real[mirror_ch][i];    // Real part (same)
                x->fft_in[ch][1] = -x->in_imag[mirror_ch][i];   // Imaginary part (negated)
            }
        }
        
        // Perform IFFT
        fftwf_execute(x->plan);
        
        // Copy IFFT result to output channels without normalization
        for (int ch = 0; ch < x->n_channels; ch++) {
            x->out[ch][i] = x->fft_out[ch];
        }
    }
    
    return (w+3);
}

static void mc_ifft_tilde_dsp(t_mc_ifft_tilde *x, t_signal **sp)
{
    x->vec_size = sp[0]->s_n;
    int new_n_channels = sp[0]->s_nchans;
    
    // Reallocate and reinitialize if necessary
    if (new_n_channels != x->n_channels) {
        x->n_channels = new_n_channels;
        
        if (x->in_real) freebytes(x->in_real, x->n_channels * sizeof(t_sample *));
        if (x->in_imag) freebytes(x->in_imag, x->n_channels * sizeof(t_sample *));
        if (x->out) freebytes(x->out, x->n_channels * sizeof(t_sample *));
        
        x->in_real = (t_sample **)getbytes(x->n_channels * sizeof(t_sample *));
        x->in_imag = (t_sample **)getbytes(x->n_channels * sizeof(t_sample *));
        x->out = (t_sample **)getbytes(x->n_channels * sizeof(t_sample *));
        
        if (x->fft_in) fftwf_free(x->fft_in);
        if (x->fft_out) fftwf_free(x->fft_out);
        if (x->plan) fftwf_destroy_plan(x->plan);
        
        x->fft_in = (fftwf_complex *)fftwf_malloc(sizeof(fftwf_complex) * x->n_channels);
        x->fft_out = (float *)fftwf_malloc(sizeof(float) * x->n_channels);
        x->plan = fftwf_plan_dft_c2r_1d(x->n_channels, x->fft_in, x->fft_out, FFTW_ESTIMATE);
    }
    
    // Assign input channels
    for (int i = 0; i < x->n_channels; i++) {
        x->in_real[i] = sp[0]->s_vec + x->vec_size * i;
        x->in_imag[i] = sp[1]->s_vec + x->vec_size * i;
    }
    
    // Set up output channels
    signal_setmultiout(&sp[2], x->n_channels);
    for (int i = 0; i < x->n_channels; i++) {
        x->out[i] = sp[2]->s_vec + sp[2]->s_n * i;
    }
    
    dsp_add(mc_ifft_tilde_perform, 2, x, sp[0]->s_n);
}

void mc_ifft_tilde_setup(void)
{
    mc_ifft_tilde_class = class_new(gensym("mc_ifft~"),
        (t_newmethod)mc_ifft_tilde_new,
        (t_method)mc_ifft_tilde_free,
        sizeof(t_mc_ifft_tilde),
        CLASS_MULTICHANNEL,
        0);
    
    class_addmethod(mc_ifft_tilde_class, (t_method)mc_ifft_tilde_dsp, gensym("dsp"), A_CANT, 0);
    CLASS_MAINSIGNALIN(mc_ifft_tilde_class, t_mc_ifft_tilde, f);
}
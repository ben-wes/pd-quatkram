#include "m_pd.h"
#include <complex.h>
#include <fftw3.h>
#include <math.h>
#include <string.h>

static t_class *frft_tilde_class;

typedef struct _frft_tilde {
    t_object x_obj;
    t_float f_dummy;
    int block_size;
    fftw_complex *buf;
    fftw_complex *work;
    fftw_plan fft;
    fftw_plan ifft;
    t_inlet *in_imag, *in_power;
    t_outlet *out_real, *out_imag;
} t_frft_tilde;

static t_int *frft_tilde_perform(t_int *w) {
    t_frft_tilde *x = (t_frft_tilde *)(w[1]);
    t_sample *in_r = (t_sample *)(w[2]), *in_i = (t_sample *)(w[3]), *power = (t_sample *)(w[4]);
    t_sample *out_r = (t_sample *)(w[5]), *out_i = (t_sample *)(w[6]);
    int N = (int)(w[7]);
    double sN = sqrt(N);
    
    // Copy input
    for(int i = 0; i < N; i++)
        x->buf[i] = in_r[i] + I * in_i[i];

    // Get power mod 4
    double a = fmod(power[0], 4.0);
    if(a < 0) a += 4.0;

    // Special cases
    if(fabs(a) < 1e-10) {
        // do nothing - keep input
    }
    else if(fabs(a - 2.0) < 1e-10) {
        // flip signal
        for(int i = 0; i < N/2; i++) {
            fftw_complex temp = x->buf[i];
            x->buf[i] = x->buf[N-1-i];
            x->buf[N-1-i] = temp;
        }
    }
    else if(fabs(a - 1.0) < 1e-10) {
        // First copy to work with shift
        for(int i = 0; i < N; i++)
            x->work[i] = x->buf[(i + N/2) % N];
            
        // FFT with scaling
        fftw_execute(x->fft);
        
        // Copy back with shift
        for(int i = 0; i < N; i++)
            x->buf[(i + N/2) % N] = x->work[i] / sN;
    }
    else if(fabs(a - 3.0) < 1e-10) {
        // First copy to work with shift
        for(int i = 0; i < N; i++)
            x->work[i] = x->buf[(i + N/2) % N];
            
        // IFFT with scaling
        fftw_execute(x->ifft);
        
        // Copy back with shift
        for(int i = 0; i < N; i++)
            x->buf[(i + N/2) % N] = x->work[i] * sN / N;
    }
    else {
        // Reduce to interval 0.5 < a < 1.5
        if(a > 2.0) {
            a -= 2.0;
            // Flip signal
            for(int i = 0; i < N/2; i++) {
                fftw_complex temp = x->buf[i];
                x->buf[i] = x->buf[N-1-i];
                x->buf[N-1-i] = temp;
            }
        }
        
        if(a > 1.5) {
            a -= 1.0;
            // FFT with shift
            for(int i = 0; i < N; i++)
                x->work[i] = x->buf[(i + N/2) % N];
            fftw_execute(x->fft);
            for(int i = 0; i < N; i++)
                x->buf[(i + N/2) % N] = x->work[i] / sN;
        }
        
        if(a < 0.5) {
            a += 1.0;
            // IFFT with shift
            for(int i = 0; i < N; i++)
                x->work[i] = x->buf[(i + N/2) % N];
            fftw_execute(x->ifft);
            for(int i = 0; i < N; i++)
                x->buf[(i + N/2) % N] = x->work[i] * sN / N;
        }
        
        // Core FRFT for 0.5 < a < 1.5
        double alpha = a * M_PI / 2;
        double tana2 = tan(alpha / 2);
        double sina = sin(alpha);
        
        // First chirp multiplication
        for(int i = 0; i < N; i++) {
            double idx = i - N/2;  // Center around 0
            fftw_complex chirp = cexp(-I * M_PI/N * tana2/4 * idx * idx);
            x->buf[i] *= chirp;
        }
        
        // Second chirp multiplication with convolution
        double c = M_PI/N/sina/4;
        for(int i = 0; i < N; i++) {
            double idx = i - N/2;
            fftw_complex chirp = cexp(I * c * idx * idx);
            x->buf[i] *= chirp;
        }
    }

    // Copy output
    for(int i = 0; i < N; i++) {
        out_r[i] = creal(x->buf[i]);
        out_i[i] = cimag(x->buf[i]);
    }

    return (w + 8);
}

static void frft_tilde_dsp(t_frft_tilde *x, t_signal **sp) {
    int vec_size = sp[0]->s_n;
    
    if(vec_size != x->block_size) {
        if(x->buf) fftw_free(x->buf);
        if(x->work) fftw_free(x->work);
        if(x->fft) fftw_destroy_plan(x->fft);
        if(x->ifft) fftw_destroy_plan(x->ifft);
        
        x->block_size = vec_size;
        x->buf = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * vec_size);
        x->work = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * vec_size);
        x->fft = fftw_plan_dft_1d(vec_size, x->work, x->work, FFTW_FORWARD, FFTW_MEASURE);
        x->ifft = fftw_plan_dft_1d(vec_size, x->work, x->work, FFTW_BACKWARD, FFTW_MEASURE);
    }
    
    dsp_add(frft_tilde_perform, 7, x,
            sp[0]->s_vec,  // real input
            sp[1]->s_vec,  // imag input
            sp[2]->s_vec,  // power input
            sp[3]->s_vec,  // real output
            sp[4]->s_vec,  // imag output
            vec_size);
}

static void *frft_tilde_new(t_floatarg f) {
    t_frft_tilde *x = (t_frft_tilde *)pd_new(frft_tilde_class);
    
    x->block_size = 0;
    x->buf = NULL;
    x->work = NULL;
    x->fft = NULL;
    x->ifft = NULL;
    
    x->in_imag = signalinlet_new(&x->x_obj, 0.0);
    x->in_power = signalinlet_new(&x->x_obj, f);
    x->out_real = outlet_new(&x->x_obj, &s_signal);
    x->out_imag = outlet_new(&x->x_obj, &s_signal);
    
    return (void *)x;
}

static void frft_tilde_free(t_frft_tilde *x) {
    if(x->buf) fftw_free(x->buf);
    if(x->work) fftw_free(x->work);
    if(x->fft) fftw_destroy_plan(x->fft);
    if(x->ifft) fftw_destroy_plan(x->ifft);
}

void frft_tilde_setup(void) {
    frft_tilde_class = class_new(gensym("frft~"),
        (t_newmethod)frft_tilde_new,
        (t_method)frft_tilde_free,
        sizeof(t_frft_tilde),
        CLASS_DEFAULT,
        A_DEFFLOAT, 0);
    
    class_addmethod(frft_tilde_class,
        (t_method)frft_tilde_dsp, gensym("dsp"), A_CANT, 0);
    
    CLASS_MAINSIGNALIN(frft_tilde_class, t_frft_tilde, f_dummy);
}

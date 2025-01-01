#include "m_pd.h"
#include <complex.h>
#include <fftw3.h>
#include <math.h>
#include <string.h>

#define FFTW_COMPLEX(r,i) ((r) + I * (i))

static t_class *frft_tilde_class;

typedef struct _frft_tilde {
    t_object x_obj;
    t_float f_dummy;
    int block_size;
    int conv_size;          // Size of main buffers
    int saved_chirp_size;   // Size of chirp buffer
    fftw_complex *buf;
    fftw_complex *work;
    fftw_plan signal_fft;    // For conv_size
    fftw_plan kernel_fft;    // For conv_size
    fftw_plan ifft_conv;     // For conv_size convolution
    fftw_plan ifft_shift;    // For size N shift operations
    fftw_plan fft_shift;     // For size N shift operations
    t_inlet *in_imag, *in_power;
    t_outlet *out_real, *out_imag;
    fftw_complex *saved_chirp;  // Store first chirp
    t_float current_power;      // Track last power value
    t_float current_alpha;      // Track last alpha value
} t_frft_tilde;

// Helper function declarations
static void flip_signal(fftw_complex *buf, int N);
static void fft_with_shift(t_frft_tilde *x, int N);
static void ifft_with_shift(t_frft_tilde *x, int N);
static void frft_tilde_free(t_frft_tilde *x);

static void sincinterp(t_frft_tilde *x, int N) {
    // Step 1: Zero-pad signal
    int work_size = 2*N-1;
    memset(x->work, 0, sizeof(fftw_complex) * work_size);
    for(int i = 0; i < N; i++) {
        x->work[2*i] = x->buf[i];
    }

    // Step 2: Create sinc kernel
    int kernel_size = 4*N-5;
    for(int i = 0; i < kernel_size; i++) {
        t_float t = (i - (2*N-3)) / 2.0;
        x->buf[i] = (fabs(t) < 1e-10) ? 1.0 : sin(M_PI * t)/(M_PI * t);
    }

    // Step 3: Convolution setup
    int conv_size = 8*N-6;
    memset(x->work + work_size, 0, sizeof(fftw_complex) * (conv_size - work_size));
    memset(x->buf + kernel_size, 0, sizeof(fftw_complex) * (conv_size - kernel_size));

    // FFT of signal and kernel
    fftw_execute(x->signal_fft);    
    fftw_execute(x->kernel_fft);
    
    // Multiply spectrums
    for(int i = 0; i < conv_size; i++) {
        x->work[i] *= x->buf[i];
    }
    
    // IFFT
    fftw_execute(x->ifft_conv);
    
    // Scale and extract result
    t_float scale = 1.0/conv_size;
    for(int i = 0; i < work_size; i++) {
        x->buf[i] = x->work[i + 2*N-3] * scale;
    }
}

static void flip_signal(fftw_complex *buf, int N) {
    for(int i = 0; i < N/2; i++) {
        fftw_complex temp = buf[i];
        buf[i] = buf[N-1-i];
        buf[N-1-i] = temp;
    }
}

static void fft_with_shift(t_frft_tilde *x, int N) {
    // Copy with shift
    for(int i = 0; i < N; i++) {
        int shifted_idx = (i + N/2) % N;
        x->work[i] = x->buf[shifted_idx];
    }
    
    fftw_execute(x->fft_shift);
    
    // Scale by 1/sqrt(N) to match numpy.fft
    t_float scale = 1.0/sqrt(N);
    for(int i = 0; i < N; i++) {
        int shifted_idx = (i + N/2) % N;
        x->buf[shifted_idx] = x->work[i] * scale;
    }
}

static void ifft_with_shift(t_frft_tilde *x, int N) {
    // First copy to work with shift
    for(int i = 0; i < N; i++) {
        x->work[i] = x->buf[(i + N/2) % N];
    }
    
    fftw_execute(x->ifft_shift);
    
    // Copy back with shift, correct scaling to match numpy's ifft
    t_float sN = sqrt(N);
    for(int i = 0; i < N; i++) {
        x->buf[(i + N/2) % N] = x->work[i] * sN / N;
    }
}

static t_int *frft_tilde_perform(t_int *w) {
    t_frft_tilde *x = (t_frft_tilde *)(w[1]);
    t_sample *in_r = (t_sample *)(w[2]), 
             *in_i = (t_sample *)(w[3]),
             *power = (t_sample *)(w[4]), 
             *out_r = (t_sample *)(w[5]),
             *out_i = (t_sample *)(w[6]);
    int N = (int)(w[7]);

    // Zero main buffers
    memset(x->buf, 0, sizeof(fftw_complex) * x->conv_size);
    memset(x->work, 0, sizeof(fftw_complex) * x->conv_size);

    // Copy input
    for(int i = 0; i < N; i++) {
        x->buf[i] = FFTW_COMPLEX(in_r[i], in_i[i]);
    }

    t_float a = fmod(power[0], 4.0);
    if(a < 0) a += 4.0;

    // Handle special cases
    if(fabs(a) < 1e-10) {
        for(int i = 0; i < N; i++) {
            out_r[i] = creal(x->buf[i]);
            out_i[i] = cimag(x->buf[i]);
        }
        return (w + 8);
    }
    if(fabs(a - 2.0) < 1e-10) {
        flip_signal(x->buf, N);
        for(int i = 0; i < N; i++) {
            out_r[i] = creal(x->buf[i]);
            out_i[i] = cimag(x->buf[i]);
        }
        return (w + 8);
    }
    if(fabs(a - 1.0) < 1e-10) {
        fft_with_shift(x, N);
        for(int i = 0; i < N; i++) {
            out_r[i] = creal(x->buf[i]);
            out_i[i] = cimag(x->buf[i]);
        }
        return (w + 8);
    }
    if(fabs(a - 3.0) < 1e-10) {
        ifft_with_shift(x, N);
        for(int i = 0; i < N; i++) {
            out_r[i] = creal(x->buf[i]);
            out_i[i] = cimag(x->buf[i]);
        }
        return (w + 8);
    }

    // Reduce to interval 0.5 < a < 1.5
    if(a > 2.0) {
        a = a - 2.0;
        flip_signal(x->buf, N);
    }
    if(a > 1.5) {
        a = a - 1.0;
        fft_with_shift(x, N);
    }
    if(a < 0.5) {
        a = a + 1.0;
        ifft_with_shift(x, N);
    }

    // Core FRFT computation
    t_float alpha = a * M_PI / 2;
    t_float tana2 = tan(alpha / 2);
    t_float sina = sin(alpha);

    // 1. Sinc interpolation
    sincinterp(x, N);

    // 2. Simple zero padding
    int sinc_size = 2*N-1;
    memcpy(x->work, x->buf, sizeof(fftw_complex) * sinc_size);
    
    int final_size = 4*N-3;
    memset(x->buf, 0, sizeof(fftw_complex) * final_size);
    
    int left_pad = N-1;
    memcpy(x->buf + left_pad, x->work, sizeof(fftw_complex) * sinc_size);

    // 3. First chirp multiplication
    int pad_size = 4*N-3;
    
    // Cache power value and only recompute chirps if power changed
    if (fabs(power[0] - x->current_power) > 1e-10) {
        x->current_power = power[0];
        x->current_alpha = alpha;
        
        // Recompute saved_chirp only when power changes
        for(int i = 0; i < pad_size; i++) {
            t_float t = -2*N + 2 + i;
            t_float phase = -M_PI/N * tana2/4 * t * t;
            x->saved_chirp[i] = cexp(I * phase);
        }
    }

    // Use pre-computed chirp
    for(int i = 0; i < pad_size; i++) {
        x->work[i] = x->buf[i] * x->saved_chirp[i];
    }

    // 4. Second chirp (convolution)
    t_float c = M_PI/N/sina/4;
    
    memcpy(x->buf, x->work, sizeof(fftw_complex) * pad_size);
    memset(x->work, 0, sizeof(fftw_complex) * x->conv_size);
    memcpy(x->work, x->buf, sizeof(fftw_complex) * pad_size);
    
    for(int i = 0; i < x->conv_size; i++) {
        t_float t = -(4*N - 4) + i;
        x->buf[i] = cexp(I * c * t * t);
    }

    fftw_execute(x->signal_fft);
    fftw_execute(x->kernel_fft);
    
    for(int i = 0; i < x->conv_size; i++) {
        x->work[i] *= x->buf[i];
    }
    
    fftw_execute(x->ifft_conv);
    
    t_float scale = sqrt(c/M_PI) / x->conv_size;
    int slice_start = 4*N-4;
    
    for(int i = 0; i < pad_size; i++) {
        x->buf[i] = x->work[i + slice_start] * scale * x->saved_chirp[i];
    }

    // 5. Final chirp
    t_float final_phase = -M_PI/4 * (1-a);
    fftw_complex final_scale = cexp(I * final_phase);
    
    for(int i = 0; i < N; i++) {
        x->work[i] = x->buf[N-1 + 2*i] * final_scale;
    }

    // Output
    for(int i = 0; i < N; i++) {
        out_r[i] = creal(x->work[i]);
        out_i[i] = cimag(x->work[i]);
    }

    return (w + 8);
}

static void frft_tilde_dsp(t_frft_tilde *x, t_signal **sp) {
    int N = sp[0]->s_n;
    int conv_size = 8 * N - 6;

    if (N != x->block_size) {
        // Clean up old resources
        frft_tilde_free(x);
        x->block_size = 0;  // Reset these since they might have been used in free
        x->conv_size = 0;
        x->saved_chirp_size = 0;

        // Allocate new resources
        x->block_size = N;
        x->conv_size = conv_size;
        x->saved_chirp_size = 4*N-3;
        x->buf = (fftw_complex *)getbytes(sizeof(fftw_complex) * conv_size);
        x->work = (fftw_complex *)getbytes(sizeof(fftw_complex) * conv_size);
        x->saved_chirp = (fftw_complex *)getbytes(sizeof(fftw_complex) * x->saved_chirp_size);

        if (!x->buf || !x->work || !x->saved_chirp) {
            pd_error(x, "frft~: memory allocation failed");
            frft_tilde_free(x);  // Clean up on failure
            x->block_size = 0;
            x->conv_size = 0;
            x->saved_chirp_size = 0;
            return;
        }

        // Initialize memory
        memset(x->buf, 0, sizeof(fftw_complex) * conv_size);
        memset(x->work, 0, sizeof(fftw_complex) * conv_size);
        memset(x->saved_chirp, 0, sizeof(fftw_complex) * x->saved_chirp_size);

        // Reset cached values
        x->current_power = -1;  // Force recomputation of chirp
        
        // Create FFT plans
        x->signal_fft = fftw_plan_dft_1d(conv_size, x->work, x->work, 
                                        FFTW_FORWARD, FFTW_PATIENT);
        x->kernel_fft = fftw_plan_dft_1d(conv_size, x->buf, x->buf, 
                                        FFTW_FORWARD, FFTW_PATIENT);
        x->ifft_conv = fftw_plan_dft_1d(conv_size, x->work, x->work, 
                                        FFTW_BACKWARD, FFTW_PATIENT);
        x->ifft_shift = fftw_plan_dft_1d(N, x->work, x->work, 
                                        FFTW_BACKWARD, FFTW_PATIENT);
        x->fft_shift = fftw_plan_dft_1d(N, x->work, x->work, 
                                       FFTW_FORWARD, FFTW_PATIENT);

        if (!x->signal_fft || !x->kernel_fft || !x->ifft_conv || !x->ifft_shift || !x->fft_shift) {
            pd_error(x, "frft~: FFT plan creation failed");
            frft_tilde_free(x);  // Clean up on failure
            x->block_size = 0;
            x->conv_size = 0;
            x->saved_chirp_size = 0;
            return;
        }
    }

    dsp_add(frft_tilde_perform, 7, x, 
            sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec, 
            sp[3]->s_vec, sp[4]->s_vec, N);
}

static void *frft_tilde_new(t_floatarg f) {
    t_frft_tilde *x = (t_frft_tilde *)pd_new(frft_tilde_class);
    
    // Initialize struct members
    x->block_size = 0;
    x->conv_size = 0;
    x->saved_chirp_size = 0;
    x->buf = NULL;
    x->work = NULL;
    x->signal_fft = NULL;
    x->kernel_fft = NULL;
    x->ifft_conv = NULL;
    x->ifft_shift = NULL;
    x->fft_shift = NULL;
    x->saved_chirp = NULL;
    x->current_power = -1;
    x->current_alpha = 0;
    
    // Create inlets
    x->in_imag = signalinlet_new(&x->x_obj, 0.0);
    x->in_power = signalinlet_new(&x->x_obj, f);
    
    // Create signal outlets
    x->out_real = outlet_new(&x->x_obj, &s_signal);
    x->out_imag = outlet_new(&x->x_obj, &s_signal);
    
    return (void *)x;
}

static void frft_tilde_free(t_frft_tilde *x) {
    // Free FFTW resources
    if (x->buf) freebytes(x->buf, sizeof(fftw_complex) * x->conv_size);
    if (x->work) freebytes(x->work, sizeof(fftw_complex) * x->conv_size);
    if (x->saved_chirp) freebytes(x->saved_chirp, sizeof(fftw_complex) * x->saved_chirp_size);
    if (x->signal_fft) fftw_destroy_plan(x->signal_fft);
    if (x->kernel_fft) fftw_destroy_plan(x->kernel_fft);
    if (x->ifft_shift) fftw_destroy_plan(x->ifft_shift);
    if (x->fft_shift) fftw_destroy_plan(x->fft_shift);
    if (x->ifft_conv) fftw_destroy_plan(x->ifft_conv);
    
    // Inlets/outlets are automatically freed by Pd
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

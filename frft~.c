#include "m_pd.h"
#include <complex.h>
#include <fftw3.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

// For debug printing
#include <stdio.h>

#define FFTW_COMPLEX(r,i) ((r) + I * (i))

static t_class *frft_tilde_class;

typedef struct _frft_tilde {
    t_object x_obj;
    t_float f_dummy;
    int block_size;
    fftw_complex *buf;
    fftw_complex *work;
    fftw_plan signal_fft;    // For conv_size
    fftw_plan kernel_fft;    // For conv_size
    fftw_plan ifft_conv;     // For conv_size convolution
    fftw_plan ifft_shift;    // For size N shift operations
    fftw_plan fft_shift;     // For size N shift operations
    t_inlet *in_imag, *in_power;
    t_outlet *out_real, *out_imag;
    t_outlet *list_out;
    fftw_complex *saved_chirp;  // Store first chirp
    int current_power;          // Track last power value
    double current_alpha;       // Track last alpha value
} t_frft_tilde;

// Helper function declarations (Corrected decimate_signal declaration)
static void flip_signal(fftw_complex *buf, int N);
static void fft_with_shift(t_frft_tilde *x, int N);
static void ifft_with_shift(t_frft_tilde *x, int N);

// static void output_stage(t_frft_tilde *x, const char *stage_name, fftw_complex *data, int size) {
//     // Create suffixed message names for real and imaginary parts
//     char stage_name_real[256], stage_name_imag[256];
//     snprintf(stage_name_real, sizeof(stage_name_real), "%s_r", stage_name);
//     snprintf(stage_name_imag, sizeof(stage_name_imag), "%s_i", stage_name);

//     // Allocate memory for real and imaginary parts
//     t_atom *real_data = (t_atom *)getbytes(sizeof(t_atom) * size);
//     t_atom *imag_data = (t_atom *)getbytes(sizeof(t_atom) * size);

//     // Fill lists with real and imaginary parts
//     for (int i = 0; i < size; i++) {
//         SETFLOAT(real_data + i, creal(data[i]));
//         SETFLOAT(imag_data + i, cimag(data[i]));
//     }

//     // Output stage-specific messages for real and imaginary parts
//     outlet_anything(x->list_out, gensym(stage_name_real), size, real_data);
//     outlet_anything(x->list_out, gensym(stage_name_imag), size, imag_data);

//     // Free allocated memory
//     freebytes(real_data, sizeof(t_atom) * size);
//     freebytes(imag_data, sizeof(t_atom) * size);
// }

static void sincinterp(t_frft_tilde *x, int N) {
    // Step 1: Zero-pad signal
    int work_size = 2*N-1;
    memset(x->work, 0, sizeof(fftw_complex) * work_size);
    for(int i = 0; i < N; i++) {
        x->work[2*i] = x->buf[i];
    }
    // output_stage(x, "padded", x->work, work_size);

    // Step 2: Create sinc kernel
    int kernel_size = 4*N-5;
    for(int i = 0; i < kernel_size; i++) {
        double t = (i - (2*N-3)) / 2.0;
        x->buf[i] = (fabs(t) < 1e-10) ? 1.0 : sin(M_PI * t)/(M_PI * t);
    }
    // output_stage(x, "kernel", x->buf, kernel_size);

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
    double scale = 1.0/conv_size;
    for(int i = 0; i < work_size; i++) {
        x->buf[i] = x->work[i + 2*N-3] * scale;
    }
    // output_stage(x, "after_sinc", x->buf, work_size);
}

// Special case helpers
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
    double scale = 1.0/sqrt(N);
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
    double sN = sqrt(N);
    for(int i = 0; i < N; i++) {
        x->buf[(i + N/2) % N] = x->work[i] * sN / N;  // Scale by sN/N matches numpy
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
    int conv_size = 8*N-6;
    memset(x->buf, 0, sizeof(fftw_complex) * conv_size);
    memset(x->work, 0, sizeof(fftw_complex) * conv_size);

    // Copy input and handle a
    for(int i = 0; i < N; i++) {
        x->buf[i] = FFTW_COMPLEX(in_r[i], in_i[i]);
    }
    // output_stage(x, "original", x->buf, N);

    double a = fmod(power[0], 4.0);
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
    double alpha = a * M_PI / 2;
    double tana2 = tan(alpha / 2);
    double sina = sin(alpha);

    // 1. Sinc interpolation gives us 2N-1 samples in x->buf
    sincinterp(x, N);
    // output_stage(x, "after_sinc", x->buf, 2*N-1);

    // 2. Simple zero padding
    // First save our interpolated signal
    int sinc_size = 2*N-1;
    memcpy(x->work, x->buf, sizeof(fftw_complex) * sinc_size);
    
    // Clear the entire buffer
    int final_size = 4*N-3;
    memset(x->buf, 0, sizeof(fftw_complex) * final_size);
    
    // Put the interpolated signal in the middle
    int left_pad = N-1;
    memcpy(x->buf + left_pad, x->work, sizeof(fftw_complex) * sinc_size);
    
    // output_stage(x, "after_pad", x->buf, final_size);

    // 3. First chirp multiplication
    int pad_size = 4*N-3;
    
    // Cache power value and only recompute chirps if power changed
    if (fabs(power[0] - x->current_power) > 1e-10) {
        x->current_power = power[0];
        x->current_alpha = a * M_PI / 2;
        
        // Recompute saved_chirp only when power changes
        for(int i = 0; i < pad_size; i++) {
            double t = -2*N + 2 + i;
            double phase = -M_PI/N * tana2/4 * t * t;
            x->saved_chirp[i] = cexp(I * phase);
        }
    }

    // Use pre-computed chirp instead of recreating it
    for(int i = 0; i < pad_size; i++) {
        x->work[i] = x->buf[i] * x->saved_chirp[i];
    }
    // output_stage(x, "after_chirp1", x->work, pad_size);

    // 4. Second chirp (convolution)
    double c = M_PI/N/sina/4;
    
    // Copy our chirped signal
    memcpy(x->buf, x->work, sizeof(fftw_complex) * pad_size);
    
    // Clear work and prepare for convolution
    memset(x->work, 0, sizeof(fftw_complex) * conv_size);
    memcpy(x->work, x->buf, sizeof(fftw_complex) * pad_size);
    
    // Create convolution chirp (this is different from first chirp!)
    for(int i = 0; i < conv_size; i++) {
        double t = -(4*N - 4) + i;
        x->buf[i] = cexp(I * c * t * t);
    }

    // Do convolution via FFT
    fftw_execute(x->signal_fft);
    fftw_execute(x->kernel_fft);
    
    for(int i = 0; i < conv_size; i++) {
        x->work[i] *= x->buf[i];
    }
    
    fftw_execute(x->ifft_conv);
    
    // After convolution, take slice and THEN multiply with first chirp again
    double scale = sqrt(c/M_PI) / conv_size;
    int slice_start = 4*N-4;
    int slice_len = 4*N-3;
    
    // Take slice and multiply with saved first chirp
    for(int i = 0; i < slice_len; i++) {
        x->buf[i] = x->work[i + slice_start] * scale * x->saved_chirp[i];
    }
    // output_stage(x, "after_conv", x->buf, slice_len);

    // 5. Final chirp
    double final_phase = -M_PI/4 * (1-a);
    fftw_complex final_scale = cexp(I * final_phase);
    
    for(int i = 0; i < N; i++) {
        x->work[i] = x->buf[N-1 + 2*i] * final_scale;
    }
    // output_stage(x, "final", x->work, N);

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
        // Clean up old plans and buffers
        if (x->buf) fftw_free(x->buf);
        if (x->work) fftw_free(x->work);
        if (x->saved_chirp) fftw_free(x->saved_chirp);  // Free old chirp buffer
        if (x->signal_fft) fftw_destroy_plan(x->signal_fft);
        if (x->kernel_fft) fftw_destroy_plan(x->kernel_fft);
        if (x->ifft_shift) fftw_destroy_plan(x->ifft_shift);
        if (x->ifft_conv) fftw_destroy_plan(x->ifft_conv);

        // Allocate new buffers
        x->block_size = N;
        x->buf = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * conv_size);
        x->work = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * conv_size);
        x->saved_chirp = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * (4*N-3));  // Allocate chirp buffer

        if (!x->buf || !x->work || !x->saved_chirp) {
            pd_error(x, "frft~: memory allocation failed");
            // Clean up on failure
            if (x->buf) fftw_free(x->buf);
            if (x->work) fftw_free(x->work);
            if (x->saved_chirp) fftw_free(x->saved_chirp);
            x->buf = x->work = x->saved_chirp = NULL;
            x->block_size = 0;
            return;
        }

        // Initialize buffers
        memset(x->buf, 0, sizeof(fftw_complex) * conv_size);
        memset(x->work, 0, sizeof(fftw_complex) * conv_size);
        memset(x->saved_chirp, 0, sizeof(fftw_complex) * (4*N-3));

        // Reset cached values
        x->current_power = -1;  // Force recomputation of chirp
        
        // In frft_tilde_dsp:
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

        if (!x->signal_fft || !x->kernel_fft || !x->ifft_conv || !x->ifft_shift) {
            pd_error(x, "frft~: FFT plan creation failed");
            if (x->buf) fftw_free(x->buf);
            if (x->work) fftw_free(x->work);
            x->buf = x->work = NULL;
            x->block_size = 0;
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
    x->buf = NULL;
    x->work = NULL;
    x->signal_fft = NULL;
    x->kernel_fft = NULL;
    x->ifft_conv = NULL;
    x->ifft_shift = NULL;
    x->fft_shift = NULL;
    x->saved_chirp = NULL;        // Initialize to NULL
    x->current_power = -1;        // Invalid initial value
    x->current_alpha = 0.0;       // Initialize alpha
    
    // Create inlets
    x->in_imag = signalinlet_new(&x->x_obj, 0.0);
    x->in_power = signalinlet_new(&x->x_obj, f);
    
    // Create signal outlets
    x->out_real = outlet_new(&x->x_obj, &s_signal);
    x->out_imag = outlet_new(&x->x_obj, &s_signal);
    
    // Create list outlet for visualization
    x->list_out = outlet_new(&x->x_obj, &s_list);
    
    return (void *)x;
}

static void frft_tilde_free(t_frft_tilde *x) {
    // Free FFTW resources
    if (x->buf) fftw_free(x->buf);
    if (x->work) fftw_free(x->work);
    if (x->signal_fft) fftw_destroy_plan(x->signal_fft);
    if (x->kernel_fft) fftw_destroy_plan(x->kernel_fft);
    if (x->ifft_shift) fftw_destroy_plan(x->ifft_shift);
    if (x->fft_shift) fftw_destroy_plan(x->fft_shift);
    if (x->ifft_conv) fftw_destroy_plan(x->ifft_conv);
    if (x->saved_chirp) fftw_free(x->saved_chirp);
    
    // Outlets are automatically freed by Pd
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

#include "m_pd.h"
#include <complex.h>
#include <fftw3.h>
#include <math.h>
#include <string.h>

static t_class *frft_class;

typedef struct _frft {
    t_object x_obj;
    t_symbol *x_array_in;     // Input array name
    t_symbol *x_array_out;    // Output array name
    t_float x_alpha;          // FRFT power (using t_float instead of double)
    
    // FFTW resources
    int last_size;            // Last processed array size
    int conv_size;            // Size of main buffers
    int saved_chirp_size;     // Size of chirp buffer
    fftw_complex *buf;        // Main buffer
    fftw_complex *work;       // Work buffer
    fftw_complex *saved_chirp; // Cached chirp for current alpha
    fftw_plan signal_fft;     // FFT plans
    fftw_plan kernel_fft;
    fftw_plan ifft_conv;
    fftw_plan ifft_shift;
    fftw_plan fft_shift;
    
    // List processing output
    t_outlet *list_out_real;
    t_outlet *list_out_imag;
} t_frft;

// Function declarations
static void frft_free(t_frft *x);
static void frft_bang(t_frft *x);
static void frft_list(t_frft *x, t_symbol *s, int argc, t_atom *argv);
static void *frft_new(t_symbol *s, int argc, t_atom *argv);
static void frft_set(t_frft *x, t_symbol *s1, t_symbol *s2);
static void frft_alpha(t_frft *x, t_floatarg f);
static void frft_process(t_frft *x, int N);

// Helper functions (from original frft~)
static void flip_signal(fftw_complex *buf, int N) {
    for(int i = 0; i < N/2; i++) {
        fftw_complex temp = buf[i];
        buf[i] = buf[N-1-i];
        buf[N-1-i] = temp;
    }
}

static void fft_with_shift(t_frft *x, int N) {
    for(int i = 0; i < N; i++) {
        int shifted_idx = (i + N/2) % N;
        x->work[i] = x->buf[shifted_idx];
    }
    
    fftw_execute(x->fft_shift);
    
    double scale = 1.0/sqrt(N);
    for(int i = 0; i < N; i++) {
        int shifted_idx = (i + N/2) % N;
        x->buf[shifted_idx] = x->work[i] * scale;
    }
}

static void ifft_with_shift(t_frft *x, int N) {
    for(int i = 0; i < N; i++) {
        x->work[i] = x->buf[(i + N/2) % N];
    }
    
    fftw_execute(x->ifft_shift);
    
    double sN = sqrt(N);
    for(int i = 0; i < N; i++) {
        x->buf[(i + N/2) % N] = x->work[i] * sN / N;
    }
}

static void sincinterp(t_frft *x, int N) {
    int work_size = 2*N-1;
    memset(x->work, 0, sizeof(fftw_complex) * work_size);
    for(int i = 0; i < N; i++) {
        x->work[2*i] = x->buf[i];
    }

    int kernel_size = 4*N-5;
    for(int i = 0; i < kernel_size; i++) {
        double t = (i - (2*N-3)) / 2.0;
        x->buf[i] = (fabs(t) < 1e-10) ? 1.0 : sin(M_PI * t)/(M_PI * t);
    }

    int conv_size = 8*N-6;
    memset(x->work + work_size, 0, sizeof(fftw_complex) * (conv_size - work_size));
    memset(x->buf + kernel_size, 0, sizeof(fftw_complex) * (conv_size - kernel_size));

    fftw_execute(x->signal_fft);    
    fftw_execute(x->kernel_fft);
    
    for(int i = 0; i < conv_size; i++) {
        x->work[i] *= x->buf[i];
    }
    
    fftw_execute(x->ifft_conv);
    
    double scale = 1.0/conv_size;
    for(int i = 0; i < work_size; i++) {
        x->buf[i] = x->work[i + 2*N-3] * scale;
    }
}

static void frft_process(t_frft *x, int N) {
    double a = fmod(x->x_alpha, 4.0);
    if(a < 0) a += 4.0;

    // Handle special cases
    if(fabs(a) < 1e-10) return;  // a = 0: identity
    if(fabs(a - 2.0) < 1e-10) {  // a = 2: signal flip
        flip_signal(x->buf, N);
        return;
    }
    if(fabs(a - 1.0) < 1e-10) {  // a = 1: regular FFT
        fft_with_shift(x, N);
        return;
    }
    if(fabs(a - 3.0) < 1e-10) {  // a = 3: inverse FFT
        ifft_with_shift(x, N);
        return;
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
    for(int i = 0; i < pad_size; i++) {
        double t = -2*N + 2 + i;
        x->work[i] = x->buf[i] * cexp(I * -M_PI/N * tana2/4 * t * t);
    }

    // 4. Second chirp (convolution)
    int conv_size = 8*N-6;
    double c = M_PI/N/sina/4;
    
    memcpy(x->buf, x->work, sizeof(fftw_complex) * pad_size);
    memset(x->work, 0, sizeof(fftw_complex) * conv_size);
    memcpy(x->work, x->buf, sizeof(fftw_complex) * pad_size);
    
    for(int i = 0; i < conv_size; i++) {
        double t = -(4*N - 4) + i;
        x->buf[i] = cexp(I * c * t * t);
    }

    fftw_execute(x->signal_fft);
    fftw_execute(x->kernel_fft);
    
    for(int i = 0; i < conv_size; i++) {
        x->work[i] *= x->buf[i];
    }
    
    fftw_execute(x->ifft_conv);
    
    double scale = sqrt(c/M_PI) / conv_size;
    int slice_start = 4*N-4;
    
    for(int i = 0; i < pad_size; i++) {
        double t = -2*N + 2 + i;
        x->buf[i] = x->work[i + slice_start] * scale * 
                    cexp(I * -M_PI/N * tana2/4 * t * t);
    }

    // 5. Final phase factor and extract result
    double final_phase = -M_PI/4 * (1-a);
    fftw_complex final_scale = cexp(I * final_phase);
    
    for(int i = 0; i < N; i++) {
        x->work[i] = x->buf[N-1 + 2*i] * final_scale;
    }
    
    memcpy(x->buf, x->work, sizeof(fftw_complex) * N);
}

// Helper function for buffer management
static int prepare_buffers(t_frft *x, int size) {
    if (size != x->last_size) {
        // Clean up old resources
        frft_free(x);
        x->last_size = 0;  // Reset since it might have been used in free
        x->conv_size = 0;
        x->saved_chirp_size = 0;

        int conv_size = 8 * size - 6;        
        // Allocate new resources
        x->conv_size = conv_size;
        x->saved_chirp_size = 4 * size - 3;
        x->buf = (fftw_complex *)getbytes(sizeof(fftw_complex) * conv_size);
        x->work = (fftw_complex *)getbytes(sizeof(fftw_complex) * conv_size);
        x->saved_chirp = (fftw_complex *)getbytes(sizeof(fftw_complex) * x->saved_chirp_size);
        
        if (!x->buf || !x->work || !x->saved_chirp) {
            pd_error(x, "out of memory");
            frft_free(x);  // Clean up on failure
            x->last_size = 0;
            x->conv_size = 0;
            x->saved_chirp_size = 0;
            return 0;
        }
        
        x->signal_fft = fftw_plan_dft_1d(conv_size, x->work, x->work, 
                                        FFTW_FORWARD, FFTW_PATIENT);
        x->kernel_fft = fftw_plan_dft_1d(conv_size, x->buf, x->buf, 
                                        FFTW_FORWARD, FFTW_PATIENT);
        x->ifft_conv = fftw_plan_dft_1d(conv_size, x->work, x->work, 
                                       FFTW_BACKWARD, FFTW_PATIENT);
        x->ifft_shift = fftw_plan_dft_1d(size, x->work, x->work, 
                                        FFTW_BACKWARD, FFTW_PATIENT);
        x->fft_shift = fftw_plan_dft_1d(size, x->work, x->work, 
                                       FFTW_FORWARD, FFTW_PATIENT);

        if (!x->signal_fft || !x->kernel_fft || !x->ifft_conv || 
            !x->ifft_shift || !x->fft_shift) {
            pd_error(x, "FFT plan creation failed");
            frft_free(x);  // Clean up on failure
            x->last_size = 0;
            x->conv_size = 0;
            x->saved_chirp_size = 0;
            return 0;
        }

        x->last_size = size;
    }
    return 1;
}

static void frft_bang(t_frft *x)
{
    t_garray *array_in, *array_out;
    t_word *vec_in, *vec_out;
    int size_in, size_out;
    
    // Get arrays and verify them
    if (!(array_in = (t_garray *)pd_findbyclass(x->x_array_in, garray_class))) {
        pd_error(x, "%s: no such array", x->x_array_in->s_name);
        return;
    }
    if (!(array_out = (t_garray *)pd_findbyclass(x->x_array_out, garray_class))) {
        pd_error(x, "%s: no such array", x->x_array_out->s_name);
        return;
    }
    if (!garray_getfloatwords(array_in, &size_in, &vec_in)) {
        pd_error(x, "%s: bad template", x->x_array_in->s_name);
        return;
    }
    if (!garray_getfloatwords(array_out, &size_out, &vec_out)) {
        pd_error(x, "%s: bad template", x->x_array_out->s_name);
        return;
    }

    // Check sizes and warn if output array is too small
    if (size_out < size_in) {
        pd_error(x, "warning: output array %s is smaller than input array (%d < %d)", 
                x->x_array_out->s_name, size_out, size_in);
    }

    // Prepare buffers
    if (!prepare_buffers(x, size_in)) return;

    // Copy input array to buffer
    for (int i = 0; i < size_in; i++) {
        x->buf[i] = vec_in[i].w_float;
    }

    // Process
    frft_process(x, size_in);

    // Copy result to output array (only up to the size of the output array)
    int write_size = (size_out < size_in) ? size_out : size_in;
    for (int i = 0; i < write_size; i++) {
        vec_out[i].w_float = creal(x->buf[i]);
    }

    garray_redraw(array_out);
}

static void frft_list(t_frft *x, t_symbol *s, int argc, t_atom *argv)
{
    (void)s;
    if (argc == 0) return;
    
    // Prepare buffers
    if (!prepare_buffers(x, argc)) return;
    
    // Copy input list to buffer
    for (int i = 0; i < argc; i++) {
        x->buf[i] = atom_getfloat(argv + i);
    }
    
    // Process
    frft_process(x, argc);
    
    // Create output lists
    t_atom *outv_real = (t_atom *)getbytes(sizeof(t_atom) * argc);
    t_atom *outv_imag = (t_atom *)getbytes(sizeof(t_atom) * argc);
    
    if (!outv_real || !outv_imag) {
        pd_error(x, "out of memory");
        if (outv_real) freebytes(outv_real, sizeof(t_atom) * argc);
        if (outv_imag) freebytes(outv_imag, sizeof(t_atom) * argc);
        return;
    }
    
    // Fill output lists
    for (int i = 0; i < argc; i++) {
        SETFLOAT(outv_real + i, creal(x->buf[i]));
        SETFLOAT(outv_imag + i, cimag(x->buf[i]));
    }
    
    // Output results
    outlet_list(x->list_out_imag, &s_list, argc, outv_imag);
    outlet_list(x->list_out_real, &s_list, argc, outv_real);
    
    // Clean up
    freebytes(outv_real, sizeof(t_atom) * argc);
    freebytes(outv_imag, sizeof(t_atom) * argc);
}

static void frft_set(t_frft *x, t_symbol *s1, t_symbol *s2)
{
    x->x_array_in = s1;
    x->x_array_out = s2;
}

static void frft_alpha(t_frft *x, t_floatarg f)
{
    x->x_alpha = f;
}

static void *frft_new(t_symbol *s, int argc, t_atom *argv)
{
    (void)s;
    t_frft *x = (t_frft *)pd_new(frft_class);
    
    // Parse arguments: [array1] [array2] [alpha]
    x->x_array_in = (argc > 0 && argv[0].a_type == A_SYMBOL) ? 
        atom_getsymbol(argv) : gensym("array1");
    x->x_array_out = (argc > 1 && argv[1].a_type == A_SYMBOL) ? 
        atom_getsymbol(argv + 1) : gensym("array2");
    x->x_alpha = (argc > 2 && argv[2].a_type == A_FLOAT) ? 
        atom_getfloat(argv + 2) : 1.0;
    
    // Create outlets for list processing
    x->list_out_real = outlet_new(&x->x_obj, &s_list);
    x->list_out_imag = outlet_new(&x->x_obj, &s_list);
    
    // Initialize buffers
    x->last_size = 0;
    x->conv_size = 0;
    x->saved_chirp_size = 0;
    x->buf = NULL;
    x->work = NULL;
    x->saved_chirp = NULL;
    x->signal_fft = NULL;
    x->kernel_fft = NULL;
    x->ifft_conv = NULL;
    x->ifft_shift = NULL;
    x->fft_shift = NULL;
    
    return (void *)x;
}

static void frft_free(t_frft *x)
{
    // Free FFTW resources
    if (x->buf) freebytes(x->buf, sizeof(fftw_complex) * x->conv_size);
    if (x->work) freebytes(x->work, sizeof(fftw_complex) * x->conv_size);
    if (x->saved_chirp) freebytes(x->saved_chirp, sizeof(fftw_complex) * x->saved_chirp_size);
    if (x->signal_fft) fftw_destroy_plan(x->signal_fft);
    if (x->kernel_fft) fftw_destroy_plan(x->kernel_fft);
    if (x->ifft_shift) fftw_destroy_plan(x->ifft_shift);
    if (x->fft_shift) fftw_destroy_plan(x->fft_shift);
    if (x->ifft_conv) fftw_destroy_plan(x->ifft_conv);
    
    // Outlets are automatically freed by Pd
}

void frft_setup(void)
{
    frft_class = class_new(gensym("frft"),
        (t_newmethod)frft_new,
        (t_method)frft_free,
        sizeof(t_frft),
        CLASS_DEFAULT,
        A_GIMME, 0);
    
    class_addbang(frft_class, frft_bang);
    class_addmethod(frft_class, (t_method)frft_set, 
        gensym("set"), A_SYMBOL, A_SYMBOL, 0);
    class_addmethod(frft_class, (t_method)frft_alpha, 
        gensym("alpha"), A_FLOAT, 0);
    class_addlist(frft_class, frft_list);
}

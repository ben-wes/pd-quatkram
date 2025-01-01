#include "m_pd.h"
#include <complex.h>
#include <fftw3.h>
#include <math.h>
#include <string.h>

#define FFTW_COMPLEX(r,i) ((r) + I * (i))
#define UNUSED(s) ((void)s)


static t_class *frft_class;

typedef struct _frft {
    t_object x_obj;
    t_symbol *x_array_in_real;    // Input array name (real)
    t_symbol *x_array_in_imag;    // Input array name (imaginary)
    t_symbol *x_array_out_real;   // Output array name (real)
    t_symbol *x_array_out_imag;   // Output array name (imaginary)
    t_float x_alpha;              // FRFT power
    
    // FFTW resources
    int last_size;            // Last processed array size
    int conv_size;            // Size of main buffers
    int saved_chirp_size;     // Size of chirp buffer
    fftw_complex *buf;        // Main buffer
    fftw_complex *work;       // Work buffer
    fftw_complex *saved_chirp;// Cached chirp
    fftw_plan signal_fft;     // FFT plans
    fftw_plan kernel_fft;
    fftw_plan ifft_conv;
    fftw_plan ifft_shift;
    fftw_plan fft_shift;
    
    // For list processing
    t_float *imag_buf;        // Store imaginary values from right inlet
    int imag_buf_size;        // Size of stored imaginary values
    t_inlet *in_imag;         // Inlet for imaginary part
    t_outlet *list_out_real;  // Real part output
    t_outlet *list_out_imag;  // Imaginary part output
} t_frft;

// Function declarations
static void frft_free(t_frft *x);
static void frft_bang(t_frft *x);
static void frft_list(t_frft *x, t_symbol *s, int argc, t_atom *argv);
static void *frft_new(t_symbol *s, int argc, t_atom *argv);
static void frft_set(t_frft *x, t_symbol *s1, t_symbol *s2, t_symbol *s3, t_symbol *s4);
static void frft_alpha(t_frft *x, t_floatarg f);
static void frft_process(t_frft *x, int N);
static void frft_imag(t_frft *x, t_symbol *s, int argc, t_atom *argv);

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
    t_garray *array_in_real, *array_in_imag = NULL;
    t_garray *array_out_real, *array_out_imag = NULL;
    t_word *vec_in_real, *vec_in_imag = NULL;
    t_word *vec_out_real, *vec_out_imag = NULL;
    int size_in_real, size_in_imag = 0;
    int size_out_real, size_out_imag = 0;
    
    // Get real input array
    if (!(array_in_real = (t_garray *)pd_findbyclass(x->x_array_in_real, garray_class))) {
        pd_error(x, "%s: no such array", x->x_array_in_real->s_name);
        return;
    }
    if (!garray_getfloatwords(array_in_real, &size_in_real, &vec_in_real)) {
        pd_error(x, "%s: bad template", x->x_array_in_real->s_name);
        return;
    }

    // Get imaginary input array if specified
    if (x->x_array_in_imag && x->x_array_in_imag != &s_) {
        if (!(array_in_imag = (t_garray *)pd_findbyclass(x->x_array_in_imag, garray_class))) {
            pd_error(x, "%s: no such array", x->x_array_in_imag->s_name);
            return;
        }
        if (!garray_getfloatwords(array_in_imag, &size_in_imag, &vec_in_imag)) {
            pd_error(x, "%s: bad template", x->x_array_in_imag->s_name);
            return;
        }
        if (size_in_imag != size_in_real) {
            pd_error(x, "imaginary input array size doesn't match real input array size");
            return;
        }
    }

    // Get real output array
    if (!(array_out_real = (t_garray *)pd_findbyclass(x->x_array_out_real, garray_class))) {
        pd_error(x, "%s: no such array", x->x_array_out_real->s_name);
        return;
    }
    if (!garray_getfloatwords(array_out_real, &size_out_real, &vec_out_real)) {
        pd_error(x, "%s: bad template", x->x_array_out_real->s_name);
        return;
    }

    // Get imaginary output array if specified
    if (x->x_array_out_imag && x->x_array_out_imag != &s_) {
        if (!(array_out_imag = (t_garray *)pd_findbyclass(x->x_array_out_imag, garray_class))) {
            pd_error(x, "%s: no such array", x->x_array_out_imag->s_name);
            return;
        }
        if (!garray_getfloatwords(array_out_imag, &size_out_imag, &vec_out_imag)) {
            pd_error(x, "%s: bad template", x->x_array_out_imag->s_name);
            return;
        }
    }

    // Check output sizes
    if (size_out_real < size_in_real) {
        pd_error(x, "warning: real output array %s is smaller than input array (%d < %d)", 
                x->x_array_out_real->s_name, size_out_real, size_in_real);
    }
    if (array_out_imag && size_out_imag < size_in_real) {
        pd_error(x, "warning: imaginary output array %s is smaller than input array (%d < %d)", 
                x->x_array_out_imag->s_name, size_out_imag, size_in_real);
    }

    // Prepare buffers
    if (!prepare_buffers(x, size_in_real)) return;

    // Copy input array to buffer
    for (int i = 0; i < size_in_real; i++) {
        t_float imag = vec_in_imag ? vec_in_imag[i].w_float : 0.0f;
        x->buf[i] = FFTW_COMPLEX(vec_in_real[i].w_float, imag);
    }

    // Process
    frft_process(x, size_in_real);

    // Copy to output arrays
    int write_size_real = (size_out_real < size_in_real) ? size_out_real : size_in_real;
    for (int i = 0; i < write_size_real; i++) {
        vec_out_real[i].w_float = creal(x->buf[i]);
    }
    
    if (array_out_imag) {
        int write_size_imag = (size_out_imag < size_in_real) ? size_out_imag : size_in_real;
        for (int i = 0; i < write_size_imag; i++) {
            vec_out_imag[i].w_float = cimag(x->buf[i]);
        }
        garray_redraw(array_out_imag);
    }

    garray_redraw(array_out_real);
}

void frft_imag(t_frft *x, t_symbol *s, int argc, t_atom *argv)
{
    UNUSED(s);
    
    // Free old buffer if it exists
    if (x->imag_buf) {
        freebytes(x->imag_buf, x->imag_buf_size * sizeof(t_float));
        x->imag_buf = NULL;
        x->imag_buf_size = 0;  // Reset size immediately after freeing
    }
    
    if (argc == 0) return;  // No new data to store
    
    // Allocate and fill new buffer
    x->imag_buf = (t_float *)getbytes(argc * sizeof(t_float));
    
    if (!x->imag_buf) {
        pd_error(x, "out of memory");
        x->imag_buf_size = 0;  // Make sure size is 0 if allocation failed
        return;
    }
    
    // Set size only after successful allocation
    x->imag_buf_size = argc;
    
    // Copy values
    for (int i = 0; i < argc; i++) {
        x->imag_buf[i] = atom_getfloat(argv + i);
    }
}

static void frft_list(t_frft *x, t_symbol *s, int argc, t_atom *argv)
{
    UNUSED(s);
    if (argc == 0) return;
    
    // Check for size mismatch with stored imaginary values
    if (x->imag_buf) {  // Only check if we have imaginary values stored
        if (x->imag_buf_size != argc) {
            pd_error(x, "size mismatch between real (%d) and imaginary (%d) parts", 
                    argc, x->imag_buf_size);
            // Clear stored imaginary values since they can't be used
            freebytes(x->imag_buf, x->imag_buf_size * sizeof(t_float));
            x->imag_buf = NULL;
            x->imag_buf_size = 0;
            return;
        }
        // Additional safety check for non-null buffer with zero size
        if (x->imag_buf_size == 0) {
            // Cleanup inconsistent state
            freebytes(x->imag_buf, sizeof(t_float));  // Free minimum size
            x->imag_buf = NULL;
        }
    }
        
    // Prepare buffers
    if (!prepare_buffers(x, argc)) return;
    
    // Copy input list to buffer
    for (int i = 0; i < argc; i++) {
        t_float imag = 0.0f;  // Default to 0
        
        // Use stored imaginary values if available
        if (x->imag_buf) {
            imag = x->imag_buf[i];
        }
        
        x->buf[i] = FFTW_COMPLEX(atom_getfloat(argv + i), imag);
    }
    
    // Clear imaginary buffer after use
    if (x->imag_buf) {
        freebytes(x->imag_buf, x->imag_buf_size * sizeof(t_float));
        x->imag_buf = NULL;
        x->imag_buf_size = 0;
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

static void frft_set(t_frft *x, t_symbol *s1, t_symbol *s2, t_symbol *s3, t_symbol *s4)
{
    x->x_array_in_real = s1;
    x->x_array_out_real = s2;
    x->x_array_in_imag = s3;
    x->x_array_out_imag = s4;
}

static void frft_alpha(t_frft *x, t_floatarg f)
{
    x->x_alpha = f;
}

static void *frft_new(t_symbol *s, int argc, t_atom *argv)
{
    UNUSED(s);
    t_frft *x = (t_frft *)pd_new(frft_class);
    
    // Parse arguments: [array_in_real] [array_out_real] [array_in_imag] [array_out_imag] [alpha]
    x->x_array_in_real = (argc > 0 && argv[0].a_type == A_SYMBOL) ? 
        atom_getsymbol(argv) : gensym("array1");
    x->x_array_out_real = (argc > 1 && argv[1].a_type == A_SYMBOL) ? 
        atom_getsymbol(argv + 1) : gensym("array2");
    x->x_array_in_imag = (argc > 2 && argv[2].a_type == A_SYMBOL) ? 
        atom_getsymbol(argv + 2) : &s_;  // Use empty symbol for no array
    x->x_array_out_imag = (argc > 3 && argv[3].a_type == A_SYMBOL) ? 
        atom_getsymbol(argv + 3) : &s_;
    x->x_alpha = (argc > 4 && argv[4].a_type == A_FLOAT) ? 
        atom_getfloat(argv + 4) : 1.0;
    
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
    
    // Initialize list input buffer
    x->imag_buf = NULL;
    x->imag_buf_size = 0;
    
    // Create inlet for imaginary part lists
    x->in_imag = inlet_new(&x->x_obj, &x->x_obj.ob_pd, &s_list, gensym("imag"));
    
    // Create outlets for list processing
    x->list_out_real = outlet_new(&x->x_obj, &s_list);
    x->list_out_imag = outlet_new(&x->x_obj, &s_list);
    
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
    
    // Free list input buffer
    if (x->imag_buf) freebytes(x->imag_buf, x->imag_buf_size * sizeof(t_float));
    
    // Inlets/outlets are automatically freed by Pd
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
        gensym("set"), A_SYMBOL, A_SYMBOL, A_SYMBOL, A_SYMBOL, 0);
    class_addmethod(frft_class, (t_method)frft_alpha, 
        gensym("alpha"), A_FLOAT, 0);
    class_addmethod(frft_class, (t_method)frft_imag,  // Add method for imaginary input
        gensym("imag"), A_GIMME, 0);
    class_addlist(frft_class, frft_list);
}

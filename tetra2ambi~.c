/*
tetra2ambi~ - Tetrahedral microphone to ambisonic converter for Pure Data

2024, Ben Wesch

Functionality:
- Converts 4-channel A-format signal from tetrahedral microphone to B-format ambisonics
- Uses configurable 4x4 transformation matrix for conversion
- Supports real-time matrix parameter adjustment
- Includes preset matrices for common tetrahedral microphones
- Optional high-pass filtering for low-frequency correction

Usage:
1. Connect 4-channel A-format signal to inlet
2. Receive 4-channel B-format signal (W, X, Y, Z) from outlet
3. Use 'matrix' message to set custom conversion matrix
4. Use 'preset' message to load microphone-specific presets

Messages:
- matrix [16 floats]: Set 4x4 conversion matrix (row-major order)
- preset [symbol]: Load preset matrix ("tetramic", "soundfield", "default")
- filter [float]: Set high-pass filter frequency in Hz (0 = off)
- reset: Reset to default matrix

Note: This code was developed with assistance from the Anthropic Claude AI language model.
*/

#include "m_pd.h"
#include <math.h>
#include <string.h>

static t_class *tetra2ambi_tilde_class;

typedef struct _tetra2ambi_tilde {
    t_object x_obj;
    t_sample f_dummy;           // Dummy float for signal inlet
    t_sample **a_format_in;    // 4-channel A-format input
    t_sample **b_format_out;   // 4-channel B-format output
    t_float conversion_matrix[16]; // 4x4 conversion matrix (row-major)
    t_float filter_freq;        // High-pass filter frequency
    t_float filter_coeff;       // Filter coefficient
    t_float filter_state[4];   // Filter state for each channel
    int vec_size;              // Vector size
    t_sample *zero_buffer;      // Zero buffer for unconnected channels
} t_tetra2ambi_tilde;

// Default conversion matrix (tetrahedral to first-order ambisonic)
// This is a standard tetrahedral to B-format conversion matrix
static const t_float default_matrix[16] = {
    0.5f, 0.5f, 0.5f, 0.5f,   // W = (A0+A1+A2+A3)/2 (omnidirectional)
    0.5f, -0.5f, -0.5f, 0.5f, // X = (A0-A1-A2+A3)/2 (front-back)
    0.5f, 0.5f, -0.5f, -0.5f, // Y = (A0+A1-A2-A3)/2 (left-right)
    0.5f, -0.5f, 0.5f, -0.5f  // Z = (A0-A1+A2-A3)/2 (up-down)
};

// Core Sound TetraMic matrix (example - would need actual calibration)
static const t_float tetramic_matrix[16] = {
    0.5f, 0.5f, 0.5f, 0.5f,   // W = (A0+A1+A2+A3)/2
    0.5f, -0.5f, -0.5f, 0.5f, // X = (A0-A1-A2+A3)/2
    0.5f, 0.5f, -0.5f, -0.5f, // Y = (A0+A1-A2-A3)/2
    0.5f, -0.5f, 0.5f, -0.5f  // Z = (A0-A1+A2-A3)/2
};

// SoundField SPS200 matrix (example - would need actual calibration)
static const t_float soundfield_matrix[16] = {
    0.5f, 0.5f, 0.5f, 0.5f,   // W = (A0+A1+A2+A3)/2
    0.5f, -0.5f, -0.5f, 0.5f, // X = (A0-A1-A2+A3)/2
    0.5f, 0.5f, -0.5f, -0.5f, // Y = (A0+A1-A2-A3)/2
    0.5f, -0.5f, 0.5f, -0.5f  // Z = (A0-A1+A2-A3)/2
};

// Calculate filter coefficient from frequency
static t_float calculate_filter_coeff(t_float freq, t_float sample_rate)
{
    if (freq <= 0.0f) return 0.0f;
    t_float omega = 2.0f * M_PI * freq / sample_rate;
    return 1.0f - exp(-omega);
}

// Apply high-pass filter to a sample
static t_float apply_filter(t_float input, t_float *state, t_float coeff)
{
    if (coeff <= 0.0f) return input;
    t_float output = input - *state;
    *state += coeff * output;
    return output;
}

// Perform matrix multiplication: B = M * A
static void matrix_multiply(const t_float *matrix, const t_float *input, t_float *output)
{
    for (int i = 0; i < 4; i++) {
        output[i] = 0.0f;
        for (int j = 0; j < 4; j++) {
            output[i] += matrix[i * 4 + j] * input[j];
        }
    }
}

// Main perform routine
static t_int *tetra2ambi_tilde_perform(t_int *w)
{
    t_tetra2ambi_tilde *x = (t_tetra2ambi_tilde *)(w[1]);
    int n = (int)(w[2]);
    int i, j;
    
    for (i = 0; i < n; i++) {
        // Collect A-format input samples
        t_float a_format[4];
        for (j = 0; j < 4; j++) {
            a_format[j] = x->a_format_in[j][i];
        }
        
        // Apply matrix conversion
        t_float b_format[4];
        matrix_multiply(x->conversion_matrix, a_format, b_format);
        
        // Apply high-pass filter if enabled
        for (j = 0; j < 4; j++) {
            b_format[j] = apply_filter(b_format[j], &x->filter_state[j], x->filter_coeff);
        }
        
        // Output B-format samples
        for (j = 0; j < 4; j++) {
            x->b_format_out[j][i] = b_format[j];
        }
    }
    
    return (w + 3);
}

// DSP setup
static void tetra2ambi_tilde_dsp(t_tetra2ambi_tilde *x, t_signal **sp)
{
    int vec_size = sp[0]->s_n;
    int a_format_chans = sp[0]->s_nchans;
    
    // Update vector size and recalculate filter coefficient
    x->vec_size = vec_size;
    t_float sample_rate = sys_getsr();
    x->filter_coeff = calculate_filter_coeff(x->filter_freq, sample_rate);
    
    // Reallocate zero buffer if necessary
    if (x->zero_buffer) {
        freebytes(x->zero_buffer, x->vec_size * sizeof(t_sample));
    }
    x->zero_buffer = (t_sample *)getbytes(vec_size * sizeof(t_sample));
    if (x->zero_buffer) {
        memset(x->zero_buffer, 0, vec_size * sizeof(t_sample));
    }
    
    // Assign A-format input channels
    for (int i = 0; i < 4; i++) {
        if (i < a_format_chans) {
            x->a_format_in[i] = sp[0]->s_vec + vec_size * i;
        } else {
            x->a_format_in[i] = x->zero_buffer;
        }
    }
    
    // Set up B-format output channels
    signal_setmultiout(&sp[1], 4);
    for (int i = 0; i < 4; i++) {
        x->b_format_out[i] = sp[1]->s_vec + sp[1]->s_n * i;
    }
    
    dsp_add(tetra2ambi_tilde_perform, 2, x, vec_size);
}

// Set conversion matrix
static void tetra2ambi_tilde_matrix(t_tetra2ambi_tilde *x, t_symbol *s, int argc, t_atom *argv)
{
    (void)s; // Unused parameter
    
    if (argc != 16) {
        pd_error(x, "tetra2ambi~: matrix message requires exactly 16 float arguments (4x4 matrix)");
        return;
    }
    
    // Copy matrix values
    for (int i = 0; i < 16; i++) {
        x->conversion_matrix[i] = atom_getfloat(argv + i);
    }
    
    post("tetra2ambi~: conversion matrix updated");
}

// Load preset matrix
static void tetra2ambi_tilde_preset(t_tetra2ambi_tilde *x, t_symbol *s)
{
    const t_float *matrix = NULL;
    const char *name = NULL;
    
    if (s == gensym("tetramic")) {
        matrix = tetramic_matrix;
        name = "TetraMic";
    } else if (s == gensym("soundfield")) {
        matrix = soundfield_matrix;
        name = "SoundField";
    } else if (s == gensym("default")) {
        matrix = default_matrix;
        name = "default";
    } else {
        pd_error(x, "tetra2ambi~: unknown preset '%s'. Available: tetramic, soundfield, default", s->s_name);
        return;
    }
    
    // Copy matrix values
    for (int i = 0; i < 16; i++) {
        x->conversion_matrix[i] = matrix[i];
    }
    
    post("tetra2ambi~: loaded %s preset matrix", name);
}

// Set high-pass filter frequency
static void tetra2ambi_tilde_filter(t_tetra2ambi_tilde *x, t_floatarg freq)
{
    x->filter_freq = (freq < 0.0f) ? 0.0f : freq;
    t_float sample_rate = sys_getsr();
    x->filter_coeff = calculate_filter_coeff(x->filter_freq, sample_rate);
    
    if (x->filter_freq > 0.0f) {
        post("tetra2ambi~: high-pass filter set to %.1f Hz", x->filter_freq);
    } else {
        post("tetra2ambi~: high-pass filter disabled");
    }
}

// Reset to default matrix
static void tetra2ambi_tilde_reset(t_tetra2ambi_tilde *x)
{
    for (int i = 0; i < 16; i++) {
        x->conversion_matrix[i] = default_matrix[i];
    }
    post("tetra2ambi~: reset to default matrix");
}

// Print current matrix
static void tetra2ambi_tilde_print(t_tetra2ambi_tilde *x)
{
    post("tetra2ambi~: current conversion matrix:");
    for (int i = 0; i < 4; i++) {
        post("  [%.3f %.3f %.3f %.3f]", 
             x->conversion_matrix[i*4], x->conversion_matrix[i*4+1],
             x->conversion_matrix[i*4+2], x->conversion_matrix[i*4+3]);
    }
    post("tetra2ambi~: filter frequency: %.1f Hz", x->filter_freq);
}

// Constructor
static void *tetra2ambi_tilde_new(t_floatarg filter_freq)
{
    t_tetra2ambi_tilde *x = (t_tetra2ambi_tilde *)pd_new(tetra2ambi_tilde_class);
    
    // Initialize with default matrix
    for (int i = 0; i < 16; i++) {
        x->conversion_matrix[i] = default_matrix[i];
    }
    
    // Set filter frequency
    x->filter_freq = (filter_freq < 0.0f) ? 0.0f : filter_freq;
    x->filter_coeff = 0.0f; // Will be calculated in dsp
    
    // Initialize filter states
    for (int i = 0; i < 4; i++) {
        x->filter_state[i] = 0.0f;
    }
    
    // Initialize pointers
    x->a_format_in = (t_sample **)getbytes(4 * sizeof(t_sample *));
    x->b_format_out = (t_sample **)getbytes(4 * sizeof(t_sample *));
    x->zero_buffer = NULL;
    x->vec_size = 0;
    
    if (!x->a_format_in || !x->b_format_out) {
        pd_error(x, "tetra2ambi~: out of memory");
        return NULL;
    }
    
    // Create signal outlet
    outlet_new(&x->x_obj, &s_signal);
    
    return (void *)x;
}

// Destructor
static void tetra2ambi_tilde_free(t_tetra2ambi_tilde *x)
{
    if (x->a_format_in) freebytes(x->a_format_in, 4 * sizeof(t_sample *));
    if (x->b_format_out) freebytes(x->b_format_out, 4 * sizeof(t_sample *));
    if (x->zero_buffer) freebytes(x->zero_buffer, x->vec_size * sizeof(t_sample));
}

// Setup function
void tetra2ambi_tilde_setup(void)
{
    tetra2ambi_tilde_class = class_new(gensym("tetra2ambi~"),
        (t_newmethod)tetra2ambi_tilde_new,
        (t_method)tetra2ambi_tilde_free,
        sizeof(t_tetra2ambi_tilde),
        CLASS_MULTICHANNEL,
        A_DEFFLOAT, 0);
    
    class_addmethod(tetra2ambi_tilde_class, (t_method)tetra2ambi_tilde_dsp, gensym("dsp"), A_CANT, 0);
    class_addmethod(tetra2ambi_tilde_class, (t_method)tetra2ambi_tilde_matrix, gensym("matrix"), A_GIMME, 0);
    class_addmethod(tetra2ambi_tilde_class, (t_method)tetra2ambi_tilde_preset, gensym("preset"), A_SYMBOL, 0);
    class_addmethod(tetra2ambi_tilde_class, (t_method)tetra2ambi_tilde_filter, gensym("filter"), A_FLOAT, 0);
    class_addmethod(tetra2ambi_tilde_class, (t_method)tetra2ambi_tilde_reset, gensym("reset"), 0);
    class_addmethod(tetra2ambi_tilde_class, (t_method)tetra2ambi_tilde_print, gensym("print"), 0);
    
    CLASS_MAINSIGNALIN(tetra2ambi_tilde_class, t_tetra2ambi_tilde, f_dummy);
}

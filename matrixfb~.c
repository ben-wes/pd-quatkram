#include "m_pd.h"
#include <math.h>
#include <stdlib.h>  // For rand() and RAND_MAX
#include <time.h>    // For time()

#define MAX_CHANNELS 32
#define CLIP_THRESHOLD 1.0f  // For signal clipping

static t_class *matrixfb_tilde_class;

typedef struct _matrixfb_tilde {
    t_object x_obj;
    t_float f;
    t_float sample_rate;
    int n_channels;
    t_float *state;        // Integrator state for each channel
    t_float *coeffs;       // Feedback matrix coefficients
    t_float leak;          // Leaky integrator coefficient
    t_float *dc_state;     // DC blocking filter state
    t_float dc_coeff;      // DC blocking filter coefficient
    t_float *lpf_state;    // Low-pass filter state
    t_float lpf_freq;      // Low-pass filter cutoff frequency
    t_float lpf_lag;       // Low-pass filter lag time
    t_float gain;          // Output gain
    t_float *prev_output;  // Previous output for rotation
    int debug_block;       // Flag for single block debug output
} t_matrixfb_tilde;

// Utility function for clipping
static inline t_float clip(t_float value, t_float threshold) {
    if (value > threshold) return threshold;
    if (value < -threshold) return -threshold;
    return value;
}

// Utility function to convert linear range -1 to 1 to exponential range min to max
static inline t_float linexp(t_float value, t_float min, t_float max) {
    // Map from -1,1 to 0,1
    t_float norm = (value + 1.0f) * 0.5f;
    // Clip to 0,1 range
    if (norm < 0.0f) norm = 0.0f;
    if (norm > 1.0f) norm = 1.0f;
    // Convert to exponential range
    return min * powf(max/min, norm);
}

// Exponential smoothing (similar to lag in SC)
static inline t_float lag(t_float current, t_float target, t_float coeff) {
    return current + (target - current) * (1.0f - coeff);
}

// Function to apply matrix multiplication and sum
static void apply_matrix_and_sum(t_matrixfb_tilde *x, t_float *channel_outputs) {
    // For each output channel
    for (int i = 0; i < x->n_channels; i++) {
        channel_outputs[i] = 0.0f;
        
        // Apply matrix coefficients from all input channels
        for (int k = 0; k < x->n_channels; k++) {
            // Scale coefficient to -1,1 range as in SC
            t_float coeff = (x->coeffs[i * x->n_channels + k] * 2.0f - 1.0f);
            
            // Multiply by state and accumulate
            channel_outputs[i] += coeff * x->state[k];
        }
        
        // Scale to prevent extreme values
        channel_outputs[i] *= 0.5f;
    }
}

// Function to apply DC blocking filter (high-pass filter)
static void apply_dc_blocking(t_matrixfb_tilde *x, t_float *signals) {
    for (int i = 0; i < x->n_channels; i++) {
        // Simple one-pole high-pass filter
        t_float dc_blocked = signals[i] - x->dc_state[i];
        x->dc_state[i] = x->dc_state[i] * x->dc_coeff + signals[i] * (1.0f - x->dc_coeff);
        signals[i] = dc_blocked;
    }
}

static t_int *matrixfb_tilde_perform(t_int *w)
{
    t_matrixfb_tilde *x = (t_matrixfb_tilde *)(w[1]);
    t_sample *in = (t_sample *)(w[2]);
    t_sample *out = (t_sample *)(w[3]);
    int n = (int)(w[4]);
    int n_inlets = (int)(w[5]);
    
    // Temporary arrays for channel processing
    t_float channel_outputs[MAX_CHANNELS];
    t_float channel_signals[MAX_CHANNELS];
    
    // Debug: Print block info if debug is requested
    if (x->debug_block) {
        post("matrixfb~: processing block of %d samples, %d input channels", n, n_inlets);
    }
    
    // Process each sample in the block
    for (int j = 0; j < n; j++) {
        // 1. Update state with input and feedback
        for (int i = 0; i < x->n_channels; i++) {
            // Get input from corresponding channel, wrapping if needed
            t_float in_sample = in[(i % n_inlets) * n + j];
            
            // Simple feedback - just add the previous output for this channel
            t_float feedback = x->prev_output[i];
            
            // Combine input and feedback
            t_float combined = in_sample + feedback;
            
            // Apply proper SC-style integrator: out(0) = in(0) + (coef * out(-1))
            // where coef is our leak parameter
            x->state[i] = combined + (x->leak * x->state[i]);
            
            // Simple limiting of state to prevent excessive buildup
            if (fabsf(x->state[i]) > 100.0f) {
                x->state[i] = (x->state[i] > 0) ? 100.0f : -100.0f;
            }
        }
        
        // 2. Apply matrix multiplication to generate output for each channel
        apply_matrix_and_sum(x, channel_outputs);
        
        // 3. Apply DC blocking filter (high-pass filter)
        apply_dc_blocking(x, channel_outputs);
        
        // 4. Apply clipping
        for (int i = 0; i < x->n_channels; i++) {
            channel_signals[i] = clip(channel_outputs[i], CLIP_THRESHOLD);
        }
        
        // 5. Apply LPF and generate output
        for (int i = 0; i < x->n_channels; i++) {
            // Simple one-pole low-pass filter
            t_float lpf_coeff = 1.0f - (x->lpf_freq / x->sample_rate) * 2.0f * M_PI;
            if (lpf_coeff < 0.0f) lpf_coeff = 0.0f;
            if (lpf_coeff > 0.999f) lpf_coeff = 0.999f;
            
            x->lpf_state[i] = x->lpf_state[i] * lpf_coeff + channel_signals[i] * (1.0f - lpf_coeff);
            
            // Apply gain
            t_float final_output = x->lpf_state[i] * x->gain;
            
            // Store for feedback
            x->prev_output[i] = final_output;
            
            // Output to corresponding channel
            out[i * n + j] = final_output;
        }
        
        // Debug output if requested
        if (x->debug_block && j < 10) {  // Limit debug to first 10 samples
            for (int i = 0; i < x->n_channels; i++) {
                post("  %d ch%d: in=%f state=%f out=%f matrixOut=%f", 
                    j, i, in[(i % n_inlets) * n + j], x->state[i], x->prev_output[i], channel_outputs[i]);
            }
        }
    }
    
    // Debug: Print final state
    if (x->debug_block) {
        for (int i = 0; i < x->n_channels; i++) {
            post("matrixfb~: channel %d final state: %f, output: %f", 
                i, x->state[i], x->prev_output[i]);
        }
        post("matrixfb~: end of block");
    }
    
    // Reset debug flag after processing the block
    x->debug_block = 0;
    
    return (w + 6);
}

static void matrixfb_tilde_dsp(t_matrixfb_tilde *x, t_signal **sp)
{
    // Store sample rate
    x->sample_rate = sp[0]->s_sr;
    
    // Set output to match our defined number of channels
    signal_setmultiout(&sp[1], x->n_channels);
    
    // Add the perform routine to the DSP chain with number of input channels
    dsp_add(matrixfb_tilde_perform, 5, x, sp[0]->s_vec, sp[1]->s_vec, sp[0]->s_n, sp[0]->s_nchans);
}

static void *matrixfb_tilde_new(t_floatarg n_channels)
{
    t_matrixfb_tilde *x = (t_matrixfb_tilde *)pd_new(matrixfb_tilde_class);
    
    // Initialize number of channels from creation argument
    x->n_channels = (int)n_channels;
    if (x->n_channels < 1) x->n_channels = 1;
    if (x->n_channels > MAX_CHANNELS) x->n_channels = MAX_CHANNELS;
    
    // Create signal outlet
    outlet_new(&x->x_obj, &s_signal);
    
    // Allocate memory for our defined number of channels
    x->state = (t_float *)getbytes(x->n_channels * sizeof(t_float));
    x->coeffs = (t_float *)getbytes(x->n_channels * x->n_channels * sizeof(t_float));
    x->dc_state = (t_float *)getbytes(x->n_channels * sizeof(t_float));
    x->lpf_state = (t_float *)getbytes(x->n_channels * sizeof(t_float));
    x->prev_output = (t_float *)getbytes(x->n_channels * sizeof(t_float));
    
    // Initialize all states to zero
    for (int i = 0; i < x->n_channels; i++) {
        x->state[i] = 0.0f;
        x->dc_state[i] = 0.0f;
        x->lpf_state[i] = 0.0f;
        x->prev_output[i] = 0.0f;  // Zero initial feedback
    }
    
    // Initialize coefficients to default values (0.5 maps to 0 in -1,1 range)
    for (int i = 0; i < x->n_channels * x->n_channels; i++) {
        x->coeffs[i] = 0.5f;
    }
    
    // Set default parameters
    x->leak = 0.99f;               // Leaky integrator coefficient - exact SC default
    x->dc_coeff = 0.995f;          // DC blocking filter coefficient
    x->lpf_freq = 12000.0f;        // Default LPF max frequency
    x->lpf_lag = 0.1f;             // Default LPF lag time
    x->gain = 1.0f;                // Output gain
    x->sample_rate = 44100.0f;     // Default sample rate, will be updated in dsp method
    
    // Initialize debug flag to off
    x->debug_block = 0;
    
    return (void *)x;
}

static void matrixfb_tilde_free(t_matrixfb_tilde *x)
{
    if (x->state) freebytes(x->state, x->n_channels * sizeof(t_float));
    if (x->coeffs) freebytes(x->coeffs, x->n_channels * x->n_channels * sizeof(t_float));
    if (x->dc_state) freebytes(x->dc_state, x->n_channels * sizeof(t_float));
    if (x->lpf_state) freebytes(x->lpf_state, x->n_channels * sizeof(t_float));
    if (x->prev_output) freebytes(x->prev_output, x->n_channels * sizeof(t_float));
}

// Message handler for 'coeffs' message
static void matrixfb_tilde_coeffs(t_matrixfb_tilde *x, t_symbol *s, int argc, t_atom *argv)
{
    (void)s; // Suppress unused parameter warning
    if (argc != x->n_channels * x->n_channels) {
        pd_error(x, "matrixfb~: coeffs message requires %d values", 
                x->n_channels * x->n_channels);
        return;
    }
    
    for (int i = 0; i < argc; i++) {
        x->coeffs[i] = atom_getfloat(&argv[i]);
        // Clamp to 0-1 range
        if (x->coeffs[i] < 0.0f) x->coeffs[i] = 0.0f;
        if (x->coeffs[i] > 1.0f) x->coeffs[i] = 1.0f;
    }
}

// Message handler for 'leak' message
static void matrixfb_tilde_leak(t_matrixfb_tilde *x, t_floatarg f)
{
    x->leak = f;
    if (x->leak < 0.0f) x->leak = 0.0f;
    if (x->leak > 1.0f) x->leak = 1.0f;
}

// Message handler for 'lpf_freq' message (max LPF frequency)
static void matrixfb_tilde_lpf_freq(t_matrixfb_tilde *x, t_floatarg f)
{
    x->lpf_freq = f;
    if (x->lpf_freq < 1.0f) x->lpf_freq = 1.0f;
    if (x->lpf_freq > x->sample_rate / 2.0f) x->lpf_freq = x->sample_rate / 2.0f;
}

// Message handler for 'lpf_lag' message
static void matrixfb_tilde_lpf_lag(t_matrixfb_tilde *x, t_floatarg f)
{
    x->lpf_lag = f;
    if (x->lpf_lag < 0.0f) x->lpf_lag = 0.0f;
    if (x->lpf_lag > 1.0f) x->lpf_lag = 1.0f;
}

// Message handler for 'gain' message
static void matrixfb_tilde_gain(t_matrixfb_tilde *x, t_floatarg f)
{
    x->gain = f;
}

// Message handler for 'reset' message
static void matrixfb_tilde_reset(t_matrixfb_tilde *x, t_floatarg f)
{
    for (int i = 0; i < x->n_channels; i++) {
        x->state[i] = f;
        x->dc_state[i] = 0.0f;
        x->lpf_state[i] = 0.0f;
        x->prev_output[i] = 0.0f;
    }
}

// Message handler for 'random' message with optional seed
static void matrixfb_tilde_random(t_matrixfb_tilde *x, t_symbol *s, int argc, t_atom *argv)
{
    (void)s; // Suppress unused parameter warning
    
    // Get seed from argument or use current time if not provided
    unsigned int seed;
    if (argc > 0) {
        seed = (unsigned int)atom_getfloat(&argv[0]);
    } else {
        seed = (unsigned int)time(NULL);
    }
    
    // Initialize random number generator with seed
    srand(seed);
    
    // Generate random coefficients between 0 and 1
    for (int i = 0; i < x->n_channels * x->n_channels; i++) {
        x->coeffs[i] = (float)rand() / RAND_MAX;
    }
    
    // Post the seed to the console so it can be reproduced
    post("matrixfb~: random coefficients generated with seed %u", seed);
}

// Message handler for 'randomize_state' message
static void matrixfb_tilde_randomize_state(t_matrixfb_tilde *x, t_floatarg scale)
{
    if (scale <= 0.0f) scale = 0.01f;  // Default small scale if not provided
    
    // Set random state values with the given scale
    for (int i = 0; i < x->n_channels; i++) {
        // Create unique random values for each channel
        t_float random_val = ((float)rand() / RAND_MAX * 2.0f - 1.0f) * scale;
        
        // Apply slightly different randomization per channel
        if (i % 2 == 0) {
            x->state[i] = random_val;
        } else {
            x->state[i] = -random_val;  // Invert for odd channels
        }
        
        // Also randomize the lpf and dc states to increase variability
        x->lpf_state[i] = ((float)rand() / RAND_MAX * 2.0f - 1.0f) * scale * 0.5f;
        x->dc_state[i] = ((float)rand() / RAND_MAX * 2.0f - 1.0f) * scale * 0.1f;
        
        // Add some small randomness to prev_output
        x->prev_output[i] = ((float)rand() / RAND_MAX * 2.0f - 1.0f) * scale * 0.2f;
    }
    
    post("matrixfb~: state randomized with scale %f", scale);
}

// Message handler for 'debug' message
static void matrixfb_tilde_debug(t_matrixfb_tilde *x)
{
    x->debug_block = 1;
    post("matrixfb~: debug output enabled for next block");
}

void matrixfb_tilde_setup(void)
{
    matrixfb_tilde_class = class_new(gensym("matrixfb~"),
        (t_newmethod)matrixfb_tilde_new,
        (t_method)matrixfb_tilde_free,
        sizeof(t_matrixfb_tilde),
        CLASS_MULTICHANNEL,
        A_DEFFLOAT, 0);
    
    class_addmethod(matrixfb_tilde_class,
        (t_method)matrixfb_tilde_dsp, gensym("dsp"), A_CANT, 0);
    class_addmethod(matrixfb_tilde_class,
        (t_method)matrixfb_tilde_coeffs, gensym("coeffs"), A_GIMME, 0);
    class_addmethod(matrixfb_tilde_class,
        (t_method)matrixfb_tilde_leak, gensym("leak"), A_FLOAT, 0);
    class_addmethod(matrixfb_tilde_class,
        (t_method)matrixfb_tilde_reset, gensym("reset"), A_FLOAT, 0);
    class_addmethod(matrixfb_tilde_class,
        (t_method)matrixfb_tilde_random, gensym("random"), A_GIMME, 0);
    class_addmethod(matrixfb_tilde_class,
        (t_method)matrixfb_tilde_randomize_state, gensym("randomize_state"), A_DEFFLOAT, 0);
    class_addmethod(matrixfb_tilde_class,
        (t_method)matrixfb_tilde_debug, gensym("debug"), A_NULL, 0);
    class_addmethod(matrixfb_tilde_class,
        (t_method)matrixfb_tilde_lpf_freq, gensym("lpf_freq"), A_FLOAT, 0);
    class_addmethod(matrixfb_tilde_class,
        (t_method)matrixfb_tilde_lpf_lag, gensym("lpf_lag"), A_FLOAT, 0);
    class_addmethod(matrixfb_tilde_class,
        (t_method)matrixfb_tilde_gain, gensym("gain"), A_FLOAT, 0);
    
    CLASS_MAINSIGNALIN(matrixfb_tilde_class, t_matrixfb_tilde, f);
}
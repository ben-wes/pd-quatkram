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
    t_float leak;          // Global leaky integrator coefficient
    t_float *leak_per_ch;  // Per-channel leaky integrator coefficients
    t_float *dc_state_x;   // DC blocking filter x[n-1] state
    t_float *dc_state_y;   // DC blocking filter y[n-1] state
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
            
            // Multiply by state and accumulate (direct from SC)
            channel_outputs[i] += coeff * x->state[k];
        }
        
        // Apply 1000x multiplier as in SC: "snd = snd * 1000 * (coeff.clump(n)*2 - 1);"
        channel_outputs[i] *= 1000.0f;
    }
}

// Function to apply DC blocking filter according to SC's LeakDC implementation
// y[n] = x[n] - x[n-1] + coef * y[n-1]
static void apply_dc_blocking(t_matrixfb_tilde *x, t_float *signals) {
    for (int i = 0; i < x->n_channels; i++) {
        // Apply the exact LeakDC formula from SC
        t_float current_x = signals[i];
        t_float dc_blocked = current_x - x->dc_state_x[i] + x->dc_coeff * x->dc_state_y[i];
        
        // Store states for next time
        x->dc_state_x[i] = current_x;
        x->dc_state_y[i] = dc_blocked;
        
        // Update signal
        signals[i] = dc_blocked;
    }
}

// Equivalent to clip2 in SC - clips signal to +/- threshold
static inline t_float clip2(t_float value, t_float threshold) {
    if (value > threshold) return threshold;
    if (value < -threshold) return -threshold;
    return value;
}

static t_int *matrixfb_tilde_perform(t_int *w)
{
    t_matrixfb_tilde *x = (t_matrixfb_tilde *)(w[1]);
    t_sample *in = (t_sample *)(w[2]);
    t_sample *out = (t_sample *)(w[3]);
    int n = (int)(w[4]);
    int n_inlets = (int)(w[5]);
    
    // Debug: Print block info if debug is requested
    if (x->debug_block) {
        post("matrixfb~: processing block of %d samples, %d input channels", n, n_inlets);
    }
    
    // Array to store the signals at each stage
    t_float inputs[MAX_CHANNELS];
    t_float combined[MAX_CHANNELS];
    t_float matrix_output[MAX_CHANNELS];
    t_float processed_signals[MAX_CHANNELS];
    
    // Process each sample in the block
    for (int j = 0; j < n; j++) {
        // Get inputs from each inlet
        for (int i = 0; i < x->n_channels; i++) {
            // Get input directly from inlet if available, otherwise wrap
            inputs[i] = in[(i % n_inlets) * n + j];
            
            // Combine with feedback
            combined[i] = inputs[i] + x->prev_output[i];
            
            // Apply leaky integrator with PER-CHANNEL coefficient
            t_float leak_coeff = x->leak_per_ch[i]; // Use individual leak coefficients
            x->state[i] = combined[i] + (leak_coeff * x->state[i]);
        }
        
        // Apply matrix multiplication
        for (int i = 0; i < x->n_channels; i++) {
            matrix_output[i] = 0.0f;
            for (int k = 0; k < x->n_channels; k++) {
                t_float coeff = x->coeffs[i * x->n_channels + k] * 2.0f - 1.0f;
                matrix_output[i] += coeff * x->state[k] * 1000.0f;
            }
            
            // Apply DC blocking (LeakDC) to each channel individually
            // This is the SuperCollider LeakDC implementation: y[n] = x[n] - x[n-1] + coef * y[n-1]
            t_float current_x = matrix_output[i];
            t_float dc_blocked = current_x - x->dc_state_x[i] + x->dc_coeff * x->dc_state_y[i];
            x->dc_state_x[i] = current_x;  // Store current input
            x->dc_state_y[i] = dc_blocked; // Store current output
            
            // Apply clipping per channel
            processed_signals[i] = clip2(dc_blocked, 1.0f);
        }
        
        // Apply rotation - shift the ARRAY, not compute indices
        t_float temp = processed_signals[x->n_channels - 1];
        for (int i = x->n_channels - 1; i > 0; i--) {
            processed_signals[i] = processed_signals[i-1];
        }
        processed_signals[0] = temp;
        
        // Apply LPF to each channel and output - match SC's implementation:
        // LPF.ar(snd.rotate(1), snd.linexp(-1,1,1,lpfMaxFreq, 1/1).lag(lpfLag));
        for (int i = 0; i < x->n_channels; i++) {
            // Calculate LPF coefficient based on all processed signals
            // In SC it's based on the signal value itself, so use each channel's own signal
            t_float lpf_freq = linexp(processed_signals[i], 1.0f, x->lpf_freq);
            
            // Apply lag to lpf_freq (simplified lag implementation)
            // SC uses .lag(lpfLag) which is an exponential lag
            t_float lag_coeff = expf(-1.0f / (x->lpf_lag * x->sample_rate));
            x->lpf_state[i] = x->lpf_state[i] * lag_coeff + lpf_freq * (1.0f - lag_coeff);
            lpf_freq = x->lpf_state[i];
            
            // Calculate LPF coefficient based on frequency
            t_float lpf_coeff = expf(-2.0f * M_PI * lpf_freq / x->sample_rate);
            if (lpf_coeff < 0.0f) lpf_coeff = 0.0f;
            if (lpf_coeff > 0.999f) lpf_coeff = 0.999f;
            
            // Apply LPF to signal
            t_float lpf_output = (lpf_coeff * x->lpf_state[i]) + ((1.0f - lpf_coeff) * processed_signals[i]);
            x->lpf_state[i] = lpf_output;
            
            // Apply per-channel gain (in SC each channel gets its own gain value)
            t_float final_output = lpf_output * x->gain;
            
            // Store for feedback next time
            x->prev_output[i] = final_output;
            
            // Output
            out[i * n + j] = final_output;
        }
        
        // Debug output
        if (x->debug_block && j < 10) {
            for (int i = 0; i < x->n_channels; i++) {
                post("  %d ch%d: in=%f combined=%f state=%f matrix=%f processed=%f rotated=%f out=%f leak=%f", 
                    j, i, inputs[i], combined[i], x->state[i], matrix_output[i], 
                    processed_signals[i], processed_signals[i], x->prev_output[i], x->leak_per_ch[i]);
            }
        }
    }
    
    // Debug: Print final state
    if (x->debug_block) {
        for (int i = 0; i < x->n_channels; i++) {
            post("matrixfb~: channel %d final state: %f, output: %f, leak: %f", 
                i, x->state[i], x->prev_output[i], x->leak_per_ch[i]);
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
    x->leak_per_ch = (t_float *)getbytes(x->n_channels * sizeof(t_float));
    x->dc_state_x = (t_float *)getbytes(x->n_channels * sizeof(t_float));
    x->dc_state_y = (t_float *)getbytes(x->n_channels * sizeof(t_float));
    x->lpf_state = (t_float *)getbytes(x->n_channels * sizeof(t_float));
    x->prev_output = (t_float *)getbytes(x->n_channels * sizeof(t_float));
    
    // Initialize states and coefficients
    srand((unsigned int)time(NULL)); // Seed random number generator
    
    for (int i = 0; i < x->n_channels; i++) {
        // Initialize each channel with slightly different values to break symmetry
        x->state[i] = 0.01f * i;
        x->dc_state_x[i] = 0.0f;
        x->dc_state_y[i] = 0.0f;
        x->lpf_state[i] = 0.0f;
        x->prev_output[i] = 0.0f;
        
        // Set different leak coefficient for each channel
        x->leak_per_ch[i] = 0.99f + (0.005f * i);
        if (x->leak_per_ch[i] > 0.999f) x->leak_per_ch[i] = 0.999f;
    }
    
    // Initialize coefficients to random values (0.5 maps to 0 in -1,1 range)
    for (int i = 0; i < x->n_channels * x->n_channels; i++) {
        x->coeffs[i] = 0.5f + ((float)rand() / RAND_MAX - 0.5f) * 0.01f;
    }
    
    // Set default parameters
    x->leak = 0.995f;
    x->dc_coeff = 0.995f;
    x->lpf_freq = 20.0f;
    x->lpf_lag = 0.1f;
    x->gain = 1.0f;
    x->sample_rate = 44100.0f;
    
    // Initialize debug flag to off
    x->debug_block = 0;
    
    return (void *)x;
}

static void matrixfb_tilde_free(t_matrixfb_tilde *x)
{
    if (x->state) freebytes(x->state, x->n_channels * sizeof(t_float));
    if (x->coeffs) freebytes(x->coeffs, x->n_channels * x->n_channels * sizeof(t_float));
    if (x->leak_per_ch) freebytes(x->leak_per_ch, x->n_channels * sizeof(t_float));
    if (x->dc_state_x) freebytes(x->dc_state_x, x->n_channels * sizeof(t_float));
    if (x->dc_state_y) freebytes(x->dc_state_y, x->n_channels * sizeof(t_float));
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
    
    // Also update all per-channel leak coefficients
    for (int i = 0; i < x->n_channels; i++) {
        x->leak_per_ch[i] = x->leak;
    }
}

// Message handler for channel-specific leak message
static void matrixfb_tilde_leak_ch(t_matrixfb_tilde *x, t_symbol *s, int argc, t_atom *argv)
{
    (void)s; // Suppress unused parameter warning
    
    if (argc < 2) {
        pd_error(x, "matrixfb~: leak_ch message requires at least 2 arguments: channel index and value");
        return;
    }
    
    int ch_idx = (int)atom_getfloat(&argv[0]);
    t_float value = atom_getfloat(&argv[1]);
    
    if (ch_idx < 0 || ch_idx >= x->n_channels) {
        pd_error(x, "matrixfb~: invalid channel index %d (valid range: 0-%d)", 
                 ch_idx, x->n_channels-1);
        return;
    }
    
    // Clamp to 0-1 range
    if (value < 0.0f) value = 0.0f;
    if (value > 1.0f) value = 1.0f;
    
    x->leak_per_ch[ch_idx] = value;
    post("matrixfb~: set leak coefficient for channel %d to %f", ch_idx, value);
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
        // Initialize each channel with a slightly different value to break symmetry
        x->state[i] = f + (i * 0.01f);
        x->dc_state_x[i] = 0.0f;
        x->dc_state_y[i] = 0.0f;
        x->lpf_state[i] = 0.0f;
        x->prev_output[i] = 0.0f;
    }
}

// Message handler for 'dc_coeff' message to set DC blocking coefficient
static void matrixfb_tilde_dc_coeff(t_matrixfb_tilde *x, t_floatarg f)
{
    x->dc_coeff = f;
    // More proper bounds for LeakDC coefficient - to remove DC but preserve low frequencies
    // 0.995 to 0.999 is the typical range for LeakDC
    if (x->dc_coeff < 0.9f) x->dc_coeff = 0.9f;
    if (x->dc_coeff > 0.9999f) x->dc_coeff = 0.9999f;
    post("matrixfb~: set DC blocking coefficient to %f", x->dc_coeff);
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
    
    // Also randomize leak coefficients for each channel (between 0.9 and 0.999)
    for (int i = 0; i < x->n_channels; i++) {
        x->leak_per_ch[i] = 0.9f + ((float)rand() / RAND_MAX) * 0.099f;
    }
    
    // Post the seed to the console so it can be reproduced
    post("matrixfb~: random coefficients and leak values generated with seed %u", seed);
}

// Message handler for 'randomize_leak' message
static void matrixfb_tilde_randomize_leak(t_matrixfb_tilde *x, t_floatarg min_val, t_floatarg max_val)
{
    // Default values if not specified
    if (min_val <= 0.0f || min_val > 1.0f) min_val = 0.9f;
    if (max_val <= min_val || max_val > 1.0f) max_val = 0.999f;
    
    // Generate random leak coefficients for each channel
    for (int i = 0; i < x->n_channels; i++) {
        x->leak_per_ch[i] = min_val + ((float)rand() / RAND_MAX) * (max_val - min_val);
    }
    
    post("matrixfb~: randomized leak coefficients between %f and %f", min_val, max_val);
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
        x->dc_state_x[i] = ((float)rand() / RAND_MAX * 2.0f - 1.0f) * scale * 0.1f;
        x->dc_state_y[i] = ((float)rand() / RAND_MAX * 2.0f - 1.0f) * scale * 0.1f;
        
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
        (t_method)matrixfb_tilde_leak_ch, gensym("leak_ch"), A_GIMME, 0);
    class_addmethod(matrixfb_tilde_class,
        (t_method)matrixfb_tilde_reset, gensym("reset"), A_FLOAT, 0);
    class_addmethod(matrixfb_tilde_class,
        (t_method)matrixfb_tilde_random, gensym("random"), A_GIMME, 0);
    class_addmethod(matrixfb_tilde_class,
        (t_method)matrixfb_tilde_randomize_leak, gensym("randomize_leak"), A_DEFFLOAT, A_DEFFLOAT, 0);
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
    class_addmethod(matrixfb_tilde_class,
        (t_method)matrixfb_tilde_dc_coeff, gensym("dc_coeff"), A_FLOAT, 0);
    
    CLASS_MAINSIGNALIN(matrixfb_tilde_class, t_matrixfb_tilde, f);
}
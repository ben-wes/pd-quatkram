/*
noisen~ - Normal distribution noise generator external for Pure Data

2024, Ben Wesch

Functionality:
- Generates normally distributed noise as a signal output
- Uses the Box-Muller transform to generate normally distributed values

Usage:
Send 'seed <value>' message to set a new seed for the random number generator

Note: This code was developed with assistance from the Anthropic Claude AI language model.
*/

#include "m_pd.h"
#include <stdlib.h>
#include <math.h>
#include <time.h>

static t_class *noisen_tilde_class;

typedef struct _noisen_tilde {
    t_object x_obj;
    uint32_t state;
    int have_spare;
    t_float spare;
    t_outlet *x_out;
} t_noisen_tilde;

// Xorshift32 PRNG
static uint32_t xorshift32(uint32_t *state) {
    uint32_t x = *state;
    x ^= x << 13;
    x ^= x >> 17;
    x ^= x << 5;
    return *state = x;
}

// Generate a random float between 0 and 1
static float rand_float(uint32_t *state) {
    return (xorshift32(state) >> 8) / 16777216.0f;
}

static t_int *noisen_tilde_perform(t_int *w) {
    t_noisen_tilde *x = (t_noisen_tilde *)(w[1]);
    t_sample *out = (t_sample *)(w[2]);
    int n = (int)(w[3]);
    
    while (n--) {
        if (x->have_spare) {
            x->have_spare = 0;
            *out++ = x->spare;
        } else {
            float u, v, s;
            do {
                u = 2.0f * rand_float(&x->state) - 1.0f;
                v = 2.0f * rand_float(&x->state) - 1.0f;
                s = u * u + v * v;
            } while (s >= 1.0f || s == 0.0f);
            
            s = sqrtf(-2.0f * logf(s) / s);
            *out++ = u * s;
            x->spare = v * s;
            x->have_spare = 1;
        }
    }
    
    return (w+4);
}

static void noisen_tilde_dsp(t_noisen_tilde *x, t_signal **sp) {
    dsp_add(noisen_tilde_perform, 3, x, sp[0]->s_vec, sp[0]->s_n);
}

static void noisen_tilde_seed(t_noisen_tilde *x, t_floatarg f) {
    uint32_t seed = (uint32_t)f;
    if (seed == 0) {
        // Handle zero seed by using current time
        seed = (uint32_t)time(NULL);
        post("noisen~: Zero seed detected. Using current time as seed: %u", seed);
    }
    x->state = seed;
    x->have_spare = 0;  // Reset the spare value when seed changes
}

static void *noisen_tilde_new(void) {
    t_noisen_tilde *x = (t_noisen_tilde *)pd_new(noisen_tilde_class);
    x->state = (uint32_t)time(NULL);  // Use current time as default seed
    x->have_spare = 0;
    x->x_out = outlet_new(&x->x_obj, &s_signal);
    post("noisen~: Initialized with seed: %u", x->state);
    return (void *)x;
}

void noisen_tilde_setup(void) {
    noisen_tilde_class = class_new(gensym("noisen~"),
        (t_newmethod)noisen_tilde_new,
        0, sizeof(t_noisen_tilde),
        CLASS_DEFAULT,
        0);
    
    class_addmethod(noisen_tilde_class, (t_method)noisen_tilde_dsp, gensym("dsp"), A_CANT, 0);
    class_addmethod(noisen_tilde_class, (t_method)noisen_tilde_seed, gensym("seed"), A_FLOAT, 0);
}
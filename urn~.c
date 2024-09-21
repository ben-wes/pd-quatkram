#include "m_pd.h"
#include <stdlib.h>
#include <stdint.h>

static t_class *urn_tilde_class;

typedef struct _urn_tilde {
    t_object  x_obj;
    t_sample f_dummy;  // dummy float for signal inlet
    int *deck;
    int count;
    int range;
    int held;
    t_outlet *x_out1;
    t_outlet *x_out2;
    uint32_t rng_state;  // State for LCG RNG
} t_urn_tilde;

// Simple Linear Congruential Generator
static uint32_t lcg_random(uint32_t *state) {
    *state = (uint32_t)(((uint64_t)*state * 1103515245 + 12345) & 0x7fffffff);
    return *state;
}

// Generate a random integer between 0 and max-1
static int random_int(t_urn_tilde *x, int max) {
    return (int)(lcg_random(&x->rng_state) % max);
}

static void urn_tilde_reset(t_urn_tilde *x) {
    for (int i = 0; i < x->count; i++) {
        x->deck[i] = i;
    }
    x->range = x->count;
}

static void urn_tilde_seed(t_urn_tilde *x, t_floatarg f) {
    x->rng_state = (uint32_t)((int)f);
}

static void urn_tilde_range(t_urn_tilde *x, t_floatarg f) {
    int new_count = (int)f;
    if (new_count > 0) {
        x->count = new_count;
        x->deck = (int *)realloc(x->deck, sizeof(int) * x->count);
        urn_tilde_reset(x);
    }
}

static t_int *urn_tilde_perform(t_int *w) {
    t_urn_tilde *x = (t_urn_tilde *)(w[1]);
    t_sample *in1 = (t_sample *)(w[2]);
    t_sample *in2 = (t_sample *)(w[3]);
    t_sample *out1 = (t_sample *)(w[4]);
    t_sample *out2 = (t_sample *)(w[5]);
    int n = (int)(w[6]);

    while (n--) {
        float input = *in1++;
        *out2 = 0;  // Reset impulse output

        if (input > 0) {
            if (x->range == 0) {
                urn_tilde_reset(x);
                *out2 = 1;  // Output impulse when deck is reset
            }
            
            int rand_index = random_int(x, x->range);
            x->held = x->deck[rand_index];
            x->deck[rand_index] = x->deck[x->range - 1];
            x->range--;

            if (*in2 > 0) {  // Only reseed when input is strictly greater than 0
                urn_tilde_seed(x, *in2);
            }
        }

        *out1++ = x->held;
        out2++;
        in2++;
    }

    return (w+7);
}

static void urn_tilde_dsp(t_urn_tilde *x, t_signal **sp) {
    dsp_add(urn_tilde_perform, 6, x,
            sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec, sp[3]->s_vec, sp[0]->s_n);
}

static void urn_tilde_free(t_urn_tilde *x) {
    free(x->deck);
    outlet_free(x->x_out1);
    outlet_free(x->x_out2);
}

static void *urn_tilde_new(t_floatarg f) {
    t_urn_tilde *x = (t_urn_tilde *)pd_new(urn_tilde_class);
    
    x->count = f > 0 ? (int)f : 1;
    x->deck = (int *)malloc(sizeof(int) * x->count);
    x->held = 0;
    x->rng_state = 1;  // Default seed
    
    urn_tilde_reset(x);

    inlet_new(&x->x_obj, &x->x_obj.ob_pd, &s_signal, &s_signal);
    x->x_out1 = outlet_new(&x->x_obj, &s_signal);
    x->x_out2 = outlet_new(&x->x_obj, &s_signal);

    return (void *)x;
}

void urn_tilde_setup(void) {
    urn_tilde_class = class_new(gensym("urn~"),
        (t_newmethod)urn_tilde_new,
        (t_method)urn_tilde_free,
        sizeof(t_urn_tilde),
        CLASS_DEFAULT,
        A_DEFFLOAT, 0);

    class_addmethod(urn_tilde_class, (t_method)urn_tilde_dsp, gensym("dsp"), A_CANT, 0);
    class_addmethod(urn_tilde_class, (t_method)urn_tilde_reset, gensym("reset"), 0);
    class_addmethod(urn_tilde_class, (t_method)urn_tilde_seed, gensym("seed"), A_FLOAT, 0);
    class_addmethod(urn_tilde_class, (t_method)urn_tilde_range, gensym("range"), A_FLOAT, 0);
    CLASS_MAINSIGNALIN(urn_tilde_class, t_urn_tilde, f_dummy);
}